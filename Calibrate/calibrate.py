#!/usr/bin/env python3
import sys
import os
import math
import numpy as np
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed

# Paths to shared libraries
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

LIB_DIR = os.path.join(SCRIPT_DIR, "..", "lib")
if os.path.isdir(LIB_DIR) and LIB_DIR not in sys.path:
    sys.path.insert(0, LIB_DIR)

COSME_DIR = os.path.join(SCRIPT_DIR, "..", "Cosme")
if os.path.isdir(COSME_DIR) and COSME_DIR not in sys.path:
    sys.path.insert(0, COSME_DIR)

import batch_utils
from cosme import read_cosme_list, apply_cosmetic_correction

# Tile size for dark optimization (fixed constant)
OPT_TILE_SIZE = 256


def print_usage():
    print(
        "Usage:\n"
        "  calibrate.py rawfiles001.fit outfiles.fit "
        "[-d dark.fit] [-b bias.fit] [-f flat.fit 1000] [-c cosme.lst] [-optimize|-o]\n\n"
        "rawfiles001.fit : input file, numbered sequence or pattern/list handled by batch_utils\n"
        "outfiles.fit    : output base name, numbered sequence or pattern handled by batch_utils\n"
        "-d dark.fit     : MASTER_DARK (required)\n"
        "-b bias.fit     : MASTER_BIAS (optional)\n"
        "-f flat.fit K   : MASTER_FLAT and MULTIPLIER_K (optional)\n"
        "-c cosme.lst    : cosmetic correction list of pixels (optional)\n"
        "-optimize|-o    : per-frame optimal dark multiplier (OPTIMIZ), else OPTIMIZ=1\n"
    )


def parse_args(argv):
    if len(argv) < 3:
        print_usage()
        sys.exit(1)

    in_spec = argv[1]
    out_spec = argv[2]

    dark_path = None
    bias_path = None
    flat_path = None
    flat_k = 1.0
    cosme_path = None
    optimize = False

    i = 3
    while i < len(argv):
        a = argv[i].lower()
        if a == "-d" and i + 1 < len(argv):
            dark_path = argv[i + 1]
            i += 2
        elif a == "-b" and i + 1 < len(argv):
            bias_path = argv[i + 1]
            i += 2
        elif a == "-f" and i + 2 < len(argv):
            flat_path = argv[i + 1]
            try:
                flat_k = float(argv[i + 2])
            except ValueError:
                print("ERROR: Invalid flat multiplier K (must be a number).")
                sys.exit(1)
            i += 3
        elif a == "-c" and i + 1 < len(argv):
            cosme_path = argv[i + 1]
            i += 2
        elif a in ("-optimize", "-o"):
            optimize = True
            i += 1
        else:
            print(f"ERROR: Unknown or incomplete argument: {argv[i]}")
            print_usage()
            sys.exit(1)

    if dark_path is None:
        print("ERROR: Master dark (-d dark.fit) is mandatory.")
        sys.exit(1)

    return in_spec, out_spec, dark_path, bias_path, flat_path, flat_k, cosme_path, optimize


def get_io_pairs(in_spec, out_spec):
    """Build list of (input, output) pairs using batch_utils."""
    try:
        io_pairs = batch_utils.build_io_file_lists(in_spec, out_spec)
    except Exception as e:
        print(f"ERROR: Failed to build IO file lists: {e}")
        sys.exit(1)

    if not io_pairs:
        print("ERROR: No input files to process.")
        sys.exit(1)

    return io_pairs


def load_master_frame(path, name):
    """Load master calibration frame as float64."""
    if not os.path.isfile(path):
        print(f"ERROR: {name} file not found: {path}")
        sys.exit(1)

    with fits.open(path, memmap=False) as hdul:
        if hdul[0].data is None:
            print(f"ERROR: {name} file has no primary image data: {path}")
            sys.exit(1)
        data = hdul[0].data.astype(np.float64, copy=True)

    return data


def load_raw_frame(path):
    """Load raw frame as float64 plus header and original dtype."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    with fits.open(path, memmap=False) as hdul:
        if hdul[0].data is None:
            raise RuntimeError(f"No image data in file: {path}")
        orig_dtype = hdul[0].data.dtype
        header = hdul[0].header.copy()
        data = hdul[0].data.astype(np.float64, copy=True)

    return data, header, orig_dtype


def clamp_to_dtype(data, dtype):
    """Clamp float64 data to full range of given dtype and cast back."""
    if np.issubdtype(dtype, np.integer):
        info = np.iinfo(dtype)
        clipped = np.clip(data, info.min, info.max)
        return clipped.astype(dtype)
    elif np.issubdtype(dtype, np.floating):
        return data.astype(dtype)
    else:
        return data.astype(np.float64)


def compute_optimal_dark_scale(raw, bias, dark, src_dtype):
    """
    Compute optimal scaling factor for master dark per frame using tiles.

    - Uses tiles OPT_TILE_SIZE x OPT_TILE_SIZE.
    - Excludes tiles containing values <= 0 or == MAX (for integer types).
    - For remaining tiles computes k that minimizes sum((R - B - kD)^2).
    - Returns median k over valid tiles.
    """
    tile = OPT_TILE_SIZE
    h, w = raw.shape

    if bias is not None:
        r = raw - bias
    else:
        r = raw.copy()

    d = dark

    max_value = None
    if np.issubdtype(src_dtype, np.integer):
        max_value = np.iinfo(src_dtype).max

    k_values = []

    for y in range(0, h, tile):
        for x in range(0, w, tile):
            sub_r = r[y:y + tile, x:x + tile]
            sub_d = d[y:y + tile, x:x + tile]

            if sub_r.size == 0 or sub_d.size == 0:
                continue

            # Reject tile if invalid pixels present
            if max_value is not None:
                if (
                    np.any(sub_r <= 0) or
                    np.any(sub_d <= 0) or
                    np.any(sub_r == max_value) or
                    np.any(sub_d == max_value)
                ):
                    continue

            denom = np.sum(sub_d * sub_d)
            if denom <= 0:
                continue

            num = np.sum(sub_d * sub_r)
            k = num / denom
            if np.isfinite(k) and k > 0:
                k_values.append(k)

    if not k_values:
        raise RuntimeError(
            "Failed to find suitable tiles for dark optimization: "
            "all tiles contain non-positive or saturated pixels."
        )

    return float(np.median(k_values))


def apply_cosmetic_if_needed(data, cosme_coords):
    """Apply cosmetic correction if cosme list is provided."""
    if cosme_coords is None:
        return data
    return apply_cosmetic_correction(data, cosme_coords)


def print_progress(done, total, last_optimiz=None):
    """Thread-safe progress printing (called only from main thread)."""
    bar_len = 40
    frac = done / total
    filled = int(bar_len * frac)
    bar = "#" * filled + "-" * (bar_len - filled)

    if last_optimiz is not None:
        msg = f"[{bar}] {done}/{total}  OPTIMIZ={last_optimiz:.6f}"
    else:
        msg = f"[{bar}] {done}/{total}"

    sys.stdout.write("\r" + msg)
    sys.stdout.flush()
    if done == total:
        sys.stdout.write("\n")


def process_one(index, in_path, out_path,
                dark, bias, flat, flat_k, cosme_coords,
                optimize,
                dark_path, bias_path, flat_path):
    """
    Worker: calibrate single file.
    No printing from here to keep output thread-safe.
    """
    raw, header, src_dtype = load_raw_frame(in_path)

    if raw.shape != dark.shape:
        raise RuntimeError(
            f"DARK shape {dark.shape} does not match RAW '{in_path}' shape {raw.shape}"
        )
    if bias is not None and bias.shape != raw.shape:
        raise RuntimeError(
            f"BIAS shape {bias.shape} does not match RAW '{in_path}' shape {raw.shape}"
        )
    if flat is not None and flat.shape != raw.shape:
        raise RuntimeError(
            f"FLAT shape {flat.shape} does not match RAW '{in_path}' shape {raw.shape}"
        )

    optimiz = 1.0
    if optimize:
        optimiz = compute_optimal_dark_scale(raw, bias, dark, src_dtype)

    # Apply calibration formula
    work = raw
    if bias is not None:
        work = work - bias

    if optimize:
        work = work - dark * optimiz
    else:
        work = work - dark

    if flat is not None:
        safe_flat = flat.copy()
        safe_flat[safe_flat <= 0] = np.nan
        work = work * flat_k / safe_flat

    # Cosmetic correction after calibration
    work = apply_cosmetic_if_needed(work, cosme_coords)

    # Replace NaNs/Infs by finite numbers (NaN, ±inf -> 0)
    work = np.nan_to_num(work, copy=False)

    # Convert back to original dtype with clamping
    out_data = clamp_to_dtype(work, src_dtype)

    # Update header
    header["HISTORY"] = "Calibrated by calibrate.py"
    header["HISTORY"] = f"  MASTER_DARK = {os.path.basename(dark_path)}"
    if bias is not None and bias_path is not None:
        header["HISTORY"] = f"  MASTER_BIAS = {os.path.basename(bias_path)}"
    if flat is not None and flat_path is not None:
        header["HISTORY"] = f"  MASTER_FLAT = {os.path.basename(flat_path)}"
        header["HISTORY"] = f"  FLAT_K = {flat_k}"
    if optimize:
        header["OPTIMIZ"] = (optimiz, "Optimal master dark multiplier for this frame")
    else:
        header["OPTIMIZ"] = (1.0, "Dark multiplier (no optimization)")
    if cosme_coords is not None:
        header["HISTORY"] = "  COSME list applied"

    # Ensure output directory exists
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fits.writeto(out_path, out_data, header, overwrite=True)

    return index, optimiz


def main(argv=None):
    if argv is None:
        argv = sys.argv

    (
        in_spec,
        out_spec,
        dark_path,
        bias_path,
        flat_path,
        flat_k,
        cosme_path,
        optimize,
    ) = parse_args(argv)

    io_pairs = get_io_pairs(in_spec, out_spec)
    total = len(io_pairs)

    print(f"Found {total} input file(s).")
    if optimize:
        print("Dark optimization enabled (-optimize).")
    else:
        print("Dark optimization disabled (OPTIMIZ=1).")
    if flat_path is not None:
        print(f"Using flat multiplier K = {flat_k}")
    if cosme_path is not None:
        print(f"Using cosmetic list: {cosme_path}")

    # Load masters
    dark = load_master_frame(dark_path, "MASTER_DARK")
    bias = load_master_frame(bias_path, "MASTER_BIAS") if bias_path is not None else None
    flat = load_master_frame(flat_path, "MASTER_FLAT") if flat_path is not None else None

    # Load cosmetic coordinates once (if provided)
    cosme_coords = None
    if cosme_path is not None:
        if not os.path.isfile(cosme_path):
            print(f"WARNING: Cosmetic list file '{cosme_path}' not found. Ignoring.")
        else:
            cosme_coords = read_cosme_list(cosme_path)

    # Determine number of workers: up to (logical_cores - 1), but not more than number of files
    try:
        cpu_count = os.cpu_count() or 1
    except Exception:
        cpu_count = 1

    max_workers = max(1, cpu_count - 1)
    max_workers = min(max_workers, total)
    if max_workers < 1:
        max_workers = 1

    print(f"Using {max_workers} worker thread(s).")

    # Run multithreaded processing
    futures = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for idx, (in_path, out_path) in enumerate(io_pairs, start=1):
            futures.append(
                executor.submit(
                    process_one,
                    idx,
                    in_path,
                    out_path,
                    dark,
                    bias,
                    flat,
                    flat_k,
                    cosme_coords,
                    optimize,
                    dark_path,
                    bias_path,
                    flat_path,
                )
            )

        done = 0
        last_opt = None

        for f in as_completed(futures):
            exc = f.exception()
            if exc is not None:
                # Newline to avoid overwriting progress bar line
                sys.stdout.write("\n")
                print(f"ERROR: {exc}")
                sys.exit(1)

            idx, optimiz = f.result()
            done += 1
            last_opt = optimiz if optimize else None
            print_progress(done, total, last_opt)

    print("Calibration completed.")


if __name__ == "__main__":
    main()
