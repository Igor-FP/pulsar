#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
from astropy.io import fits

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils  # noqa: E402

# Try to import cosmetic correction helpers from ..\Cosme\cosme.py
COSME_MODULE = None
try:
    cosme_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "Cosme")
    )
    if os.path.isdir(cosme_dir):
        sys.path.insert(0, cosme_dir)
        import cosme as _cosme  # type: ignore  # noqa: E402

        # Ensure required functions exist
        if hasattr(_cosme, "read_cosme_list") and hasattr(
            _cosme, "apply_cosmetic_correction"
        ):
            COSME_MODULE = _cosme
except Exception:
    COSME_MODULE = None

# Tile size (pixels) for K estimation
TILE_SIZE = 64


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  darkopt.py input.fit output.fit master_dark.fit [cosme.lst]\n"
        "\n"
        "    input.fit       - single file OR numbered pattern (e.g. light0001.fit)\n"
        "    output.fit      - single file OR numbered pattern for results\n"
        "    master_dark.fit - master dark frame (single FITS file)\n"
        "    cosme.lst       - optional cosmetic correction list (P x y),\n"
        "                      applied AFTER optimized dark subtraction\n"
        "\n"
        "For each input frame X, this script computes an optimal scale factor K\n"
        "using the darkest suitable tile, then applies:\n"
        "    Y = X - K * DARK\n"
        "and optionally cosmetic correction of hot pixels.\n"
        "\n"
        "Tile selection rules for K estimation:\n"
        f"  - Frame is split into tiles of size {TILE_SIZE}x{TILE_SIZE} "
        "(edge tiles may be smaller).\n"
        "  - Tiles are excluded if they contain ANY pixel <= 0 or == +MAX\n"
        "    (checked for both light and dark, using their native dtypes).\n"
        "  - Among remaining tiles, the one with the lowest mean light level is used.\n"
        "  - If no suitable tiles remain, an error is raised.\n"
        "\n"
        "Result is written in the same data type as input.\n"
        "For integer types, values are clamped to the full dtype range\n"
        "(signed types may have negative values).\n"
    )
    sys.exit(1)


def parse_args(argv):
    if len(argv) not in (4, 5):
        usage()

    input_pattern = argv[1]
    output_pattern = argv[2]
    dark_path = argv[3]
    cosme_path = argv[4] if len(argv) == 5 else None

    return input_pattern, output_pattern, dark_path, cosme_path


def load_master_dark(dark_path):
    """Load master dark frame as float64 array and return (data_f64, original_dtype)."""
    if not os.path.exists(dark_path):
        raise FileNotFoundError(f"Master dark file '{dark_path}' not found.")

    with fits.open(dark_path, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"Master dark '{dark_path}' has no primary image data.")
        if hdul[0].data.ndim != 2:
            raise ValueError(f"Master dark '{dark_path}' is not a 2D image.")
        dark_raw = hdul[0].data
        dark_dtype = dark_raw.dtype
        dark_f = dark_raw.astype(np.float64)

    return dark_f, dark_dtype


def build_invalid_mask(data_f64, orig_dtype):
    """
    Build boolean mask of invalid pixels for tile selection.

    Rules:
      - Integer types: invalid if value <= 0 or value == max(dtype)
      - Float types:   invalid if value <= 0
    """
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        maxv = float(info.max)
        return (data_f64 <= 0.0) | (data_f64 == maxv)
    else:
        return data_f64 <= 0.0


def find_darkest_tile(light_f, dark_f, light_dtype, dark_dtype, tile_size):
    """
    Find coordinates of the darkest valid tile.

    A tile is valid iff:
      - It has no invalid pixels in either light or dark:
          * invalid if <= 0 or == max (for integer types),
          * invalid if <= 0 (for float types).
    Among valid tiles, choose the one with minimal mean(light).

    Returns:
        (y0, y1, x0, x1) or None if no valid tile exists.
    """
    h, w = light_f.shape
    ts = tile_size

    invalid_light = build_invalid_mask(light_f, light_dtype)
    invalid_dark = build_invalid_mask(dark_f, dark_dtype)
    invalid_any = invalid_light | invalid_dark

    tiles_y = (h + ts - 1) // ts
    tiles_x = (w + ts - 1) // ts

    best_mean = None
    best_coords = None

    for ty in range(tiles_y):
        y0 = ty * ts
        y1 = min(y0 + ts, h)
        tile_h = y1 - y0
        if tile_h <= 1:
            continue

        for tx in range(tiles_x):
            x0 = tx * ts
            x1 = min(x0 + ts, w)
            tile_w = x1 - x0
            if tile_w <= 1:
                continue

            tile_invalid = invalid_any[y0:y1, x0:x1]
            if np.any(tile_invalid):
                continue

            lt = light_f[y0:y1, x0:x1]
            mean_l = float(np.mean(lt))

            if (best_mean is None) or (mean_l < best_mean):
                best_mean = mean_l
                best_coords = (y0, y1, x0, x1)

    return best_coords


def estimate_scale_factor(light_f, dark_f, light_dtype, dark_dtype, tile_size):
    """
    Estimate optimal scale factor K for Y = X - K * DARK.

    - Select darkest valid tile.
    - On that tile, compute least squares solution:

        K = sum(L * D) / sum(D^2)

      which minimizes RMS of (L - K*D) on that tile.

    Returns:
        float K

    Raises:
        ValueError if no suitable tile is found or degenerate.
    """
    if light_f.shape != dark_f.shape:
        raise ValueError(
            f"Shape mismatch: light {light_f.shape} vs dark {dark_f.shape}."
        )

    coords = find_darkest_tile(light_f, dark_f, light_dtype, dark_dtype, tile_size)
    if coords is None:
        raise ValueError(
            "No suitable tile found for K estimation "
            "(all candidate tiles contain non-positive or saturated pixels)."
        )

    y0, y1, x0, x1 = coords

    lt = light_f[y0:y1, x0:x1].ravel()
    dt = dark_f[y0:y1, x0:x1].ravel()

    mask = np.isfinite(lt) & np.isfinite(dt) & (dt != 0.0)
    if mask.sum() < 16:
        raise ValueError(
            "Not enough valid pixels in selected tile for K estimation."
        )

    lt = lt[mask]
    dt = dt[mask]

    sum_ld = float(np.dot(lt, dt))
    sum_d2 = float(np.dot(dt, dt))

    if sum_d2 <= 0.0:
        raise ValueError(
            "Degenerate dark signal in selected tile (sum(D^2) <= 0)."
        )

    k = sum_ld / sum_d2
    if k < 0.0:
        k = 0.0

    return k


def apply_dark_optimized(light_data, dark_f, dark_dtype, tile_size):
    """
    Apply optimized dark subtraction:

        result = light - K * dark

    where K is estimated from the darkest valid tile.

    Returns:
        (result_array, K)
      - result has same dtype as light_data
      - integer: clamped to full dtype range (signed may be negative)
      - float: no artificial clamping
    """
    if light_data is None or light_data.ndim != 2:
        raise ValueError("Expected 2D primary image in light frame.")

    light_dtype = light_data.dtype
    light_f = light_data.astype(np.float64, copy=False)

    if light_f.shape != dark_f.shape:
        raise ValueError(
            f"Shape mismatch: light {light_f.shape} vs dark {dark_f.shape}."
        )

    k = estimate_scale_factor(light_f, dark_f, light_dtype, dark_dtype, tile_size)
    work = light_f - k * dark_f

    if np.issubdtype(light_dtype, np.floating):
        if light_dtype == np.float32:
            return work.astype(np.float32), k
        if light_dtype == np.float64:
            return work.astype(np.float64), k
        return work.astype(light_dtype), k

    info = np.iinfo(light_dtype)
    np.clip(work, info.min, info.max, out=work)
    work = np.rint(work)
    return work.astype(light_dtype), k


def process_file(infile, outfile, dark_f, dark_dtype, tile_size, cosme_coords=None):
    """
    Process a single light frame with optimized dark subtraction and
    optional cosmetic correction.

    Returns:
        K used for this frame.
    """
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None or hdul[0].data.ndim != 2:
            raise ValueError(f"File '{infile}' has invalid or missing 2D image data.")

        light = hdul[0].data
        header = hdul[0].header

        if light.shape != dark_f.shape:
            raise ValueError(
                f"Shape mismatch: light {light.shape} vs dark {dark_f.shape} "
                f"for file '{infile}'."
            )

        new_data, k = apply_dark_optimized(light, dark_f, dark_dtype, tile_size)

        # Optional cosmetic correction (after dark subtraction)
        if cosme_coords is not None and COSME_MODULE is not None:
            new_data = COSME_MODULE.apply_cosmetic_correction(new_data, cosme_coords)

        header["HIERARCH DARKSCALE"] = (
            float(k),
            "Dark optimization scale factor K for this frame",
        )

        hdul[0].data = new_data
        hdul[0].header = header
        hdul.writeto(outfile, overwrite=True)

        return k


def main():
    input_pattern, output_pattern, dark_path, cosme_path = parse_args(sys.argv)

    # Build IO pairs
    try:
        io_pairs = batch_utils.build_io_file_lists(input_pattern, output_pattern)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    # Load master dark
    try:
        dark_f, dark_dtype = load_master_dark(dark_path)
    except Exception as e:
        sys.stderr.write(f"Error loading master dark: {e}\n")
        sys.exit(1)

    # Optional cosmetic correction
    cosme_coords = None
    if cosme_path is not None:
        if COSME_MODULE is None:
            sys.stderr.write(
                "Error: cosme.py not found or invalid in '../Cosme', "
                "cannot apply cosmetic correction.\n"
            )
            sys.exit(1)
        try:
            cosme_coords = COSME_MODULE.read_cosme_list(cosme_path)
        except SystemExit:
            # Propagate exits from cosme.read_cosme_list if it calls sys.exit
            raise
        except Exception as e:
            sys.stderr.write(
                f"Error reading cosmetic list '{cosme_path}': {e}\n"
            )
            sys.exit(1)

    total = len(io_pairs)

    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            outdir = os.path.dirname(os.path.abspath(outfile))
            if outdir and not os.path.isdir(outdir):
                os.makedirs(outdir, exist_ok=True)

            k = process_file(
                infile, outfile, dark_f, dark_dtype, TILE_SIZE, cosme_coords
            )

            msg = (
                f"\rProcessed {i} / {total} files | "
                f"K={k:.4f}"
            )
            if cosme_coords is not None:
                msg += " | COSME"
            msg += f" | {os.path.basename(infile)}"
            sys.stderr.write(msg)
            sys.stderr.flush()

        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")
            # Uncomment to abort on first error:
            # sys.exit(1)

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
