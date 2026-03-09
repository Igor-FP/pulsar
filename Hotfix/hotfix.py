#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
hotfix - Remove single hot pixels from FITS images.

For each pixel, computes mean and std of its 8 neighbors.
If pixel > mean + sigma_threshold * std, replaces it with the neighbor mean.

This catches isolated hot/warm pixels but preserves stars (which always
span multiple pixels due to PSF, seeing, guiding, focus).

Fully vectorized: neighbor stats computed via shifted padded array
using only 2 accumulator buffers (sum, sum_of_squares) for memory efficiency.

Supports 2D (mono) and 3D (color, N×H×W) FITS images.
"""

import sys
import os
import numpy as np
from astropy.io import fits

# Import shared utilities (two levels up from Personal/Hotfix/)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "hotfix - Remove single hot pixels from FITS images.\n"
        "\n"
        "For each pixel, checks if it exceeds mean + N*sigma of its 8 neighbors.\n"
        "If so, replaces with neighbor mean. Stars are preserved (PSF > 1 pixel).\n"
        "\n"
        "Usage:\n"
        "  hotfix.py input_spec output_spec [options]\n"
        "\n"
        "Positional arguments:\n"
        "  input_spec       Input file(s): single file, wildcard (*.fit),\n"
        "                   numbered sequence (img0001.fit), or list (@list.txt)\n"
        "  output_spec      Output: single file, numbered pattern (out0001.fit),\n"
        "                   or directory\n"
        "\n"
        "Options:\n"
        "  --sigma N        Detection threshold in std deviations (default: 5)\n"
        "  --floor N        Minimum noise level in ADU (default: auto = median std)\n"
        "                   In uniform regions neighbor std can be 0 due to integer\n"
        "                   quantization; floor ensures meaningful threshold.\n"
        "  --cold           Also fix cold (dead) pixels below mean - N*sigma\n"
        "  --debug          Print diagnostic dump for the first file\n"
        "\n"
        "Examples:\n"
        "  hotfix.py img0001.fit out0001.fit\n"
        "  hotfix.py *.fit fixed/ --sigma 7\n"
        "  hotfix.py img0001.fit out0001.fit --sigma 4 --cold\n"
        "  hotfix.py img.fit out.fit --floor 30 --sigma 5\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    sigma_threshold = 5.0
    noise_floor = None  # auto
    fix_cold = False
    debug = False

    i = 0
    positional = []
    while i < len(args):
        if args[i] == "--sigma":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --sigma requires a value.\n")
                sys.exit(1)
            try:
                sigma_threshold = float(args[i+1])
            except ValueError:
                sys.stderr.write("Error: --sigma value must be a number.\n")
                sys.exit(1)
            i += 2
        elif args[i] == "--floor":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --floor requires a value.\n")
                sys.exit(1)
            try:
                noise_floor = float(args[i+1])
            except ValueError:
                sys.stderr.write("Error: --floor value must be a number.\n")
                sys.exit(1)
            i += 2
        elif args[i] == "--cold":
            fix_cold = True
            i += 1
        elif args[i] == "--debug":
            debug = True
            i += 1
        else:
            positional.append(args[i])
            i += 1

    if len(positional) < 2:
        usage()

    return positional[0], positional[1], sigma_threshold, noise_floor, fix_cold, debug


# ---------------------------------------------------------------------------
# Hot pixel detection & fix (vectorized)
# ---------------------------------------------------------------------------

# 8-neighbor offsets (dy, dx)
_NEIGHBOR_OFFSETS = [
    (-1, -1), (-1, 0), (-1, 1),
    ( 0, -1),          ( 0, 1),
    ( 1, -1), ( 1, 0), ( 1, 1),
]


def fix_hot_pixels_2d(channel_2d, sigma_threshold, noise_floor, fix_cold,
                      debug_label=None):
    """Fix hot (and optionally cold) pixels in a single 2D plane.

    Returns (fixed_channel, num_fixed).
    """
    h, w = channel_2d.shape
    work = channel_2d.astype(np.float64)

    # Pad by 1 pixel using edge reflection
    padded = np.pad(work, 1, mode='reflect')

    # Accumulate neighbor sum and sum-of-squares (memory-efficient)
    neighbor_sum = np.zeros((h, w), dtype=np.float64)
    neighbor_sq = np.zeros((h, w), dtype=np.float64)

    for dy, dx in _NEIGHBOR_OFFSETS:
        n = padded[1 + dy: 1 + dy + h, 1 + dx: 1 + dx + w]
        neighbor_sum += n
        neighbor_sq += n * n

    mean = neighbor_sum * 0.125
    variance = neighbor_sq * 0.125 - mean * mean
    std = np.sqrt(np.maximum(variance, 0.0))

    # --- DEBUG DUMP ---
    if debug_label:
        zero_std_count = int(np.count_nonzero(std == 0))
        deviation = np.abs(work - mean)
        sys.stderr.write(
            f"\n  [{debug_label}] shape={h}x{w} dtype={channel_2d.dtype}\n"
            f"  pixel  : min={work.min():.1f}  max={work.max():.1f}"
            f"  median={np.median(work[::10, ::10]):.1f}\n"
            f"  std    : min={std.min():.4f}  max={std.max():.1f}"
            f"  median={np.median(std):.4f}  zero_count={zero_std_count}\n"
            f"  |dev|  : min={deviation.min():.4f}  max={deviation.max():.1f}"
            f"  median={np.median(deviation):.4f}\n"
            f"  threshold at sigma={sigma_threshold}:"
            f"  median_std*sigma={np.median(std)*sigma_threshold:.1f}\n"
        )
        # Show top-10 outliers
        dev_flat = deviation.ravel()
        std_flat = std.ravel()
        mean_flat = mean.ravel()
        work_flat = work.ravel()
        top_idx = np.argpartition(dev_flat, -10)[-10:]
        top_idx = top_idx[np.argsort(dev_flat[top_idx])[::-1]]
        sys.stderr.write("  Top-10 outliers:\n")
        for idx in top_idx:
            y, x = divmod(int(idx), w)
            sys.stderr.write(
                f"    [{y},{x}] val={work_flat[idx]:.1f}"
                f"  mean={mean_flat[idx]:.1f}  std={std_flat[idx]:.4f}"
                f"  dev={dev_flat[idx]:.1f}"
                f"  dev/std={dev_flat[idx]/std_flat[idx] if std_flat[idx]>0 else 'inf'}\n"
            )
    # --- END DEBUG ---

    # Noise floor: minimum std to avoid false detections in uniform regions
    # where integer quantization makes neighbor std = 0.
    if noise_floor is None:
        # Auto: median of per-pixel std across the image
        floor = float(np.median(std))
        if floor < 1e-10:
            floor = 1.0
    else:
        floor = noise_floor
    std = np.maximum(std, floor)

    # Detect hot pixels: value > mean + sigma * std
    hot_mask = work > (mean + sigma_threshold * std)

    # Optionally detect cold pixels: value < mean - sigma * std
    if fix_cold:
        cold_mask = work < (mean - sigma_threshold * std)
        bad_mask = hot_mask | cold_mask
    else:
        bad_mask = hot_mask

    if debug_label:
        n_hot = int(np.count_nonzero(hot_mask))
        n_cold = int(np.count_nonzero(cold_mask)) if fix_cold else 0
        sys.stderr.write(
            f"  noise_floor={floor:.4f}"
            f"  hot={n_hot}  cold={n_cold}  total={n_hot+n_cold}\n"
        )

    # Replace bad pixels with neighbor mean
    result = np.where(bad_mask, mean, work)
    num_fixed = int(np.count_nonzero(bad_mask))

    return result, num_fixed


# ---------------------------------------------------------------------------
# dtype conversion
# ---------------------------------------------------------------------------

def convert_to_original_dtype(data, orig_dtype):
    """Convert float64 work array back to original FITS dtype."""
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        arr = np.clip(data, info.min, info.max)
        arr = np.rint(arr)
        return arr.astype(orig_dtype)
    if np.issubdtype(orig_dtype, np.floating):
        arr = data.astype(orig_dtype)
        bad = ~np.isfinite(arr)
        if np.any(bad):
            arr[bad] = 0
        return arr
    arr = data.astype(np.float32)
    bad = ~np.isfinite(arr)
    if np.any(bad):
        arr[bad] = 0
    return arr


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------

def process_file(infile, outfile, sigma_threshold, noise_floor, fix_cold,
                  debug=False):
    """Process a single FITS file. Returns total number of fixed pixels."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")

        data = hdul[0].data
        header = hdul[0].header.copy()
        orig_dtype = data.dtype

    if debug:
        sys.stderr.write(
            f"\nDEBUG: {os.path.basename(infile)}"
            f"  ndim={data.ndim}  shape={data.shape}  dtype={orig_dtype}\n"
        )

    total_fixed = 0
    ch_labels = {0: "R", 1: "G", 2: "B"}

    if data.ndim == 2:
        label = "MONO" if debug else None
        result, n_fixed = fix_hot_pixels_2d(
            data, sigma_threshold, noise_floor, fix_cold, label)
        total_fixed = n_fixed
        out_data = convert_to_original_dtype(result, orig_dtype)
    elif data.ndim == 3:
        planes = []
        for ch in range(data.shape[0]):
            label = ch_labels.get(ch, str(ch)) if debug else None
            fixed_ch, n_fixed = fix_hot_pixels_2d(
                data[ch], sigma_threshold, noise_floor, fix_cold, label)
            planes.append(fixed_ch)
            total_fixed += n_fixed
        result = np.stack(planes, axis=0)
        out_data = convert_to_original_dtype(result, orig_dtype)
    else:
        raise ValueError(f"Unsupported NAXIS={data.ndim}, expected 2 or 3.")

    # Update header
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    cold_str = "+cold" if fix_cold else ""
    header["HISTORY"] = (
        f"Hot pixel fix by hotfix.py: sigma={sigma_threshold}{cold_str}, "
        f"fixed={total_fixed}")

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fits.PrimaryHDU(out_data, header=header).writeto(outfile, overwrite=True)

    return total_fixed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    input_spec, output_spec, sigma_threshold, noise_floor, fix_cold, debug = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    cold_str = "+cold" if fix_cold else ""
    floor_str = f", floor={noise_floor}" if noise_floor is not None else ", floor=auto"
    print(f"Hot pixel fix{cold_str}: {total} file(s), sigma={sigma_threshold}{floor_str}")

    grand_total = 0
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            do_debug = debug and (i == 1)  # debug dump for first file only
            n_fixed = process_file(
                infile, outfile, sigma_threshold, noise_floor, fix_cold, do_debug)
            grand_total += n_fixed
            sys.stderr.write(f"\r{i}/{total}  fixed={n_fixed}")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write(f"\nDone. Total pixels fixed: {grand_total}\n")


if __name__ == "__main__":
    main()
