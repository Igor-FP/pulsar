#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ngain.py - Normalize by gain (multiplication)

Analogous to NGAIN in Iris. Multiplies each image by a coefficient
so that the median becomes equal to the target value.

Formula: result = input * (target_median / current_median)

Supports all FITS data types:
  - Integer: signed/unsigned 8, 16, 32, 64 bit
  - Float: 32, 64 bit
All calculations performed in float64, then converted back to original dtype.
"""

import sys
import os
import numpy as np
from astropy.io import fits

# Add path to shared utilities
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  ngain.py input_spec output_spec target_median\n"
        "\n"
        "    input_spec    - single file OR numbered pattern (e.g. light0001.fit)\n"
        "                    OR wildcard mask (e.g. *.fit) OR @list.txt\n"
        "    output_spec   - single file OR numbered pattern for results\n"
        "    target_median - target median value for all output images\n"
        "\n"
        "Normalizes each image by multiplication:\n"
        "    result = input * (target_median / current_median)\n"
        "\n"
        "All calculations in float64, result converted back to original dtype.\n"
    )
    sys.exit(1)


def parse_args(argv):
    if len(argv) != 4:
        usage()

    input_pattern = argv[1]
    output_pattern = argv[2]

    try:
        target_median = float(argv[3])
    except ValueError:
        sys.stderr.write("Error: target_median must be a number.\n")
        sys.exit(1)

    return input_pattern, output_pattern, target_median


def convert_to_original_dtype(data, orig_dtype):
    """
    Convert float64 data back to original dtype with proper clamping.

    For integer types:
      - Clamp to full dtype range
      - Round to nearest integer (banker's rounding avoided, use floor+0.5)
    For float types:
      - Cast directly, replace NaN/Inf with 0
    """
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        # Round: add 0.5 and floor for positive, subtract 0.5 and ceil for negative
        # Using np.rint which rounds to nearest even - standard for scientific
        # But for minimal loss, we use np.round which is same as rint
        arr = np.clip(data, info.min, info.max)
        arr = np.rint(arr)
        return arr.astype(orig_dtype)

    if np.issubdtype(orig_dtype, np.floating):
        arr = data.astype(orig_dtype)
        # Handle NaN/Inf
        bad = ~np.isfinite(arr)
        if np.any(bad):
            arr[bad] = 0
        return arr

    # Fallback: float32
    arr = data.astype(np.float32)
    bad = ~np.isfinite(arr)
    if np.any(bad):
        arr[bad] = 0
    return arr


def process_file(infile, outfile, target_median):
    """
    Process a single FITS file:
      1. Read image and header
      2. Compute current median
      3. Compute gain coefficient
      4. Multiply image by coefficient
      5. Convert back to original dtype
      6. Write output

    Returns: (current_median, gain_coefficient)
    """
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"File '{infile}' has no primary image data.")
        if hdul[0].data.ndim != 2:
            raise ValueError(f"File '{infile}' is not a 2D image.")

        data = hdul[0].data
        header = hdul[0].header
        orig_dtype = data.dtype

        # Convert to float64 for calculations
        work = data.astype(np.float64)

        # Compute current median (using fast sampled median for large images)
        current_median = batch_utils.fast_median(work)

        # Avoid division by zero
        if current_median == 0.0:
            # If median is zero, cannot normalize by multiplication
            # Just copy the image as-is
            sys.stderr.write(
                f"\nWarning: '{infile}' has zero median, copying unchanged.\n"
            )
            gain = 1.0
        else:
            gain = target_median / current_median
            work = work * gain

        # Convert back to original dtype
        out_data = convert_to_original_dtype(work, orig_dtype)

        # Update header
        header["HISTORY"] = f"Normalized by ngain.py: target_median={target_median}"
        header["HISTORY"] = f"  Original median={current_median:.6f}, gain={gain:.6f}"

        # Write output
        out_dir = os.path.dirname(outfile)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        hdul[0].data = out_data
        hdul[0].header = header
        hdul.writeto(outfile, overwrite=True)

    return current_median, gain


def print_progress(current, total, median_val, gain_val):
    """Print progress bar with current stats."""
    bar_len = 30
    frac = current / total if total > 0 else 1.0
    filled = int(bar_len * frac)
    bar = "#" * filled + "-" * (bar_len - filled)
    msg = f"\r[{bar}] {current}/{total}  median={median_val:.2f} gain={gain_val:.6f}"
    sys.stdout.write(msg)
    sys.stdout.flush()
    if current == total:
        sys.stdout.write("\n")


def main():
    input_pattern, output_pattern, target_median = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_pattern, output_pattern)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    print(f"Normalizing {total} file(s) to target median = {target_median}")

    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            median_val, gain_val = process_file(infile, outfile, target_median)
            print_progress(i, total, median_val, gain_val)
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    print("Done.")


if __name__ == "__main__":
    main()
