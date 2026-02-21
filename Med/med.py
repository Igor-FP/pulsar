#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
med.py - Median combine FITS images.

Usage:
  med.py input_spec output.fits [--memory N]

input_spec:
  firstNNNN.fits   - numbered sequence
  list.txt / .lst  - text file with list of FITS paths
  @list.txt        - same as above
  mask*.fit        - wildcard mask

Options:
  --memory N       - max memory in GB for stack (default: 4.0)

Uses parallel I/O with ThreadPoolExecutor for fast loading.
"""

import argparse
import sys
import os
import time
import warnings
import numpy as np
from astropy.io import fits

# Try to import VerifyWarning from Astropy
try:
    from astropy.io.fits.verify import VerifyWarning
except Exception:
    class VerifyWarning(Warning):
        pass

# Add path to shared utilities
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Median-combine FITS images into a single output FITS.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("input_spec", help="Input: file, list, @list, or wildcard mask")
    parser.add_argument("output", help="Output FITS file")
    parser.add_argument(
        "--memory", "-m",
        type=float,
        default=16.0,
        help="Max memory in GB for image stack (default: 16.0)",
    )
    # Keep --tile for backward compatibility, but it now maps to memory
    parser.add_argument(
        "--tile", "-t",
        type=int,
        default=None,
        help=argparse.SUPPRESS,  # Hidden, for backward compatibility
    )
    return parser.parse_args()


def resolve_inputs(input_spec):
    """Resolve input specification to list of files."""
    try:
        return batch_utils.expand_input_spec(input_spec)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def load_reference(filename):
    """Load header, shape, and dtype from reference file."""
    if not os.path.isfile(filename):
        print(f"Reference file not found: {filename}", file=sys.stderr)
        sys.exit(1)

    try:
        with fits.open(filename, memmap=False) as hdul:
            data = hdul[0].data
            if data is None:
                raise ValueError("No data in primary HDU.")
            if data.ndim < 2:
                raise ValueError("Primary HDU does not contain a 2D image.")
            return hdul[0].header.copy(), data.shape, data.dtype
    except Exception as e:
        print(f"Error reading '{filename}': {e}", file=sys.stderr)
        sys.exit(1)


def print_progress(current, total, message):
    """Simple progress callback."""
    bar_len = 40
    fraction = current / total if total > 0 else 1.0
    filled = int(bar_len * fraction)
    bar = "#" * filled + "-" * (bar_len - filled)
    sys.stdout.write(f"\r{message}: [{bar}] {current}/{total}")
    sys.stdout.flush()
    if current >= total:
        sys.stdout.write("\n")


def safe_add_history(header, text):
    """Append HISTORY lines safely."""
    max_len = 70
    for line in text.split("\n"):
        line = line.strip()
        while len(line) > max_len:
            header.add_history(line[:max_len])
            line = line[max_len:]
        if line:
            header.add_history(line)


def main():
    args = parse_args()

    input_files = resolve_inputs(args.input_spec)
    n_files = len(input_files)

    print(f"Found {n_files} input file(s).")
    if n_files < 2:
        print("Need at least 2 input files for median combine.", file=sys.stderr)
        sys.exit(1)

    ref_header, ref_shape, ref_dtype = load_reference(input_files[0])
    height, width = ref_shape

    is_int = np.issubdtype(ref_dtype, np.integer)
    is_float = np.issubdtype(ref_dtype, np.floating)
    if not (is_int or is_float):
        print(f"Unsupported data type: {ref_dtype}", file=sys.stderr)
        sys.exit(1)

    # Determine memory limit
    max_memory = args.memory
    if args.tile is not None:
        # Backward compatibility: --tile 0 means full-frame (high memory)
        if args.tile == 0:
            max_memory = 64.0  # Use lots of memory
        # --tile N is ignored, we use memory-based decisions now

    print(f"Image size: {width}x{height}, combining {n_files} files")
    print(f"Max memory: {max_memory:.1f} GB")

    start_time = time.time()

    # Use optimized median combine with parallel I/O
    try:
        median = batch_utils.fast_median_combine(
            input_files, ref_shape, ref_dtype,
            max_memory_gb=max_memory,
            progress_callback=print_progress
        )
    except Exception as e:
        print(f"\nError during median combine: {e}", file=sys.stderr)
        sys.exit(1)

    elapsed = time.time() - start_time
    print(f"Median computed in {elapsed:.2f}s")

    # Convert to original dtype
    if is_int:
        median = np.floor(median)
        info = np.iinfo(ref_dtype)
        median = np.clip(median, info.min, info.max)
        out_data = median.astype(ref_dtype)
    else:
        out_data = median.astype(ref_dtype, copy=False)

    # Update header
    safe_add_history(ref_header, f"MEDIAN combine from {n_files} input file(s).")
    safe_add_history(ref_header, "Exact pixel-wise median using np.partition.")
    if is_int:
        safe_add_history(ref_header, "Integer: median = floor((v_low + v_high) / 2).")

    # Write output
    output_path = args.output
    print(f"Writing output: {output_path}")

    hdu = fits.PrimaryHDU(data=out_data, header=ref_header)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=VerifyWarning)
        hdu.writeto(output_path, overwrite=True)

    print(f"Done. Processed {n_files} files into {output_path}")


if __name__ == "__main__":
    main()
