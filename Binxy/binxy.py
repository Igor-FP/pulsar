#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
binxy - Software 2x2 or 4x4 pixel binning for FITS images.

Reduces image dimensions by 2 (or 4) in both axes.
  - Integer types: sums 2x2 blocks; promotes output dtype only if actual
    values exceed the original range (uint8->uint16->uint32->uint64,
    int8->int16->int32->int64).
  - Integer + --keep_bitness: sums then integer-divides by 4, keeping
    original dtype (trades quantization noise for smaller file size).
  - Float types: averages 2x2 blocks (sum * 0.25) to preserve value range.

Supports 2D (H, W) and 3D (N, H, W) FITS images.
4x4 binning is implemented as two passes of 2x2.
"""

import sys
import os
import numpy as np
from astropy.io import fits

# Import shared utilities (two levels up from Personal/Binxy/)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "binxy - Software pixel binning for FITS images.\n"
        "\n"
        "Reduces image dimensions by summing adjacent 2x2 pixel blocks.\n"
        "Supports 2D (H x W) and 3D (N x H x W) FITS images.\n"
        "\n"
        "Usage:\n"
        "  binxy.py input_spec output_spec --2|--4 [--keep_bitness]\n"
        "\n"
        "Positional arguments:\n"
        "  input_spec       Input file(s): single file, wildcard (*.fit),\n"
        "                   numbered sequence (img0001.fit), or list (@list.txt)\n"
        "  output_spec      Output: single file, numbered pattern (out0001.fit),\n"
        "                   or directory\n"
        "\n"
        "Required (one of):\n"
        "  --2              2x2 binning (reduce each axis by 2)\n"
        "  --4              4x4 binning (reduce each axis by 4, two 2x2 passes)\n"
        "\n"
        "Options:\n"
        "  --keep_bitness   Keep original dtype by integer-dividing sum by 4.\n"
        "                   Without this flag, integer sums may promote dtype\n"
        "                   (e.g. uint16 -> uint32) if values exceed range.\n"
        "                   Float images always average (x0.25) regardless.\n"
        "\n"
        "Behavior:\n"
        "  Integer images:  sum 2x2 block -> promote dtype if needed\n"
        "  + keep_bitness:  sum 2x2 block -> //4 -> keep original dtype\n"
        "  Float images:    sum 2x2 block -> x0.25 (average), dtype unchanged\n"
        "\n"
        "  Odd dimensions are cropped to even before binning.\n"
        "  Header fields XPIXSZ, YPIXSZ, XBINNING, YBINNING are updated.\n"
        "\n"
        "Examples:\n"
        "  binxy.py img0001.fit out0001.fit --2\n"
        "  binxy.py *.fit binned/ --4\n"
        "  binxy.py img0001.fit out0001.fit --2 --keep_bitness\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    factor = None
    keep_bitness = False
    positional = []
    for a in args:
        if a == "--2":
            factor = 2
        elif a == "--4":
            factor = 4
        elif a == "--keep_bitness":
            keep_bitness = True
        else:
            positional.append(a)

    if factor is None or len(positional) < 2:
        usage()

    if len(positional) < 2:
        usage()

    return positional[0], positional[1], factor, keep_bitness


# ---------------------------------------------------------------------------
# dtype promotion
# ---------------------------------------------------------------------------

# Promotion chains: each dtype maps to the next wider type in the same family
_UNSIGNED_CHAIN = [np.uint8, np.uint16, np.uint32, np.uint64]
_SIGNED_CHAIN = [np.int8, np.int16, np.int32, np.int64]


def _find_min_dtype(orig_dtype, actual_min, actual_max):
    """Find the smallest integer dtype in the same family that fits actual values."""
    if np.issubdtype(orig_dtype, np.unsignedinteger):
        chain = _UNSIGNED_CHAIN
    else:
        chain = _SIGNED_CHAIN

    # Start from original dtype position in the chain
    start = 0
    for i, dt in enumerate(chain):
        if np.dtype(dt) == np.dtype(orig_dtype):
            start = i
            break

    for dt in chain[start:]:
        info = np.iinfo(dt)
        if actual_min >= info.min and actual_max <= info.max:
            return dt

    # Fallback: largest in chain
    return chain[-1]


# ---------------------------------------------------------------------------
# 2x2 binning core
# ---------------------------------------------------------------------------

def bin2x_2d(data_2d, keep_bitness=False):
    """Bin a single 2D plane by 2x2. Returns (binned, out_dtype)."""
    h, w = data_2d.shape
    # Crop to even dimensions
    h_even = h - (h % 2)
    w_even = w - (w % 2)
    cropped = data_2d[:h_even, :w_even]

    orig_dtype = data_2d.dtype

    if np.issubdtype(orig_dtype, np.floating):
        # Float: always average
        work = cropped.astype(np.float64)
        work = work.reshape(h_even // 2, 2, w_even // 2, 2)
        binned = work.sum(axis=(1, 3)) * 0.25
        return binned.astype(orig_dtype), orig_dtype
    else:
        # Integer: sum in 64-bit
        if np.issubdtype(orig_dtype, np.unsignedinteger):
            work_dt = np.uint64
        else:
            work_dt = np.int64
        work = cropped.astype(work_dt)
        work = work.reshape(h_even // 2, 2, w_even // 2, 2)
        binned = work.sum(axis=(1, 3))

        if keep_bitness:
            # Integer divide by 4 to stay in original dtype range
            binned = binned // 4
            return binned.astype(orig_dtype), orig_dtype
        else:
            # Find minimum dtype that fits actual values
            out_dtype = _find_min_dtype(orig_dtype, int(binned.min()), int(binned.max()))
            return binned.astype(out_dtype), out_dtype


def bin2x(data, keep_bitness=False):
    """Bin 2D or 3D array by 2x2. Returns (binned_data, out_dtype)."""
    if data.ndim == 2:
        return bin2x_2d(data, keep_bitness)
    elif data.ndim == 3:
        planes = []
        out_dtype = None
        for ch in range(data.shape[0]):
            binned_ch, dt = bin2x_2d(data[ch], keep_bitness)
            planes.append(binned_ch)
            # Use widest dtype across channels
            if out_dtype is None:
                out_dtype = dt
            elif np.dtype(dt).itemsize > np.dtype(out_dtype).itemsize:
                out_dtype = dt
        # Ensure all channels have the same dtype
        result = np.stack([p.astype(out_dtype) for p in planes], axis=0)
        return result, out_dtype
    else:
        raise ValueError(f"Unsupported NAXIS={data.ndim}, expected 2 or 3.")


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------

def process_file(infile, outfile, factor, keep_bitness):
    """Process a single FITS file. Returns (orig_shape, out_shape, orig_dtype, out_dtype)."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")

        data = hdul[0].data
        header = hdul[0].header.copy()
        orig_dtype = data.dtype
        orig_shape = data.shape

    # First pass: 2x2
    result, out_dtype = bin2x(data, keep_bitness)

    # Second pass for 4x4
    if factor == 4:
        result, out_dtype = bin2x(result, keep_bitness)

    out_shape = result.shape

    # Update header
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    bin_factor = factor
    if "XPIXSZ" in header:
        header["XPIXSZ"] = header["XPIXSZ"] * bin_factor
    if "YPIXSZ" in header:
        header["YPIXSZ"] = header["YPIXSZ"] * bin_factor
    if "XBINNING" in header:
        header["XBINNING"] = header["XBINNING"] * bin_factor
    if "YBINNING" in header:
        header["YBINNING"] = header["YBINNING"] * bin_factor

    dtype_str = str(np.dtype(out_dtype))
    header["HISTORY"] = f"Binned {bin_factor}x{bin_factor} by binxy.py"
    header["HISTORY"] = (
        f"  {orig_shape[-1]}x{orig_shape[-2]} -> {out_shape[-1]}x{out_shape[-2]}"
        f", dtype {orig_dtype} -> {dtype_str}")

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fits.PrimaryHDU(result, header=header).writeto(outfile, overwrite=True)

    return orig_shape, out_shape, orig_dtype, np.dtype(out_dtype)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    input_spec, output_spec, factor, keep_bitness = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    kb_str = ", keep_bitness" if keep_bitness else ""
    print(f"Binning {total} file(s), factor={factor}x{factor}{kb_str}")

    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            orig_shape, out_shape, orig_dt, out_dt = process_file(infile, outfile, factor, keep_bitness)
            promoted = " PROMOTED" if out_dt != orig_dt else ""
            sys.stderr.write(
                f"\r{i}/{total}  "
                f"{orig_shape[-1]}x{orig_shape[-2]} -> {out_shape[-1]}x{out_shape[-2]}"
                f"  {orig_dt}->{out_dt}{promoted}")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\nDone.\n")


if __name__ == "__main__":
    main()
