#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from astropy.io import fits
import numpy as np

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


def parse_args(argv):
    """Parse command line arguments."""
    test_mode = False

    # argv: [script, in, out, cosme, (-t)?]
    if len(argv) == 5 and argv[4] == "-t":
        test_mode = True
    elif len(argv) != 4:
        sys.stderr.write(
            "Usage:\n"
            "  cosme.py input.fit output.fit cosme.lst [-t]\n"
            "\n"
            "    input.fit  - single file OR numbered pattern (e.g. light0001.fit)\n"
            "    output.fit - single file OR numbered pattern for results\n"
            "    cosme.lst  - list of hot pixels: lines 'P x y'\n"
            "    -t         - test mode: create mask image instead of correction\n"
        )
        sys.exit(1)

    input_pattern = argv[1]
    output_pattern = argv[2]
    cosme_path = argv[3]
    return input_pattern, output_pattern, cosme_path, test_mode


def read_cosme_list(path):
    """Read cosmetic pixel list: lines 'P x y'."""
    if not os.path.isfile(path):
        sys.stderr.write(f"Error: cosme list file '{path}' not found.\n")
        sys.exit(1)

    coords = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            if parts[0].upper() != "P":
                continue
            try:
                x = int(parts[1])
                y = int(parts[2])
            except ValueError:
                continue
            coords.append((x, y))

    if not coords:
        sys.stderr.write("Warning: no valid hot pixels found in cosme list.\n")

    return coords


def generate_hot_mask(shape, coords):
    """Build boolean mask of hot pixels from coordinate list."""
    h, w = shape
    mask = np.zeros((h, w), dtype=bool)
    for x, y in coords:
        if 0 <= x < w and 0 <= y < h:
            mask[y, x] = True
    return mask


def apply_cosmetic_correction(data, coords):
    """
    Replace listed pixels by the mean of valid neighboring pixels.

    - For floating types: use mean directly.
    - For integer types: round to nearest integer and clamp to dtype range.
    - If no valid neighbors found for a pixel, keep original value.
    """
    if data is None or data.ndim != 2:
        raise ValueError("Expected 2D primary image for cosmetic correction.")

    h, w = data.shape
    hot_mask = generate_hot_mask((h, w), coords)
    out = data.copy()

    is_float = np.issubdtype(out.dtype, np.floating)
    if not is_float:
        dtype_info = np.iinfo(out.dtype)

    src = data  # original data for neighbor sampling

    for x, y in coords:
        if not (0 <= x < w and 0 <= y < h):
            continue
        if not hot_mask[y, x]:
            # Coordinate outside mask (e.g. duplicates or invalid) -> skip
            continue

        s = 0.0
        cnt = 0

        # 3x3 neighborhood around (x, y), excluding the hot pixel itself
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dx == 0 and dy == 0:
                    continue
                ny = y + dy
                nx = x + dx
                if 0 <= nx < w and 0 <= ny < h and not hot_mask[ny, nx]:
                    v = src[ny, nx]
                    if is_float and not np.isfinite(v):
                        continue
                    s += float(v)
                    cnt += 1

        if cnt > 0:
            mean_val = s / cnt
            if is_float:
                out[y, x] = mean_val
            else:
                val = int(round(mean_val))
                if val < dtype_info.min:
                    val = dtype_info.min
                elif val > dtype_info.max:
                    val = dtype_info.max
                out[y, x] = val
        # If cnt == 0: keep the original (potentially hot) value

    return out


def generate_test_image(shape, dtype, coords):
    """
    Generate mask image for test mode:
    - All zeros, hot pixels set to:
        1.0 for float types,
        max for integer types.
    """
    if len(shape) != 2:
        raise ValueError("Expected 2D primary image for test mode.")

    h, w = shape
    img = np.zeros((h, w), dtype=dtype)

    is_float = np.issubdtype(dtype, np.floating)
    if is_float:
        white = 1.0
    else:
        dtype_info = np.iinfo(dtype)
        white = dtype_info.max

    for x, y in coords:
        if 0 <= x < w and 0 <= y < h:
            img[y, x] = white

    return img


def process_file(infile, outfile, coords, test_mode):
    """Apply cosmetic correction or generate test mask for a single file."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"File '{infile}' has no primary image data.")
        if hdul[0].data.ndim != 2:
            raise ValueError(f"File '{infile}' is not a 2D image.")

        data = hdul[0].data
        header = hdul[0].header

        if test_mode:
            new_data = generate_test_image(data.shape, data.dtype, coords)
        else:
            new_data = apply_cosmetic_correction(data, coords)

        hdul[0].data = new_data
        hdul[0].header = header

        hdul.writeto(outfile, overwrite=True)


def main():
    input_pattern, output_pattern, cosme_path, test_mode = parse_args(sys.argv)
    coords = read_cosme_list(cosme_path)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_pattern, output_pattern)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)

    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            process_file(infile, outfile, coords, test_mode)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
