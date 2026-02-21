#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from astropy.io import fits
import numpy as np

# Number of hottest pixels to record
TOP_N_PIXELS = 10000


def find_hottest_pixels(data, n):
    # Ensure 2D (no cubes etc.)
    if data is None:
        raise ValueError("No image data in primary HDU.")
    if data.ndim != 2:
        raise ValueError(f"Expected 2D image, got {data.ndim}D.")

    # Work on a flat view of the data
    flat = data.reshape(-1)

    # Handle NaN / inf for floating-point images: treat them as very low values
    if np.issubdtype(flat.dtype, np.floating):
        flat = flat.copy()
        bad = ~np.isfinite(flat)
        if np.any(bad):
            flat[bad] = -np.inf

    total_pixels = flat.size
    if total_pixels == 0:
        return []

    n = min(n, total_pixels)
    if n <= 0:
        return []

    # Use argpartition for efficiency, then sort selected indices by value desc
    idx_part = np.argpartition(flat, -n)[-n:]
    idx_sorted = idx_part[np.argsort(flat[idx_part])[::-1]]

    width = data.shape[1]

    coords = []
    for idx in idx_sorted:
        y, x = divmod(int(idx), width)
        # Coordinates are zero-based (Iris/ISIS COSME.LST style)
        coords.append((x, y))

    return coords


def write_cosme_lst(coords, output_path):
    # Write CRLF line endings explicitly
    with open(output_path, "w", newline="") as f:
        for x, y in coords:
            f.write(f"P {x} {y}\r\n")


def main():
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage:\n"
            "  make_cosmetic.py input.fits cosme.lst\n"
        )
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    try:
        # memmap=False to avoid issues with BZERO/BSCALE/BLANK when mmap'd
        with fits.open(input_path, memmap=False) as hdul:
            data = hdul[0].data

        coords = find_hottest_pixels(data, TOP_N_PIXELS)
        write_cosme_lst(coords, output_path)

    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
