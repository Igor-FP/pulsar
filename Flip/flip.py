#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
flip - Mirror FITS images along the X and/or Y axis.

Terminology (READ THIS - "flip along an axis" is ambiguous in everyday speech
and people guess the wrong axis about half the time). Here it means REVERSING
THAT COORDINATE, defined strictly by the pixel formula:

  --y  flip along Y = reverse the Y (row) coordinate:   Ynew = H - Yold - 1
         => TOP row <-> BOTTOM row (a vertical flip; image turned upside down).
         It does NOT swap left/right.
  --x  flip along X = reverse the X (column) coordinate: Xnew = W - Xold - 1
         => LEFT column <-> RIGHT column (a horizontal / left-right mirror).
         It does NOT swap top/bottom.

So "flip along Y" reverses the Y coordinate (top<->bottom); it is NOT a
reflection across the Y-axis line (that would be left<->right, which is --x).
When unsure, trust the coordinate formula, not the words.

Both may be combined (--x --y = a 180-degree rotation). Pixel values and dtype
are preserved exactly (a flip is a lossless geometric remap). WCS keywords are
NOT adjusted - the flip is applied to pixels only.

Supports 2D (mono) and 3-channel colour ((3,H,W) or (H,W,3)); only the spatial
axes are flipped, the channel axis is left alone.
"""

import sys
import os
import numpy as np
from astropy.io import fits

# Import shared utilities
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "flip - Mirror FITS images along the X and/or Y axis.\n"
        "\n"
        "Usage:\n"
        "  flip.py input_spec output_spec (--x | --y | --x --y)\n"
        "\n"
        "  input_spec   - single file, wildcard (*.fit), numbered (img0001.fit),\n"
        "                 or @list.txt\n"
        "  output_spec  - single file, numbered pattern, or directory\n"
        "\n"
"Axis (at least one required; may be combined). \"Flip along an axis\" is\n"
        "ambiguous in plain speech, so it is defined here strictly by the pixel\n"
        "coordinate formula - reversing that coordinate:\n"
        "  --y   flip along Y = reverse the Y (row) coordinate:  Ynew = H - Yold - 1\n"
        "        -> TOP <-> BOTTOM (vertical flip, image upside down). NOT left/right.\n"
        "  --x   flip along X = reverse the X (column) coordinate: Xnew = W - Xold - 1\n"
        "        -> LEFT <-> RIGHT (horizontal / left-right mirror). NOT top/bottom.\n"
        "\n"
        "Note: \"flip along Y\" reverses the Y coordinate (top<->bottom); it is NOT a\n"
        "reflection across the Y-axis line (that is left<->right, i.e. --x). When in\n"
        "doubt, trust the coordinate formula, not the words. This is a pixel-index\n"
        "reversal only: pixel values and dtype are preserved, and WCS keywords are\n"
        "NOT adjusted. Colour images (3xHxW or HxWx3) flip on the spatial axes only.\n"
        "\n"
        "Examples:\n"
        "  flip.py in.fit out.fit --y          (vertical flip: top <-> bottom)\n"
        "  flip.py *.fit flipped/ --x          (horizontal flip, batch)\n"
        "  flip.py in.fit out.fit --x --y      (both axes = 180-degree rotation)\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]
    flip_x = False
    flip_y = False
    positional = []

    i = 0
    while i < len(args):
        a = args[i]
        if a in ("-h", "--help"):
            usage()
        if a == "--x":
            flip_x = True
            i += 1
            continue
        if a == "--y":
            flip_y = True
            i += 1
            continue
        if a.startswith("--"):
            sys.stderr.write(f"Error: unknown option: {a}\n")
            usage()
        positional.append(a)
        i += 1

    if len(positional) < 2:
        usage()
    if not (flip_x or flip_y):
        sys.stderr.write("Error: specify at least one axis to flip: --x and/or --y.\n")
        usage()

    return positional[0], positional[1], flip_x, flip_y


# ---------------------------------------------------------------------------
# Processing
# ---------------------------------------------------------------------------

def _spatial_axes(data):
    """Return (y_axis, x_axis) numpy axis indices for the image layout.

    2D -> (0, 1); (3,H,W) -> (1, 2); (H,W,3) -> (0, 1)."""
    if data.ndim == 2:                 # (H, W)
        return 0, 1
    if data.ndim == 3:
        if data.shape[0] == 3:         # (C, H, W) - PULSAR convention
            return 1, 2
        if data.shape[2] == 3:         # (H, W, C)
            return 0, 1
    raise ValueError(f"unsupported image shape {data.shape}")


def flip_data(data, flip_x, flip_y):
    """Reverse the requested spatial axes. Values and dtype are preserved.

    flip_y: Ynew = H - Yold - 1 (reverse the row axis).
    flip_x: Xnew = W - Xold - 1 (reverse the column axis)."""
    y_ax, x_ax = _spatial_axes(data)
    out = data
    if flip_y:
        out = np.flip(out, axis=y_ax)
    if flip_x:
        out = np.flip(out, axis=x_ax)
    # np.flip returns a reversed view; realize it so writeto stores it plainly.
    return np.ascontiguousarray(out)


def process_file(infile, outfile, flip_x, flip_y):
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("no image data")
        data = hdul[0].data
        header = hdul[0].header.copy()

    # The input file is closed here, BEFORE writing, so an in-place output
    # (output == input, e.g. `flip *.fit .`) does not overwrite a file whose
    # handle is still open (which fails on Windows with a PermissionError).
    out = flip_data(data, flip_x, flip_y)

    # Drop stale scaling; astropy re-derives BZERO/BSCALE from the output dtype
    # (e.g. uint16 -> BZERO=32768), so every integer/float type round-trips.
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    axes = []
    if flip_x:
        axes.append("X")
    if flip_y:
        axes.append("Y")
    header["HISTORY"] = (
        "flip.py: flipped " + "+".join(axes)
        + " (pixel-index reversal; WCS not adjusted)")

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    fits.PrimaryHDU(out, header=header).writeto(outfile, overwrite=True)


def main():
    input_spec, output_spec, flip_x, flip_y = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no input files.\n")
        sys.exit(1)

    total = len(io_pairs)
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            process_file(infile, outfile, flip_x, flip_y)
            sys.stderr.write(f"\rProcessed {i} / {total}")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
