#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import numpy as np
from astropy.io import fits

try:
    from PIL import Image
except ImportError:
    sys.stderr.write(
        "Error: Pillow is required but not installed.\n"
        "Install it with: pip install Pillow\n"
    )
    sys.exit(1)

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  fits2tiff.py [--bits 8|16|32] [--flip] input_spec output_spec\n"
        "\n"
        "    input_spec   - single file OR numbered pattern (e.g. light0001.fit)\n"
        "                   OR wildcard mask (e.g. *.fit) OR @list.txt\n"
        "    output_spec  - single file (e.g. image.tif) OR numbered pattern\n"
        "                   (e.g. out0001.tif) OR wildcard (e.g. *.tif)\n"
        "                   Wildcard: * is replaced with input filename stem\n"
        "\n"
        "    --bits 8     - linear stretch [min,max] -> [0,255], uint8 output\n"
        "    --bits 16    - clamp to [0,65535], uint16 output\n"
        "    --bits 32    - float32 output, no scaling\n"
        "    (default)    - auto: uint8/int8->8, uint16/int16->16, else->32\n"
        "\n"
        "    --flip       - flip image vertically (FITS bottom-left to TIFF top-left)\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]
    bits = None
    flip = False

    # Extract --flip flag
    if "--flip" in args:
        flip = True
        args.remove("--flip")

    # Extract --bits option
    if "--bits" in args:
        idx = args.index("--bits")
        if idx + 1 >= len(args):
            sys.stderr.write("Error: --bits requires a value (8, 16, or 32).\n")
            sys.exit(1)
        bits_str = args[idx + 1]
        if bits_str not in ("8", "16", "32"):
            sys.stderr.write("Error: --bits must be 8, 16, or 32.\n")
            sys.exit(1)
        bits = int(bits_str)
        args = args[:idx] + args[idx + 2:]

    if len(args) != 2:
        usage()

    input_spec = args[0]
    output_spec = args[1]

    return input_spec, output_spec, bits, flip


# Numbered pattern for TIFF output: prefix + digits + .tif(f)
_TIFF_SEQ_RE = re.compile(r"^(.*?)(\d+)(\.tiff?)$", re.IGNORECASE)


def _apply_wildcard_output(inputs, output_spec):
    """Replace * in output_spec with each input file's stem (name without ext)."""
    out_dir = os.path.dirname(output_spec)
    out_pattern = os.path.basename(output_spec)

    pairs = []
    for inp in inputs:
        stem = os.path.splitext(os.path.basename(inp))[0]
        outname = out_pattern.replace("*", stem)
        if out_dir:
            outpath = os.path.join(out_dir, outname)
        else:
            outpath = outname
        pairs.append((inp, outpath))
    return pairs


def build_tiff_io_pairs(input_spec, output_spec):
    """Build (input, output) pairs supporting .tif/.tiff output patterns."""
    inputs = batch_utils.expand_input_spec(input_spec)
    if not inputs:
        raise ValueError("No input files found.")

    # Wildcard output: *.tif, prefix_*.tiff, etc.
    if "*" in output_spec:
        return _apply_wildcard_output(inputs, output_spec)

    # Single input -> single output
    if len(inputs) == 1:
        return [(inputs[0], output_spec)]

    # Multiple inputs => output must be numbered pattern
    base = os.path.basename(output_spec)
    m = _TIFF_SEQ_RE.match(base)
    if not m:
        raise ValueError(
            "Output pattern must contain a numeric field when multiple "
            "input files are provided (e.g. out0001.tif), or use "
            "wildcard (e.g. *.tif) to preserve input names."
        )

    prefix, digits, ext = m.group(1), m.group(2), m.group(3)
    width = len(digits)
    start_index = int(digits)
    out_dir = os.path.dirname(os.path.abspath(output_spec)) or "."

    pairs = []
    idx = start_index
    for inp in inputs:
        fname = f"{prefix}{str(idx).zfill(width)}{ext}"
        pairs.append((inp, os.path.join(out_dir, fname)))
        idx += 1

    return pairs


def auto_bits(dtype):
    """Choose output bit depth based on FITS data type."""
    if dtype in (np.uint8, np.int8):
        return 8
    if dtype in (np.uint16, np.int16):
        return 16
    return 32


def convert_file(infile, outfile, bits, flip):
    """Convert a single FITS file to TIFF."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"File '{infile}' has no primary image data.")
        if hdul[0].data.ndim != 2:
            raise ValueError(f"File '{infile}' is not a 2D image.")

        data = hdul[0].data
        orig_dtype = data.dtype

    # Determine output bits
    out_bits = bits if bits is not None else auto_bits(orig_dtype)

    # Sanitize NaN/Inf
    work = data.astype(np.float64)
    work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)

    if flip:
        work = np.flipud(work)

    if out_bits == 8:
        vmin = work.min()
        vmax = work.max()
        if vmax > vmin:
            scaled = (work - vmin) / (vmax - vmin) * 255.0
        else:
            scaled = np.zeros_like(work)
        scaled = np.clip(scaled, 0, 255)
        result = np.rint(scaled).astype(np.uint8)
        img = Image.fromarray(result)

    elif out_bits == 16:
        clipped = np.clip(work, 0, 65535)
        result = np.rint(clipped).astype(np.uint16)
        img = Image.fromarray(result)

    else:  # 32-bit float, normalize to [0, 1] for Photoshop compatibility
        vmin = work.min()
        vmax = work.max()
        if vmax > vmin:
            work = (work - vmin) / (vmax - vmin)
        else:
            work = np.zeros_like(work)
        result = work.astype(np.float32)
        img = Image.fromarray(result)

    img.save(outfile)


def main():
    input_spec, output_spec, bits, flip = parse_args(sys.argv)

    try:
        io_pairs = build_tiff_io_pairs(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            convert_file(infile, outfile, bits, flip)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
