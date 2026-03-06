#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import glob
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
        "  tiff2fits.py [--flip] [--header source.fit] input_spec output_spec\n"
        "\n"
        "    input_spec   - single TIFF file, wildcard (*.tif), numbered pattern\n"
        "                   (img0001.tif), or @list.txt\n"
        "    output_spec  - single file (e.g. image.fit) OR numbered pattern\n"
        "                   (e.g. out0001.fit) OR wildcard (e.g. *.fit)\n"
        "                   Wildcard: * is replaced with input filename stem\n"
        "\n"
        "    --header F   - take FITS header from file F (applied to all outputs)\n"
        "                   Without --header: if output .fit already exists,\n"
        "                   its header is preserved before overwriting\n"
        "\n"
        "    --flip       - flip image vertically (reverse fits2tiff --flip)\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]
    flip = False
    header_file = None

    # Extract --flip flag
    if "--flip" in args:
        flip = True
        args.remove("--flip")

    # Extract --header option
    if "--header" in args:
        idx = args.index("--header")
        if idx + 1 >= len(args):
            sys.stderr.write("Error: --header requires a FITS filename.\n")
            sys.exit(1)
        header_file = args[idx + 1]
        if not os.path.isfile(header_file):
            sys.stderr.write(f"Error: header file not found: {header_file}\n")
            sys.exit(1)
        args = args[:idx] + args[idx + 2:]

    if len(args) != 2:
        usage()

    input_spec = args[0]
    output_spec = args[1]

    return input_spec, output_spec, flip, header_file


# ---------------------------------------------------------------
# TIFF input expansion (batch_utils only handles .fit/.fits sequences)
# ---------------------------------------------------------------

_TIFF_SEQ_RE = re.compile(r"^(.*?)(\d+)(\.tiff?)$", re.IGNORECASE)
_FITS_SEQ_RE = re.compile(r"^(.*?)(\d+)(\.(?:fit|fits))$", re.IGNORECASE)


def _find_tiff_sequence(first_path):
    """Discover contiguous numbered TIFF sequence starting from first_path."""
    base = os.path.basename(first_path)
    m = _TIFF_SEQ_RE.match(base)
    if not m:
        return []

    prefix, digits, ext = m.group(1), m.group(2), m.group(3)
    width = len(digits)
    start_index = int(digits)
    dir_ = os.path.dirname(os.path.abspath(first_path)) or "."

    seq = []
    idx = start_index
    while True:
        name = f"{prefix}{str(idx).zfill(width)}{ext}"
        full = os.path.join(dir_, name)
        if not os.path.isfile(full):
            break
        seq.append(full)
        idx += 1

    return seq


def expand_tiff_input_spec(spec):
    """Expand input spec supporting .tif/.tiff numbered sequences."""
    spec = spec.strip()

    # Wildcard
    if "*" in spec or "?" in spec:
        matches = glob.glob(spec)
        files = [os.path.abspath(m) for m in matches if os.path.isfile(m)]
        files.sort()
        if not files:
            raise FileNotFoundError(f"No files match wildcard pattern: {spec}")
        return files

    # @list.txt
    if spec.startswith("@"):
        return batch_utils._expand_list_file(spec[1:])

    # list.txt / list.lst
    if spec.lower().endswith((".txt", ".lst")) and os.path.isfile(spec):
        return batch_utils._expand_list_file(spec)

    # Single file or numbered sequence
    if os.path.isfile(spec):
        seq = _find_tiff_sequence(spec)
        if seq:
            return seq
        return [os.path.abspath(spec)]

    raise FileNotFoundError(f"Input not recognized or file not found: {spec}")


# ---------------------------------------------------------------
# Output pairing
# ---------------------------------------------------------------

def _apply_wildcard_output(inputs, output_spec):
    """Replace * in output_spec with each input file's stem."""
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


def build_io_pairs(input_spec, output_spec):
    """Build (input, output) pairs for TIFF→FITS conversion."""
    inputs = expand_tiff_input_spec(input_spec)
    if not inputs:
        raise ValueError("No input files found.")

    # Wildcard output: *.fit, prefix_*.fits, etc.
    if "*" in output_spec:
        return _apply_wildcard_output(inputs, output_spec)

    # Single input -> single output
    if len(inputs) == 1:
        return [(inputs[0], output_spec)]

    # Multiple inputs => try numbered FITS pattern
    base = os.path.basename(output_spec)
    m = _FITS_SEQ_RE.match(base)
    if m:
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

    raise ValueError(
        "Output pattern must contain a numeric field when multiple "
        "input files are provided (e.g. out0001.fit), or use "
        "wildcard (e.g. *.fit) to preserve input names."
    )


# ---------------------------------------------------------------
# Header handling
# ---------------------------------------------------------------

# Keywords managed by astropy — must not be copied from source header
_STRUCTURAL_KEYS = frozenset([
    "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2",
    "BZERO", "BSCALE", "EXTEND",
])


def read_fits_header(path):
    """Read FITS header, stripping structural keywords."""
    with fits.open(path, memmap=False) as hdul:
        hdr = hdul[0].header.copy()
    for kw in _STRUCTURAL_KEYS:
        hdr.pop(kw, None)
    return hdr


def get_header_for_file(outfile, header_file):
    """
    Get FITS header for output file.

    Priority:
    1. Explicit --header file (if provided)
    2. Existing output FITS file (read before overwriting)
    3. None (minimal header created by astropy)
    """
    if header_file:
        return read_fits_header(header_file)

    if os.path.isfile(outfile):
        try:
            return read_fits_header(outfile)
        except Exception:
            return None

    return None


# ---------------------------------------------------------------
# Conversion
# ---------------------------------------------------------------

def convert_file(infile, outfile, flip, header_file):
    """Convert a single TIFF file to FITS."""
    # Read header BEFORE overwriting the output file
    header = get_header_for_file(outfile, header_file)

    im = Image.open(infile)

    # Convert to numpy based on TIFF mode
    if im.mode == "I;16":
        data = np.array(im, dtype=np.uint16)
    elif im.mode == "I":
        data = np.array(im, dtype=np.int32)
    elif im.mode == "F":
        data = np.array(im, dtype=np.float32)
    elif im.mode == "L":
        data = np.array(im, dtype=np.uint8)
    else:
        # RGB or other — convert to grayscale
        im = im.convert("L")
        data = np.array(im, dtype=np.uint8)

    im.close()

    if data.ndim != 2:
        raise ValueError(f"File '{infile}' is not a 2D grayscale image (mode: {im.mode}).")

    if flip:
        data = np.flipud(data)

    # Sanitize NaN/Inf for float data
    if np.issubdtype(data.dtype, np.floating):
        data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)

    hdu = fits.PrimaryHDU(data, header=header)
    hdu.writeto(outfile, overwrite=True)


def main():
    input_spec, output_spec, flip, header_file = parse_args(sys.argv)

    try:
        io_pairs = build_io_pairs(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            convert_file(infile, outfile, flip, header_file)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
