#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
xisf2fits - Convert XISF (PixInsight) files to FITS format.

Preserves pixel data as-is (no type conversion by default) and restores
FITS header keywords from XISF FITSKeyword metadata.

Requires: xisf (pip install xisf)
"""

import sys
import os
import numpy as np
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils

try:
    from xisf import XISF
except ImportError:
    print("ERROR: 'xisf' package is required for xisf2fits.")
    print("Install it with:  pip install xisf")
    sys.exit(1)


# FITS structural keywords that astropy generates automatically
_SKIP_KEYWORDS = frozenset({
    'SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'NAXIS3',
    'EXTEND', 'BZERO', 'BSCALE', 'END', '',
})


def usage():
    prog = os.path.basename(sys.argv[0])
    sys.stderr.write(
        f"xisf2fits - Convert XISF to FITS\n"
        f"\n"
        f"Usage:\n"
        f"  {prog} input_spec output_spec\n"
        f"\n"
        f"  input_spec  - single .xisf file, wildcard (*.xisf), or @list.txt\n"
        f"  output_spec - single .fit file, numbered pattern (out0001.fit),\n"
        f"                or directory (outdir/)\n"
        f"\n"
        f"Pixel data is preserved as-is (float32/float64/uint16 etc.).\n"
        f"FITS header keywords are restored from XISF metadata.\n"
        f"\n"
        f"Examples:\n"
        f"  {prog} image.xisf image.fit\n"
        f"  {prog} *.xisf converted/\n"
        f"  {prog} @list.txt out0001.fit\n"
    )
    sys.exit(1)


def expand_xisf_input(spec):
    """Expand input spec for .xisf files (not .fit, so custom expansion)."""
    import glob

    if spec.startswith("@"):
        listfile = spec[1:]
        if not os.path.isfile(listfile):
            raise FileNotFoundError(f"List file not found: {listfile}")
        files = []
        with open(listfile, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    if os.path.isfile(line):
                        files.append(os.path.abspath(line))
                    else:
                        raise FileNotFoundError(f"File not found: {line}")
        if not files:
            raise FileNotFoundError(f"No files in list: {listfile}")
        return files

    if "*" in spec or "?" in spec:
        files = sorted(f for f in glob.glob(spec) if os.path.isfile(f))
        if not files:
            raise FileNotFoundError(f"No files match: {spec}")
        return [os.path.abspath(f) for f in files]

    if os.path.isfile(spec):
        return [os.path.abspath(spec)]

    raise FileNotFoundError(f"Not found: {spec}")


def build_output_list(input_files, output_spec):
    """Build output file paths from input list and output spec."""
    n = len(input_files)

    # Directory: preserve base name, change extension
    if output_spec.endswith(("/", "\\")) or os.path.isdir(output_spec):
        outdir = output_spec.rstrip("/\\")
        result = []
        for f in input_files:
            base = os.path.basename(f)
            root, _ = os.path.splitext(base)
            result.append(os.path.join(outdir, root + ".fit"))
        return list(zip(input_files, result))

    # Single file
    if n == 1:
        return [(input_files[0], output_spec)]

    # Numbered pattern: use batch_utils
    return batch_utils.build_io_file_lists_from_list(input_files, output_spec)


def cast_fits_value(key, value_str):
    """Try to cast a FITSKeyword string value to the appropriate Python type."""
    s = value_str.strip()

    # Boolean
    if s == 'T':
        return True
    if s == 'F':
        return False

    # Remove surrounding quotes for string values
    if len(s) >= 2 and s[0] == "'" and s[-1] == "'":
        return s[1:-1].rstrip()

    # Integer
    try:
        return int(s)
    except ValueError:
        pass

    # Float
    try:
        return float(s)
    except ValueError:
        pass

    # Return as string
    return s


def build_fits_header(xisf_meta):
    """Build astropy FITS header from XISF image metadata."""
    header = fits.Header()

    if not xisf_meta or 'FITSKeywords' not in xisf_meta:
        return header

    fits_kw = xisf_meta['FITSKeywords']

    for key, entries in fits_kw.items():
        if key in _SKIP_KEYWORDS:
            continue

        for entry in entries:
            value_str = entry.get('value', '')
            comment = entry.get('comment', '')

            if key == 'COMMENT':
                header['COMMENT'] = comment or value_str
                continue
            if key == 'HISTORY':
                header['HISTORY'] = comment or value_str
                continue

            value = cast_fits_value(key, value_str)
            try:
                header[key] = (value, comment)
            except (ValueError, KeyError):
                # Skip malformed keywords
                pass

    return header


def convert_file(infile, outfile):
    """Convert a single XISF file to FITS."""
    xisf = XISF(infile)
    im_data = xisf.read_image(0)
    meta_list = xisf.get_images_metadata()
    meta = meta_list[0] if meta_list else {}

    # Handle channel axis: xisf returns HWC, FITS expects CHW or 2D
    if im_data.ndim == 3:
        if im_data.shape[2] == 1:
            # Mono: squeeze to 2D
            im_data = im_data[:, :, 0]
        else:
            # RGB: HWC -> CHW for FITS
            im_data = np.moveaxis(im_data, 2, 0)

    # Build FITS header from XISF metadata
    header = build_fits_header(meta)
    header['HISTORY'] = f'Converted from XISF by xisf2fits'

    # Ensure C-contiguous
    im_data = np.ascontiguousarray(im_data)

    # Write FITS
    outdir = os.path.dirname(outfile)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    hdu = fits.PrimaryHDU(im_data, header=header)
    hdu.writeto(outfile, overwrite=True)


def parse_args(argv):
    args = argv[1:]
    if len(args) < 2:
        usage()
    return args[0], args[1]


def main():
    input_spec, output_spec = parse_args(sys.argv)

    input_files = expand_xisf_input(input_spec)
    if not input_files:
        print("ERROR: No input files found.")
        sys.exit(1)

    io_pairs = build_output_list(input_files, output_spec)

    total = len(io_pairs)
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            convert_file(infile, outfile)
            sys.stderr.write(f"\r  [{i}/{total}] {os.path.basename(infile)} -> {os.path.basename(outfile)}")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\n  ERROR: {os.path.basename(infile)}: {e}\n")

    sys.stderr.write("\n")
    print(f"Done. Converted {total} file(s).")


if __name__ == "__main__":
    main()
