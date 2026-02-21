#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from astropy.io import fits
import numpy as np

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  add.py input_spec output_spec operand [offset]\n"
        "\n"
        "    input_spec   - single file OR numbered pattern (e.g. light0001.fit)\n"
        "                   OR wildcard mask (e.g. *.fit, light_*.fit)\n"
        "    output_spec  - single file OR numbered pattern for results;\n"
        "                   when input_spec is a mask, outputs are numbered\n"
        "                   according to batch_utils.build_io_file_lists rules\n"
        "    operand      - number OR FITS file OR numbered FITS pattern\n"
        "                   (value/image to add to input)\n"
        "    offset       - optional numeric value added to the result\n"
    )
    sys.exit(1)


def parse_args(argv):
    if len(argv) not in (4, 5):
        usage()

    input_pattern = argv[1]
    output_pattern = argv[2]
    operand_str = argv[3]

    offset = 0.0
    if len(argv) == 5:
        try:
            offset = float(argv[4])
        except ValueError:
            sys.stderr.write("Error: offset must be a number.\n")
            sys.exit(1)

    return input_pattern, output_pattern, operand_str, offset


def apply_add_operation(base_data, operand, offset):
    """
    Core arithmetic: result = base + operand + offset

    - base_data: ndarray from input (2D)
    - operand: scalar or ndarray (2D same shape)
    - offset: scalar
    - For integer types: compute in float64, clamp to dtype range, round, cast back.
    - For floats: compute in float64, cast back to original float dtype.
    """
    if base_data is None or base_data.ndim != 2:
        raise ValueError("Expected 2D primary image in input.")

    if isinstance(operand, np.ndarray):
        if operand.shape != base_data.shape:
            raise ValueError(
                f"Operand image shape {operand.shape} does not match input {base_data.shape}."
            )
        op = operand
    else:
        op = float(operand)

    if np.issubdtype(base_data.dtype, np.floating):
        work = base_data.astype(np.float64)
        work = work + op + offset

        if base_data.dtype == np.float32:
            return work.astype(np.float32)
        if base_data.dtype == np.float64:
            return work.astype(np.float64)
        return work.astype(base_data.dtype)

    info = np.iinfo(base_data.dtype)
    work = base_data.astype(np.float64)
    work = work + op + offset

    np.clip(work, info.min, info.max, out=work)
    work = np.rint(work)
    return work.astype(base_data.dtype)


def process_file(infile, outfile, operand_spec, file_index, offset):
    """Load, process, and save single FITS file."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"File '{infile}' has no primary image data.")
        if hdul[0].data.ndim != 2:
            raise ValueError(f"File '{infile}' is not a 2D image.")

        data = hdul[0].data
        header = hdul[0].header

        # Get operand (scalar or file path)
        operand_raw = batch_utils.get_operand_for_file(operand_spec, file_index)

        # Resolve to numpy array (handles both scalar constants and FITS files)
        operand = batch_utils.resolve_operand_value(
            operand_raw, data.shape, data.dtype
        )

        new_data = apply_add_operation(data, operand, offset)

        hdul[0].data = new_data
        hdul[0].header = header
        hdul.writeto(outfile, overwrite=True)


def main():
    input_pattern, output_pattern, operand_str, offset = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_pattern, output_pattern)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    try:
        operand_spec = batch_utils.build_operand_spec(operand_str, len(io_pairs))
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    total = len(io_pairs)
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            process_file(infile, outfile, operand_spec, i - 1, offset)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
