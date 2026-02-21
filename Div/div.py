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
        "  div.py input.fit output.fit operand [scale]\n"
        "\n"
        "    input.fit   - single file OR numbered pattern (e.g. light0001.fit)\n"
        "    output.fit  - single file OR numbered pattern for results\n"
        "    operand     - number OR FITS file OR numbered FITS pattern\n"
        "                  (value to divide input by)\n"
        "    scale       - optional numeric multiplier applied to the result\n"
    )
    sys.exit(1)


def parse_args(argv):
    if len(argv) not in (4, 5):
        usage()

    input_pattern = argv[1]
    output_pattern = argv[2]
    operand_str = argv[3]

    scale = 1.0
    if len(argv) == 5:
        try:
            scale = float(argv[4])
        except ValueError:
            sys.stderr.write("Error: scale must be a number.\n")
            sys.exit(1)

    return input_pattern, output_pattern, operand_str, scale


def apply_div_operation(base_data, operand, scale):
    """
    Core math: result = (base / operand) * scale

    Division by zero:
      - For any zero operand value, result is set to 0 at that position.
    Type handling:
      - Float types: computed in float64, then cast back.
      - Integer types: computed in float64, clamped to dtype range, rounded, cast back.
    """
    if base_data is None or base_data.ndim != 2:
        raise ValueError("Expected 2D primary image in input.")

    base = base_data.astype(np.float64)

    if isinstance(operand, np.ndarray):
        if operand.shape != base_data.shape:
            raise ValueError(
                f"Operand image shape {operand.shape} does not match input {base_data.shape}."
            )
        op = operand.astype(np.float64)
        work = np.zeros_like(base, dtype=np.float64)
        mask = (op != 0.0)
        work[mask] = (base[mask] / op[mask]) * scale
    else:
        op = float(operand)
        work = np.zeros_like(base, dtype=np.float64)
        if op != 0.0:
            work = (base / op) * scale

    if np.issubdtype(base_data.dtype, np.floating):
        if base_data.dtype == np.float32:
            return work.astype(np.float32)
        if base_data.dtype == np.float64:
            return work.astype(np.float64)
        return work.astype(base_data.dtype)
    else:
        info = np.iinfo(base_data.dtype)
        np.clip(work, info.min, info.max, out=work)
        work = np.rint(work)
        return work.astype(base_data.dtype)


def process_file(infile, outfile, operand_spec, file_index, scale):
    """Load, process, and save single FITS file."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None or hdul[0].data.ndim != 2:
            raise ValueError(f"File '{infile}' has invalid image data.")

        data = hdul[0].data
        header = hdul[0].header

        # Get operand (scalar or file path)
        operand_raw = batch_utils.get_operand_for_file(operand_spec, file_index)

        # Resolve to numpy array (handles both scalar constants and FITS files)
        operand = batch_utils.resolve_operand_value(
            operand_raw, data.shape, data.dtype
        )

        new_data = apply_div_operation(data, operand, scale)

        hdul[0].data = new_data
        hdul[0].header = header
        hdul.writeto(outfile, overwrite=True)


def main():
    input_pattern, output_pattern, operand_str, scale = parse_args(sys.argv)

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
            process_file(infile, outfile, operand_spec, i - 1, scale)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
