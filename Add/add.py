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
        "  add.py input_spec output_spec operand --screen\n"
        "\n"
        "    input_spec   - single file OR numbered pattern (e.g. light0001.fit)\n"
        "                   OR wildcard mask (e.g. *.fit, light_*.fit)\n"
        "    output_spec  - single file OR numbered pattern for results;\n"
        "                   when input_spec is a mask, outputs are numbered\n"
        "                   according to batch_utils.build_io_file_lists rules\n"
        "    operand      - number OR FITS file OR numbered FITS pattern\n"
        "                   (value/image to add to, or screen with, input)\n"
        "    offset       - optional numeric value added to the result (add mode only)\n"
        "    --screen     - screen blend instead of add:\n"
        "                     result = 1 - (1 - input) * (1 - operand)\n"
        "                   the counterpart of unscreen, used to recombine a starless\n"
        "                   layer with a stars layer. Works on 2D and 3D (RGB) images.\n"
        "                   Float data is treated as [0,1] (white = 1.0); integer data\n"
        "                   uses its dtype range (white = dtype max). Takes no offset.\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    screen = False
    if "--screen" in args:
        screen = True
        args = [a for a in args if a != "--screen"]

    if len(args) not in (3, 4):
        usage()

    input_pattern = args[0]
    output_pattern = args[1]
    operand_str = args[2]

    offset = 0.0
    if len(args) == 4:
        if screen:
            sys.stderr.write("Error: --screen does not take an offset.\n")
            sys.exit(1)
        try:
            offset = float(args[3])
        except ValueError:
            sys.stderr.write("Error: offset must be a number.\n")
            sys.exit(1)

    return input_pattern, output_pattern, operand_str, offset, screen


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


def apply_screen_operation(base_data, operand):
    """
    Screen blend: result = 1 - (1 - base) * (1 - operand).

    Screen is defined on [0, 1] (the counterpart of unscreen). Elementwise, so it
    handles 2D and 3D (RGB) images alike.
    - Float types: data is treated as already in [0, 1] (white = 1.0).
    - Integer types: normalized by the dtype range (white = dtype max), screened,
      scaled back, clamped, and rounded.
    """
    if base_data is None or base_data.ndim not in (2, 3):
        raise ValueError("Expected a 2D or 3D primary image in input.")

    if isinstance(operand, np.ndarray):
        if operand.shape != base_data.shape:
            raise ValueError(
                f"Operand image shape {operand.shape} does not match input {base_data.shape}."
            )
        op = operand.astype(np.float64)
    else:
        op = float(operand)

    if np.issubdtype(base_data.dtype, np.floating):
        a = base_data.astype(np.float64)
        work = 1.0 - (1.0 - a) * (1.0 - op)
        return work.astype(base_data.dtype)

    info = np.iinfo(base_data.dtype)
    white = float(info.max)
    if white <= 0:
        white = 1.0
    a = base_data.astype(np.float64) / white
    b = op / white
    work = (1.0 - (1.0 - a) * (1.0 - b)) * white
    np.clip(work, info.min, info.max, out=work)
    work = np.rint(work)
    return work.astype(base_data.dtype)


def process_file(infile, outfile, operand_spec, file_index, offset, screen):
    """Load, process, and save single FITS file."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"File '{infile}' has no primary image data.")
        allowed = (2, 3) if screen else (2,)
        if hdul[0].data.ndim not in allowed:
            kind = "2D or 3D" if screen else "2D"
            raise ValueError(f"File '{infile}' is not a {kind} image.")

        data = hdul[0].data
        header = hdul[0].header

        # Get operand (scalar or file path)
        operand_raw = batch_utils.get_operand_for_file(operand_spec, file_index)

        # Resolve to numpy array (handles both scalar constants and FITS files)
        operand = batch_utils.resolve_operand_value(
            operand_raw, data.shape, data.dtype
        )

        if screen:
            new_data = apply_screen_operation(data, operand)
        else:
            new_data = apply_add_operation(data, operand, offset)

        hdul[0].data = new_data
        hdul[0].header = header
        hdul.writeto(outfile, overwrite=True)


def main():
    input_pattern, output_pattern, operand_str, offset, screen = parse_args(sys.argv)

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
            process_file(infile, outfile, operand_spec, i - 1, offset, screen)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
