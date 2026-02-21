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
        "  arith.py (add|sub|mul|div) input.fit output.fit operand [param]\n"
        "\n"
        "    operation : add  -> result = base + operand + offset\n"
        "                sub  -> result = base - operand + offset\n"
        "                mul  -> result = base * operand * scale\n"
        "                div  -> result = (base / operand) * scale\n"
        "\n"
        "    input.fit : single file OR numbered pattern (e.g. light0001.fit)\n"
        "    output.fit: single file OR numbered pattern for results\n"
        "\n"
        "    operand   : number OR FITS file OR numbered FITS pattern\n"
        "\n"
        "    param     : for add/sub -> offset (default = 0.0)\n"
        "                for mul/div -> scale  (default = 1.0)\n"
    )
    sys.exit(1)


def parse_args(argv):
    if len(argv) < 5 or len(argv) > 6:
        usage()

    op = argv[1].lower()
    if op not in ("add", "sub", "mul", "div"):
        sys.stderr.write("Error: first argument must be one of: add, sub, mul, div.\n")
        usage()

    input_pattern = argv[2]
    output_pattern = argv[3]
    operand_str = argv[4]

    if op in ("add", "sub"):
        default_param = 0.0
    else:
        default_param = 1.0

    param = default_param
    if len(argv) == 6:
        try:
            param = float(argv[5])
        except ValueError:
            sys.stderr.write("Error: param must be a number.\n")
            sys.exit(1)

    return op, input_pattern, output_pattern, operand_str, param


def apply_add(base_data, operand, offset):
    """result = base + operand + offset"""
    if base_data is None or base_data.ndim != 2:
        raise ValueError("Expected 2D primary image in input.")

    if isinstance(operand, np.ndarray):
        if operand.shape != base_data.shape:
            raise ValueError(
                f"Operand shape {operand.shape} does not match input {base_data.shape}."
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


def apply_sub(base_data, operand, offset):
    """result = base - operand + offset"""
    if base_data is None or base_data.ndim != 2:
        raise ValueError("Expected 2D primary image in input.")

    if isinstance(operand, np.ndarray):
        if operand.shape != base_data.shape:
            raise ValueError(
                f"Operand shape {operand.shape} does not match input {base_data.shape}."
            )
        op = operand
    else:
        op = float(operand)

    if np.issubdtype(base_data.dtype, np.floating):
        work = base_data.astype(np.float64)
        work = work - op + offset
        if base_data.dtype == np.float32:
            return work.astype(np.float32)
        if base_data.dtype == np.float64:
            return work.astype(np.float64)
        return work.astype(base_data.dtype)

    info = np.iinfo(base_data.dtype)
    work = base_data.astype(np.float64)
    work = work - op + offset
    np.clip(work, info.min, info.max, out=work)
    work = np.rint(work)
    return work.astype(base_data.dtype)


def apply_mul(base_data, operand, scale):
    """result = base * operand * scale"""
    if base_data is None or base_data.ndim != 2:
        raise ValueError("Expected 2D primary image in input.")

    if isinstance(operand, np.ndarray):
        if operand.shape != base_data.shape:
            raise ValueError(
                f"Operand shape {operand.shape} does not match input {base_data.shape}."
            )
        op = operand
    else:
        op = float(operand)

    if np.issubdtype(base_data.dtype, np.floating):
        work = base_data.astype(np.float64)
        work = work * op * scale
        if base_data.dtype == np.float32:
            return work.astype(np.float32)
        if base_data.dtype == np.float64:
            return work.astype(np.float64)
        return work.astype(base_data.dtype)

    info = np.iinfo(base_data.dtype)
    work = base_data.astype(np.float64)
    work = work * op * scale
    np.clip(work, info.min, info.max, out=work)
    work = np.rint(work)
    return work.astype(base_data.dtype)


def apply_div(base_data, operand, scale):
    """result = (base / operand) * scale, division by zero -> 0"""
    if base_data is None or base_data.ndim != 2:
        raise ValueError("Expected 2D primary image in input.")

    base = base_data.astype(np.float64)

    if isinstance(operand, np.ndarray):
        if operand.shape != base_data.shape:
            raise ValueError(
                f"Operand shape {operand.shape} does not match input {base_data.shape}."
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

    info = np.iinfo(base_data.dtype)
    np.clip(work, info.min, info.max, out=work)
    work = np.rint(work)
    return work.astype(base_data.dtype)


def process_file(op, infile, outfile, operand_spec, file_index, param):
    """Load input, apply selected operation, write output."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None or hdul[0].data.ndim != 2:
            raise ValueError(f"File '{infile}' has invalid or missing 2D image data.")

        data = hdul[0].data
        header = hdul[0].header

        # Get operand (scalar or file path)
        operand_raw = batch_utils.get_operand_for_file(operand_spec, file_index)

        # Resolve to numpy array (handles both scalar constants and FITS files)
        operand = batch_utils.resolve_operand_value(
            operand_raw, data.shape, data.dtype
        )

        if op == "add":
            new_data = apply_add(data, operand, param)
        elif op == "sub":
            new_data = apply_sub(data, operand, param)
        elif op == "mul":
            new_data = apply_mul(data, operand, param)
        elif op == "div":
            new_data = apply_div(data, operand, param)
        else:
            raise RuntimeError("Unsupported operation.")

        hdul[0].data = new_data
        hdul[0].header = header
        hdul.writeto(outfile, overwrite=True)


def main():
    op, input_pattern, output_pattern, operand_str, param = parse_args(sys.argv)

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
            process_file(op, infile, outfile, operand_spec, i - 1, param)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
