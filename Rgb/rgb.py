#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rgb - Merge/split RGB channels in FITS files.

Merge: combine 3 monochrome FITS (R, G, B) into one RGB FITS (3, H, W).
Split: extract 3 monochrome FITS from one RGB FITS (3, H, W).

Header is taken from the first input (R channel for merge, RGB for split).
"""

import sys
import os
import numpy as np
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


def usage():
    prog = os.path.basename(sys.argv[0])
    sys.stderr.write(
        f"rgb - Merge/split RGB channels in FITS files.\n"
        f"\n"
        f"Usage:\n"
        f"  {prog} --merge inR inG inB --out outRGB\n"
        f"  {prog} --split inRGB --out outR outG outB\n"
        f"\n"
        f"Merge mode:\n"
        f"  Combines 3 monochrome 2D FITS into one 3-channel RGB FITS (3, H, W).\n"
        f"  Each input can be: single file, wildcard (*.fit), numbered (img0001.fit),\n"
        f"  or @list.txt. All three inputs must expand to the same number of files.\n"
        f"\n"
        f"Split mode:\n"
        f"  Extracts R, G, B channels from RGB FITS (3, H, W) into 3 mono FITS.\n"
        f"  Input and each output can be file specs as above.\n"
        f"\n"
        f"Examples:\n"
        f"  {prog} --merge r.fit g.fit b.fit --out rgb.fit\n"
        f"  {prog} --merge r0001.fit g0001.fit b0001.fit --out rgb0001.fit\n"
        f"  {prog} --split rgb.fit --out r.fit g.fit b.fit\n"
        f"  {prog} --split rgb0001.fit --out r0001.fit g0001.fit b0001.fit\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    if not args or args[0] in ("-h", "--help"):
        usage()

    mode = None
    merge_inputs = None   # (r_spec, g_spec, b_spec)
    split_input = None
    out_specs = None

    i = 0
    while i < len(args):
        if args[i] in ("-h", "--help"):
            usage()
        elif args[i] == "--merge":
            mode = "merge"
            if i + 3 >= len(args):
                sys.stderr.write("Error: --merge requires 3 arguments (R G B).\n")
                sys.exit(1)
            merge_inputs = (args[i+1], args[i+2], args[i+3])
            i += 4
        elif args[i] == "--split":
            mode = "split"
            if i + 1 >= len(args):
                sys.stderr.write("Error: --split requires 1 argument (input).\n")
                sys.exit(1)
            split_input = args[i+1]
            i += 2
        elif args[i] == "--out":
            if mode == "merge":
                if i + 1 >= len(args):
                    sys.stderr.write("Error: --out requires 1 argument for merge.\n")
                    sys.exit(1)
                out_specs = (args[i+1],)
                i += 2
            elif mode == "split":
                if i + 3 >= len(args):
                    sys.stderr.write("Error: --out requires 3 arguments (R G B) for split.\n")
                    sys.exit(1)
                out_specs = (args[i+1], args[i+2], args[i+3])
                i += 4
            else:
                sys.stderr.write("Error: --out must come after --merge or --split.\n")
                sys.exit(1)
        elif args[i].startswith("--"):
            sys.stderr.write(f"Error: unknown option '{args[i]}'.\n\n")
            usage()
        else:
            sys.stderr.write(f"Error: unexpected argument '{args[i]}'.\n\n")
            usage()
        # end while

    if mode is None:
        sys.stderr.write("Error: specify --merge or --split.\n\n")
        usage()

    if out_specs is None:
        sys.stderr.write("Error: --out is required.\n\n")
        usage()

    return mode, merge_inputs, split_input, out_specs


def do_merge(r_spec, g_spec, b_spec, out_spec):
    """Merge R, G, B mono FITS into RGB FITS."""
    r_files = batch_utils.expand_input_spec(r_spec)
    g_files = batch_utils.expand_input_spec(g_spec)
    b_files = batch_utils.expand_input_spec(b_spec)

    n = len(r_files)
    if len(g_files) != n or len(b_files) != n:
        sys.stderr.write(
            f"Error: channel file counts differ: R={len(r_files)}, "
            f"G={len(g_files)}, B={len(b_files)}.\n")
        sys.exit(1)

    # Build output list from R files (use R as reference for naming)
    out_pairs = batch_utils.build_io_file_lists_from_list(r_files, out_spec)
    out_files = [p[1] for p in out_pairs]

    total = n
    print(f"RGB merge: {total} file(s)")

    for i in range(total):
        try:
            with fits.open(r_files[i], memmap=False) as hdul_r:
                data_r = hdul_r[0].data
                header = hdul_r[0].header.copy()
                if data_r is None or data_r.ndim != 2:
                    raise ValueError(f"R channel not 2D: {r_files[i]}")

            with fits.open(g_files[i], memmap=False) as hdul_g:
                data_g = hdul_g[0].data
                if data_g is None or data_g.ndim != 2:
                    raise ValueError(f"G channel not 2D: {g_files[i]}")

            with fits.open(b_files[i], memmap=False) as hdul_b:
                data_b = hdul_b[0].data
                if data_b is None or data_b.ndim != 2:
                    raise ValueError(f"B channel not 2D: {b_files[i]}")

            if data_r.shape != data_g.shape or data_r.shape != data_b.shape:
                raise ValueError(
                    f"Shape mismatch: R={data_r.shape}, G={data_g.shape}, "
                    f"B={data_b.shape}")

            # Stack into (3, H, W)
            rgb = np.stack([data_r, data_g, data_b], axis=0)

            # Update header
            header["NAXIS"] = 3
            header["NAXIS3"] = 3
            header["HISTORY"] = "Merged RGB by rgb.py"

            out_dir = os.path.dirname(out_files[i])
            if out_dir and not os.path.exists(out_dir):
                os.makedirs(out_dir, exist_ok=True)

            fits.PrimaryHDU(rgb, header=header).writeto(
                out_files[i], overwrite=True)

            sys.stderr.write(f"\r  [{i+1}/{total}] {os.path.basename(out_files[i])}")
            sys.stderr.flush()

        except Exception as e:
            sys.stderr.write(f"\n  Error: {e}\n")

    sys.stderr.write("\n")
    print(f"Done. {total} file(s) merged.")


def do_split(in_spec, out_r_spec, out_g_spec, out_b_spec):
    """Split RGB FITS into R, G, B mono FITS."""
    in_files = batch_utils.expand_input_spec(in_spec)
    n = len(in_files)

    # Build output lists
    r_pairs = batch_utils.build_io_file_lists_from_list(in_files, out_r_spec)
    g_pairs = batch_utils.build_io_file_lists_from_list(in_files, out_g_spec)
    b_pairs = batch_utils.build_io_file_lists_from_list(in_files, out_b_spec)

    r_files = [p[1] for p in r_pairs]
    g_files = [p[1] for p in g_pairs]
    b_files = [p[1] for p in b_pairs]

    total = n
    print(f"RGB split: {total} file(s)")

    for i in range(total):
        try:
            with fits.open(in_files[i], memmap=False) as hdul:
                data = hdul[0].data
                header = hdul[0].header.copy()

                if data is None or data.ndim != 3 or data.shape[0] != 3:
                    raise ValueError(
                        f"Expected RGB (3, H, W), got shape "
                        f"{data.shape if data is not None else 'None'}")

            # Remove 3rd axis info from header for 2D output
            h2d = header.copy()
            if "NAXIS3" in h2d:
                del h2d["NAXIS3"]
            h2d["NAXIS"] = 2

            for ch, outfile, label in [
                (0, r_files[i], "R"),
                (1, g_files[i], "G"),
                (2, b_files[i], "B"),
            ]:
                ch_header = h2d.copy()
                ch_header["FILTER"] = label
                ch_header["HISTORY"] = f"Split {label} channel by rgb.py"

                out_dir = os.path.dirname(outfile)
                if out_dir and not os.path.exists(out_dir):
                    os.makedirs(out_dir, exist_ok=True)

                fits.PrimaryHDU(data[ch], header=ch_header).writeto(
                    outfile, overwrite=True)

            sys.stderr.write(
                f"\r  [{i+1}/{total}] {os.path.basename(in_files[i])}")
            sys.stderr.flush()

        except Exception as e:
            sys.stderr.write(f"\n  Error: {e}\n")

    sys.stderr.write("\n")
    print(f"Done. {total} file(s) split.")


def main():
    mode, merge_inputs, split_input, out_specs = parse_args(sys.argv)

    if mode == "merge":
        do_merge(merge_inputs[0], merge_inputs[1], merge_inputs[2],
                 out_specs[0])
    elif mode == "split":
        do_split(split_input, out_specs[0], out_specs[1], out_specs[2])


if __name__ == "__main__":
    main()
