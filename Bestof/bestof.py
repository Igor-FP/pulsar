#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
bestof - Select best frames by FWHM (seeing quality).

Analyses star FWHM across a set of FITS images and either:
  - Copies the best N% (lowest FWHM) to output
  - Copies frames below a FWHM threshold to output
  - Writes a CSV report (fwhm, filename) sorted by FWHM

Uses SEP-based star detection from lib/star_utils.py.
"""

import sys
import os
import shutil
import time
import threading
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils

try:
    import sep
except ImportError:
    print("ERROR: 'sep' package is required for bestof.")
    print("Install it with:  pip install sep")
    sys.exit(1)

import star_utils


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    prog = os.path.basename(sys.argv[0])
    sys.stderr.write(
        f"bestof - Select best frames by FWHM\n"
        f"\n"
        f"Usage:\n"
        f"  {prog} input_spec [output_spec] [options]\n"
        f"\n"
        f"Selection modes (at least one required unless --csv only):\n"
        f"  --best P         Copy best P percent of frames (lowest FWHM), 1-99\n"
        f"  --threshold T    Copy frames with FWHM <= T pixels\n"
        f"\n"
        f"Report:\n"
        f"  --csv FILE       Write CSV report: fwhm,filename (sorted by FWHM)\n"
        f"\n"
        f"Options:\n"
        f"  --snr N          SNR threshold for star detection (default: 10)\n"
        f"\n"
        f"Examples:\n"
        f"  {prog} light0001.fit best0001.fit --best 70\n"
        f"  {prog} *.fit selected/ --threshold 3.5\n"
        f"  {prog} *.fit --csv report.csv\n"
        f"  {prog} *.fit best/ --best 80 --csv report.csv\n"
        f"  {prog} @list.txt best0001.fit --best 50 --snr 15\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    best_pct = None
    threshold = None
    csv_file = None
    snr = 10.0

    positional = []
    i = 0
    while i < len(args):
        if args[i] == "--best":
            if i + 1 >= len(args):
                print("ERROR: --best requires a percentage value")
                sys.exit(1)
            try:
                best_pct = float(args[i + 1])
            except ValueError:
                print(f"ERROR: --best must be a number, got: {args[i + 1]}")
                sys.exit(1)
            if best_pct <= 0 or best_pct >= 100:
                print(f"ERROR: --best must be between 1 and 99, got: {best_pct}")
                sys.exit(1)
            i += 2
        elif args[i] == "--threshold":
            if i + 1 >= len(args):
                print("ERROR: --threshold requires a FWHM value")
                sys.exit(1)
            try:
                threshold = float(args[i + 1])
            except ValueError:
                print(f"ERROR: --threshold must be a number, got: {args[i + 1]}")
                sys.exit(1)
            if threshold <= 0:
                print(f"ERROR: --threshold must be positive, got: {threshold}")
                sys.exit(1)
            i += 2
        elif args[i] == "--csv":
            if i + 1 >= len(args):
                print("ERROR: --csv requires a file path")
                sys.exit(1)
            csv_file = args[i + 1]
            i += 2
        elif args[i] == "--snr":
            if i + 1 >= len(args):
                print("ERROR: --snr requires a value")
                sys.exit(1)
            try:
                snr = float(args[i + 1])
            except ValueError:
                print(f"ERROR: --snr must be a number, got: {args[i + 1]}")
                sys.exit(1)
            i += 2
        elif args[i].startswith("--"):
            print(f"ERROR: Unknown option: {args[i]}")
            usage()
        else:
            positional.append(args[i])
            i += 1

    # Validate arguments
    has_copy = best_pct is not None or threshold is not None

    if not has_copy and csv_file is None:
        print("ERROR: Specify at least --best, --threshold, or --csv")
        usage()

    if best_pct is not None and threshold is not None:
        print("ERROR: --best and --threshold are mutually exclusive")
        sys.exit(1)

    if has_copy:
        if len(positional) < 2:
            print("ERROR: Need input_spec and output_spec for copy mode")
            usage()
        input_spec = positional[0]
        output_spec = positional[1]
    else:
        # CSV-only mode
        if len(positional) < 1:
            print("ERROR: Need input_spec")
            usage()
        input_spec = positional[0]
        output_spec = None

    return input_spec, output_spec, best_pct, threshold, csv_file, snr


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    if argv is None:
        argv = sys.argv

    input_spec, output_spec, best_pct, threshold, csv_file, snr = parse_args(argv)

    # Expand input files
    input_files = batch_utils.expand_input_spec(input_spec)
    if not input_files:
        print("ERROR: No input files found.")
        sys.exit(1)

    n_files = len(input_files)
    ncpu = os.cpu_count() or 1
    workers = min(max(1, ncpu - 1), n_files)
    print(f"Analysing {n_files} file(s) using {workers} thread(s)...")

    # Worker function: read FITS and measure FWHM
    def measure_fwhm(fpath):
        with fits.open(fpath, memmap=False) as hdul:
            data = hdul[0].data
            if data is None:
                return None
            if data.ndim == 3:
                data = data[0]
        fwhm = star_utils.estimate_fwhm(data, snr=snr)
        if np.isnan(fwhm):
            return None
        return fwhm

    # Measure FWHM in parallel
    results = []  # [(fwhm, filename), ...]
    t0 = time.time()
    done = 0
    lock = threading.Lock()

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(measure_fwhm, fpath): fpath
            for fpath in input_files
        }

        for future in as_completed(futures):
            fpath = futures[future]
            fname = os.path.basename(fpath)
            with lock:
                done += 1
                counter = done
            try:
                fwhm = future.result()
                if fwhm is None:
                    sys.stderr.write(f"  [{counter}/{n_files}] {fname} — skipped\n")
                else:
                    results.append((fwhm, fpath))
                    sys.stderr.write(f"  [{counter}/{n_files}] {fname}  FWHM = {fwhm:.2f} px\n")
            except Exception as e:
                sys.stderr.write(f"  [{counter}/{n_files}] {fname} — error: {e}\n")

    elapsed = time.time() - t0

    if not results:
        print("ERROR: No files with measurable FWHM.")
        sys.exit(1)

    # Sort by FWHM ascending
    results.sort(key=lambda r: r[0])

    fwhm_values = [r[0] for r in results]
    print(f"\n{len(results)} files analysed in {elapsed:.1f}s")
    print(f"  FWHM range: {min(fwhm_values):.2f} — {max(fwhm_values):.2f} px")
    print(f"  FWHM median: {np.median(fwhm_values):.2f} px")

    # Write CSV report
    if csv_file is not None:
        with open(csv_file, "w") as f:
            f.write("fwhm,filename\n")
            for fwhm, fpath in results:
                f.write(f"{fwhm:.3f},{os.path.basename(fpath)}\n")
        print(f"\nCSV report: {csv_file}")

    # Select files for copying
    if best_pct is not None:
        n_keep = max(1, int(len(results) * best_pct / 100.0 + 0.5))
        selected = results[:n_keep]
        fwhm_cutoff = selected[-1][0]
        print(f"\nSelecting best {best_pct:.0f}%: {len(selected)} of {len(results)} files (FWHM <= {fwhm_cutoff:.2f} px)")

    elif threshold is not None:
        selected = [(f, p) for f, p in results if f <= threshold]
        if not selected:
            print(f"\nNo files with FWHM <= {threshold:.2f} px")
            sys.exit(0)
        print(f"\nSelecting FWHM <= {threshold:.2f}: {len(selected)} of {len(results)} files")

    else:
        # CSV-only mode, no copying
        return

    # Copy selected files to output
    selected_paths = [p for _, p in selected]
    io_pairs = batch_utils.build_io_file_lists_from_list(selected_paths, output_spec)

    for src, dst in io_pairs:
        dst_dir = os.path.dirname(dst)
        if dst_dir and not os.path.exists(dst_dir):
            os.makedirs(dst_dir, exist_ok=True)
        shutil.copy2(src, dst)

    print(f"Copied {len(io_pairs)} file(s).")


if __name__ == "__main__":
    main()
