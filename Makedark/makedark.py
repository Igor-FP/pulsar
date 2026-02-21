#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
makedark.py - Create master dark frames from raw dark files.

Scans input for files with IMAGETYP='Dark Frame', groups by exposure time,
subtracts bias, median-combines, and generates cosmetic correction lists.

Output files are written to the current working directory.

Uses direct function calls and threading for speed (no subprocess overhead).
"""

import sys
import os
import tempfile
import shutil
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from astropy.io import fits
import numpy as np

# Add path to shared utilities
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils

# Thread-safe progress output
_progress_lock = threading.Lock()


def print_progress(current, total, prefix="  Progress"):
    """Thread-safe progress output with line clearing."""
    with _progress_lock:
        msg = f"\r{prefix}: {current}/{total}".ljust(60)
        sys.stderr.write(msg)
        sys.stderr.flush()
        if current >= total:
            sys.stderr.write("\n")


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  makedark.py input_spec [bias_spec]\n"
        "\n"
        "    input_spec  - directory OR file mask OR numbered sequence OR @list.txt\n"
        "                  containing dark frames (IMAGETYP='Dark Frame')\n"
        "    bias_spec   - optional: numeric constant OR FITS file OR file list/mask\n"
        "                  (default: 0)\n"
        "\n"
        "Output (to current directory):\n"
        "    dark<exp>.fit   - master dark for each exposure time\n"
        "    cosme<exp>.lst  - hot pixel list for each master dark\n"
        "    bias.fit        - master bias (if bias_spec was a file list)\n"
        "\n"
        "Exposure format: 300s, 500ms (ms if < 1 second)\n"
    )
    sys.exit(1)


def parse_args(argv):
    if len(argv) < 2 or len(argv) > 3:
        usage()

    input_spec = argv[1]
    bias_spec = "0" if len(argv) < 3 else argv[2]

    return input_spec, bias_spec


def expand_input_to_files(spec):
    """Expand input specification to list of files."""
    spec = spec.strip()

    if os.path.isdir(spec):
        import glob
        pattern = os.path.join(spec, "*.fit")
        files = glob.glob(pattern)
        files.extend(glob.glob(os.path.join(spec, "*.fits")))
        files = sorted(set(os.path.abspath(f) for f in files))
        if not files:
            raise FileNotFoundError(f"No FITS files found in directory: {spec}")
        return files

    return batch_utils.expand_input_spec(spec)


def get_header_value(filepath, *keys):
    """Get first available header value from list of keys."""
    with fits.open(filepath, memmap=False) as hdul:
        header = hdul[0].header
        for key in keys:
            if key in header:
                return header[key]
    return None


def filter_dark_frames(files):
    """Filter files to only those with IMAGETYP='Dark Frame'."""
    darks = []
    for f in files:
        try:
            imgtype = get_header_value(f, "IMAGETYP", "IMAGETYPE")
            if imgtype and imgtype.strip().lower() == "dark frame":
                darks.append(f)
        except Exception as e:
            sys.stderr.write(f"Warning: cannot read '{f}': {e}\n")
    return darks


def group_by_exposure(files):
    """Group files by exposure time."""
    groups = {}
    for f in files:
        exp = get_header_value(f, "EXPTIME", "EXPOSURE")
        if exp is None:
            sys.stderr.write(f"Warning: no EXPTIME/EXPOSURE in '{f}', skipping.\n")
            continue
        exp = float(exp)
        if exp not in groups:
            groups[exp] = []
        groups[exp].append(f)
    return groups


def format_exposure_suffix(exp_seconds):
    """Format exposure time for filename."""
    if exp_seconds >= 1.0:
        if exp_seconds == int(exp_seconds):
            return f"{int(exp_seconds)}s"
        else:
            ms = int(round(exp_seconds * 1000))
            return f"{ms}ms"
    else:
        ms = int(round(exp_seconds * 1000))
        return f"{ms}ms"


# ---------------------------------------------------------
# Core processing functions (no subprocess)
# ---------------------------------------------------------

def subtract_bias(infile, outfile, bias_value):
    """Subtract bias from a single file. bias_value can be scalar or ndarray."""
    with fits.open(infile, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        orig_dtype = data.dtype

        work = data.astype(np.float64) - bias_value

        if np.issubdtype(orig_dtype, np.integer):
            info = np.iinfo(orig_dtype)
            work = np.clip(work, info.min, info.max)
            work = np.rint(work)

        out_data = work.astype(orig_dtype)
        hdul[0].data = out_data
        hdul.writeto(outfile, overwrite=True)


def median_combine_files(file_list, output_path):
    """Median combine multiple FITS files using optimized parallel I/O."""
    if not file_list:
        raise ValueError("No files to combine")

    with fits.open(file_list[0], memmap=False) as hdul:
        ref_header = hdul[0].header.copy()
        ref_shape = hdul[0].data.shape
        ref_dtype = hdul[0].data.dtype

    # Use fast median combine with parallel I/O
    median = batch_utils.fast_median_combine(
        file_list, ref_shape, ref_dtype,
        max_memory_gb=16.0,
        progress_callback=None
    )

    if np.issubdtype(ref_dtype, np.integer):
        info = np.iinfo(ref_dtype)
        median = np.clip(median, info.min, info.max)
        median = np.floor(median)
        out_data = median.astype(ref_dtype)
    else:
        out_data = median.astype(ref_dtype)

    hdu = fits.PrimaryHDU(data=out_data, header=ref_header)
    hdu.writeto(output_path, overwrite=True)


def find_hottest_pixels(data, n=10000):
    """Find N hottest pixels in image."""
    if data is None or data.ndim != 2:
        raise ValueError("Expected 2D image")

    flat = data.ravel()

    if np.issubdtype(flat.dtype, np.floating):
        flat = flat.copy()
        bad = ~np.isfinite(flat)
        if np.any(bad):
            flat[bad] = -np.inf

    total = flat.size
    if total == 0:
        return []

    n = min(n, total)
    if n <= 0:
        return []

    idx_part = np.argpartition(flat, -n)[-n:]
    idx_sorted = idx_part[np.argsort(flat[idx_part])[::-1]]

    width = data.shape[1]
    coords = []
    for idx in idx_sorted:
        y, x = divmod(int(idx), width)
        coords.append((x, y))

    return coords


def write_cosme_list(coords, output_path):
    """Write cosmetic pixel list."""
    with open(output_path, "w", newline="") as f:
        for x, y in coords:
            f.write(f"P {x} {y}\r\n")


def prepare_bias(bias_spec):
    """
    Prepare bias value/file.
    Returns: (bias_value, master_bias_path)
    bias_value: float scalar or ndarray
    master_bias_path: path if created, None otherwise
    """
    # Try parsing as numeric constant
    const = batch_utils.parse_numeric_constant(bias_spec)
    if const is not None:
        return float(const), None

    # Try as file specification
    try:
        files = batch_utils.expand_input_spec(bias_spec)
    except FileNotFoundError:
        sys.stderr.write(f"Error: bias_spec '{bias_spec}' is not a valid constant or file.\n")
        sys.exit(1)

    if len(files) == 1:
        # Single file - load as array
        with fits.open(files[0], memmap=False) as hdul:
            bias_data = hdul[0].data.astype(np.float64)
        return bias_data, None

    # Multiple files - create master bias via median
    sys.stderr.write(f"Creating master bias from {len(files)} files...\n")
    master_bias_path = os.path.join(os.getcwd(), "bias.fit")
    median_combine_files(files, master_bias_path)

    with fits.open(master_bias_path, memmap=False) as hdul:
        bias_data = hdul[0].data.astype(np.float64)

    return bias_data, master_bias_path


def process_exposure_group(dark_files, exp_seconds, bias_value, temp_dir, max_workers=None):
    """Process a group of dark files with the same exposure time using threading."""
    suffix = format_exposure_suffix(exp_seconds)
    output_dark = os.path.join(os.getcwd(), f"dark{suffix}.fit")
    output_cosme = os.path.join(os.getcwd(), f"cosme{suffix}.lst")
    n_files = len(dark_files)

    sys.stderr.write(f"\n=== Processing {n_files} darks with exposure {exp_seconds}s ({suffix}) ===\n")

    # Step 1: Subtract bias (parallel)
    sys.stderr.write(f"Subtracting bias from {n_files} files...\n")
    temp_files = []

    def do_subtract(args):
        i, dark_file = args
        temp_out = os.path.join(temp_dir, f"temp_dark_{suffix}_{i:04d}.fit")
        subtract_bias(dark_file, temp_out, bias_value)
        return temp_out

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(do_subtract, (i, f)) for i, f in enumerate(dark_files)]
        for i, future in enumerate(as_completed(futures)):
            temp_files.append(future.result())
            print_progress(i + 1, n_files, "  Bias subtracted")

    temp_files.sort()

    # Step 2: Median combine
    sys.stderr.write(f"Median combining {n_files} bias-subtracted darks...\n")
    median_combine_files(temp_files, output_dark)

    # Step 3: Generate cosmetic correction list
    sys.stderr.write(f"Generating cosmetic list...\n")
    with fits.open(output_dark, memmap=False) as hdul:
        coords = find_hottest_pixels(hdul[0].data, n=10000)
    write_cosme_list(coords, output_cosme)

    # Step 4: Clean up temp files
    for tf in temp_files:
        if os.path.exists(tf):
            os.unlink(tf)

    sys.stderr.write(f"Created: {output_dark}\n")
    sys.stderr.write(f"Created: {output_cosme}\n")


def main():
    input_spec, bias_spec = parse_args(sys.argv)

    try:
        all_files = expand_input_to_files(input_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    sys.stderr.write(f"Found {len(all_files)} FITS files in input.\n")

    dark_files = filter_dark_frames(all_files)
    if not dark_files:
        sys.stderr.write("Error: no files with IMAGETYP='Dark Frame' found.\n")
        sys.exit(1)

    sys.stderr.write(f"Found {len(dark_files)} dark frames.\n")

    groups = group_by_exposure(dark_files)
    if not groups:
        sys.stderr.write("Error: no dark frames with valid EXPTIME/EXPOSURE found.\n")
        sys.exit(1)

    sys.stderr.write(f"Found {len(groups)} exposure group(s): {sorted(groups.keys())}\n")

    # Prepare bias
    bias_value, master_bias_path = prepare_bias(bias_spec)
    if isinstance(bias_value, np.ndarray):
        sys.stderr.write(f"Bias: using image data\n")
    else:
        sys.stderr.write(f"Bias: {bias_value}\n")

    temp_dir = tempfile.mkdtemp(prefix="makedark_")

    try:
        for exp_seconds in sorted(groups.keys()):
            process_exposure_group(groups[exp_seconds], exp_seconds, bias_value, temp_dir)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    sys.stderr.write("\nDone.\n")


if __name__ == "__main__":
    main()
