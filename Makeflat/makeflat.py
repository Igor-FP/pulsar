#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
makeflat.py - Create master flat frames from raw flat files.

Scans input for files with IMAGETYP='Flat Frame', groups by FILTER,
validates exposure times, subtracts darks, normalizes, median-combines,
and applies cosmetic correction.

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
        # Clear line and print progress
        msg = f"\r{prefix}: {current}/{total}".ljust(60)
        sys.stderr.write(msg)
        sys.stderr.flush()
        if current >= total:
            sys.stderr.write("\n")


# Standard filter code mapping (case-insensitive)
FILTER_CODES = {
    "oiii": "o",
    "sii": "s",
    "ha": "h",
    "l": "l",
    "r": "r",
    "g": "g",
    "b": "b",
}


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  makeflat.py input_spec [target_median]\n"
        "\n"
        "    input_spec     - directory OR file mask OR numbered sequence OR @list.txt\n"
        "                     containing flat frames (IMAGETYP='Flat Frame')\n"
        "    target_median  - optional: normalization target (default: 5000)\n"
        "\n"
        "Output (to current directory):\n"
        "    flat_<filter>.fit - master flat for each filter\n"
        "                        (e.g., flat_r.fit, flat_g.fit, flat_Clear.fit)\n"
        "\n"
        "Requires dark<exp>.fit and cosme<exp>.lst for each exposure time.\n"
        "If not found, runs makedark.py automatically.\n"
    )
    sys.exit(1)


def parse_args(argv):
    if len(argv) < 2 or len(argv) > 3:
        usage()

    input_spec = argv[1]
    target_median = 5000.0

    if len(argv) == 3:
        try:
            target_median = float(argv[2])
        except ValueError:
            sys.stderr.write("Error: target_median must be a number.\n")
            sys.exit(1)

    return input_spec, target_median


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
        return files, os.path.abspath(spec)

    files = batch_utils.expand_input_spec(spec)
    base_dir = os.path.dirname(os.path.abspath(files[0])) if files else os.getcwd()
    return files, base_dir


def get_header_value(filepath, *keys):
    """Get first available header value from list of keys."""
    with fits.open(filepath, memmap=False) as hdul:
        header = hdul[0].header
        for key in keys:
            if key in header:
                return header[key]
    return None


def filter_flat_frames(files):
    """Filter files to only those with IMAGETYP='Flat Frame'."""
    flats = []
    for f in files:
        try:
            imgtype = get_header_value(f, "IMAGETYP", "IMAGETYPE")
            if imgtype and imgtype.strip().lower() == "flat frame":
                flats.append(f)
        except Exception as e:
            sys.stderr.write(f"Warning: cannot read '{f}': {e}\n")
    return flats


def get_filter_code(filter_name):
    """Get short filter code from filter name."""
    if not filter_name:
        return "unknown"
    trimmed = filter_name.strip()
    lower = trimmed.lower()
    if lower in FILTER_CODES:
        return FILTER_CODES[lower]
    return trimmed


def group_by_filter(files):
    """Group files by FILTER header (trimmed)."""
    groups = {}
    for f in files:
        try:
            filter_val = get_header_value(f, "FILTER")
            if filter_val is None:
                filter_val = "NoFilter"
            else:
                filter_val = filter_val.strip()

            exp = get_header_value(f, "EXPTIME", "EXPOSURE")
            if exp is None:
                sys.stderr.write(f"Warning: no EXPTIME/EXPOSURE in '{f}', skipping.\n")
                continue
            exp = float(exp)

            if filter_val not in groups:
                groups[filter_val] = []
            groups[filter_val].append((f, exp))
        except Exception as e:
            sys.stderr.write(f"Warning: error reading '{f}': {e}\n")
    return groups


def validate_exposure_consistency(groups):
    """Validate that all files in each filter group have the same exposure time."""
    valid_groups = {}
    for filter_name, file_exp_list in groups.items():
        exposures = set(exp for _, exp in file_exp_list)
        if len(exposures) > 1:
            sys.stderr.write(
                f"Warning: filter '{filter_name}' has inconsistent exposures: {sorted(exposures)}\n"
            )
            valid_groups[filter_name] = None
        else:
            valid_groups[filter_name] = list(exposures)[0]
    return valid_groups


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


def find_dark_and_cosme(exp_seconds, search_dirs):
    """Find dark and cosme files for given exposure time."""
    suffix = format_exposure_suffix(exp_seconds)
    dark_name = f"dark{suffix}.fit"
    cosme_name = f"cosme{suffix}.lst"

    dark_path = None
    cosme_path = None

    for dir_path in search_dirs:
        if not os.path.isdir(dir_path):
            continue
        candidate_dark = os.path.join(dir_path, dark_name)
        candidate_cosme = os.path.join(dir_path, cosme_name)

        if dark_path is None and os.path.isfile(candidate_dark):
            dark_path = candidate_dark
        if cosme_path is None and os.path.isfile(candidate_cosme):
            cosme_path = candidate_cosme
        if dark_path and cosme_path:
            break

    return dark_path, cosme_path


# ---------------------------------------------------------
# Core processing functions (no subprocess)
# ---------------------------------------------------------

def subtract_dark(infile, outfile, dark_data):
    """Subtract dark from a single file."""
    with fits.open(infile, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        orig_dtype = data.dtype

        # Subtraction in float64
        work = data.astype(np.float64) - dark_data

        # Convert back to original dtype
        if np.issubdtype(orig_dtype, np.integer):
            info = np.iinfo(orig_dtype)
            work = np.clip(work, info.min, info.max)
            work = np.rint(work)
        out_data = work.astype(orig_dtype)

        hdul[0].data = out_data
        hdul.writeto(outfile, overwrite=True)


def normalize_ngain(infile, outfile, target_median):
    """Normalize a single file by multiplication (ngain)."""
    with fits.open(infile, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        orig_dtype = data.dtype

        work = data.astype(np.float64)
        current_median = batch_utils.fast_median(work)

        if current_median == 0.0:
            gain = 1.0
        else:
            gain = target_median / current_median
            work = work * gain

        # Convert back
        if np.issubdtype(orig_dtype, np.integer):
            info = np.iinfo(orig_dtype)
            work = np.clip(work, info.min, info.max)
            work = np.rint(work)
        elif np.issubdtype(orig_dtype, np.floating):
            bad = ~np.isfinite(work)
            if np.any(bad):
                work[bad] = 0

        out_data = work.astype(orig_dtype)
        hdul[0].data = out_data
        hdul.writeto(outfile, overwrite=True)

    return current_median, gain


def apply_cosme(infile, outfile, cosme_coords):
    """Apply cosmetic correction to a file."""
    with fits.open(infile, memmap=False) as hdul:
        data = hdul[0].data.copy()
        header = hdul[0].header
        height, width = data.shape

        for x, y in cosme_coords:
            if 0 <= x < width and 0 <= y < height:
                # Get 3x3 neighborhood, excluding center
                y0, y1 = max(0, y - 1), min(height, y + 2)
                x0, x1 = max(0, x - 1), min(width, x + 2)
                neighborhood = data[y0:y1, x0:x1].flatten()
                # Exclude center pixel
                center_idx = (y - y0) * (x1 - x0) + (x - x0)
                if 0 <= center_idx < len(neighborhood):
                    neighborhood = np.delete(neighborhood, center_idx)
                if len(neighborhood) > 0:
                    data[y, x] = np.mean(neighborhood)

        hdul[0].data = data
        hdul.writeto(outfile, overwrite=True)


def read_cosme_list(path):
    """Read cosmetic pixel list."""
    coords = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[0].upper() == "P":
                try:
                    x, y = int(parts[1]), int(parts[2])
                    coords.append((x, y))
                except ValueError:
                    pass
    return coords


def median_combine_files(file_list, output_path):
    """Median combine multiple FITS files using optimized parallel I/O."""
    if not file_list:
        raise ValueError("No files to combine")

    # Load reference header and shape
    with fits.open(file_list[0], memmap=False) as hdul:
        ref_header = hdul[0].header.copy()
        ref_shape = hdul[0].data.shape
        ref_dtype = hdul[0].data.dtype

    # Use fast median combine with parallel I/O
    median = batch_utils.fast_median_combine(
        file_list, ref_shape, ref_dtype,
        max_memory_gb=16.0,
        progress_callback=None  # Silent, we report at higher level
    )

    # Convert to original dtype
    if np.issubdtype(ref_dtype, np.integer):
        info = np.iinfo(ref_dtype)
        median = np.clip(median, info.min, info.max)
        median = np.floor(median)
        out_data = median.astype(ref_dtype)
    else:
        out_data = median.astype(ref_dtype)

    # Write output
    hdu = fits.PrimaryHDU(data=out_data, header=ref_header)
    hdu.writeto(output_path, overwrite=True)


def run_makedark(input_dir):
    """Run makedark.py on input directory."""
    import subprocess
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    MAKEDARK_PY = os.path.join(SCRIPT_DIR, "..", "Makedark", "makedark.py")

    sys.stderr.write(f"\n=== Running makedark on {input_dir} ===\n")
    result = subprocess.run([sys.executable, MAKEDARK_PY, input_dir, "0"])
    if result.returncode != 0:
        raise RuntimeError("makedark failed")


def process_filter_group(filter_name, files, exp_seconds, dark_path, cosme_path,
                         target_median, temp_dir, max_workers=None):
    """Process a group of flat files for one filter using threading."""
    filter_code = get_filter_code(filter_name)
    output_flat = os.path.join(os.getcwd(), f"flat_{filter_code}.fit")
    n_files = len(files)

    sys.stderr.write(f"\n=== Processing filter '{filter_name}' ({n_files} files, exp={exp_seconds}s) -> flat_{filter_code}.fit ===\n")

    # Load dark data once
    with fits.open(dark_path, memmap=False) as hdul:
        dark_data = hdul[0].data.astype(np.float64)

    # Step 1: Subtract dark (parallel)
    sys.stderr.write(f"Subtracting dark from {n_files} files...\n")
    dark_sub_files = []

    def do_subtract(args):
        i, flat_file = args
        temp_out = os.path.join(temp_dir, f"darksub_{filter_code}_{i:04d}.fit")
        subtract_dark(flat_file, temp_out, dark_data)
        return temp_out

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(do_subtract, (i, f)) for i, f in enumerate(files)]
        for i, future in enumerate(as_completed(futures)):
            dark_sub_files.append(future.result())
            print_progress(i + 1, n_files, "  Dark subtracted")

    # Sort by index to maintain order
    dark_sub_files.sort()

    # Step 2: Normalize with ngain (parallel)
    sys.stderr.write(f"Normalizing {n_files} files to median={target_median}...\n")
    norm_files = []

    def do_normalize(args):
        i, ds_file = args
        temp_out = os.path.join(temp_dir, f"norm_{filter_code}_{i:04d}.fit")
        normalize_ngain(ds_file, temp_out, target_median)
        return temp_out

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(do_normalize, (i, f)) for i, f in enumerate(dark_sub_files)]
        for i, future in enumerate(as_completed(futures)):
            norm_files.append(future.result())
            print_progress(i + 1, n_files, "  Normalized")

    norm_files.sort()

    # Step 3: Median combine
    sys.stderr.write(f"Median combining {n_files} normalized flats...\n")
    temp_median = os.path.join(temp_dir, f"median_{filter_code}.fit")
    median_combine_files(norm_files, temp_median)

    # Step 4: Apply cosmetic correction
    sys.stderr.write(f"Applying cosmetic correction...\n")
    cosme_coords = read_cosme_list(cosme_path)
    apply_cosme(temp_median, output_flat, cosme_coords)

    # Clean up temp files
    for f in dark_sub_files + norm_files + [temp_median]:
        if os.path.exists(f):
            os.unlink(f)

    sys.stderr.write(f"Created: {output_flat}\n")


def main():
    input_spec, target_median = parse_args(sys.argv)

    try:
        all_files, input_dir = expand_input_to_files(input_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    sys.stderr.write(f"Found {len(all_files)} FITS files in input.\n")
    sys.stderr.write(f"Input directory: {input_dir}\n")

    flat_files = filter_flat_frames(all_files)
    if not flat_files:
        sys.stderr.write("Error: no files with IMAGETYP='Flat Frame' found.\n")
        sys.exit(1)

    sys.stderr.write(f"Found {len(flat_files)} flat frames.\n")

    groups = group_by_filter(flat_files)
    if not groups:
        sys.stderr.write("Error: no valid flat frames found.\n")
        sys.exit(1)

    sys.stderr.write(f"Found {len(groups)} filter group(s): {list(groups.keys())}\n")

    filter_exposures = validate_exposure_consistency(groups)

    unique_exposures = set()
    for filter_name, exp in filter_exposures.items():
        if exp is not None:
            unique_exposures.add(exp)

    if not unique_exposures:
        sys.stderr.write("Error: no valid filter groups.\n")
        sys.exit(1)

    sys.stderr.write(f"Unique exposure times: {sorted(unique_exposures)}\n")

    search_dirs = [os.getcwd()]
    if input_dir != os.getcwd():
        search_dirs.append(input_dir)

    dark_cosme_map = {}
    missing_exposures = []

    for exp in unique_exposures:
        dark_path, cosme_path = find_dark_and_cosme(exp, search_dirs)
        if dark_path and cosme_path:
            dark_cosme_map[exp] = (dark_path, cosme_path)
            sys.stderr.write(f"Found dark/cosme for {format_exposure_suffix(exp)}\n")
        else:
            missing_exposures.append(exp)

    if missing_exposures:
        sys.stderr.write(f"\nMissing darks for: {[format_exposure_suffix(e) for e in missing_exposures]}\n")
        try:
            run_makedark(input_dir)
        except RuntimeError as e:
            sys.stderr.write(f"Warning: makedark failed: {e}\n")

        for exp in missing_exposures[:]:
            dark_path, cosme_path = find_dark_and_cosme(exp, search_dirs)
            if dark_path and cosme_path:
                dark_cosme_map[exp] = (dark_path, cosme_path)
                missing_exposures.remove(exp)

    if missing_exposures:
        sys.stderr.write(f"\n*** WARNING: Still missing darks for: ***\n")
        for exp in missing_exposures:
            sys.stderr.write(f"  - {format_exposure_suffix(exp)}\n")

    temp_dir = tempfile.mkdtemp(prefix="makeflat_")

    try:
        processed = 0
        skipped = []

        for filter_name, file_exp_list in groups.items():
            exp = filter_exposures.get(filter_name)
            if exp is None:
                skipped.append((filter_name, "inconsistent exposures"))
                continue

            if exp not in dark_cosme_map:
                skipped.append((filter_name, f"no dark for {format_exposure_suffix(exp)}"))
                continue

            dark_path, cosme_path = dark_cosme_map[exp]
            files = [f for f, _ in file_exp_list]

            process_filter_group(filter_name, files, exp, dark_path, cosme_path,
                                 target_median, temp_dir)
            processed += 1

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    sys.stderr.write(f"\n=== Summary ===\n")
    sys.stderr.write(f"Processed {processed} filter group(s).\n")
    if skipped:
        sys.stderr.write(f"Skipped {len(skipped)} group(s):\n")
        for name, reason in skipped:
            sys.stderr.write(f"  - {name}: {reason}\n")

    sys.stderr.write("\nDone.\n")


if __name__ == "__main__":
    main()
