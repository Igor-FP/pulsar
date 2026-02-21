#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils as bu

USAGE = """Usage:
  normalize.py input_spec output_spec [basefile.fit] [method]

Where:
  input_spec   - input FITS: numbered sequence (e.g. img0001.fit),
                 wildcard (e.g. *.fit) or single file.
  output_spec  - output FITS name or pattern (same semantics as in other tools).
  basefile.fit - optional reference frame:
                   - if omitted: use first input file
                     (first in numbered sequence or sorted wildcard list).
  method       - optional method code:
                   1 = linear regression vs reference (per-frame)
                   2 = robust regression vs reference (sigma-clipped)
                   3 = global iterative normalization of all frames
                 default = 1

Notes:
  - Model: I = B * R + C
    Normalized output is: (I - C) / B
  - Method 3 ignores basefile if provided.
"""


# ---------- Progress bar ----------

def print_progress(current, total, prefix=""):
    """Print console progress bar on a single updating line."""
    width = 40
    if total <= 0:
        total = 1
    ratio = float(current) / float(total)
    ratio = max(0.0, min(1.0, ratio))
    filled = int(width * ratio)
    bar = "#" * filled + "-" * (width - filled)
    line = f"{prefix}[{bar}] {current}/{total}"
    if current < total:
        sys.stdout.write("\r" + line)
        sys.stdout.flush()
    else:
        sys.stdout.write("\r" + line + "\n")
        sys.stdout.flush()


# ---------- FITS I/O helpers ----------

def read_fits_data(fname):
    """Read 2D FITS primary HDU as float64 array, with original header and dtype."""
    with fits.open(fname, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    if data is None or data.ndim != 2:
        raise RuntimeError(f"No 2D primary image data in file: {fname}")
    orig_dtype = data.dtype
    return np.array(data, dtype=np.float64), header, orig_dtype


def convert_to_original_dtype(data, orig_dtype):
    """
    Convert float64 normalized data back to original dtype with clamping.

    For integer types: clamp to full dtype range.
    For float types: cast back, replace NaN/Inf with 0.
    """
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        arr = np.clip(data, info.min, info.max)
        arr = np.rint(arr)
        return arr.astype(orig_dtype)

    if np.issubdtype(orig_dtype, np.floating):
        arr = np.array(data, dtype=orig_dtype)
        bad = ~np.isfinite(arr)
        if np.any(bad):
            arr[bad] = 0
        return arr

    # Fallback: use float32
    arr = np.array(data, dtype=np.float32)
    bad = ~np.isfinite(arr)
    if np.any(bad):
        arr[bad] = 0
    return arr


def write_fits_data(fname, data, header_like, orig_dtype):
    """Write FITS file preserving header, using original dtype and safe clamping."""
    out_data = convert_to_original_dtype(data, orig_dtype)
    header = header_like.copy()

    # Remove BSCALE/BZERO so data are stored directly
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    header.add_history("Normalized by normalize.py")

    out_dir = os.path.dirname(os.path.abspath(fname))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    hdu = fits.PrimaryHDU(data=out_data, header=header)
    hdu.writeto(fname, overwrite=True)


# ---------- Regression helpers ----------

def sample_pixels(ref, img, max_samples=2_000_000):
    """Collect paired finite pixels from ref and img, with optional subsampling."""
    mask = np.isfinite(ref) & np.isfinite(img)
    x = ref[mask].ravel()
    y = img[mask].ravel()
    n = x.size
    if n == 0:
        raise RuntimeError("No valid pixels for regression.")
    if n > max_samples:
        idx = np.random.choice(n, size=max_samples, replace=False)
        x = x[idx]
        y = y[idx]
    return x, y


def linear_fit(ref, img):
    """Ordinary least squares fit: img = B * ref + C."""
    x, y = sample_pixels(ref, img)
    xm = x.mean()
    ym = y.mean()
    dx = x - xm
    denom = np.sum(dx * dx)
    if denom <= 0:
        return 1.0, ym - xm
    B = np.sum(dx * (y - ym)) / denom
    C = ym - B * xm
    return B, C


def robust_fit_sigma_clipped(ref, img, max_iter=5, clip_sigma=3.0):
    """Robust linear fit with iterative sigma-clipping."""
    x, y = sample_pixels(ref, img)

    # Initial OLS
    xm = x.mean()
    ym = y.mean()
    dx = x - xm
    denom = np.sum(dx * dx)
    if denom <= 0:
        return 1.0, ym - xm
    B = np.sum(dx * (y - ym)) / denom
    C = ym - B * xm

    for _ in range(max_iter):
        model = B * x + C
        resid = y - model
        med = np.median(resid)
        mad = np.median(np.abs(resid - med))
        if mad == 0:
            break
        sigma = 1.4826 * mad
        good = np.abs(resid - med) <= clip_sigma * sigma
        if np.count_nonzero(good) < 10:
            break

        xg = x[good]
        yg = y[good]
        xm = xg.mean()
        ym = yg.mean()
        dx = xg - xm
        denom = np.sum(dx * dx)
        if denom <= 0:
            break
        B = np.sum(dx * (yg - ym)) / denom
        C = ym - B * xm

    return B, C


# ---------- Global iterative normalization (method 3) ----------

def global_iterative_normalization(data_list, max_iter=5):
    """
    Global normalization across all frames.

    Model: I_i = B_i * R + C_i

    Steps:
      1) Initialize (B_i, C_i) using frame 0 as temporary reference.
      2) Build common reference R = mean of normalized frames.
      3) Refit each frame to R.
      4) Repeat.

    Gauge: frame 0 is fixed to B_0 = 1, C_0 = 0.
    Returns:
        B_list, C_list for all frames.
    """
    n = len(data_list)
    if n == 0:
        raise RuntimeError("No frames for global normalization.")

    B = [1.0] + [1.0] * (n - 1)
    C = [0.0] + [0.0] * (n - 1)

    ref0 = data_list[0]
    # Initial guess vs frame 0
    for i in range(1, n):
        bi, ci = linear_fit(ref0, data_list[i])
        if bi == 0:
            bi = 1.0
        B[i] = bi
        C[i] = ci

    for _ in range(max_iter):
        # Build common reference from current normalization
        norm_stack = []
        for i in range(n):
            bi = B[i] if B[i] != 0 else 1.0
            norm_stack.append((data_list[i] - C[i]) / bi)
        R = np.mean(norm_stack, axis=0)

        # Refit, keeping frame 0 as gauge
        B[0] = 1.0
        C[0] = 0.0
        for i in range(1, n):
            bi, ci = linear_fit(R, data_list[i])
            if bi == 0:
                bi = 1.0
            B[i] = bi
            C[i] = ci

    return B, C


# ---------- Argument parsing ----------

def parse_args(argv):
    if len(argv) < 3:
        print(USAGE)
        sys.exit(1)

    input_spec = argv[1]
    output_spec = argv[2]

    basefile = None
    method = 1

    if len(argv) >= 4:
        arg3 = argv[3]
        if arg3.isdigit():
            method = int(arg3)
        else:
            basefile = arg3

    if len(argv) >= 5:
        arg4 = argv[4]
        if not arg4.isdigit():
            print("Method argument must be an integer (1, 2 or 3).")
            sys.exit(1)
        method = int(arg4)

    if method not in (1, 2, 3):
        print("Invalid method code. Must be 1, 2 or 3.")
        sys.exit(1)

    return input_spec, output_spec, basefile, method


# ---------- Main ----------

def main():
    input_spec, output_spec, basefile, method = parse_args(sys.argv)

    # Build (input, output) pairs using shared helper.
    try:
        io_pairs = bu.build_io_file_lists(input_spec, output_spec)
    except Exception as e:
        print(f"Error building IO file lists: {e}")
        sys.exit(1)

    input_files = [p[0] for p in io_pairs]
    output_files = [p[1] for p in io_pairs]
    n_files = len(input_files)

    if n_files == 0:
        print("No input files found.")
        sys.exit(1)

    # Load all input frames
    data_list = []
    headers = []
    orig_dtypes = []
    for fname in input_files:
        data, hdr, dt = read_fits_data(fname)
        data_list.append(data)
        headers.append(hdr)
        orig_dtypes.append(dt)

    # Prepare reference for methods 1 and 2
    ref_data = None
    if method in (1, 2):
        if basefile is not None:
            if basefile in input_files:
                idx = input_files.index(basefile)
                ref_data = data_list[idx]
            else:
                if not os.path.exists(basefile):
                    print(f"Base file not found: {basefile}")
                    sys.exit(1)
                ref_data, _, _ = read_fits_data(basefile)
        else:
            ref_data = data_list[0]

    # Compute normalization coefficients
    if method == 1:
        B_list = []
        C_list = []
        for img in data_list:
            if ref_data is img:
                B_list.append(1.0)
                C_list.append(0.0)
            else:
                B, C = linear_fit(ref_data, img)
                if B == 0:
                    B = 1.0
                B_list.append(B)
                C_list.append(C)

    elif method == 2:
        B_list = []
        C_list = []
        for img in data_list:
            if ref_data is img:
                B_list.append(1.0)
                C_list.append(0.0)
            else:
                B, C = robust_fit_sigma_clipped(ref_data, img)
                if B == 0:
                    B = 1.0
                B_list.append(B)
                C_list.append(C)

    else:
        # Method 3: global iterative normalization (basefile ignored)
        B_list, C_list = global_iterative_normalization(data_list)

    # Apply normalization and write outputs
    print(f"Normalizing {n_files} file(s) using method {method}...")
    for i in range(n_files):
        B = B_list[i]
        C = C_list[i]
        if B == 0:
            B = 1.0
        norm = (data_list[i] - C) / B
        write_fits_data(output_files[i], norm, headers[i], orig_dtypes[i])
        print_progress(i + 1, n_files, prefix="Writing: ")

    print("Done.")


if __name__ == "__main__":
    main()
