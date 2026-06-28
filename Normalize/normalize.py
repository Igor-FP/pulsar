#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils as bu

USAGE = """Usage:
  normalize.py input_spec output_spec [basefile.fit] [method] [--sat F] [--zero]

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
  --sat F      - optional pre-cut: before fitting, ignore pixels brighter
                 than F * P99.5 (per frame), where P99.5 is the 99.5th
                 percentile value (a robust "max"). F in (0, 1]. Off by
                 default. A coarse guard on top of the residual clipping.
  --zero, -z   - preserve no-data zeros: after the fit (which already ignores
                 zero pixels in all methods) and the affine transform, restore
                 the input frame's exact-zero pixels to 0 in the output, so
                 normalization never shifts a 0 (registration border / no-data)
                 to a non-zero value. Off by default.

Notes:
  - Model: I = B * R + C
    Normalized output is: (I - C) / B
  - Zero pixels (no-data borders) are excluded from the fit in all methods.
    With --zero they are additionally restored to 0 in the output (otherwise
    the affine maps 0 -> -C / B, a small non-zero value).
  - Saturation is handled by residual sigma-clipping (methods 2 and 3),
    not a fixed threshold: after calibration the clip level varies across
    the field (raw saturation scaled by 1/flat), so it is not a constant.
  - Method 1 (plain OLS) excludes zeros but does not reject saturation;
    use method 2 or 3 for frames with saturated stars.
  - Method 3 ignores basefile if provided.

Examples:
  normalize.py light0001.fit norm0001.fit
  normalize.py light0001.fit norm0001.fit reference.fit 2
  normalize.py *.fit norm0001.fit 3
  normalize.py *.fit norm0001.fit 2 --sat 0.8
  normalize.py *.fit norm0001.fit 2 --zero
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
        arr = np.clip(data, info.min, info.max)   # also bounds +/-Inf to max/min
        arr = np.rint(arr)
        # NaN passes through clip/rint unchanged and casts to garbage on integers
        # (e.g. INT_MIN); never write NaN to a FITS file.
        arr[~np.isfinite(arr)] = 0
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


def write_fits_data(fname, data, header_like, orig_dtype, restore_zeros=False):
    """Write FITS file preserving header, using original dtype and safe clamping."""
    out_data = convert_to_original_dtype(data, orig_dtype)
    header = header_like.copy()

    # Remove BSCALE/BZERO so data are stored directly
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    header.add_history("Normalized by normalize.py")
    if restore_zeros:
        header.add_history("  input zero pixels preserved (--zero)")

    out_dir = os.path.dirname(os.path.abspath(fname))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    hdu = fits.PrimaryHDU(data=out_data, header=header)
    hdu.writeto(fname, overwrite=True)


# ---------- Regression helpers ----------

# Robust "max" for the optional --sat cut. A high percentile resists a single
# hot/cosmic pixel that would otherwise set the scale; 99.9 can still be
# dominated by hot pixels on some CCDs, so 99.5 is used.
_SAT_PERCENTILE = 99.5


def sample_pixels(ref, img, max_samples=2_000_000, sat_frac=None):
    """Collect paired valid pixels from ref and img, with optional subsampling.

    Excluded from the fit:
      - non-finite values (NaN/Inf);
      - zero pixels in EITHER frame. Zero is the no-data sentinel in this
        project (registration borders, masked/unexposed areas). A matched
        (0, 0) population otherwise anchors the regression at the origin via
        leverage and biases the offset C toward 0; mismatched (0, y)/(x, 0)
        border pixels act as pure outliers.
      - if sat_frac is given: pixels brighter than sat_frac * P99.5 in EITHER
        frame, where P99.5 is the 99.5th-percentile value of that frame's
        valid pixels (a robust "max"). An optional, blunt value pre-cut.

    Saturation by value is OFF by default: after calibration the clip level is
    (raw_sat - bias - dark) / flat, which varies across the field and is not a
    constant, so no fixed threshold tracks it. The always-on saturation
    handling is the residual sigma-clipping in robust_fit_sigma_clipped
    (methods 2 and 3); --sat is only a coarse extra guard.
    """
    mask = np.isfinite(ref) & np.isfinite(img) & (ref != 0) & (img != 0)
    if sat_frac is not None and np.any(mask):
        ref_hi = np.percentile(ref[mask], _SAT_PERCENTILE)
        img_hi = np.percentile(img[mask], _SAT_PERCENTILE)
        mask &= (ref <= sat_frac * ref_hi) & (img <= sat_frac * img_hi)
    x = ref[mask].ravel()
    y = img[mask].ravel()
    n = x.size
    if n == 0:
        raise RuntimeError("No overlapping non-zero pixels for regression.")
    if n > max_samples:
        idx = np.random.choice(n, size=max_samples, replace=False)
        x = x[idx]
        y = y[idx]
    return x, y


def linear_fit(ref, img, sat_frac=None):
    """Ordinary least squares fit: img = B * ref + C."""
    x, y = sample_pixels(ref, img, sat_frac=sat_frac)
    xm = x.mean()
    ym = y.mean()
    dx = x - xm
    denom = np.sum(dx * dx)
    if denom <= 0:
        return 1.0, ym - xm
    B = np.sum(dx * (y - ym)) / denom
    C = ym - B * xm
    return B, C


def robust_fit_sigma_clipped(ref, img, max_iter=5, clip_sigma=3.0, sat_frac=None):
    """Robust linear fit with iterative sigma-clipping."""
    x, y = sample_pixels(ref, img, sat_frac=sat_frac)

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

def global_iterative_normalization(data_list, max_iter=5, sat_frac=None):
    """
    Global normalization across all frames.

    Model: I_i = B_i * R + C_i

    Steps:
      1) Initialize (B_i, C_i) using frame 0 as temporary reference.
      2) Build common reference R = mean of normalized frames.
      3) Refit each frame to R.
      4) Repeat.

    All fits use robust_fit_sigma_clipped, so zero pixels (excluded in
    sample_pixels) and saturated pixels do not bias the coefficients.

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
    # Initial guess vs frame 0 (robust: reject saturation/outliers by residual)
    for i in range(1, n):
        bi, ci = robust_fit_sigma_clipped(ref0, data_list[i], sat_frac=sat_frac)
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
            bi, ci = robust_fit_sigma_clipped(R, data_list[i], sat_frac=sat_frac)
            if bi == 0:
                bi = 1.0
            B[i] = bi
            C[i] = ci

    return B, C


# ---------- Argument parsing ----------

def parse_args(argv):
    argv = list(argv)

    if "-h" in argv or "--help" in argv:
        print(USAGE)
        sys.exit(0)

    sat_frac = None
    if "--sat" in argv:
        i = argv.index("--sat")
        if i + 1 >= len(argv):
            print("Error: --sat requires a fraction in (0, 1].")
            sys.exit(1)
        try:
            sat_frac = float(argv[i + 1])
        except ValueError:
            print("Error: --sat value must be a number in (0, 1].")
            sys.exit(1)
        if not (0.0 < sat_frac <= 1.0):
            print("Error: --sat must be in (0, 1].")
            sys.exit(1)
        del argv[i:i + 2]

    restore_zeros = ("--zero" in argv) or ("-z" in argv)
    argv = [a for a in argv if a not in ("--zero", "-z")]

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

    return input_spec, output_spec, basefile, method, sat_frac, restore_zeros


# ---------- Main ----------

def main():
    input_spec, output_spec, basefile, method, sat_frac, restore_zeros = parse_args(sys.argv)

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
                B, C = linear_fit(ref_data, img, sat_frac=sat_frac)
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
                B, C = robust_fit_sigma_clipped(ref_data, img, sat_frac=sat_frac)
                if B == 0:
                    B = 1.0
                B_list.append(B)
                C_list.append(C)

    else:
        # Method 3: global iterative normalization (basefile ignored)
        B_list, C_list = global_iterative_normalization(data_list, sat_frac=sat_frac)

    # Apply normalization and write outputs
    print(f"Normalizing {n_files} file(s) using method {method}...")
    for i in range(n_files):
        B = B_list[i]
        C = C_list[i]
        if B == 0:
            B = 1.0
        norm = (data_list[i] - C) / B
        if restore_zeros:
            # --zero: keep the input frame's no-data zeros at 0 instead of
            # letting the affine (I - C) / B shift them to -C / B.
            norm[data_list[i] == 0] = 0.0
        write_fits_data(output_files[i], norm, headers[i], orig_dtypes[i],
                        restore_zeros=restore_zeros)
        print_progress(i + 1, n_files, prefix="Writing: ")

    print("Done.")


if __name__ == "__main__":
    main()
