#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
from astropy.io import fits

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils

# Global debug flag
DEBUG = False


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  autoflat.py [-d] input.fit output.fit [poly_order]\n"
        "\n"
        "    -d           - enable debug outputs (save intermediate FITS images)\n"
        "    input.fit    - single file OR numbered pattern OR list file (via batch_utils)\n"
        "    output.fit   - single file OR numbered pattern for results\n"
        "    poly_order   - optional polynomial order for background surface (default: 1)\n"
        "                   0 = constant, 1 = plane, 2+ = higher order 2D polynomial\n"
    )
    sys.exit(1)


def parse_args(argv):
    debug = False
    args = argv[1:]

    if args and args[0] == "-d":
        debug = True
        args = args[1:]

    if len(args) not in (2, 3):
        usage()

    input_pattern = args[0]
    output_pattern = args[1]

    if len(args) == 3:
        try:
            poly_order = int(args[2])
            if poly_order < 0:
                raise ValueError
        except ValueError:
            sys.stderr.write("Error: poly_order must be a non-negative integer.\n")
            sys.exit(1)
    else:
        poly_order = 1

    return input_pattern, output_pattern, poly_order, debug


def log(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def log_stats(tag, arr):
    arr = np.asarray(arr)
    if arr.size == 0:
        log(f"  [{tag}] empty")
        return
    vmin = float(np.nanmin(arr))
    vmax = float(np.nanmax(arr))
    mean = float(np.nanmean(arr))
    log(f"  [{tag}] min={vmin:.3f} max={vmax:.3f} mean={mean:.3f}")


# =========================
# FITS helpers
# =========================

def to_float64(data):
    return np.asarray(data, dtype=np.float64)


def from_float64(data64, orig_dtype):
    if np.issubdtype(orig_dtype, np.floating):
        return data64.astype(orig_dtype)

    if np.issubdtype(orig_dtype, np.signedinteger) or np.issubdtype(orig_dtype, np.unsignedinteger):
        info = np.iinfo(orig_dtype)
        arr = np.rint(data64)
        arr = np.clip(arr, info.min, info.max)
        return arr.astype(orig_dtype)

    return data64.astype(np.float32)


def make_header_for_shape(template_header, shape):
    h = template_header.copy() if template_header is not None else fits.Header()
    if len(shape) == 2:
        h["NAXIS"] = 2
        h["NAXIS1"] = int(shape[1])
        h["NAXIS2"] = int(shape[0])
    return h


def write_debug_image(base_path, tag, data, header=None):
    if not DEBUG or base_path is None:
        return
    try:
        root, _ = os.path.splitext(base_path)
        path = f"{root}_{tag}.fit"
        h = make_header_for_shape(header, data.shape)
        fits.PrimaryHDU(data=data.astype(np.float32), header=h).writeto(path, overwrite=True)
        log(f"  [debug] Saved {tag} -> '{path}'")
    except Exception as e:
        log(f"  [debug] Failed to save {tag}: {e}")


# =========================
# Step 1: Zero expansion
# =========================

def expand_zero_mask(data, iterations=2):
    log(f"  [step1] Expanding zero mask, iterations={iterations}")
    mask = (data == 0.0)
    if not np.any(mask):
        log("  [step1] No zero pixels found, skipping")
        log_stats("step1_out", data)
        return data

    h, w = data.shape
    for it in range(iterations):
        log(f"    [step1] Iteration {it + 1}, size={w}x{h}")
        padded = np.pad(mask, 1, mode="edge")
        neigh = (
            padded[0:h,     0:w]     | padded[0:h,     1:w+1] |
            padded[0:h,     2:w+2]   | padded[1:h+1,   0:w]   |
            padded[1:h+1,   1:w+1]   | padded[1:h+1,   2:w+2] |
            padded[2:h+2,   0:w]     | padded[2:h+2,   1:w+1] |
            padded[2:h+2,   2:w+2]
        )
        mask = neigh

    out = data.copy()
    out[mask] = 0.0
    log("  [step1] Completed")
    log_stats("step1_out", out)
    return out


# =========================
# Step 2: Median (fast, zero-aware)
# =========================

try:
    from scipy.ndimage import median_filter, distance_transform_edt
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False
    log("  [init] SciPy not available, will use slow fallbacks")


def fill_zeros_with_nearest(data):
    """Fill zeros with nearest non-zero using EDT (if available)."""
    mask = (data != 0.0)
    if np.all(mask):
        return data

    if HAVE_SCIPY:
        # Build array where non-zero pixels are features (0), others 1
        feat = np.where(mask, 0, 1)
        _, indices = distance_transform_edt(feat, return_indices=True)
        filled = data.copy()
        # For zero pixels, copy value from nearest non-zero
        yy = indices[0][~mask]
        xx = indices[1][~mask]
        filled[~mask] = data[yy, xx]
        return filled

    # Fallback: simple brute force (should be rare)
    filled = data.copy()
    nz_y, nz_x = np.nonzero(mask)
    if nz_y.size == 0:
        return filled
    nz_points = np.vstack((nz_y, nz_x)).T
    zero_y, zero_x = np.nonzero(~mask)
    for y, x in zip(zero_y, zero_x):
        dy = nz_points[:, 0] - y
        dx = nz_points[:, 1] - x
        j = np.argmin(dy * dy + dx * dx)
        yy, xx = nz_points[j]
        filled[y, x] = data[yy, xx]
    return filled


def median_aperture_ignore_zeros(data, radius=2):
    """
    Fast median filter:
      - Uses SciPy median_filter (C implementation) when available.
      - Zeros are ignored via sentinel trick:
          zeros -> sentinel >> real values
          median over window
          if median >= sentinel -> result 0 (no signal)
      - Then zeros are filled from nearest non-zero.
    """
    size = 2 * radius + 1
    h, w = data.shape
    log(f"  [step2] Median(ignore zeros) radius={radius}, size={size}, image={w}x{h}")

    if not HAVE_SCIPY:
        # Slow exact fallback
        log("  [step2] SciPy missing, using slow pure-numpy median(ignore zeros)")
        pad = radius
        padded = np.pad(data, pad, mode="edge")
        out = np.zeros_like(data, dtype=np.float64)
        for y in range(h):
            if y % 64 == 0:
                log(f"    [step2] row {y + 1}/{h}")
            ys = y
            ye = y + 2 * radius + 1
            for x in range(w):
                xs = x
                xe = x + 2 * radius + 1
                window = padded[ys:ye, xs:xe]
                vals = window[window != 0.0]
                out[y, x] = np.median(vals) if vals.size > 0 else 0.0
        out = fill_zeros_with_nearest(out)
        log("  [step2] Completed (pure numpy)")
        log_stats("step2_out", out)
        return out

    # SciPy fast path
    zero_mask = (data == 0.0)
    if not np.any(zero_mask):
        # No zeros -> plain median
        out = median_filter(data, size=size, mode="nearest")
        log("  [step2] Completed (no zeros, plain median_filter)")
        log_stats("step2_out", out)
        return out

    nz = data[~zero_mask]
    if nz.size == 0:
        log("  [step2] All pixels zero, returning zeros")
        return data.copy()

    vmax = float(np.max(nz))
    sentinel = vmax + (abs(vmax) + 1.0)

    work = np.where(zero_mask, sentinel, data)
    med = median_filter(work, size=size, mode="nearest")

    # Windows без сигнала -> медиана >= sentinel -> считаем 0
    med[med >= sentinel] = 0.0

    out = fill_zeros_with_nearest(med)
    log("  [step2] Completed (median_filter + sentinel)")
    log_stats("step2_out", out)
    return out


# =========================
# Step 3: Min-binning 2x2 (darkest)
# =========================

def min_pool_2x_ignore_zeros(data):
    """
    2x2 min pooling ignoring zeros:
      - if есть не-нулевые -> min по ним;
      - если все нули -> 0.
    """
    h, w = data.shape
    if h == 1 and w == 1:
        return data.copy()

    pad_h = h % 2
    pad_w = w % 2
    if pad_h or pad_w:
        padded = np.pad(data, ((0, pad_h), (0, pad_w)), mode="edge")
    else:
        padded = data

    a = padded[0::2, 0::2]
    b = padded[0::2, 1::2]
    c = padded[1::2, 0::2]
    d = padded[1::2, 1::2]

    stack = np.stack([a, b, c, d], axis=0)
    nonzero = (stack != 0.0)
    any_nz = np.any(nonzero, axis=0)

    # Replace zeros with +inf only where there is at least one non-zero
    large = np.finfo(np.float64).max
    stack_adj = np.where(nonzero, stack, large)
    out = np.min(stack_adj, axis=0)
    out[~any_nz] = 0.0

    return out


def iterative_min_binning(data, target_min_size=8):
    out = data.copy()
    it = 0
    while min(out.shape) > target_min_size:
        it += 1
        h, w = out.shape
        log(f"  [step3] Min-binning iter {it}, size={w}x{h}")
        out = min_pool_2x_ignore_zeros(out)
        log_stats(f"step3_iter{it}", out)
    log("  [step3] Completed")
    log_stats("step3_out", out)
    return out


# =========================
# Step 4: Polynomial fit & render
# =========================

def build_poly_terms(order):
    return [(i, j) for i in range(order + 1) for j in range(order + 1 - i)]


def fit_poly2d(data, order):
    """
    Fit 2D polynomial surface of given order to data != 0.
    If requested order is too high for the number of points, it is reduced
    so that number of terms <= number of non-zero pixels.
    Coordinates are normalized to [-1, 1] for stability.
    """
    h, w = data.shape
    log(f"  [step4] Requested poly order={order} on reduced image {w}x{h}")

    yy, xx = np.mgrid[0:h, 0:w]
    mask = (data != 0.0)
    N = int(mask.sum())
    if N == 0:
        log("  [step4] No non-zero pixels for fit, aborting")
        return None, None

    # Choose maximal allowed order so that num_terms <= N
    def num_terms(k):
        return (k + 1) * (k + 2) // 2

    used_order = min(order, 20)  # hard safety cap, optional
    while used_order > 0 and num_terms(used_order) > N:
        used_order -= 1

    if used_order < order:
        log(f"  [step4] Reducing poly order from {order} to {used_order} "
            f"(points={N}, terms={num_terms(used_order)})")
    else:
        log(f"  [step4] Using poly order={used_order}, points={N}, terms={num_terms(used_order)}")

    order = used_order

    z = data[mask]
    x = xx[mask].astype(np.float64)
    y = yy[mask].astype(np.float64)

    if w > 1:
        x = 2.0 * (x / (w - 1.0)) - 1.0
    else:
        x = np.zeros_like(x)
    if h > 1:
        y = 2.0 * (y / (h - 1.0)) - 1.0
    else:
        y = np.zeros_like(y)

    # Build terms for chosen order
    terms = [(i, j) for i in range(order + 1) for j in range(order + 1 - i)]
    cols = []
    for (i, j) in terms:
        cols.append((x ** i) * (y ** j))
    A = np.vstack(cols).T

    coeffs, _, _, _ = np.linalg.lstsq(A, z, rcond=None)
    log("  [step4] Polynomial fit completed")
    return coeffs, terms


def fit_poly2d_OLD(data, order):
    h, w = data.shape
    log(f"  [step4] Fitting poly order={order} on {w}x{h}")
    yy, xx = np.mgrid[0:h, 0:w]
    mask = (data != 0.0)
    if not np.any(mask):
        log("  [step4] No non-zero pixels, abort fit")
        return None, None

    z = data[mask]
    x = xx[mask].astype(np.float64)
    y = yy[mask].astype(np.float64)

    if w > 1:
        x = 2.0 * (x / (w - 1.0)) - 1.0
    else:
        x = np.zeros_like(x)
    if h > 1:
        y = 2.0 * (y / (h - 1.0)) - 1.0
    else:
        y = np.zeros_like(y)

    terms = build_poly_terms(order)
    log(f"  [step4] Terms={len(terms)}, points={z.size}")

    cols = []
    for i, j in terms:
        cols.append((x ** i) * (y ** j))
    A = np.vstack(cols).T

    coeffs, _, _, _ = np.linalg.lstsq(A, z, rcond=None)
    log("  [step4] Fit done")
    return coeffs, terms


def render_poly2d(height, width, coeffs, terms):
    log(f"  [step6] Rendering model {width}x{height}")
    yy, xx = np.mgrid[0:height, 0:width]

    if width > 1:
        x = 2.0 * (xx / (width - 1.0)) - 1.0
    else:
        x = np.zeros_like(xx, dtype=np.float64)
    if height > 1:
        y = 2.0 * (yy / (height - 1.0)) - 1.0
    else:
        y = np.zeros_like(yy, dtype=np.float64)

    model = np.zeros((height, width), dtype=np.float64)
    for idx, (c, (i, j)) in enumerate(zip(coeffs, terms), start=1):
        if c != 0.0:
            model += c * (x ** i) * (y ** j)
        if idx % 10 == 0:
            log(f"    [step6] applied {idx} terms")

    log_stats("step6_model", model)
    log("  [step6] Model done")
    return model


# =========================
# Pipeline
# =========================

def autoflat_frame(data64, poly_order, debug_prefix=None, header=None):
    log("  [pipeline] Start")

    step1 = expand_zero_mask(data64, iterations=2)
    write_debug_image(debug_prefix, "step1_zeroexp", step1, header)

    bg = median_aperture_ignore_zeros(step1, radius=2)
    write_debug_image(debug_prefix, "step2_median", bg, header)

    small = iterative_min_binning(bg, target_min_size=8)
    write_debug_image(debug_prefix, "step3_binning", small, header)

    coeffs, terms = fit_poly2d(small, poly_order)
    if coeffs is None:
        log("  [pipeline] Fit failed, return original")
        return data64.copy()

    nz = small[small != 0.0]
    offset = float(np.median(nz)) if nz.size > 0 else 0.0
    log(f"  [step5] OFFSET={offset:.6f}")

    h, w = data64.shape
    model_full = render_poly2d(h, w, coeffs, terms)
    write_debug_image(debug_prefix, "step6_model", model_full, header)

    log("  [step7] Applying correction: orig - model + OFFSET")
    corrected = data64 - model_full + offset
    log_stats("step7_out", corrected)

    log("  [pipeline] Done")
    return corrected


# =========================
# Per-file
# =========================

def process_file(infile, outfile, poly_order):
    log(f"[file] {infile} -> {outfile}")
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None or hdul[0].data.ndim != 2:
            raise ValueError(f"'{infile}' is not a 2D FITS image")

        orig = hdul[0].data
        header = hdul[0].header

        data64 = to_float64(orig)
        log_stats("input", data64)

        debug_prefix = outfile if DEBUG else None
        corrected64 = autoflat_frame(data64, poly_order,
                                     debug_prefix=debug_prefix,
                                     header=header)

        out_data = from_float64(corrected64, orig.dtype)
        hdul[0].data = out_data
        hdul[0].header = header
        hdul.writeto(outfile, overwrite=True)

    log(f"[file] Done {infile}")


# =========================
# Main
# =========================

def main():
    global DEBUG

    input_pattern, output_pattern, poly_order, debug = parse_args(sys.argv)
    DEBUG = debug

    try:
        io_pairs = batch_utils.build_io_file_lists(input_pattern, output_pattern)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    log(f"Found {len(io_pairs)} file(s). Poly order = {poly_order}, debug = {DEBUG}")

    total = len(io_pairs)
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            process_file(infile, outfile, poly_order)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")
    sys.stderr.flush()


if __name__ == "__main__":
    main()
