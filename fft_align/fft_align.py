#!/usr/bin/env python3
import argparse
import math
import os
import sys
import shutil
import multiprocessing as mp

import numpy as np
from astropy.io import fits
from scipy import ndimage

# ----------------------------------------------------------
#  Locate and import batch_utils.py (../lib relative to this script)
# ----------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "lib"))
if LIB_DIR not in sys.path:
    sys.path.insert(0, LIB_DIR)

import batch_utils  # type: ignore


# ==========================================================
#                     FITS I/O
# ==========================================================

def read_fits(path):
    """Read FITS as float32 ndarray and return (array, header)."""
    hdul = fits.open(path, memmap=True)
    data = hdul[0].data.astype(np.float32)
    header = hdul[0].header
    hdul.close()
    return data, header


def write_fits(path, data, header):
    """Write float32 ndarray to FITS with given header."""
    hdr = header.copy()
    hdu = fits.PrimaryHDU(data=data.astype(np.float32), header=hdr)
    hdu.writeto(path, overwrite=True)


# ==========================================================
#               NORMALIZATION (MAD-based)
# ==========================================================

def robust_normalize(a):
    """Fallback percentile normalization (1..99%)."""
    a = a.astype(np.float32)
    lo, hi = np.percentile(a, [1, 99])
    if hi <= lo:
        hi = lo + 1.0
    a = np.clip(a, lo, hi)
    return (a - lo) / (hi - lo)


def robust_normalize_mad(a):
    """
    Normalize image using median and MAD:
        norm = (a - median) / (1.4826 * MAD)
    """
    a = a.astype(np.float32)
    med = np.median(a)
    mad = np.median(np.abs(a - med))
    if mad < 1e-12:
        return robust_normalize(a)
    return (a - med) / (1.4826 * mad)


# ==========================================================
#                 DOWNSAMPLING
# ==========================================================

def downsample_to_max(a, target_max=256):
    """
    Downsample image so that max(width, height) ~= target_max.
    Returns: (downsampled_image, scale_inverse).
    """
    h, w = a.shape
    max_dim = max(h, w)
    if max_dim <= target_max:
        return a.copy(), 1.0

    scale = target_max / float(max_dim)
    a_small = ndimage.zoom(a, zoom=scale, order=1)
    return a_small, 1.0 / scale


def make_hann_window(shape):
    """2D Hann window for FFT edge suppression."""
    h, w = shape
    return np.outer(np.hanning(h), np.hanning(w)).astype(np.float32)


# ==========================================================
#                PHASE CORRELATION
# ==========================================================

def phase_correlate(a, b):
    """
    Compute subpixel shift (dx, dy) between images a and b using phase correlation.
    Returns (dx, dy, peak).

    NOTE:
      If b is a shifted version of a by (+sy, +sx), the returned (dx,dy)
      will be approximately (-sx, -sy).
    """
    a0 = a - np.mean(a)
    b0 = b - np.mean(b)

    fa = np.fft.fft2(a0)
    fb = np.fft.fft2(b0)

    denom = np.abs(fa * np.conj(fb))
    denom[denom < 1e-8] = 1e-8
    r = fa * np.conj(fb) / denom

    cc = np.abs(np.fft.ifft2(r))
    maxpos = np.unravel_index(np.argmax(cc), cc.shape)
    peak = cc[maxpos]

    shifts = np.array(maxpos, dtype=np.float32)
    h, w = a.shape

    if shifts[0] > h / 2:
        shifts[0] -= h
    if shifts[1] > w / 2:
        shifts[1] -= w

    dy, dx = shifts
    return float(dx), float(dy), float(peak)


# ==========================================================
#                  GEOMETRIC WARPS
# ==========================================================

def warp_rotate_scale(image, angle_deg, scale, order=1):
    """
    Rotate + isotropic scale around center (inverse mapping).
    order=1 for pyramid, order>=3 for final high-quality resampling.
    """
    h, w = image.shape
    center = np.array([(h - 1) / 2.0, (w - 1) / 2.0], dtype=np.float32)

    theta = np.deg2rad(angle_deg)
    c = math.cos(theta)
    s = math.sin(theta)

    R = np.array([[c, -s],
                  [s,  c]], dtype=np.float32)

    M = (1.0 / scale) * R.T
    offset = center - M @ center

    return ndimage.affine_transform(
        image,
        matrix=M,
        offset=offset,
        order=order,
        mode="constant",
        cval=0.0,
        prefilter=(order > 1),
    )


def warp_shift(image, dx, dy, order=3):
    """Subpixel translation dx,dy."""
    offset = [-dy, -dx]
    return ndimage.affine_transform(
        image,
        matrix=np.eye(2, dtype=np.float32),
        offset=offset,
        order=order,
        mode="constant",
        cval=0.0,
        prefilter=(order > 1),
    )


# ==========================================================
#                 GRID SEARCH (ANGLE+SCALE)
# ==========================================================

def grid_search_angle_scale(fixed_small, moving_small,
                            angle_range, angle_step,
                            scale_range, scale_step):
    """
    Brute-force search over angle + scale on small images.

    Returns best (angle, scale, dx, dy, peak).
    """
    h0, w0 = fixed_small.shape
    if moving_small.shape != fixed_small.shape:
        zoom_y = h0 / moving_small.shape[0]
        zoom_x = w0 / moving_small.shape[1]
        zoom = 0.5 * (zoom_y + zoom_x)
        moving_small = ndimage.zoom(moving_small, zoom=zoom, order=1)
        if moving_small.shape != fixed_small.shape:
            tmp = np.zeros_like(fixed_small)
            yy = min(tmp.shape[0], moving_small.shape[0])
            xx = min(tmp.shape[1], moving_small.shape[1])
            tmp[:yy, :xx] = moving_small[:yy, :xx]
            moving_small = tmp

    win = make_hann_window(fixed_small.shape)
    fixed_w = fixed_small * win
    moving_base = moving_small * win

    best_peak = -1.0
    best_params = (0.0, 1.0, 0.0, 0.0)

    angles = np.arange(angle_range[0], angle_range[1] + 1e-9, angle_step)
    scales = np.arange(scale_range[0], scale_range[1] + 1e-9, scale_step)

    for angle in angles:
        for scale in scales:
            warped = warp_rotate_scale(moving_base, angle, scale, order=1)
            dx, dy, peak = phase_correlate(fixed_w, warped)
            if peak > best_peak:
                best_peak = peak
                best_params = (float(angle), float(scale), float(dx), float(dy))

    return (*best_params, best_peak)


# ==========================================================
#            SUPERFINE HIERARCHICAL REGISTRATION
# ==========================================================

def make_pyramid_targets(max_dim, top_size=128):
    """Build list of target sizes (top_size, 2*top_size, ..., full)."""
    sizes = []
    s = top_size
    while s < max_dim:
        sizes.append(s)
        s *= 2
    sizes.append(max_dim)
    return sizes


def register_superfine(fixed, moving,
                       max_angle, angle_step,
                       scale_delta, scale_step,
                       top_size=128,
                       normalize=True):
    """
    Multi-level pyramid registration:
        128 → 256 → 512 → ... → full resolution.

    At each level:
      - search in a narrow band around previous best
      - steps are halved on each level
    """
    if normalize:
        fixed_reg = robust_normalize_mad(fixed)
        moving_reg = robust_normalize_mad(moving)
    else:
        fixed_reg = robust_normalize(fixed)
        moving_reg = robust_normalize(moving)

    h, w = fixed_reg.shape
    max_dim = max(h, w)
    levels = make_pyramid_targets(max_dim, top_size)
    peaks = []

    angle_center = 0.0
    scale_center = 1.0
    cur_angle_step = angle_step
    cur_scale_step = scale_step
    first = True

    for tgt in levels:
        fixed_small, _ = downsample_to_max(fixed_reg, tgt)
        moving_small, _ = downsample_to_max(moving_reg, tgt)

        if first:
            ang_range = (-max_angle, max_angle)
            scl_range = (1.0 - scale_delta, 1.0 + scale_delta)
        else:
            ang_range = (angle_center - 2 * cur_angle_step,
                         angle_center + 2 * cur_angle_step)
            scl_range = (scale_center - 2 * cur_scale_step,
                         scale_center + 2 * cur_scale_step)

        angle_best, scale_best, _, _, peak = grid_search_angle_scale(
            fixed_small, moving_small,
            angle_range=ang_range, angle_step=cur_angle_step,
            scale_range=scl_range, scale_step=cur_scale_step,
        )

        peaks.append((tgt, peak))
        angle_center = angle_best
        scale_center = scale_best

        cur_angle_step *= 0.5
        cur_scale_step *= 0.5
        first = False

    angle_final = angle_center
    scale_final = scale_center

    # Final shift estimation at full (normalized) resolution
    moving_rs = warp_rotate_scale(moving_reg, angle_final, scale_final, order=1)
    dx_full, dy_full, peak_full = phase_correlate(fixed_reg, moving_rs)

    return angle_final, scale_final, dx_full, dy_full, (peaks, peak_full)


# ==========================================================
#            COARSE+FINE REGISTRATION
# ==========================================================

def register_coarse_fine(fixed, moving,
                         angle_range, angle_step,
                         scale_range, scale_step,
                         ds_max=256,
                         normalize=True):
    """
    Two-level registration:
      - coarse search at ds_max
      - fine search at ~2*ds_max (up to 512)
      - final subpixel shift estimation at full resolution
    """
    if normalize:
        fixed_reg = robust_normalize_mad(fixed)
        moving_reg = robust_normalize_mad(moving)
    else:
        fixed_reg = robust_normalize(fixed)
        moving_reg = robust_normalize(moving)

    # Coarse
    fixed_c, _ = downsample_to_max(fixed_reg, ds_max)
    moving_c, _ = downsample_to_max(moving_reg, ds_max)

    ang_c, scl_c, _, _, peak_c = grid_search_angle_scale(
        fixed_c, moving_c,
        angle_range, angle_step,
        scale_range, scale_step,
    )

    # Fine
    ds_fine = min(ds_max * 2, 512)
    fixed_f, _ = downsample_to_max(fixed_reg, ds_fine)
    moving_f, _ = downsample_to_max(moving_reg, ds_fine)

    ang_span = 2 * angle_step
    scl_span = 2 * scale_step

    ang_range_f = (ang_c - ang_span, ang_c + ang_span)
    scl_range_f = (scl_c - scl_span, scl_c + scl_span)

    ang_step_f = angle_step * 0.25
    scl_step_f = scale_step * 0.25

    ang_f, scl_f, _, _, peak_f = grid_search_angle_scale(
        fixed_f, moving_f,
        angle_range=ang_range_f, angle_step=ang_step_f,
        scale_range=scl_range_f, scale_step=scl_step_f,
    )

    moving_rs = warp_rotate_scale(moving_reg, ang_f, scl_f, order=1)
    dx_full, dy_full, peak_full = phase_correlate(fixed_reg, moving_rs)

    return ang_f, scl_f, dx_full, dy_full, (peak_c, peak_f, peak_full)


# ==========================================================
#                APPLY FINAL TRANSFORMATION
# ==========================================================

def apply_full_transform(moving, angle, scale, dx, dy, order=3):
    """Apply final (angle,scale,dx,dy) to moving image."""
    img = warp_rotate_scale(moving, angle, scale, order=order)
    return warp_shift(img, dx, dy, order=order)


# ==========================================================
#         POST-CORRECTION: LOCAL AFFINE REFINEMENT
# ==========================================================

def compute_saturation_threshold(img):
    """
    Compute saturation threshold using:
        SAT = median + 0.8 * (max - median)
    Works for any numeric range (float/int).
    """
    med = np.median(img)
    maxv = np.max(img)
    return med + 0.8 * (maxv - med)


def pick_window_centers_auto(fixed_img, win_size, search_radius):
    """
    Automatically pick 4 good window centers for post-correction,
    around nominal positions:
        (25%,25%), (75%,25%), (25%,75%), (75%,75%)

    For each nominal position we scan a 5x5 grid within search_radius,
    maximizing a score combining:
      - local MAD (structure)
      - fraction of saturated pixels
      - gradient magnitude
    """
    h, w = fixed_img.shape
    half = win_size // 2

    gy, gx = np.gradient(fixed_img)
    grad_mag = np.sqrt(gx * gx + gy * gy)

    sat_thr = compute_saturation_threshold(fixed_img)

    centers = []
    frac_positions = [
        (0.25, 0.25),
        (0.25, 0.75),
        (0.75, 0.25),
        (0.75, 0.75),
    ]

    if search_radius <= 0:
        offsets = [0]
    else:
        step = search_radius // 2 if search_radius >= 4 else 1
        offsets = [-2 * step, -step, 0, step, 2 * step]

    for fy, fx in frac_positions:
        cy0 = fy * (h - 1)
        cx0 = fx * (w - 1)

        best_score = -1.0
        best_center = None

        for dy_off in offsets:
            for dx_off in offsets:
                cy = int(round(cy0 + dy_off))
                cx = int(round(cx0 + dx_off))

                y0 = cy - half
                y1 = cy + half
                x0 = cx - half
                x1 = cx + half

                if y0 < 0 or x0 < 0 or y1 > h or x1 > w:
                    continue

                win = fixed_img[y0:y1, x0:x1]
                if win.shape != (win_size, win_size):
                    continue

                med_loc = np.median(win)
                mad_loc = np.median(np.abs(win - med_loc))
                if mad_loc < 1e-6:
                    continue

                sat_frac = float((win >= sat_thr).mean())
                grad_win = grad_mag[y0:y1, x0:x1]
                grad_mean = float(np.mean(grad_win))

                score = mad_loc * (1.0 - sat_frac) / (1.0 + grad_mean)

                if score > best_score:
                    best_score = score
                    best_center = (cy, cx)

        if best_center is not None:
            centers.append(best_center)

    return centers


def fit_affine_from_points(fixed_points, aligned_points):
    """
    Fit affine transform that maps fixed_points (output coords)
    -> aligned_points (input coords).

    We solve:
      [y_in]   [a b c] [y_out]
      [x_in] = [d e f] [x_out]
                          [1]
    """
    assert len(fixed_points) == len(aligned_points)
    n = len(fixed_points)
    if n < 3:
        return np.eye(2, dtype=np.float32), np.zeros(2, dtype=np.float32)

    A = np.zeros((2 * n, 6), dtype=np.float64)
    b = np.zeros((2 * n,), dtype=np.float64)

    for i, (fpt, apt) in enumerate(zip(fixed_points, aligned_points)):
        y_out, x_out = fpt
        y_in, x_in = apt

        # y_in row
        A[2 * i, 0] = y_out
        A[2 * i, 1] = x_out
        A[2 * i, 2] = 1.0
        A[2 * i, 3] = 0.0
        A[2 * i, 4] = 0.0
        A[2 * i, 5] = 0.0
        b[2 * i] = y_in

        # x_in row
        A[2 * i + 1, 0] = 0.0
        A[2 * i + 1, 1] = 0.0
        A[2 * i + 1, 2] = 0.0
        A[2 * i + 1, 3] = y_out
        A[2 * i + 1, 4] = x_out
        A[2 * i + 1, 5] = 1.0
        b[2 * i + 1] = x_in

    params, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    a, b1, c, d, e, f1 = params

    matrix = np.array([[a, b1],
                       [d, e]], dtype=np.float32)
    offset = np.array([c, f1], dtype=np.float32)
    return matrix, offset


def post_correction_affine(fixed, aligned,
                           normalize=True,
                           win_size=64,
                           search_radius=32):
    """
    Post-correction step:
      - choose 4 good windows automatically on the fixed image
      - for each window, compute local dx,dy between fixed and aligned
      - fit a global affine mapping (output->input) from fixed coords to aligned coords
      - apply this affine to 'aligned' and return refined image
    """
    h, w = fixed.shape
    if aligned.shape != fixed.shape:
        return aligned, None

    if normalize:
        fixed_reg = robust_normalize_mad(fixed)
        aligned_reg = robust_normalize_mad(aligned)
    else:
        fixed_reg = robust_normalize(fixed)
        aligned_reg = robust_normalize(aligned)

    centers = pick_window_centers_auto(fixed_reg, win_size, search_radius)
    if len(centers) < 3:
        return aligned, {"used_windows": len(centers), "status": "insufficient_windows"}

    half = win_size // 2

    fixed_points = []
    aligned_points = []

    for (cy, cx) in centers:
        y0 = cy - half
        y1 = cy + half
        x0 = cx - half
        x1 = cx + half

        if y0 < 0 or x0 < 0 or y1 > h or x1 > w:
            continue

        f_win = fixed_reg[y0:y1, x0:x1]
        a_win = aligned_reg[y0:y1, x0:x1]
        if f_win.shape != (win_size, win_size) or a_win.shape != (win_size, win_size):
            continue

        dx, dy, peak = phase_correlate(f_win, a_win)

        # IMPORTANT: phase_correlate returns (dx,dy) ≈ -(true_shift)
        # If features in aligned are at (cy+u_y, cx+u_x) vs fixed at (cy,cx),
        # we get dx ≈ -u_x, dy ≈ -u_y.
        # So input (aligned) point is at (cy - dy, cx - dx).
        fixed_points.append((float(cy), float(cx)))
        aligned_points.append((float(cy - dy), float(cx - dx)))

    if len(fixed_points) < 3:
        return aligned, {"used_windows": len(fixed_points), "status": "insufficient_windows"}

    matrix, offset = fit_affine_from_points(fixed_points, aligned_points)

    refined = ndimage.affine_transform(
        aligned,
        matrix=matrix,
        offset=offset,
        order=1,           # linear here to avoid extra ringing around stars
        mode="constant",
        cval=0.0,
        prefilter=False,   # for order=1 no prefilter needed
    )

    info = {
        "used_windows": len(fixed_points),
        "status": "ok",
        "matrix": matrix,
        "offset": offset,
    }
    return refined, info


# ==========================================================
#                     WORKER FUNCTION
# ==========================================================

def _worker_job(job):
    """
    Single worker job:
      job = (ref_path, in_path, out_path, opts_dict)
    """
    ref_path, in_path, out_path, opts = job
    superfine = opts["superfine"]
    normalize = opts["normalize"]
    max_angle = opts["max_angle"]
    angle_step = opts["angle_step"]
    scale_delta = opts["scale_delta"]
    scale_step = opts["scale_step"]
    ds_max = opts["ds_max"]
    top_size = opts["top_size"]
    post_corr = opts["post_correction"]
    pc_win = opts["pcorr_win_size"]
    pc_radius = opts["pcorr_search_radius"]
    flux_mode = opts.get("flux_mode", False)

    # Main interpolation order:
    #   in flux mode we use linear kernel (order=1) to avoid ringing around
    #   hot pixels and preserve local flux.
    main_order = 1 if flux_mode else 3

    ref_abs = os.path.abspath(ref_path)
    in_abs = os.path.abspath(in_path)
    out_abs = os.path.abspath(out_path)

    # Reference file itself: copy or skip
    if in_abs == ref_abs:
        if in_abs != out_abs:
            os.makedirs(os.path.dirname(out_abs), exist_ok=True)
            shutil.copy2(in_abs, out_abs)
            status = "copied(ref)"
        else:
            status = "skipped(ref_inplace)"
        return in_path, out_path, status, None

    fixed, fixed_header = read_fits(ref_abs)
    moving, _ = read_fits(in_abs)

    if superfine:
        angle, scale, dx, dy, extra = register_superfine(
            fixed, moving,
            max_angle=max_angle,
            angle_step=angle_step,
            scale_delta=scale_delta,
            scale_step=scale_step,
            top_size=top_size,
            normalize=normalize,
        )
        peaks, peak_full = extra
    else:
        angle_range = (-max_angle, max_angle)
        scale_range = (1.0 - scale_delta, 1.0 + scale_delta)
        angle, scale, dx, dy, peaks = register_coarse_fine(
            fixed, moving,
            angle_range=angle_range,
            angle_step=angle_step,
            scale_range=scale_range,
            scale_step=scale_step,
            ds_max=ds_max,
            normalize=normalize,
        )
        peak_full = peaks[-1]

    aligned = apply_full_transform(moving, angle, scale, dx, dy, order=main_order)

    pc_info = None
    if post_corr:
        aligned, pc_info = post_correction_affine(
            fixed, aligned,
            normalize=normalize,
            win_size=pc_win,
            search_radius=pc_radius,
        )

    os.makedirs(os.path.dirname(out_abs), exist_ok=True)
    write_fits(out_abs, aligned, fixed_header)

    status = "aligned"
    info = {
        "angle": angle,
        "scale": scale,
        "dx": dx,
        "dy": dy,
        "peak_full": float(peak_full),
        "post_correction": pc_info,
    }
    return in_path, out_path, status, info


# ==========================================================
#                        PROGRESS BAR
# ==========================================================

def format_progress(current, total, width=40):
    """Return a string like: [####------] 3/10."""
    if total <= 0:
        return "[{}] 0/0".format("-" * width)
    frac = current / float(total)
    frac = max(0.0, min(1.0, frac))
    filled = int(width * frac + 0.5)
    filled = min(filled, width)
    bar = "#" * filled + "-" * (width - filled)
    return f"[{bar}] {current}/{total}"


# ==========================================================
#                          CLI
# ==========================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "FFT-based registration of FITS frames to a reference frame.\n\n"
            "Features:\n"
            "  • rotation/scale/shift search via phase correlation\n"
            "  • SUPERFINE pyramidal refinement mode\n"
            "  • local affine post-correction using 4 windows\n"
            "  • batch processing via batch_utils (masks, lists, numbered templates)\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "reference",
        help="Reference FITS file (.fit/.fits) to which all images are aligned.",
    )
    parser.add_argument(
        "input_spec",
        help=(
            "Input frame list (batch_utils syntax):\n"
            "  single filename, wildcard mask, @list, numbered template, etc."
        ),
    )
    parser.add_argument(
        "output_spec",
        help=(
            "Output list (batch_utils syntax):\n"
            "  single filename or numbered template for the whole set."
        ),
    )

    parser.add_argument(
        "--superfine",
        action="store_true",
        help=(
            "Enable SUPERFINE pyramidal mode:\n"
            "  coarse search at small size, then progressive refinement\n"
            "  on larger sizes up to the full-resolution image."
        ),
    )

    parser.add_argument(
        "--post-correction",
        action="store_true",
        help=(
            "Enable post-correction (affine refinement):\n"
            "  using 4 local windows, estimate residual affine offset\n"
            "  and apply one more global affine transform to the result."
        ),
    )

    parser.add_argument(
        "--max-angle",
        type=float,
        default=5.0,
        help="Maximum absolute rotation angle (degrees) for the search. Default: 5.",
    )
    parser.add_argument(
        "--angle-step",
        type=float,
        default=0.25,
        help="Angle step (coarse level), in degrees. Default: 0.25.",
    )
    parser.add_argument(
        "--scale-delta",
        type=float,
        default=0.01,
        help=(
            "Allowed relative scale deviation around 1.0.\n"
            "Default: 0.01 (±1%%)."
        ),
    )
    parser.add_argument(
        "--scale-step",
        type=float,
        default=0.002,
        help="Scale step on coarse search. Default: 0.002.",
    )
    parser.add_argument(
        "--ds-max",
        type=int,
        default=256,
        help=(
            "Maximum downsample size (longest side) in coarse/fine mode\n"
            "  (when SUPERFINE is disabled). Default: 256."
        ),
    )
    parser.add_argument(
        "--top-size",
        type=int,
        default=128,
        help=(
            "Top pyramid size (longest side) in SUPERFINE mode.\n"
            "Default: 128."
        ),
    )
    parser.add_argument(
        "--no-normalize",
        action="store_true",
        help=(
            "Disable robust normalization (median + MAD) for similarity search.\n"
            "By default normalization is ENABLED."
        ),
    )

    parser.add_argument(
        "--pcorr-win-size",
        type=int,
        default=64,
        help=(
            "Window size (N×N) for post-correction local registration.\n"
            "Default: 64."
        ),
    )
    parser.add_argument(
        "--pcorr-search-radius",
        type=int,
        default=32,
        help=(
            "Search radius for window centers around nominal 25/75%% positions in X/Y.\n"
            "A 5×5 grid is used inside this radius. Default: 32."
        ),
    )

    parser.add_argument(
        "--flux",
        action="store_true",
        help=(
            "Enable FLUX mode: the main interpolation (rotation+scale+shift)\n"
            "  uses linear kernel (order=1), without negative lobes,\n"
            "  preserving local flux. Recommended for real astro frames,\n"
            "  especially in presence of cosmic rays and hot pixels."
        ),
    )

    # If no arguments supplied – show help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return

    args = parser.parse_args()

    ref_path = args.reference
    if not os.path.isfile(ref_path):
        raise FileNotFoundError(f"Reference file not found: {ref_path}")

    pairs = batch_utils.build_io_file_lists(args.input_spec, args.output_spec)

    normalize = not args.no_normalize
    flux_mode = args.flux

    opts = {
        "superfine": args.superfine,
        "normalize": normalize,
        "max_angle": args.max_angle,
        "angle_step": args.angle_step,
        "scale_delta": args.scale_delta,
        "scale_step": args.scale_step,
        "ds_max": args.ds_max,
        "top_size": args.top_size,
        "post_correction": args.post_correction,
        "pcorr_win_size": args.pcorr_win_size,
        "pcorr_search_radius": args.pcorr_search_radius,
        "flux_mode": flux_mode,
    }

    jobs = [(ref_path, inp, out, opts) for (inp, out) in pairs]

    cpu_count = mp.cpu_count()
    workers = max(1, cpu_count - 1)

    total = len(jobs)
    print(f"Reference : {ref_path}")
    print(f"Total files to process: {total}")
    print(f"Using {workers} worker process(es)")
    if args.superfine:
        print("Mode      : SUPERFINE pyramid")
    else:
        print("Mode      : coarse+fine")
    if args.post_correction:
        print("Post-corr : enabled (affine refinement)")
    else:
        print("Post-corr : disabled")
    if flux_mode:
        print("Interp    : FLUX mode (linear, order=1)")
    else:
        print("Interp    : default (cubic main, linear post-corr)")

    if total == 0:
        return

    bar_width = 40
    done = 0

    with mp.Pool(processes=workers) as pool:
        for in_path, out_path, status, info in pool.imap_unordered(_worker_job, jobs):
            done += 1
            prog = format_progress(done, total, bar_width)

            # Print progress bar on one line, then details on the next line
            if status.startswith("aligned"):
                if info is not None:
                    pc = info.get("post_correction", None)
                    pc_str = ""
                    if pc is not None and isinstance(pc, dict):
                        pc_status = pc.get("status", "n/a")
                        used_w = pc.get("used_windows", 0)
                        pc_str = f", pcorr={pc_status}, windows={used_w}"
                    print(prog)
                    print(
                        f"  aligned: {in_path} -> {out_path} | "
                        f"angle={info['angle']:.6f} deg, "
                        f"scale={info['scale']:.8f}, "
                        f"dx={info['dx']:.3f}, dy={info['dy']:.3f}, "
                        f"peak={info['peak_full']:.3g}{pc_str}"
                    )
                else:
                    print(prog)
                    print(f"  aligned: {in_path} -> {out_path}")
            elif status.startswith("copied"):
                print(prog)
                print(f"  copied ref: {in_path} -> {out_path}")
            elif status.startswith("skipped"):
                print(prog)
                print(f"  skipped (ref in-place): {in_path}")
            else:
                print(prog)
                print(f"  {status}: {in_path} -> {out_path}")


if __name__ == "__main__":
    main()
