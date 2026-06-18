#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mtf - Midtone Transfer Function (as in PixInsight PixelMath).

Applies a nonlinear tone curve that passes through (0,0), (m, 0.5), (1,1)
where m is the midtones balance parameter.

Formula:  mtf(x, m) = (1-m)*x / (m + x*(1-2*m))

  m < 0.5  -> brightens midtones
  m = 0.5  -> identity (linear)
  m > 0.5  -> darkens midtones

Pixel values are normalized to [0,1] using the dtype's full range
(e.g., uint16: 0..65535) before applying the function, then scaled back.
This matches PixInsight's internal [0,1] normalization.

Supports 2D (mono) and 3D (color, 3xH×W) FITS images.

Auto black/white level detection:
  --autoblack   Compute black level from darkest background cells (median - 5*MAD)
  --autowhite   Compute white level from brightest 0.01% of pixels
  --auto        Both --autoblack and --autowhite

Color preservation (3D only):
  --keepcolor [K]  Apply MTF to luminance (R+2G+B)/4, restore color ratios.
                   K=1: full color preservation. K=0: none (blend with independent).
  --equal          Use equal-weight luminance (R+G+B)/3 instead of green-weighted.
"""

import sys
import os
import warnings
import numpy as np
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils

# Suppress warnings from nanmedian/nanstd on all-NaN slices
warnings.filterwarnings('ignore', message='.*All-NaN slice.*', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*Degrees of freedom <= 0.*', category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CLIP_MARGIN_INT = 100     # ADU offset from 0 and MAX for integer types
CLIP_MARGIN_FLOAT = 0.01  # offset from 0.0 and 1.0 for float types


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "mtf - Midtone Transfer Function (PixInsight-compatible).\n"
        "\n"
        "Usage:\n"
        "  mtf.py input_spec output_spec K [options]\n"
        "  mtf.py input_spec output_spec --median T [options]\n"
        "\n"
        "Positional arguments:\n"
        "  input_spec       Input file(s): single file, wildcard (*.fit),\n"
        "                   numbered sequence (img0001.fit), or list (@list.txt)\n"
        "  output_spec      Output: single file, numbered pattern (out0001.fit),\n"
        "                   or directory\n"
        "  K [K2 K3 ...]    One or more midtones balance values (0..1 exclusive).\n"
        "                     K < 0.5 -> brighten midtones\n"
        "                     K = 0.5 -> no change (identity)\n"
        "                     K > 0.5 -> darken midtones\n"
        "                   Multiple K: result = average of mtf(x,Ki) for each Ki.\n"
        "                   Duplicate values act as weights (counted multiple times).\n"
        "                   Example: 0.1 0.1 0.4 -> 2/3 weight on K=0.1, 1/3 on K=0.4\n"
        "\n"
        "Options:\n"
        "  --median T, -m T  Auto-pick K so the image median maps to T, with\n"
        "                    T in (0,1), instead of giving K directly. For color\n"
        "                    the median is taken from luminance. Mutually\n"
        "                    exclusive with positional K; pairs well with --auto.\n"
        "  --zero, -z        Zero handling for aligned/invalid borders: ignore\n"
        "                    zeros when computing the --median, and restore the\n"
        "                    input's zero pixels in the output so MTF never\n"
        "                    shifts a 0 to a non-zero value.\n"
        "  --clip [P]        After MTF, clip black/white by percentile P\n"
        "                    (default P=0.001). Ignores pixels <= 0.\n"
        "                    Remaps [black_P, white_P] to [0, 1].\n"
        "\n"
        "Auto black/white (linked for color: shared boundaries):\n"
        "  --autoblack       Auto black level from darkest background cells.\n"
        "  --autowhite       Auto white level from brightest 0.01% pixels.\n"
        "  --auto            Both --autoblack and --autowhite.\n"
        "\n"
        "Color options (3-channel images only):\n"
        "  --keepcolor [K]   Preserve color ratios via luminance scaling.\n"
        "                    K=1 (default): full color preservation.\n"
        "                    K=0: no preservation (same as per-channel).\n"
        "                    0<K<1: blend between per-channel and preserved.\n"
        "                    Luminance = (R + 2*G + B) / 4 (green-weighted).\n"
        "  --keepcolor-hsl [K]  Preserve color via HSL: apply MTF to L only,\n"
        "                    H and S unchanged. Better saturation control\n"
        "                    for aggressive stretch (K < 0.15).\n"
        "  --equal           Use equal-weight luminance (R+G+B)/3.\n"
        "                    Works with both --keepcolor and --keepcolor-hsl.\n"
        "\n"
        "Examples:\n"
        "  mtf.py img.fit out.fit 0.2\n"
        "  mtf.py img.fit out.fit --median 0.25 --auto\n"
        "  mtf.py img.fit out.fit -m 0.2 --keepcolor --auto\n"
        "  mtf.py img.fit out.fit --median 0.25 --auto --zero\n"
        "  mtf.py img.fit out.fit 0.1 0.2 0.4\n"
        "  mtf.py img.fit out.fit 0.1 0.1 0.4         (weight K=0.1 x2)\n"
        "  mtf.py img.fit out.fit 0.15 --auto\n"
        "  mtf.py img.fit out.fit 0.25 --keepcolor --auto\n"
        "  mtf.py img.fit out.fit 0.01 --keepcolor-hsl --auto\n"
        "  mtf.py img.fit out.fit 0.25 --keepcolor 0.7 --equal\n"
    )
    sys.exit(1)


def _parse_k(s, name):
    """Parse and validate a midtones balance value."""
    try:
        v = float(s)
    except ValueError:
        sys.stderr.write(f"Error: {name} must be a number.\n")
        sys.exit(1)
    if v <= 0.0 or v >= 1.0:
        sys.stderr.write(f"Error: {name} must be in range (0, 1) exclusive.\n")
        sys.exit(1)
    return v


def parse_args(argv):
    args = argv[1:]

    clip = None
    keepcolor_k = None      # None = disabled, float = ratio mode
    keepcolor_hsl_k = None  # None = disabled, float = HSL mode
    equal_lum = False
    autoblack = False
    autowhite = False
    median_target = None    # None = disabled; float in (0,1) = target-median mode
    zero_mode = False       # --zero: ignore zeros in stats, restore them on output

    i = 0
    positional = []
    while i < len(args):
        if args[i] in ('--median', '-m'):
            i += 1
            if i >= len(args):
                sys.stderr.write("Error: --median requires a value T in (0, 1).\n")
                sys.exit(1)
            try:
                median_target = float(args[i])
            except ValueError:
                sys.stderr.write("Error: --median value must be a number in (0, 1).\n")
                sys.exit(1)
            if median_target <= 0.0 or median_target >= 1.0:
                sys.stderr.write("Error: --median T must be in range (0, 1) exclusive.\n")
                sys.exit(1)
            i += 1
        elif args[i] == '--keepcolor-hsl':
            keepcolor_hsl_k = 1.0  # default
            i += 1
            if i < len(args) and not args[i].startswith('--'):
                try:
                    keepcolor_hsl_k = float(args[i])
                    i += 1
                except ValueError:
                    pass
        elif args[i] == '--keepcolor':
            keepcolor_k = 1.0  # default
            i += 1
            # Optional K value
            if i < len(args) and not args[i].startswith('--'):
                try:
                    keepcolor_k = float(args[i])
                    i += 1
                except ValueError:
                    pass
        elif args[i] == '--equal':
            equal_lum = True
            i += 1
        elif args[i] in ('--zero', '-z'):
            zero_mode = True
            i += 1
        elif args[i] == '--auto':
            autoblack = True
            autowhite = True
            i += 1
        elif args[i] == '--autoblack':
            autoblack = True
            i += 1
        elif args[i] == '--autowhite':
            autowhite = True
            i += 1
        elif args[i] == '--preserve_color':
            # Backward compat: treat as --keepcolor 1
            keepcolor_k = 1.0
            i += 1
        elif args[i] == '--clip':
            clip = 0.001
            i += 1
            if i < len(args) and not args[i].startswith('--'):
                try:
                    clip = float(args[i])
                    i += 1
                except ValueError:
                    pass
        else:
            positional.append(args[i])
            i += 1

    if median_target is not None:
        # Target-median mode: K is computed per file, not given.
        if len(positional) < 2:
            usage()
        if len(positional) > 2:
            sys.stderr.write(
                "Error: --median cannot be combined with explicit K value(s).\n")
            sys.exit(1)
        input_spec = positional[0]
        output_spec = positional[1]
        ks = []
    else:
        if len(positional) < 3:
            usage()
        input_spec = positional[0]
        output_spec = positional[1]
        # All remaining positional args are K values
        ks = []
        for s in positional[2:]:
            ks.append(_parse_k(s, "K"))
        if not ks:
            usage()

    return {
        'input_spec': input_spec,
        'output_spec': output_spec,
        'ks': ks,
        'median_target': median_target,
        'zero': zero_mode,
        'clip': clip,
        'keepcolor_k': keepcolor_k,
        'keepcolor_hsl_k': keepcolor_hsl_k,
        'equal_lum': equal_lum,
        'autoblack': autoblack,
        'autowhite': autowhite,
    }


# ---------------------------------------------------------------------------
# MTF core
# ---------------------------------------------------------------------------

def apply_mtf(x, m):
    """Apply midtone transfer function to normalized [0,1] array.

    mtf(x, m) = (1-m)*x / (m + x*(1-2*m))
    """
    if abs(m - 0.5) < 1e-10:
        return x.copy()

    num = (1.0 - m) * x
    den = m + x * (1.0 - 2.0 * m)

    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(den != 0.0, num / den, np.where(x <= 0.0, 0.0, 1.0))

    return np.clip(result, 0.0, 1.0)


def apply_multi_mtf(x, ks):
    """Apply multiple MTFs and average: result = mean([mtf(x, ki) for ki in ks]).

    Duplicate K values act as weights (counted multiple times in the average).
    Single K: equivalent to apply_mtf(x, K).
    """
    if len(ks) == 1:
        return apply_mtf(x, ks[0])

    acc = apply_mtf(x, ks[0])
    for k in ks[1:]:
        acc = acc + apply_mtf(x, k)
    acc /= len(ks)
    return np.clip(acc, 0.0, 1.0)


def _find_mtf_k(median_norm, target):
    """Find the MTF parameter K that maps median_norm to target.

    Analytical inverse of the MTF (copied from fits2tiff.py to keep the tool
    self-contained, per the project's tool-independence rule):
        target = (1-K)*med / (K + med*(1-2K))
      => K = med*(1 - target) / (target - 2*target*med + med)
    """
    if median_norm <= 0 or median_norm >= 1:
        return 0.5
    denom = target - 2.0 * target * median_norm + median_norm
    if abs(denom) < 1e-15:
        return 0.5
    k = median_norm * (1.0 - target) / denom
    return float(np.clip(k, 0.001, 0.999))


def representative_median(work, equal_lum, ignore_zeros=False):
    """Median of the image in [0,1] space.

    For 3-channel color the median is taken from luminance, matching the
    --median target and the keepcolor luminance choice (equal_lum).
    ignore_zeros (--zero): drop zero pixels (aligned-frame borders / no-data)
    so they do not drag the median down.
    """
    if work.ndim == 3 and work.shape[0] == 3:
        if equal_lum:
            lum = np.mean(work, axis=0)
        else:
            lum = (work[0] + 2.0 * work[1] + work[2]) / 4.0
    else:
        lum = work
    valid = lum[lum > 0] if ignore_zeros else lum.ravel()
    if valid.size == 0:
        return 0.5
    return float(np.median(valid))


# ---------------------------------------------------------------------------
# HSL color space conversion (vectorized numpy)
# ---------------------------------------------------------------------------

def _rgb_to_hsl(rgb):
    """Convert RGB (3, H, W) in [0,1] to HSL (3, H, W) in [0,1]."""
    r, g, b = rgb[0], rgb[1], rgb[2]

    cmax = np.maximum(np.maximum(r, g), b)
    cmin = np.minimum(np.minimum(r, g), b)
    delta = cmax - cmin

    l = (cmax + cmin) / 2.0

    s = np.zeros_like(l)
    mask = delta > 0
    low = mask & (l <= 0.5)
    high = mask & (l > 0.5)
    denom_low = cmax[low] + cmin[low]
    denom_high = 2.0 - cmax[high] - cmin[high]
    s[low] = np.where(denom_low > 0, delta[low] / denom_low, 0.0)
    s[high] = np.where(denom_high > 0, delta[high] / denom_high, 0.0)

    h = np.zeros_like(l)
    with np.errstate(divide='ignore', invalid='ignore'):
        inv_delta = np.where(delta > 0, 1.0 / delta, 0.0)

    rm = mask & (cmax == r)
    h[rm] = ((g[rm] - b[rm]) * inv_delta[rm]) % 6.0

    gm = mask & (cmax == g) & ~rm
    h[gm] = (b[gm] - r[gm]) * inv_delta[gm] + 2.0

    bm = mask & ~rm & ~gm
    h[bm] = (r[bm] - g[bm]) * inv_delta[bm] + 4.0

    h /= 6.0
    h = h % 1.0

    return np.stack([h, s, l], axis=0)


def _hsl_to_rgb(hsl):
    """Convert HSL (3, H, W) in [0,1] to RGB (3, H, W) in [0,1]."""
    h, s, l = hsl[0], hsl[1], hsl[2]

    c = (1.0 - np.abs(2.0 * l - 1.0)) * s
    hp = h * 6.0
    x = c * (1.0 - np.abs(hp % 2.0 - 1.0))
    m = l - c / 2.0

    r = np.zeros_like(h)
    g = np.zeros_like(h)
    b = np.zeros_like(h)

    for lo, hi, rv, gv, bv in [
        (0, 1, 'c', 'x', '0'), (1, 2, 'x', 'c', '0'),
        (2, 3, '0', 'c', 'x'), (3, 4, '0', 'x', 'c'),
        (4, 5, 'x', '0', 'c'), (5, 6, 'c', '0', 'x'),
    ]:
        mask = (hp >= lo) & (hp < hi)
        for ch, val in [(r, rv), (g, gv), (b, bv)]:
            if val == 'c':
                ch[mask] = c[mask]
            elif val == 'x':
                ch[mask] = x[mask]

    return np.clip(np.stack([r + m, g + m, b + m], axis=0), 0.0, 1.0)


# ---------------------------------------------------------------------------
# Auto black/white detection
# ---------------------------------------------------------------------------

def _compute_autoblack(plane_01):
    """Compute auto black level for a single 2D plane in [0,1] space.

    Uses darkest 3% of background cells (min 1 cell).
    Black = median - 5.0 * MAD of pixels in those cells.
    """
    from background import _build_cell_grid

    cell_size = 64
    h, w = plane_01.shape

    # Build cell grid (sigma-clipped medians per cell)
    grid = _build_cell_grid(plane_01, cell_size, clip_k=1.7)
    ny, nx = grid.shape

    # Find darkest 3% of non-zero cells
    valid_mask = grid > 0
    valid_values = grid[valid_mask]
    if valid_values.size == 0:
        return 0.0

    n_dark = max(1, int(valid_values.size * 0.03))
    threshold = np.partition(valid_values, n_dark)[n_dark - 1]

    # Get cell indices where grid value <= threshold and > 0
    dark_cells = np.argwhere((grid <= threshold) & (grid > 0))

    # Collect original pixels from those cells
    pixels = []
    for cy, cx in dark_cells:
        y0 = cy * cell_size
        y1 = min(y0 + cell_size, h)
        x0 = cx * cell_size
        x1 = min(x0 + cell_size, w)
        cell = plane_01[y0:y1, x0:x1].ravel()
        cell = cell[cell > 0]  # exclude zeros
        if cell.size > 0:
            pixels.append(cell)

    if not pixels:
        return 0.0

    all_pixels = np.concatenate(pixels)
    med = float(np.median(all_pixels))
    mad = float(np.median(np.abs(all_pixels - med)))

    black = med - 5.0 * mad
    return max(0.0, black)


def _compute_autowhite(plane_01):
    """Compute auto white level for a single 2D plane in [0,1] space.

    White = median of top 0.01% of non-zero pixels.
    """
    valid = plane_01[plane_01 > 0].ravel()
    if valid.size == 0:
        return 1.0

    n_bright = max(1, int(valid.size * 0.0001))
    idx = valid.size - n_bright
    partitioned = np.partition(valid, idx)
    white = float(np.median(partitioned[idx:]))

    return min(1.0, max(white, 0.01))


def compute_auto_levels(work, autoblack, autowhite):
    """Compute auto black/white levels. Returns (black, white).

    For 3D color images: linked mode.
      black = min(per-channel blacks)
      white = max(per-channel whites)
    For 2D mono: computed directly.
    """
    black = 0.0
    white = 1.0

    if work.ndim == 2:
        if autoblack:
            black = _compute_autoblack(work)
        if autowhite:
            white = _compute_autowhite(work)
    elif work.ndim == 3 and work.shape[0] == 3:
        if autoblack:
            blacks = [_compute_autoblack(work[ch]) for ch in range(3)]
            black = min(blacks)
        if autowhite:
            whites = [_compute_autowhite(work[ch]) for ch in range(3)]
            white = max(whites)

    return black, white


def apply_black_white(work, black, white):
    """Remap data from [black, white] to [0, 1]. Preserves zeros."""
    if black == 0.0 and white >= 1.0:
        return work  # no-op

    bw_range = white - black
    if bw_range <= 0:
        return work

    zero_mask = work <= 0
    result = (work - black) / bw_range
    result = np.clip(result, 0.0, 1.0)
    result[zero_mask] = 0.0
    return result


# ---------------------------------------------------------------------------
# Color preservation
# ---------------------------------------------------------------------------

def apply_mtf_keepcolor(work_3d, ks, keepcolor_k, equal_lum):
    """Apply MTF with color preservation.

    1. Compute luminance L_old
    2. Apply MTF to luminance: L_new = multi_MTF(L_old)
    3. Restore colors: ch_new = (ch_old / L_old) * L_new
    4. Blend with per-channel MTF result by keepcolor_k
    """
    # Compute luminance
    if equal_lum:
        lum_old = np.mean(work_3d, axis=0)  # (R+G+B)/3
    else:
        lum_old = (work_3d[0] + 2.0 * work_3d[1] + work_3d[2]) / 4.0  # (R+2G+B)/4

    # Apply MTF to luminance
    lum_new = apply_multi_mtf(lum_old, ks)

    # Scale channels by luminance ratio
    with np.errstate(divide='ignore', invalid='ignore'):
        scale = np.where(lum_old > 0.0, lum_new / lum_old, 0.0)

    preserved = work_3d * scale[np.newaxis, :, :]
    preserved = np.clip(preserved, 0.0, 1.0)

    # If keepcolor_k == 1, return fully preserved
    if keepcolor_k >= 1.0:
        return preserved

    # If keepcolor_k == 0, return per-channel MTF
    if keepcolor_k <= 0.0:
        return apply_multi_mtf(work_3d, ks)

    # Blend: result = per_channel * (1-K) + preserved * K
    per_channel = apply_multi_mtf(work_3d, ks)
    return per_channel * (1.0 - keepcolor_k) + preserved * keepcolor_k


def apply_mtf_keepcolor_hsl(work_3d, ks, keepcolor_k):
    """Apply MTF with color preservation via HSL color space.

    1. Convert RGB to HSL
    2. Apply multi-MTF to L channel only (H and S unchanged)
    3. Convert back to RGB
    4. Blend with per-channel MTF result by keepcolor_k

    Better saturation preservation than ratio method for aggressive stretch.
    """
    hsl = _rgb_to_hsl(work_3d)

    # Apply MTF only to L channel
    hsl[2] = apply_multi_mtf(hsl[2], ks)

    preserved = _hsl_to_rgb(hsl)

    if keepcolor_k >= 1.0:
        return preserved

    if keepcolor_k <= 0.0:
        return apply_multi_mtf(work_3d, ks)

    per_channel = apply_multi_mtf(work_3d, ks)
    return per_channel * (1.0 - keepcolor_k) + preserved * keepcolor_k


# ---------------------------------------------------------------------------
# Normalization helpers
# ---------------------------------------------------------------------------

_norm_scale = {'min': 0.0, 'max': 1.0}  # set by normalize_to_01, used by denormalize


def normalize_to_01(data, orig_dtype):
    """Normalize data to [0, 1] range.

    Integer types: uses dtype range (0..65535 for uint16, etc).
    Float types: uses actual data range (0..max_value).
    """
    global _norm_scale
    work = data.astype(np.float64)
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        _norm_scale = {'min': float(info.min), 'max': float(info.max)}
        work = (work - info.min) / float(info.max - info.min)
    else:
        # Float: normalize by actual data range
        valid = work[work > 0] if np.any(work > 0) else work.ravel()
        vmax = float(np.max(valid)) if valid.size > 0 else 1.0
        if vmax <= 0:
            vmax = 1.0
        _norm_scale = {'min': 0.0, 'max': vmax}
        work = work / vmax
        np.clip(work, 0.0, 1.0, out=work)
    return work


def denormalize_from_01(data, orig_dtype):
    """Convert [0, 1] data back to original range."""
    vmin = _norm_scale['min']
    vmax = _norm_scale['max']
    work = data * (vmax - vmin) + vmin
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        work = np.clip(work, info.min, info.max)
        work = np.rint(work)
        return work.astype(orig_dtype)
    if np.issubdtype(orig_dtype, np.floating):
        arr = work.astype(orig_dtype)
        bad = ~np.isfinite(arr)
        if np.any(bad):
            arr[bad] = 0
        return arr
    return work.astype(np.float32)


# ---------------------------------------------------------------------------
# Clip (legacy)
# ---------------------------------------------------------------------------

def clip_black_white(work, clip_pct, orig_dtype):
    """Clip black/white levels by percentile and remap with margin.

    Ignores pixels <= 0 (zero-padded regions from crop).
    Works on normalized [0, 1] data. Applied per-channel for 3D.

    Remaps [black_pct, white_pct] -> [margin_lo, 1 - margin_hi]
    """
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        full_range = float(info.max - info.min)
        margin = CLIP_MARGIN_INT / full_range
    else:
        margin = CLIP_MARGIN_FLOAT

    dst_lo = margin
    dst_hi = 1.0 - margin

    def _clip_plane(plane):
        valid = plane[plane > 0]
        if len(valid) == 0:
            return plane
        n = len(valid)
        n_lo = max(1, int(n * clip_pct / 100.0))
        n_hi = max(1, int(n * clip_pct / 100.0))

        partitioned = np.partition(valid, n_lo)
        black = float(np.median(partitioned[:n_lo]))

        idx = n - n_hi
        partitioned = np.partition(valid, idx)
        white = float(np.median(partitioned[idx:]))

        if white <= black:
            return plane

        zero_mask = plane <= 0
        result = dst_lo + (plane - black) / (white - black) * (dst_hi - dst_lo)
        result = np.clip(result, 0.0, 1.0)
        result[zero_mask] = 0.0
        return result

    if work.ndim == 2:
        return _clip_plane(work)

    avg = np.mean(work, axis=0)
    valid = avg[avg > 0]
    if len(valid) == 0:
        return work
    n = len(valid)
    n_lo = max(1, int(n * clip_pct / 100.0))
    n_hi = max(1, int(n * clip_pct / 100.0))

    partitioned = np.partition(valid, n_lo)
    black = float(np.median(partitioned[:n_lo]))

    idx = n - n_hi
    partitioned = np.partition(valid, idx)
    white = float(np.median(partitioned[idx:]))

    if white <= black:
        return work

    zero_mask = work <= 0
    result = dst_lo + (work - black) / (white - black) * (dst_hi - dst_lo)
    result = np.clip(result, 0.0, 1.0)
    result[zero_mask] = 0.0
    return result


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------

def process_file(infile, outfile, config):
    """Process a single FITS file."""
    ks = config['ks']
    median_target = config['median_target']
    restore_zeros = config['zero']
    clip = config['clip']
    keepcolor_k = config['keepcolor_k']
    keepcolor_hsl_k = config['keepcolor_hsl_k']
    equal_lum = config['equal_lum']
    do_autoblack = config['autoblack']
    do_autowhite = config['autowhite']

    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")

        data = hdul[0].data
        header = hdul[0].header.copy()
        orig_dtype = data.dtype

    # --zero: remember where the input is exactly 0 (aligned borders / no-data)
    zero_mask = (data == 0) if restore_zeros else None

    # Normalize to [0, 1]
    work = normalize_to_01(data, orig_dtype)

    # Auto black/white detection
    black, white = 0.0, 1.0
    if do_autoblack or do_autowhite:
        black, white = compute_auto_levels(work, do_autoblack, do_autowhite)

    # Apply black/white remapping
    work = apply_black_white(work, black, white)

    # Target-median mode: derive K so the (luminance) median maps to T
    if median_target is not None:
        med = representative_median(work, equal_lum, ignore_zeros=restore_zeros)
        ks = [_find_mtf_k(med, median_target)]

    # Apply MTF
    is_color = (work.ndim == 3 and work.shape[0] == 3)
    if is_color and keepcolor_hsl_k is not None and keepcolor_hsl_k > 0:
        work = apply_mtf_keepcolor_hsl(work, ks, keepcolor_hsl_k)
    elif is_color and keepcolor_k is not None and keepcolor_k > 0:
        work = apply_mtf_keepcolor(work, ks, keepcolor_k, equal_lum)
    else:
        work = apply_multi_mtf(work, ks)

    # Optional percentile clip (legacy)
    if clip is not None:
        work = clip_black_white(work, clip, orig_dtype)

    # Convert back
    out_data = denormalize_from_01(work, orig_dtype)

    # --zero: restore the input's zeros so MTF/normalization never shifts a 0
    if zero_mask is not None:
        out_data[zero_mask] = 0

    # Update header
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    # Build HISTORY
    ks_str = " ".join(f"{k:g}" for k in ks)
    if median_target is not None:
        parts = [f"median->{median_target:g} (K={ks_str})"]
    else:
        parts = [f"K=[{ks_str}]"]
    if do_autoblack or do_autowhite:
        parts.append(f"black={black:.6f}, white={white:.6f}")
    if is_color and keepcolor_hsl_k is not None:
        parts.append(f"keepcolor-hsl={keepcolor_hsl_k}")
    elif is_color and keepcolor_k is not None:
        lum_str = "equal" if equal_lum else "green"
        parts.append(f"keepcolor={keepcolor_k}, lum={lum_str}")
    if clip is not None:
        parts.append(f"clip={clip}%")

    header["HISTORY"] = f"MTF applied by mtf.py: {', '.join(parts)}"

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fits.PrimaryHDU(out_data, header=header).writeto(outfile, overwrite=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    config = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(
            config['input_spec'], config['output_spec'])
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    if config['median_target'] is not None:
        parts = [f"median->{config['median_target']:g}"]
    else:
        ks_str = " ".join(f"{k:g}" for k in config['ks'])
        parts = [f"K=[{ks_str}]" if len(config['ks']) > 1 else f"K={config['ks'][0]:g}"]
    if config['autoblack']:
        parts.append("autoblack")
    if config['autowhite']:
        parts.append("autowhite")
    if config['keepcolor_hsl_k'] is not None:
        parts.append(f"keepcolor-hsl={config['keepcolor_hsl_k']}")
    elif config['keepcolor_k'] is not None:
        lum = "equal" if config['equal_lum'] else "green"
        parts.append(f"keepcolor={config['keepcolor_k']}, lum={lum}")
    if config['clip'] is not None:
        parts.append(f"clip={config['clip']}%")

    ncpu = os.cpu_count() or 1
    workers = min(max(1, ncpu - 1), total)
    print(f"Applying MTF ({', '.join(parts)}) to {total} file(s), threads={workers}")

    done = 0
    errors = 0

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(process_file, infile, outfile, config):
                os.path.basename(infile)
            for infile, outfile in io_pairs
        }

        for future in as_completed(futures):
            fname = futures[future]
            done += 1
            try:
                future.result()
                sys.stdout.write(f"\r{done}/{total}")
                sys.stdout.flush()
            except Exception as e:
                errors += 1
                sys.stderr.write(f"\n  Error '{fname}': {e}\n")

    print(f"\nDone. {done - errors} OK, {errors} errors.")


if __name__ == "__main__":
    main()
