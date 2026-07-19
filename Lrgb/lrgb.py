#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
lrgb - LRGB composition for astronomical images.

Combines a high-SNR luminance channel with RGB color data to produce
a detailed color image.

Two combination methods:
  hsl   - Convert RGB to HSL, replace L, convert back (classic approach)
  ratio - Scale RGB by L_new/L_old ratio (preserves color ratios precisely)

Two operation modes:
  Without --auto: bare LRGB combine only
  With --auto:    full pipeline (RGB balance + background flatten + combine)

Supports two input variants:
  4 mono files: L R G B -> output
  L + RGB file: L RGB  -> output
"""

import sys
import os
import numpy as np
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils
import background

# Import rgbbalance functions for --auto mode (lazy, on demand)
_rgbbalance = None
def _get_rgbbalance():
    global _rgbbalance
    if _rgbbalance is None:
        sys.path.append(os.path.abspath(os.path.join(
            os.path.dirname(__file__), "../Rgbbalance")))
        import rgbbalance as _rb
        _rgbbalance = _rb
    return _rgbbalance


# =========================================================================
# CLI
# =========================================================================

def usage():
    sys.stderr.write(
        "lrgb - LRGB composition for astronomical images.\n"
        "\n"
        "Usage:\n"
        "  lrgb.py L.fit R.fit G.fit B.fit output.fit [options]\n"
        "  lrgb.py L.fit RGB.fit output.fit [options]\n"
        "\n"
        "Input modes:\n"
        "  5 args: L R G B files (mono 2D FITS) -> output RGB (3,H,W)\n"
        "  3 args: L file + RGB file (3,H,W)    -> output RGB (3,H,W)\n"
        "\n"
        "Options:\n"
        "  --auto          Full auto pipeline before LRGB combine:\n"
        "                    1. RGB balance by star photometry\n"
        "                    2. Background flattening per channel\n"
        "                    3. LRGB combination\n"
        "  --slum          Build a super-luminance master before combining:\n"
        "                  inverse-variance (SNR^2) blend of L with the synthetic\n"
        "                  luminance R+G+B (scale-matched). Opt-in: without it L\n"
        "                  is used as-is (the input L may already be a super-lum).\n"
        "  --dry-run       Measure and print the L:syn noise ratio and expected\n"
        "                  SNR gain, then exit without writing output.\n"
        "  --method M      Combination method: ratio (default) or hsl.\n"
        "                    ratio: C_out = C * L_fit / (R+G+B + eps)\n"
        "                    hsl:   RGB->HSL, replace L, HSL->RGB\n"
        "  --saturation S  Saturation boost after combination (default: 1.0).\n"
        "                  L channel tends to desaturate; try 1.1-1.5.\n"
        "  --lightness K   MTF midtones balance on lightness (default: 0.5 = no change).\n"
        "                  K < 0.5 brightens, K > 0.5 darkens.\n"
        "  --mask-center [D]  Exclude central ellipse from background estimation\n"
        "                  (default D=0.6). Used with --auto for bright central objects.\n"
        "  --warmth W     Color temperature shift (-1..+1, default: 0).\n"
        "                  W>0: warmer, W<0: cooler. Used with --auto.\n"
        "  --eps E         Denominator floor for the ratio combine (default: 0.002,\n"
        "                  normalized [0,1]). Numerical stability in shadows only.\n"
        "  --bg-desat K    Pull low-signal pixels toward neutral grey: colour is\n"
        "                  kept for signal >= K sigma above background, faded to\n"
        "                  grey below (default: 0 = off). Try 2-4. Cleans background\n"
        "                  colour noise (canonical saturation-mask approach).\n"
        "\n"
        "Examples:\n"
        "  lrgb.py L.fit R.fit G.fit B.fit result.fit\n"
        "  lrgb.py L.fit R.fit G.fit B.fit result.fit --auto\n"
        "  lrgb.py L.fit RGB.fit result.fit --auto --method ratio\n"
        "  lrgb.py L.fit R.fit G.fit B.fit result.fit --saturation 1.3\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    value_opts = {'--method', '--saturation', '--lightness', '--warmth',
                  '--eps', '--bg-desat'}
    flag_set = {'--auto', '--slum', '--dry-run'}
    opts = {}
    flags = {}
    positional = []
    mask_center_d = None

    i = 0
    while i < len(args):
        if args[i] == '--mask-center':
            if i + 1 < len(args) and not args[i + 1].startswith('-'):
                try:
                    mask_center_d = float(args[i + 1])
                    i += 1
                except ValueError:
                    mask_center_d = 0.6
            else:
                mask_center_d = 0.6
            i += 1
        elif args[i] in flag_set:
            flags[args[i]] = True
            i += 1
        elif args[i] in value_opts:
            if i + 1 >= len(args):
                sys.stderr.write(f"Error: {args[i]} requires a value\n")
                sys.exit(1)
            opts[args[i]] = args[i + 1]
            i += 2
        elif args[i].startswith('-'):
            sys.stderr.write(f"Error: unknown option '{args[i]}'\n")
            usage()
        else:
            positional.append(args[i])
            i += 1

    if len(positional) == 5:
        mode = '4mono'
        inputs = {'L': positional[0], 'R': positional[1],
                  'G': positional[2], 'B': positional[3]}
        output = positional[4]
    elif len(positional) == 3:
        mode = 'lrgb'
        inputs = {'L': positional[0], 'RGB': positional[1]}
        output = positional[2]
    else:
        usage()

    method = opts.get('--method', 'ratio')
    if method not in ('hsl', 'ratio'):
        sys.stderr.write("Error: --method must be 'hsl' or 'ratio'\n")
        sys.exit(1)

    try:
        saturation = float(opts.get('--saturation', '1.0'))
    except ValueError:
        sys.stderr.write("Error: --saturation must be a number\n")
        sys.exit(1)

    try:
        lightness = float(opts.get('--lightness', '0.5'))
    except ValueError:
        sys.stderr.write("Error: --lightness must be a number\n")
        sys.exit(1)

    try:
        warmth = float(opts.get('--warmth', '0.0'))
    except ValueError:
        sys.stderr.write("Error: --warmth must be a number\n")
        sys.exit(1)

    try:
        eps = float(opts.get('--eps', '0.002'))
    except ValueError:
        sys.stderr.write("Error: --eps must be a number\n")
        sys.exit(1)
    if eps < 0:
        sys.stderr.write("Error: --eps must be >= 0\n")
        sys.exit(1)

    try:
        bg_desat = float(opts.get('--bg-desat', '0.0'))
    except ValueError:
        sys.stderr.write("Error: --bg-desat must be a number\n")
        sys.exit(1)
    if bg_desat < 0:
        sys.stderr.write("Error: --bg-desat must be >= 0\n")
        sys.exit(1)

    return {
        'mode': mode,
        'inputs': inputs,
        'output': output,
        'method': method,
        'saturation': saturation,
        'lightness': lightness,
        'auto': flags.get('--auto', False),
        'slum': flags.get('--slum', False),
        'dry_run': flags.get('--dry-run', False),
        'mask_center_d': mask_center_d,
        'warmth': warmth,
        'eps': eps,
        'bg_desat': bg_desat,
    }


# =========================================================================
# Color space conversion (vectorized numpy)
# =========================================================================

def rgb_to_hsl(rgb):
    """Convert RGB array (3, H, W) in [0,1] to HSL array (3, H, W) in [0,1]."""
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


def hsl_to_rgb(hsl):
    """Convert HSL array (3, H, W) in [0,1] to RGB array (3, H, W) in [0,1]."""
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


# =========================================================================
# MTF for lightness adjustment
# =========================================================================

def apply_mtf(x, m):
    """Midtone transfer function: (1-m)*x / (m + x*(1-2*m))."""
    if abs(m - 0.5) < 1e-10:
        return x
    num = (1.0 - m) * x
    den = m + x * (1.0 - 2.0 * m)
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(den != 0.0, num / den, np.where(x <= 0.0, 0.0, 1.0))
    return np.clip(result, 0.0, 1.0)


# =========================================================================
# LinearFit: match L brightness to RGB luminance
# =========================================================================

def linear_fit(target, reference):
    """Scale target to match reference by median and robust scale (MAD)."""
    t_valid = target[target > 0]
    r_valid = reference[reference > 0]

    if t_valid.size == 0 or r_valid.size == 0:
        return target.copy()

    med_t = np.median(t_valid)
    med_r = np.median(r_valid)
    mad_t = np.median(np.abs(t_valid - med_t))
    mad_r = np.median(np.abs(r_valid - med_r))

    if mad_t == 0:
        return target.copy()

    scale = mad_r / mad_t
    result = (target - med_t) * scale + med_r
    return np.clip(result, 0.0, 1.0)


# =========================================================================
# Super-luminance: max-SNR master luminance from L + (R+G+B)
# =========================================================================

def _estimate_noise(channel, bg_mask):
    """Robust per-pixel noise sigma from the high-pass residual over background.
    Only the L:syn RATIO matters for the weights, so a consistent estimator (the
    same for every channel) is sufficient."""
    from scipy.ndimage import gaussian_filter
    hp = channel - gaussian_filter(channel, 1.5)
    v = hp[bg_mask] if np.count_nonzero(bg_mask) > 100 else hp.ravel()
    m = np.median(v)
    s = 1.4826 * np.median(np.abs(v - m))
    if s <= 0:
        return float(np.std(v)) or 1.0
    keep = np.abs(v - m) < 4.0 * s          # drop stars / cosmic rays
    vk = v[keep]
    return float(1.4826 * np.median(np.abs(vk - np.median(vk))))


def _match_scale_raw(target, reference, bg_mask):
    """Rescale target to reference's flux scale (robust median + MAD), no clip."""
    t = target[bg_mask]
    r = reference[bg_mask]
    if t.size == 0 or r.size == 0:
        return target.copy()
    mt, mr = np.median(t), np.median(r)
    at = 1.4826 * np.median(np.abs(t - mt))
    ar = 1.4826 * np.median(np.abs(r - mr))
    if at <= 0:
        return target.copy()
    return (target - mt) * (ar / at) + mr


def super_luminance(lum, rgb, bg_mask, verbose=True):
    """Max-SNR master luminance: inverse-variance (SNR^2) blend of the measured
    L and the synthetic luminance L_syn = R+G+B, after scale matching.

    This is the luminance of the joint BLUE estimate under the L-band == R+G+B
    model: RGB noise enters only the SUM (luminance) projection here and the
    DIFFERENCES (chrominance) in the colour, so it is used once, not double
    counted. PSF/resolution is out of scope (channel FWHMs are ~equal), so a
    flat (not frequency-dependent) blend is used. Returns L in raw units.
    """
    lsyn = rgb[0] + rgb[1] + rgb[2]
    l_scaled = _match_scale_raw(lum, lsyn, bg_mask)

    sL = _estimate_noise(l_scaled, bg_mask)
    sR = _estimate_noise(rgb[0], bg_mask)
    sG = _estimate_noise(rgb[1], bg_mask)
    sB = _estimate_noise(rgb[2], bg_mask)
    S = sR * sR + sG * sG + sB * sB          # var(L_syn)
    if sL <= 0 or S <= 0:
        if verbose:
            print("  [slum] degenerate noise estimate; using L as-is")
        return lum.copy()

    inv_L, inv_S = 1.0 / (sL * sL), 1.0 / S
    wL = inv_L / (inv_L + inv_S)
    I = (l_scaled * inv_L + lsyn * inv_S) / (inv_L + inv_S)
    gain = np.sqrt(1.0 + (sL * sL) / S)

    if verbose:
        print(f"  [slum] sigma (L_syn units): L={sL:.3g} R={sR:.3g} G={sG:.3g} "
              f"B={sB:.3g} -> syn={np.sqrt(S):.3g}")
        print(f"  [slum] weights: L={wL:.3f} syn={1.0 - wL:.3f}  |  "
              f"SNR gain over L-only = x{gain:.3f} (+{(gain - 1) * 100:.0f}%)")
    return I


# =========================================================================
# LRGB combination methods
# =========================================================================

def combine_hsl(lum, rgb, saturation=1.0, lightness_k=0.5, bg_desat=0.0):
    """LRGB via HSL: replace L channel, preserve H and S.

    (No `eps`: HSL replaces lightness directly, there is no ratio division.
    `bg_desat` attenuates saturation in low-signal areas -> neutral background.)
    """
    hsl = rgb_to_hsl(rgb)
    l_old = hsl[2].copy()

    lum_fitted = linear_fit(lum, l_old)
    if abs(lightness_k - 0.5) > 1e-10:
        lum_fitted = apply_mtf(lum_fitted, lightness_k)

    hsl[2] = lum_fitted
    sat = hsl[1]
    if bg_desat > 0:
        sat = sat * _bg_desat_weight(rgb[0] + rgb[1] + rgb[2], bg_desat)
    if saturation != 1.0:
        sat = sat * saturation
    hsl[1] = np.clip(sat, 0.0, 1.0)

    return hsl_to_rgb(hsl)


def _bg_desat_weight(l_syn, strength):
    """Signal-significance mask for background desaturation: 0 in the noise
    background, ramping to 1 for signal >= `strength` sigma above background.
    `l_syn` = R+G+B. Returns w in [0,1]; multiply colour by it, blend the rest
    toward neutral grey."""
    from scipy.ndimage import gaussian_filter
    m = l_syn > 0
    if not np.any(m):
        return np.ones_like(l_syn)
    bg = np.median(l_syn[m])
    hp = (l_syn - gaussian_filter(l_syn, 1.5))[m]
    sigma = 1.4826 * np.median(np.abs(hp - np.median(hp)))
    if sigma <= 0:
        return np.ones_like(l_syn)
    return np.clip((l_syn - bg) / (strength * sigma), 0.0, 1.0)


def combine_ratio(lum, rgb, saturation=1.0, lightness_k=0.5, eps=0.002,
                  bg_desat=0.0):
    """LRGB via ratio: C_out = C * L_fit / (L_syn + eps), with L_syn = R+G+B.

    Using L_syn = R+G+B (not the perceptual (R+2G+B)/4) keeps the pipeline
    self-consistent with the super-luminance and gives the exact identity
    sum(C_out) = L_fit (colour ratios preserved). `eps` floors the denominator
    for numerical stability in the shadows; `bg_desat` (>0) pulls low-signal
    pixels toward neutral grey (canonical saturation-mask / background clean-up).
    """
    l_syn = rgb[0] + rgb[1] + rgb[2]

    lum_fitted = linear_fit(lum, l_syn)
    if abs(lightness_k - 0.5) > 1e-10:
        lum_fitted = apply_mtf(lum_fitted, lightness_k)

    scale = lum_fitted / (l_syn + eps)
    result = np.clip(rgb * scale[np.newaxis, :, :], 0.0, 1.0)

    if bg_desat > 0:
        w = _bg_desat_weight(l_syn, bg_desat)
        grey = (lum_fitted / 3.0)[np.newaxis, :, :]    # neutral, same sum-lum
        result = w[np.newaxis, :, :] * result + (1.0 - w[np.newaxis, :, :]) * grey
        result = np.clip(result, 0.0, 1.0)

    if saturation != 1.0:
        hsl = rgb_to_hsl(result)
        hsl[1] = np.clip(hsl[1] * saturation, 0.0, 1.0)
        result = hsl_to_rgb(hsl)

    return result


# =========================================================================
# I/O helpers
# =========================================================================

def load_mono(filepath):
    """Load a 2D mono FITS as float64. Returns (data, header, orig_dtype)."""
    with fits.open(filepath, memmap=False) as hdul:
        orig_dtype = hdul[0].data.dtype
        data = hdul[0].data.astype(np.float64)
        header = hdul[0].header
    if data.ndim != 2:
        raise ValueError(f"'{filepath}' is not a 2D mono FITS")
    return data, header, orig_dtype


def load_rgb(filepath):
    """Load a 3-channel RGB FITS as float64. Returns (data_3hw, header, orig_dtype)."""
    with fits.open(filepath, memmap=False) as hdul:
        orig_dtype = hdul[0].data.dtype
        data = hdul[0].data.astype(np.float64)
        header = hdul[0].header
    if data.ndim != 3 or data.shape[0] != 3:
        raise ValueError(f"'{filepath}' is not an RGB FITS (3,H,W)")
    return data, header, orig_dtype


def to_01(data, vmax):
    """Normalize to [0,1] using given max value."""
    if vmax <= 0:
        return np.zeros_like(data)
    return np.clip(data / vmax, 0.0, 1.0)


def from_01(data, vmax, orig_dtype):
    """Convert [0,1] back to original range and dtype."""
    work = data * vmax
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        work = np.clip(work, info.min, info.max)
        return np.rint(work).astype(orig_dtype)
    if np.issubdtype(orig_dtype, np.floating):
        result = work.astype(orig_dtype)
        result[~np.isfinite(result)] = 0
        return result
    return work.astype(np.float32)


# =========================================================================
# Auto pipeline: RGB balance + background flatten
# =========================================================================

def auto_preprocess(rgb_data, zero_mask, mask_center_d=None, warmth=0.0):
    """Run automatic RGB preprocessing pipeline on raw data.

    1. RGB balance by star photometry (autostar) + warmth shift
    2. Background flattening per channel

    Works on raw float64 data (not normalized to [0,1]).
    Modifies rgb_data in-place.
    """
    rb = _get_rgbbalance()

    # Step 1: RGB balance by star photometry
    print("  [auto] RGB balance by star photometry...")
    K, n_stars = rb.compute_star_balance(rgb_data)
    K = rb.apply_warmth(K, warmth)
    print(f"  [auto] Balance coefficients: R={K[0]:.4f} G={K[1]:.4f} B={K[2]:.4f} "
          f"({n_stars} stars)"
          + (f", warmth={warmth:+.2f}" if abs(warmth) > 1e-6 else ""))

    # Apply balance: shift to common black, scale by K
    stats = rb.compute_frame_stats(rgb_data, 5, 5)
    black, range_ref, K_combined = rb.compute_reference(stats, K, False)
    rb._apply_balance(rgb_data, 5, 5, black, range_ref, K_combined)

    # Restore zeros
    rgb_data[:, zero_mask] = 0.0

    # Step 2: Background flatten each channel
    print("  [auto] Background flattening...")
    for ch, name in enumerate(['R', 'G', 'B']):
        model = background.estimate_background(rgb_data[ch],
                                                mask_center_d=mask_center_d)
        offset = float(np.median(model[model > 0])) if np.any(model > 0) else 0.0
        rgb_data[ch] = rgb_data[ch] - model + offset
        print(f"    {name}: bg offset={offset:.1f}")

    # Restore zeros again (flatten may have shifted them)
    rgb_data[:, zero_mask] = 0.0
    # Clamp negative to zero
    rgb_data[rgb_data < 0] = 0.0


# =========================================================================
# Main processing
# =========================================================================

def process(config):
    """Main LRGB processing pipeline."""
    method = config['method']
    saturation = config['saturation']
    lightness_k = config['lightness']
    eps = config['eps']
    bg_desat = config['bg_desat']
    do_auto = config['auto']

    # Load inputs
    if config['mode'] == '4mono':
        print(f"LRGB: 4 mono files, method={method}"
              + (", auto" if do_auto else ""))
        lum_data, header, orig_dtype = load_mono(config['inputs']['L'])
        r_data, _, _ = load_mono(config['inputs']['R'])
        g_data, _, _ = load_mono(config['inputs']['G'])
        b_data, _, _ = load_mono(config['inputs']['B'])

        if not (lum_data.shape == r_data.shape == g_data.shape == b_data.shape):
            sys.stderr.write("Error: all input files must have the same dimensions\n")
            sys.exit(1)

        rgb_data = np.stack([r_data, g_data, b_data], axis=0)

    else:
        print(f"LRGB: L + RGB, method={method}"
              + (", auto" if do_auto else ""))
        lum_data, _, orig_dtype = load_mono(config['inputs']['L'])
        rgb_data, header, orig_dtype_rgb = load_rgb(config['inputs']['RGB'])

        if lum_data.shape != rgb_data.shape[1:]:
            sys.stderr.write(
                f"Error: L shape {lum_data.shape} != RGB shape {rgb_data.shape[1:]}\n")
            sys.exit(1)

    h, w = lum_data.shape
    print(f"  Image: {w}x{h}, dtype={orig_dtype}")

    # Zero mask (border from alignment)
    zero_mask_l = lum_data == 0
    zero_mask_rgb = np.any(rgb_data == 0, axis=0)
    zero_mask = zero_mask_l | zero_mask_rgb

    # Auto preprocessing pipeline
    if do_auto:
        auto_preprocess(rgb_data, zero_mask_rgb, config['mask_center_d'],
                        config['warmth'])

    # Optional super-luminance (opt-in): max-SNR master luminance from L + (R+G+B).
    # Default is L as-is, since the input L may already be a super-luminance.
    if config['slum'] or config['dry_run']:
        master = super_luminance(lum_data, rgb_data, ~zero_mask, verbose=True)
        if config['dry_run']:
            print("  [dry-run] super-luminance measured; no output written.")
            return
        lum_data = master

    # Compute normalization scale (max across all channels including L)
    all_valid = np.concatenate([
        lum_data[lum_data > 0].ravel(),
        rgb_data[0][rgb_data[0] > 0].ravel(),
        rgb_data[1][rgb_data[1] > 0].ravel(),
        rgb_data[2][rgb_data[2] > 0].ravel(),
    ])
    vmax = float(np.max(all_valid)) if all_valid.size > 0 else 1.0

    # Normalize to [0, 1]
    lum = to_01(lum_data, vmax)
    rgb = np.stack([to_01(rgb_data[ch], vmax) for ch in range(3)], axis=0)

    print(f"  L median: {np.median(lum[lum > 0]):.6f}")
    l_syn = rgb[0] + rgb[1] + rgb[2]
    print(f"  L_syn (R+G+B) median: {np.median(l_syn[l_syn > 0]):.6f}")

    # Combine
    if method == 'hsl':
        result = combine_hsl(lum, rgb, saturation, lightness_k, bg_desat=bg_desat)
    else:
        result = combine_ratio(lum, rgb, saturation, lightness_k,
                               eps=eps, bg_desat=bg_desat)

    # Restore zeros
    result[:, zero_mask] = 0.0

    # Denormalize
    out_data = np.stack([from_01(result[ch], vmax, orig_dtype) for ch in range(3)], axis=0)

    # Write output
    for key in ('BSCALE', 'BZERO'):
        if key in header:
            del header[key]

    header['NAXIS'] = 3
    header['NAXIS3'] = 3

    parts = [f"method={method}"]
    if config['slum']:
        parts.append("super-luminance")
    if do_auto:
        parts.append("auto")
    if saturation != 1.0:
        parts.append(f"saturation={saturation}")
    if abs(lightness_k - 0.5) > 1e-10:
        parts.append(f"lightness={lightness_k}")
    if method != 'hsl' and abs(eps - 0.002) > 1e-9:
        parts.append(f"eps={eps}")
    if bg_desat > 0:
        parts.append(f"bg-desat={bg_desat}")
    header['HISTORY'] = f"LRGB combined by lrgb.py: {', '.join(parts)}"

    out_dir = os.path.dirname(config['output'])
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fits.PrimaryHDU(out_data, header=header).writeto(
        config['output'], overwrite=True)

    print(f"  Result saved to {config['output']}")


# =========================================================================
# Entry point
# =========================================================================

def main():
    config = parse_args(sys.argv)
    process(config)


if __name__ == "__main__":
    main()
