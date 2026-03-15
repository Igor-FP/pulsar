#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
optstack - Optimal weighted stacking with sigma-fade clipping.

Multi-pass pipeline:
  Pass 0: Analyze frames (detect stars, measure SNR/FWHM, cache bg coefficients)
  Pass 1: Weighted average of normalized frames
  Pass 2: RMS deviation map (bg-subtracted comparison)
  Pass 3: Sigma-clip with fade (accumulate original normalized frames)
  Iterations 2+: repeat passes 2-3 with refined average

Weight modes: SNR^2, QFWHM (MTF-ranked FWHM), SNR^2 * QFWHM.
Background subtraction is used ONLY for sigma comparison, not in the sum.
"""

import sys
import os
import time
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils
import background
import star_utils


# =========================================================================
# CLI
# =========================================================================

def usage():
    sys.stderr.write(
        "Usage:\n"
        "  stack.py input_spec output.fit [options]\n"
        "\n"
        "Positional:\n"
        "  input_spec    Input files: wildcard, numbered, @list.txt\n"
        "  output.fit    Output stacked image\n"
        "\n"
        "Sigma clipping (optional, without --sigma: plain weighted sum):\n"
        "  --sigma [N]   Enable sigma clipping, threshold (default: 2.3)\n"
        "  --fade N      Sigma fade width in sigmas (default: 0.7)\n"
        "                Pixels beyond sigma are faded, beyond sigma+fade excluded.\n"
        "  --iter [N]    Sigma-clip iterations (default: 2)\n"
        "\n"
        "Weight options:\n"
        "  Default weight = SNR^2 (frames with better signal-to-noise get more weight).\n"
        "  --fwhm [K]    Enable FWHM weighting: weight *= MTF(FWHM_rank, K).\n"
        "                K is MTF midtones parameter (default: 0.5).\n"
        "                Higher K = stricter preference for sharp frames.\n"
        "  --nosnr       Disable SNR^2 weighting (weight by FWHM only).\n"
        "                Requires --fwhm, otherwise all weights = 1.\n"
        "\n"
        "Border options:\n"
        "  --border N    Zero-mask expansion in pixels (default: 1)\n"
        "\n"
        "Debug:\n"
        "  --debug [FILE]  Enable verbose output with timings.\n"
        "                  Without FILE: print to stdout.\n"
        "                  With FILE: write to log file.\n"
        "\n"
        "Performance:\n"
        "  --threads N   Max worker threads (default: cpu_count-1).\n"
        "                Actual count adapts: new threads launch only\n"
        "                while RAM usage < 80%.\n"
        "\n"
        "Background options (from autoflat mode 1):\n"
        "  --cell N      Background cell size in pixels (default: 64)\n"
        "  --clip K      Sigma clipping for background cells (default: 1.7)\n"
        "  --poly N      Background polynomial order (default: 3)\n"
        "  --mask-center [D]  Exclude central ellipse D*W x D*H from\n"
        "                 background estimation (default D=0.6).\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    value_opts = {'--fade', '--border', '--cell', '--clip', '--poly', '--threads'}
    flag_set = {'--nosnr'}
    opts = {}
    flags = {}
    positional = []
    sigma_val = None   # None = disabled
    iter_val = None    # None = use default
    fwhm_mtf_k = None  # None = FWHM weighting disabled
    debug_file = None  # None = no debug; '' = stderr; 'path' = file
    debug_enabled = False
    mask_center_d = None  # None = disabled

    i = 0
    while i < len(args):
        a = args[i]
        if a == '--mask-center':
            # --mask-center [D]: optional float argument, default 0.6
            if i + 1 < len(args) and not args[i + 1].startswith('-'):
                try:
                    mask_center_d = float(args[i + 1])
                    i += 1
                except ValueError:
                    mask_center_d = 0.6
            else:
                mask_center_d = 0.6
        elif a == '--debug':
            debug_enabled = True
            # --debug [FILE]: optional file argument
            if i + 1 < len(args) and not args[i + 1].startswith('-'):
                debug_file = args[i + 1]
                i += 1
            else:
                debug_file = ''  # stderr
        elif a == '--sigma':
            # --sigma [N]: optional float argument, default 2.3
            if i + 1 < len(args) and not args[i + 1].startswith('-'):
                try:
                    sigma_val = float(args[i + 1])
                    i += 1
                except ValueError:
                    sigma_val = 2.3
            else:
                sigma_val = 2.3
        elif a == '--iter':
            # --iter [N]: optional int argument, default 2
            if i + 1 < len(args) and not args[i + 1].startswith('-'):
                try:
                    iter_val = int(args[i + 1])
                    i += 1
                except ValueError:
                    iter_val = 2
            else:
                iter_val = 2
        elif a == '--fwhm':
            # --fwhm [K]: optional float argument, default 0.5
            if i + 1 < len(args) and not args[i + 1].startswith('-'):
                try:
                    fwhm_mtf_k = float(args[i + 1])
                    i += 1
                except ValueError:
                    fwhm_mtf_k = 0.5
            else:
                fwhm_mtf_k = 0.5
        elif a in flag_set:
            flags[a] = True
        elif a in value_opts:
            if i + 1 >= len(args):
                sys.stderr.write(f"Error: {a} requires a value\n")
                sys.exit(1)
            opts[a] = args[i + 1]
            i += 1
        elif a.startswith('-'):
            sys.stderr.write(f"Error: unknown option '{a}'\n")
            usage()
        else:
            positional.append(a)
        i += 1

    if len(positional) != 2:
        usage()

    def parse_float(key, default):
        val = opts.get(key)
        if val is None:
            return default
        try:
            return float(val)
        except ValueError:
            sys.stderr.write(f"Error: {key} must be a number\n")
            sys.exit(1)

    def parse_int(key, default):
        val = opts.get(key)
        if val is None:
            return default
        try:
            v = int(val)
            if v < 0:
                raise ValueError
            return v
        except ValueError:
            sys.stderr.write(f"Error: {key} must be a non-negative integer\n")
            sys.exit(1)

    use_sigma = sigma_val is not None
    use_fwhm = fwhm_mtf_k is not None
    use_snr = not flags.get('--nosnr', False)

    return {
        'input_spec': positional[0],
        'output_file': positional[1],
        # Sigma clipping
        'use_sigma': use_sigma,
        'sigma': sigma_val if use_sigma else 0,
        'fade': parse_float('--fade', 0.7) if use_sigma else 0,
        'iterations': (iter_val if iter_val is not None else 2) if use_sigma else 1,
        # Weights
        'use_snr': use_snr,
        'use_fwhm': use_fwhm,
        'fwhm_mtf_k': fwhm_mtf_k if use_fwhm else 0.5,
        # Debug
        'debug': debug_enabled,
        'debug_file': debug_file,
        # Border
        'border': parse_int('--border', 1),
        # Background (autoflat mode 1 defaults)
        'cell_size': parse_int('--cell', 64),
        'clip_k': parse_float('--clip', 1.7),
        'poly_order': parse_int('--poly', 3),
        'mask_center_d': mask_center_d,
        # Performance
        'threads': parse_int('--threads', max(1, (os.cpu_count() or 4) - 1)),
    }


# =========================================================================
# Logging
# =========================================================================

_debug_stream = None  # set in main()
_debug_enabled = False
_t0 = 0.0  # global start time


def _init_logging(config):
    global _debug_stream, _debug_enabled, _t0
    _t0 = time.time()
    _debug_enabled = config['debug']
    if _debug_enabled and config['debug_file']:
        _debug_stream = open(config['debug_file'], 'w')
    else:
        _debug_stream = None  # use stderr


def _close_logging():
    global _debug_stream
    if _debug_stream is not None and _debug_stream is not sys.stderr:
        _debug_stream.close()
        _debug_stream = None


def _elapsed():
    return time.time() - _t0


def log(msg):
    """Always printed (to stdout)."""
    print(msg, flush=True)


def dbg(msg):
    """Printed only in debug mode (to debug stream)."""
    if not _debug_enabled:
        return
    line = f"[{_elapsed():8.2f}s] {msg}"
    if _debug_stream is not None:
        _debug_stream.write(line + "\n")
        _debug_stream.flush()
    else:
        print(line, flush=True)


def progress(msg):
    sys.stdout.write(f"\r{msg}")
    sys.stdout.flush()


class Timer:
    """Context manager for timing stages."""
    def __init__(self, label):
        self.label = label
        self.start = 0
    def __enter__(self):
        self.start = time.time()
        dbg(f">> {self.label}")
        return self
    def __exit__(self, *args):
        dt = time.time() - self.start
        dbg(f"<< {self.label}: {dt:.2f}s")


# =========================================================================
# Helpers
# =========================================================================


def _mem_percent():
    """Return system memory usage percent, or -1 if psutil not available."""
    try:
        import psutil
        return psutil.virtual_memory().percent
    except ImportError:
        return -1


def _mem_ok(threshold=80):
    """Check if system memory usage is below threshold percent."""
    pct = _mem_percent()
    return pct < 0 or pct < threshold  # -1 = no psutil, assume OK


def mtf_func(x, m):
    """Midtone transfer function: (1-m)*x / (m + x*(1-2*m))."""
    if abs(m - 0.5) < 1e-10:
        return x
    return (1.0 - m) * x / (m + x * (1.0 - 2.0 * m))


def load_channel(filepath, channel):
    """Load a single channel from FITS as float64. Returns (data2d, header)."""
    with fits.open(filepath, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    if data is None:
        raise ValueError(f"'{filepath}' has no image data")
    if data.ndim == 3:
        return data[channel].astype(np.float64), header
    return data.astype(np.float64), header


def expand_zero_mask(data, border=1):
    """Create zero mask expanded by border pixels."""
    mask = (data == 0)
    if not np.any(mask) or border <= 0:
        return mask
    h, w = mask.shape
    padded = np.pad(mask, border, mode='constant', constant_values=False)
    expanded = np.zeros_like(mask)
    size = 2 * border + 1
    for dy in range(size):
        for dx in range(size):
            expanded |= padded[dy:dy+h, dx:dx+w]
    return expanded


def _run_parallel(worker_fn, items, max_threads, label, merge_fn):
    """
    Run worker_fn on items with adaptive thread count.
    New workers are launched while RAM usage < 80%.
    When >= 80%, wait for running workers to finish before launching more.

    worker_fn(item) -> result
    merge_fn(accumulated, result) -> updated accumulated
    Returns final accumulated result.
    """
    n = len(items)
    max_threads = min(max_threads, n)
    accumulated = None
    done = 0
    pending = {}  # future -> item

    with ThreadPoolExecutor(max_workers=max_threads) as pool:
        next_idx = 0

        while done < n:
            # Submit new tasks while memory is OK and we have capacity
            while next_idx < n and len(pending) < max_threads and _mem_ok():
                future = pool.submit(worker_fn, items[next_idx])
                pending[future] = next_idx
                next_idx += 1

            # If nothing pending and nothing submitted (memory pressure from start)
            # force at least one task
            if not pending and next_idx < n:
                future = pool.submit(worker_fn, items[next_idx])
                pending[future] = next_idx
                next_idx += 1

            # Wait for at least one to complete
            completed = set()
            for future in as_completed(pending):
                result = future.result()
                accumulated = merge_fn(accumulated, result)
                completed.add(future)
                done += 1
                active = len(pending) - len(completed)
                mem = _mem_percent()
                mem_str = f"  RAM:{mem:.0f}%" if mem >= 0 else ""
                progress(f"  [{label}] {done}/{n}  threads:{active}{mem_str}   ")
                # After collecting one result, check if we can submit more
                break

            for f in completed:
                del pending[f]

    return accumulated


def measure_noise(data, bg_model):
    """Measure per-pixel noise from residual (data - bg_model), sigma-clipped."""
    residual = data - bg_model
    valid = data > 0
    r = residual[valid]
    if r.size == 0:
        return 1.0
    for _ in range(3):
        med = np.median(r)
        s = np.std(r)
        if s == 0:
            break
        mask = np.abs(r - med) < 2.5 * s
        if np.all(mask):
            break
        r = r[mask]
    return float(np.std(r)) if r.size > 0 else 1.0


# =========================================================================
# Pass 0: Analyze frames
# =========================================================================

def analyze_frames(files, channel, config):
    """
    Analyze all frames: detect reference stars, measure SNR/FWHM,
    cache background polynomial coefficients.

    First frame is processed sequentially (reference star detection).
    Remaining frames are processed in parallel (ThreadPoolExecutor).

    Returns:
        ref: dict with reference params
        frames: list of per-frame param dicts (ordered same as files)
        bg_cache: list of (coeffs, terms) tuples for deferred bg rendering
    """
    n = len(files)
    need_fwhm = config['use_fwhm']

    log(f"[pass 0] Analyzing {n} frames...")

    # --- First frame: sequential (detect reference stars) ---
    filepath0 = files[0]
    fname0 = os.path.basename(filepath0)
    progress(f"  [analyze] 1/{n}: {fname0} (reference)")
    t_frame = time.time()

    data0, _ = load_channel(filepath0, channel)
    h, w = data0.shape

    t1 = time.time()
    coeffs0, terms0 = background.estimate_background_poly(
        data0, cell_size=config['cell_size'], clip_k=config['clip_k'],
        poly_order=config['poly_order'],
        mask_center_d=config['mask_center_d'])
    dbg(f"  bg_poly {fname0}: {time.time()-t1:.2f}s")

    bg_model0 = background.render_background(h, w, coeffs0, terms0)
    dark0 = float(np.median(bg_model0))
    noise0 = measure_noise(data0, bg_model0)

    t1 = time.time()
    catalog = star_utils.detect_stars(data0, snr=10, photometry=True)

    max_val = float(np.max(data0[data0 > 0])) if np.any(data0 > 0) else 1.0
    peak_mask = catalog.peak < 0.5 * max_val
    star_x = catalog.x[peak_mask]
    star_y = catalog.y[peak_mask]
    star_flux = catalog.flux[peak_mask]

    flux_mask = star_flux > 0
    star_x = star_x[flux_mask]
    star_y = star_y[flux_mask]
    star_flux = star_flux[flux_mask]

    dbg(f"  detect_stars {fname0}: {time.time()-t1:.2f}s, "
        f"total={len(catalog.x)}, filtered={len(star_x)}")

    if len(star_x) == 0:
        sys.stderr.write("\nError: no suitable reference stars found\n")
        sys.exit(1)

    fwhm0 = star_utils.estimate_fwhm(data0)
    if not np.isfinite(fwhm0) or fwhm0 < 1:
        fwhm0 = 3.0
    r_aperture = max(1.25 * fwhm0, 2.0)
    r_inner = max(2.0 * fwhm0, r_aperture + 2.0)
    r_outer = max(3.0 * fwhm0, r_inner + 3.0)

    light0 = float(np.median(star_flux))
    snr0 = light0 / noise0 if noise0 > 0 else 0

    ref = {
        'dark': dark0, 'light': light0, 'noise': noise0, 'fwhm': fwhm0,
        'star_x': star_x, 'star_y': star_y,
        'r_aperture': r_aperture, 'r_inner': r_inner, 'r_outer': r_outer,
        'h': h, 'w': w,
    }

    frame0 = {'dark': dark0, 'light': light0, 'noise': noise0,
              'fwhm': fwhm0, 'snr': snr0}
    dbg(f"  REF {fname0}: dark={dark0:.1f} light={light0:.1f} noise={noise0:.3f} "
        f"FWHM={fwhm0:.2f} SNR={snr0:.1f} [{time.time()-t_frame:.2f}s]")

    # --- Remaining frames: parallel (ThreadPoolExecutor) ---
    def _analyze_one(idx, filepath):
        """Analyze a single non-reference frame."""
        fname = os.path.basename(filepath)
        t0 = time.time()

        data, _ = load_channel(filepath, channel)

        coeffs, terms = background.estimate_background_poly(
            data, cell_size=config['cell_size'], clip_k=config['clip_k'],
            poly_order=config['poly_order'],
            mask_center_d=config['mask_center_d'])

        bg_model = background.render_background(h, w, coeffs, terms)
        dark = float(np.median(bg_model))
        noise = measure_noise(data, bg_model)

        net_flux = star_utils.measure_flux_at(
            data, ref['star_x'], ref['star_y'],
            ref['r_aperture'], ref['r_inner'], ref['r_outer'])

        positive = net_flux[net_flux > 0]
        light = float(np.median(positive)) if positive.size > 0 else ref['light']

        if need_fwhm:
            fwhm = star_utils.estimate_fwhm(data)
            if not np.isfinite(fwhm):
                fwhm = 99.0
        else:
            fwhm = 0.0

        snr = light / noise if noise > 0 else 0
        frame = {'dark': dark, 'light': light, 'noise': noise,
                 'fwhm': fwhm, 'snr': snr}

        dbg(f"  {fname}: dark={dark:.1f} light={light:.1f} noise={noise:.3f} "
            f"FWHM={fwhm:.2f} SNR={snr:.1f} [{time.time()-t0:.2f}s]")

        return idx, frame, coeffs, terms

    # Pre-allocate result arrays (indexed by file position)
    frames = [None] * n
    bg_cache = [None] * n
    frames[0] = frame0
    bg_cache[0] = (coeffs0, terms0)

    if n > 1:
        workers = min(config['threads'], n - 1)
        dbg(f"  Parallel analysis: {workers} threads for {n-1} frames")
        done = 1  # first frame already done

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_analyze_one, i, files[i]): i
                for i in range(1, n)
            }
            for future in as_completed(futures):
                idx, frame, coeffs, terms = future.result()
                frames[idx] = frame
                bg_cache[idx] = (coeffs, terms)
                done += 1
                progress(f"  [analyze] {done}/{n}")

    log(f"\n  Reference: dark={ref['dark']:.1f}, light={ref['light']:.1f}, "
        f"noise={ref['noise']:.3f}, FWHM={ref['fwhm']:.2f}, stars={len(ref['star_x'])}")

    # Log frame statistics
    snr_arr = np.array([f['snr'] for f in frames])
    log(f"  SNR range: {np.min(snr_arr):.1f} .. {np.max(snr_arr):.1f}, "
        f"median={np.median(snr_arr):.1f}")
    if need_fwhm:
        fwhm_arr = np.array([f['fwhm'] for f in frames])
        log(f"  FWHM range: {np.min(fwhm_arr):.2f} .. {np.max(fwhm_arr):.2f}, "
            f"median={np.median(fwhm_arr):.2f}")

    return ref, frames, bg_cache


# =========================================================================
# Weight computation
# =========================================================================

def compute_weights(frames, config):
    """Compute per-frame weights: SNR^2 * QFWHM (each factor optional)."""
    n = len(frames)

    # SNR^2 factor
    if config['use_snr']:
        snr_values = np.array([f['snr'] for f in frames])
        w_snr = snr_values ** 2
    else:
        w_snr = np.ones(n)

    # QFWHM factor
    if config['use_fwhm']:
        fwhm_values = np.array([f['fwhm'] for f in frames])
        mtf_k = config['fwhm_mtf_k']

        # Rank by FWHM (ascending = best first)
        sorted_indices = np.argsort(fwhm_values)
        ranks = np.empty(n, dtype=np.float64)
        for rank, idx in enumerate(sorted_indices):
            ranks[idx] = rank

        # Normalize: best (lowest FWHM) = 1.0, worst = 0.0
        if n > 1:
            nfwhm = (n - 1 - ranks) / (n - 1)
        else:
            nfwhm = np.ones(n)

        # Apply MTF
        w_fwhm = np.array([mtf_func(v, mtf_k) for v in nfwhm])

        # Log FWHM ranking
        best_idx = sorted_indices[0]
        worst_idx = sorted_indices[-1]
        log(f"  FWHM ranking: best={fwhm_values[best_idx]:.2f}px (w={w_fwhm[best_idx]:.3f}), "
            f"worst={fwhm_values[worst_idx]:.2f}px (w={w_fwhm[worst_idx]:.3f}), "
            f"MTF K={mtf_k}")
    else:
        w_fwhm = np.ones(n)

    return w_snr * w_fwhm


# =========================================================================
# Pass 1: Weighted average
# =========================================================================

def weighted_average(files, channel, ref, frames, weights, config):
    """Compute weighted average of normalized frames (adaptive parallel)."""
    n = len(files)
    h, w = ref['h'], ref['w']
    darkRef = np.float32(ref['dark'])
    lightRef = np.float32(ref['light'])
    border = config['border']

    log(f"[pass 1] Computing weighted average ({n} files)...")

    def _process_one(idx):
        data, _ = load_channel(files[idx], channel)
        data = data.astype(np.float32)
        zmask = expand_zero_mask(data, border)
        valid = ~zmask
        lightK = lightRef / np.float32(frames[idx]['light'])
        scaled = (data - np.float32(frames[idx]['dark'])) * lightK + darkRef
        w_i = weights[idx]
        # Return sparse: only valid pixels contribute
        acc = np.zeros((h, w), dtype=np.float64)
        wsum = np.zeros((h, w), dtype=np.float64)
        acc[valid] = (w_i * scaled[valid]).astype(np.float64)
        wsum[valid] = w_i
        return acc, wsum

    def _merge(accumulated, result):
        local_acc, local_wsum = result
        if accumulated is None:
            return [local_acc, local_wsum]
        accumulated[0] += local_acc
        accumulated[1] += local_wsum
        return accumulated

    result = _run_parallel(_process_one, list(range(n)),
                           config['threads'], "average", _merge)
    acc, wsum = result

    nonzero = wsum > 0
    acc[nonzero] /= wsum[nonzero]
    acc[~nonzero] = 0

    log(f"\n  [average] Done")
    return acc


# =========================================================================
# Pass 2: Deviation map
# =========================================================================

def compute_deviation(files, channel, ref, frames, weights, avg, bg_cache,
                      config, old_dev=None):
    """Compute RMS deviation map (adaptive parallel)."""
    n = len(files)
    h, w = ref['h'], ref['w']
    lightRef = np.float32(ref['light'])
    noiseRef = np.float32(ref['noise'])
    sigma_val = np.float32(config['sigma'])
    sigma_fade = np.float32(config['fade'])
    border = config['border']

    # Background model of the average
    avrBg = background.estimate_background(
        avg, cell_size=config['cell_size'], clip_k=config['clip_k'],
        poly_order=config['poly_order'],
        mask_center_d=config['mask_center_d'])

    log(f"[pass 2] Computing deviation map ({n} files)...")

    # Shared float32 arrays (read-only across threads)
    avg32 = avg.astype(np.float32)
    avrBg32 = avrBg.astype(np.float32)
    old_dev32 = old_dev.astype(np.float32) if old_dev is not None else None

    def _process_one(idx):
        data, _ = load_channel(files[idx], channel)
        data = data.astype(np.float32)
        zmask = expand_zero_mask(data, border)
        valid = ~zmask

        lightK = lightRef / np.float32(frames[idx]['light'])
        noiseK = np.float32(frames[idx]['noise']) * lightK
        Kdev = noiseRef / noiseK if noiseK > 0 else np.float32(1.0)

        coeffs, terms = bg_cache[idx]
        bg_model = background.render_background(h, w, coeffs, terms).astype(np.float32)

        scaledV2 = (data - bg_model) * lightK + avrBg32
        dev = (scaledV2 - avg32) * Kdev

        if old_dev32 is not None:
            abs_dev = np.abs(dev)
            sig = np.zeros_like(dev)
            nz = old_dev32 > 0
            sig[nz] = abs_dev[nz] / old_dev32[nz]

            fade = np.ones_like(dev)
            in_fade = (sig > sigma_val) & (sig < sigma_val + sigma_fade)
            fade[in_fade] = (sigma_val + sigma_fade - sig[in_fade]) / sigma_fade
            fade[sig >= sigma_val + sigma_fade] = 0.0
        else:
            fade = np.ones_like(dev)

        fade[~valid] = 0.0

        sum_dev2 = (dev ** 2 * fade).astype(np.float64)
        sum_w = fade.astype(np.float64)
        return sum_dev2, sum_w

    def _merge(accumulated, result):
        local_dev2, local_w = result
        if accumulated is None:
            return [local_dev2, local_w]
        accumulated[0] += local_dev2
        accumulated[1] += local_w
        return accumulated

    result = _run_parallel(_process_one, list(range(n)),
                           config['threads'], "deviation", _merge)
    sum_dev2, sum_w = result

    # RMS
    dev_image = np.zeros((h, w), dtype=np.float64)
    nonzero = sum_w > 0
    dev_image[nonzero] = np.sqrt(sum_dev2[nonzero] / sum_w[nonzero])

    log(f"\n  [deviation] Done, median RMS={np.median(dev_image[nonzero]):.3f}")
    return dev_image


# =========================================================================
# Pass 3: Sigma-clip with fade
# =========================================================================

def sigma_clip_stack(files, channel, ref, frames, weights, avg, dev_image,
                     bg_cache, config):
    """Sigma-clip with fade, accumulate normalized original frames (adaptive parallel)."""
    n = len(files)
    h, w = ref['h'], ref['w']
    darkRef = np.float32(ref['dark'])
    lightRef = np.float32(ref['light'])
    noiseRef = np.float32(ref['noise'])
    sigma_val = np.float32(config['sigma'])
    sigma_fade_width = np.float32(config['fade'])
    border = config['border']

    # Background model of the average
    avrBg = background.estimate_background(
        avg, cell_size=config['cell_size'], clip_k=config['clip_k'],
        poly_order=config['poly_order'],
        mask_center_d=config['mask_center_d'])

    log(f"[pass 3] Sigma-clipping ({n} files, sigma={config['sigma']:g}, fade={config['fade']:g})...")

    # Shared float32 arrays (read-only across threads)
    avg32 = avg.astype(np.float32)
    avrBg32 = avrBg.astype(np.float32)
    dev32 = dev_image.astype(np.float32)

    def _process_one(idx):
        data, _ = load_channel(files[idx], channel)
        data = data.astype(np.float32)
        zmask = expand_zero_mask(data, border)
        valid = ~zmask

        lightK = lightRef / np.float32(frames[idx]['light'])
        w_i = weights[idx]

        scaledV = (data - np.float32(frames[idx]['dark'])) * lightK + darkRef

        coeffs, terms = bg_cache[idx]
        bg_model = background.render_background(h, w, coeffs, terms).astype(np.float32)
        scaledV2 = (data - bg_model) * lightK + avrBg32

        noiseK = np.float32(frames[idx]['noise']) * lightK
        Kdev = noiseRef / noiseK if noiseK > 0 else np.float32(1.0)

        dev = (scaledV2 - avg32) * Kdev

        abs_dev = np.abs(dev)
        sig = np.zeros_like(dev)
        nz = dev32 > 0
        sig[nz] = abs_dev[nz] / dev32[nz]

        fade = np.ones_like(dev)
        in_fade = (sig > sigma_val) & (sig < sigma_val + sigma_fade_width)
        fade[in_fade] = (sigma_val + sigma_fade_width - sig[in_fade]) / sigma_fade_width
        fade[sig >= sigma_val + sigma_fade_width] = 0.0
        fade[~valid] = 0.0

        pixel_weight = w_i * fade
        acc = (pixel_weight * scaledV).astype(np.float64)
        wsum = pixel_weight.astype(np.float64)
        rejected = float(np.sum(1.0 - fade[valid]))
        pixels = int(np.sum(valid))
        return acc, wsum, rejected, pixels

    def _merge(accumulated, result):
        local_acc, local_wsum, local_rej, local_pix = result
        if accumulated is None:
            return [local_acc, local_wsum, local_rej, local_pix]
        accumulated[0] += local_acc
        accumulated[1] += local_wsum
        accumulated[2] += local_rej
        accumulated[3] += local_pix
        return accumulated

    result = _run_parallel(_process_one, list(range(n)),
                           config['threads'], "sigma-clip", _merge)
    acc, wsum, total_rejected, total_pixels = result

    # Divide
    out = np.zeros((h, w), dtype=np.float64)
    nonzero = wsum > 0
    out[nonzero] = acc[nonzero] / wsum[nonzero]

    reject_pct = total_rejected / total_pixels * 100 if total_pixels > 0 else 0
    log(f"\n  [sigma-clip] Done, rejected {reject_pct:.2f}% of pixel contributions")

    return out


# =========================================================================
# Channel stacking orchestrator
# =========================================================================

def stack_channel(files, channel, config):
    """Main stacking function for a single channel."""
    n = len(files)
    t0 = time.time()

    # Pass 0: Analyze
    with Timer("Pass 0: Analyze frames"):
        ref, frames, bg_cache = analyze_frames(files, channel, config)

    log(f"  Threads: max {config['threads']} (adaptive by RAM)")

    # Compute weights
    weights = compute_weights(frames, config)

    # Log weight distribution
    w_sorted = np.sort(weights)[::-1]
    wsum = np.sum(weights)
    log(f"  Weights: min={w_sorted[-1]:.4f}, median={np.median(weights):.4f}, "
        f"max={w_sorted[0]:.4f}")
    # Top contributors
    cumw = np.cumsum(w_sorted) / wsum * 100
    for pct in (50, 80, 90):
        idx = int(np.searchsorted(cumw, pct))
        log(f"    {pct}% of total weight: top {idx+1}/{n} frames")

    # Debug: per-frame weight table
    if _debug_enabled:
        sorted_by_weight = np.argsort(weights)[::-1]
        dbg("  Per-frame weights (sorted):")
        for rank, idx in enumerate(sorted_by_weight):
            f = frames[idx]
            dbg(f"    #{rank+1:3d} w={weights[idx]:10.4f}  SNR={f['snr']:7.1f}  "
                f"FWHM={f['fwhm']:5.2f}  dark={f['dark']:7.1f}  "
                f"light={f['light']:9.1f}  noise={f['noise']:6.3f}  "
                f"{os.path.basename(files[idx])}")

    # Pass 1: Weighted average
    with Timer("Pass 1: Weighted average"):
        avg = weighted_average(files, channel, ref, frames, weights, config)

    if not config['use_sigma']:
        # No sigma clipping: plain weighted average is the result
        elapsed = time.time() - t0
        log(f"\n[stack] Channel done in {elapsed:.1f}s (no sigma clipping)")
        return avg

    # Iterative sigma clipping
    old_dev = None
    result = None

    for iteration in range(config['iterations']):
        log(f"\n[iteration {iteration+1}/{config['iterations']}]")

        # Pass 2: Deviation map
        with Timer(f"Pass 2: Deviation map (iter {iteration+1})"):
            dev = compute_deviation(files, channel, ref, frames, weights, avg,
                                    bg_cache, config, old_dev=old_dev)

        # Pass 3: Sigma clip
        with Timer(f"Pass 3: Sigma-clip (iter {iteration+1})"):
            result = sigma_clip_stack(files, channel, ref, frames, weights, avg,
                                      dev, bg_cache, config)

        # Update for next iteration
        avg = result
        old_dev = dev

    elapsed = time.time() - t0
    log(f"\n[stack] Channel done in {elapsed:.1f}s")
    dbg(f"Channel total: {elapsed:.2f}s")

    return result


# =========================================================================
# Main
# =========================================================================

def main():
    config = parse_args(sys.argv)
    _init_logging(config)

    input_spec = config['input_spec']
    output_file = config['output_file']

    # Expand input files
    files = batch_utils.expand_input_spec(input_spec)
    if len(files) < 2:
        sys.stderr.write("Error: need at least 2 files to stack\n")
        sys.exit(1)

    log(f"Weighted stacking: {len(files)} files")
    if config['use_sigma']:
        log(f"  sigma={config['sigma']}, fade={config['fade']}, "
            f"iter={config['iterations']}")
    else:
        log(f"  sigma clipping: OFF (plain weighted sum)")
    parts = []
    if config['use_snr']:
        parts.append("SNR^2")
    if config['use_fwhm']:
        parts.append(f"FWHM(MTF K={config['fwhm_mtf_k']})")
    log(f"  weight: {' x '.join(parts) if parts else '1.0 (equal)'}")
    log(f"  border={config['border']}")
    log(f"  bg: cell={config['cell_size']}, clip={config['clip_k']}, "
        f"poly={config['poly_order']}")

    # Determine mono/RGB from first file
    with fits.open(files[0], memmap=False) as hdul:
        data0 = hdul[0].data
        header = hdul[0].header.copy()
        orig_dtype = data0.dtype

    if data0.ndim == 3:
        nchannels = data0.shape[0]
        log(f"  RGB mode: {nchannels} channels\n")
        channels = []
        for ch in range(nchannels):
            log(f"{'='*60}")
            log(f"Channel {ch+1}/{nchannels}")
            log(f"{'='*60}")
            result = stack_channel(files, ch, config)
            channels.append(result)
        result_data = np.stack(channels, axis=0)
    else:
        log(f"  Mono mode\n")
        result_data = stack_channel(files, 0, config)

    # Convert back to original dtype
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        result_data = np.clip(result_data, info.min, info.max)
        result_data = np.rint(result_data).astype(orig_dtype)
    elif np.issubdtype(orig_dtype, np.floating):
        result_data = result_data.astype(orig_dtype)

    # Accumulate total exposure
    total_exposure = 0
    for f in files:
        try:
            with fits.open(f, memmap=False) as h:
                total_exposure += float(h[0].header.get('EXPTIME', 0))
        except Exception:
            pass

    # Write output
    # Remove BSCALE/BZERO to avoid astropy conflicts with unsigned int types
    for key in ('BSCALE', 'BZERO'):
        if key in header:
            del header[key]

    header['NCOMBINE'] = len(files)
    if total_exposure > 0:
        header['EXPTIME'] = total_exposure
    if config['use_sigma']:
        header['HISTORY'] = (f'stack: {len(files)} files, sigma={config["sigma"]}, '
                             f'fade={config["fade"]}, iter={config["iterations"]}')
    else:
        header['HISTORY'] = f'stack: {len(files)} files, weighted sum (no sigma clip)'
    w_parts = []
    if config['use_snr']:
        w_parts.append('SNR^2')
    if config['use_fwhm']:
        w_parts.append(f'FWHM(K={config["fwhm_mtf_k"]})')
    header['HISTORY'] = f'stack: weight={" x ".join(w_parts) if w_parts else "equal"}'

    hdu = fits.PrimaryHDU(data=result_data, header=header)
    hdu.writeto(output_file, overwrite=True)

    total_elapsed = time.time() - _t0
    log(f"\nResult saved to {output_file}")
    if total_exposure > 0:
        log(f"Total exposure: {total_exposure:.1f}s ({total_exposure/3600:.2f}h)")
    log(f"Total time: {total_elapsed:.1f}s ({total_elapsed/60:.1f}m)")
    dbg(f"Total time: {total_elapsed:.2f}s")

    _close_logging()


if __name__ == "__main__":
    main()
