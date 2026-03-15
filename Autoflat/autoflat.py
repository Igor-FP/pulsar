#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
from astropy.io import fits

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils
import background

# Global debug flag
DEBUG = False


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  autoflat.py [-d] [--mode {1,2}] input_spec [output_spec] [options]\n"
        "\n"
        "  input_spec   - single file, wildcard, numbered, @list.txt\n"
        "  output_spec  - output for subtracted result (required with --subtract)\n"
        "\n"
        "Output options (at least one required):\n"
        "  --save FILE    Save background model as FITS (pattern for batch)\n"
        "  --subtract     Subtract model from input, write to output_spec\n"
        "  (default: --subtract if neither specified)\n"
        "\n"
        "Common options:\n"
        "  -d             Enable debug outputs (save intermediate FITS)\n"
        "  --mode N       Algorithm: 1 = cell-based (default), 2 = min-binning\n"
        "  --poly N       Polynomial order (default: 3 for mode 1, 1 for mode 2)\n"
        "\n"
        "Mode 1 options (cell-based background estimation):\n"
        "  --cell N       Cell size in pixels (default: 64)\n"
        "  --clip K       Sigma clipping threshold (default: 1.7)\n"
        "  --median1 N    First median filter size on grid (default: 3)\n"
        "  --median2 N    Second median filter size on extended grid (default: 5)\n"
        "  --border N     Border extension in cells (default: 2)\n"
        "  --mask-center [D]   Exclude central ellipse D*W x D*H from background\n"
        "                 estimation (default D=0.6). Useful when a bright\n"
        "                 diffuse object (galaxy/nebula) dominates the center.\n"
        "                 The polynomial fit interpolates smoothly across the hole.\n"
        "\n"
        "Mode 2 options (min-binning + polynomial):\n"
        "  --radius N     Median filter radius (default: 2)\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    # Named options with values
    value_opts = {
        '--mode', '--poly', '--cell', '--clip',
        '--median1', '--median2', '--border', '--radius', '--save',
    }
    # Flags
    flag_set = {'-d', '--subtract'}

    flags = {}
    opts = {}
    positional = []

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
        elif a in flag_set:
            flags[a] = True
        elif a in value_opts:
            if i + 1 >= len(args):
                sys.stderr.write(f"Error: {a} requires a value\n")
                sys.exit(1)
            opts[a] = args[i + 1]
            i += 1
        elif a.startswith('-') and not a[1:].replace('.', '', 1).replace('+', '', 1).replace('-', '', 1).isdigit():
            sys.stderr.write(f"Error: unknown option '{a}'\n")
            usage()
        else:
            positional.append(a)
        i += 1

    if len(positional) < 1:
        usage()

    input_spec = positional[0]
    output_spec = positional[1] if len(positional) >= 2 else None

    if len(positional) > 2:
        usage()

    # Parse mode
    mode = int(opts.get('--mode', '1'))
    if mode not in (1, 2):
        sys.stderr.write("Error: --mode must be 1 or 2\n")
        sys.exit(1)

    # Determine operations
    subtract = flags.get('--subtract', False)
    save_spec = opts.get('--save')

    # Default: if neither --save nor --subtract, assume --subtract
    if not subtract and save_spec is None:
        subtract = True

    if subtract and output_spec is None:
        sys.stderr.write("Error: output_spec required for subtraction\n")
        usage()

    # Parse numeric options
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

    def parse_float(key, default):
        val = opts.get(key)
        if val is None:
            return default
        try:
            return float(val)
        except ValueError:
            sys.stderr.write(f"Error: {key} must be a number\n")
            sys.exit(1)

    poly_default = 3 if mode == 1 else 1

    config = {
        'mode': mode,
        'debug': flags.get('-d', False),
        'input_spec': input_spec,
        'output_spec': output_spec,
        'subtract': subtract,
        'save_spec': save_spec,
        'poly_order': parse_int('--poly', poly_default),
        # Mode 1 options
        'cell_size': parse_int('--cell', 64),
        'clip_k': parse_float('--clip', 1.7),
        'median1': parse_int('--median1', 3),
        'median2': parse_int('--median2', 5),
        'border': parse_int('--border', 2),
        'mask_center_d': mask_center_d,
        # Mode 2 options
        'radius': parse_int('--radius', 2),
    }

    return config


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
# Mode 2: Zero expansion
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
# Mode 2: Median (fast, zero-aware)
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
        feat = np.where(mask, 0, 1)
        _, indices = distance_transform_edt(feat, return_indices=True)
        filled = data.copy()
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
      - Zeros are ignored via sentinel trick.
    """
    size = 2 * radius + 1
    h, w = data.shape
    log(f"  [step2] Median(ignore zeros) radius={radius}, size={size}, image={w}x{h}")

    if not HAVE_SCIPY:
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

    med[med >= sentinel] = 0.0

    out = fill_zeros_with_nearest(med)
    log("  [step2] Completed (median_filter + sentinel)")
    log_stats("step2_out", out)
    return out


# =========================
# Mode 2: Min-binning 2x2
# =========================

def min_pool_2x_ignore_zeros(data):
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
# Shared: Polynomial fit & render
# =========================

def fit_poly2d(data, order, coords=None):
    """
    Fit 2D polynomial surface of given order to data != 0.
    coords: optional (norm_x, norm_y) arrays for grid points.
            If None, uses pixel coordinates normalized to [-1, 1].
    """
    h, w = data.shape
    log(f"  [polyfit] Requested poly order={order} on grid {w}x{h}")

    if coords is not None:
        norm_x, norm_y = coords
        mask = (data != 0.0)
        N = int(mask.sum())
        if N == 0:
            log("  [polyfit] No non-zero points for fit, aborting")
            return None, None
        z = data[mask]
        x_flat = norm_x[mask]
        y_flat = norm_y[mask]
    else:
        yy, xx = np.mgrid[0:h, 0:w]
        mask = (data != 0.0)
        N = int(mask.sum())
        if N == 0:
            log("  [polyfit] No non-zero pixels for fit, aborting")
            return None, None

        z = data[mask]
        x = xx[mask].astype(np.float64)
        y = yy[mask].astype(np.float64)

        x_flat = (2.0 * (x / (w - 1.0)) - 1.0) if w > 1 else np.zeros_like(x)
        y_flat = (2.0 * (y / (h - 1.0)) - 1.0) if h > 1 else np.zeros_like(y)

    # Choose maximal allowed order so that num_terms <= N
    def num_terms(k):
        return (k + 1) * (k + 2) // 2

    used_order = min(order, 20)
    while used_order > 0 and num_terms(used_order) > N:
        used_order -= 1

    if used_order < order:
        log(f"  [polyfit] Reducing poly order from {order} to {used_order} "
            f"(points={N}, terms={num_terms(used_order)})")
    else:
        log(f"  [polyfit] Using poly order={used_order}, points={N}, terms={num_terms(used_order)}")

    order = used_order

    terms = [(i, j) for i in range(order + 1) for j in range(order + 1 - i)]
    cols = [(x_flat ** i) * (y_flat ** j) for (i, j) in terms]
    A = np.vstack(cols).T

    coeffs, _, _, _ = np.linalg.lstsq(A, z, rcond=None)
    log("  [polyfit] Polynomial fit completed")
    return coeffs, terms


def render_poly2d(height, width, coeffs, terms):
    log(f"  [render] Rendering model {width}x{height}")
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
    for c, (i, j) in zip(coeffs, terms):
        if c != 0.0:
            model += c * (x ** i) * (y ** j)

    log_stats("render_model", model)
    log("  [render] Model done")
    return model


# =========================
# Mode 1: Cell-based background estimation (delegates to lib/background.py)
# =========================

def mode1_background(data64, cell_size, clip_k, poly_order, median1, median2, border,
                     mask_center_d=None, debug_prefix=None, header=None):
    """Mode 1: Cell-based background estimation with polynomial fit."""
    h, w = data64.shape
    log(f"  [mode1] Image {w}x{h}, cell={cell_size}, clip={clip_k}"
        + (f", mask-center={mask_center_d}" if mask_center_d else ""))

    model = background.estimate_background(
        data64, cell_size=cell_size, clip_k=clip_k, poly_order=poly_order,
        median1=median1, median2=median2, border=border,
        mask_center_d=mask_center_d)

    log_stats("mode1_model", model)
    write_debug_image(debug_prefix, "mode1_model", model, header)
    return model


# =========================
# Mode 2 pipeline (original algorithm)
# =========================

def mode2_background(data64, poly_order, radius, debug_prefix=None, header=None):
    """Mode 2: Zero-expand -> median -> min-binning -> polynomial fit."""
    log("  [mode2] Start")

    step1 = expand_zero_mask(data64, iterations=2)
    write_debug_image(debug_prefix, "step1_zeroexp", step1, header)

    bg = median_aperture_ignore_zeros(step1, radius=radius)
    write_debug_image(debug_prefix, "step2_median", bg, header)

    small = iterative_min_binning(bg, target_min_size=8)
    write_debug_image(debug_prefix, "step3_binning", small, header)

    coeffs, terms = fit_poly2d(small, poly_order)
    if coeffs is None:
        log("  [mode2] Fit failed, returning zeros")
        return np.zeros(data64.shape, dtype=np.float64)

    h, w = data64.shape
    model = render_poly2d(h, w, coeffs, terms)
    write_debug_image(debug_prefix, "step6_model", model, header)

    log("  [mode2] Done")
    return model


# =========================
# Per-file processing
# =========================

def _build_model_2d(data64, config, debug_prefix=None, header=None):
    """Build background model for a single 2D plane."""
    if config['mode'] == 1:
        return mode1_background(
            data64, config['cell_size'], config['clip_k'],
            config['poly_order'], config['median1'], config['median2'],
            config['border'], mask_center_d=config['mask_center_d'],
            debug_prefix=debug_prefix, header=header)
    else:
        return mode2_background(
            data64, config['poly_order'], config['radius'],
            debug_prefix=debug_prefix, header=header)


def process_file(infile, outfile, save_file, config):
    label = outfile or save_file
    log(f"[file] {infile} -> {label}")
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"'{infile}' has no image data")

        orig = hdul[0].data
        header = hdul[0].header
        ndim = orig.ndim

        if ndim not in (2, 3):
            raise ValueError(f"'{infile}' has unsupported shape {orig.shape} (need 2D or 3D)")

        debug_prefix = (outfile or save_file) if DEBUG else None

        if ndim == 2:
            # Mono image
            nplanes = 1
            planes = [orig]
        else:
            # RGB/multi-channel: shape = (C, H, W)
            nplanes = orig.shape[0]
            planes = [orig[c] for c in range(nplanes)]
            log(f"  [rgb] {nplanes} channels, processing independently")

        models = []
        corrected_planes = []

        for c, plane in enumerate(planes):
            if nplanes > 1:
                log(f"  [channel {c+1}/{nplanes}]")
                ch_debug = f"{debug_prefix}_ch{c+1}" if debug_prefix else None
            else:
                ch_debug = debug_prefix

            plane64 = to_float64(plane)
            log_stats(f"input_ch{c+1}" if nplanes > 1 else "input", plane64)

            model = _build_model_2d(plane64, config, debug_prefix=ch_debug, header=header)
            models.append(model)

            if outfile is not None:
                zero_mask = (plane == 0)
                offset = float(np.median(model))
                log(f"  [subtract] Correction: orig - model + offset({offset:.3f})")
                corrected = plane64 - model + offset
                corrected[zero_mask] = 0.0
                log_stats(f"corrected_ch{c+1}" if nplanes > 1 else "corrected", corrected)
                corrected_planes.append(from_float64(corrected, plane.dtype))

        # Save model if requested
        if save_file is not None:
            if ndim == 2:
                model_data = models[0].astype(np.float32)
            else:
                model_data = np.stack([m.astype(np.float32) for m in models], axis=0)
            model_hdu = fits.PrimaryHDU(data=model_data, header=header.copy())
            model_hdu.header['HISTORY'] = f'autoflat: background model (mode {config["mode"]})'
            model_hdu.writeto(save_file, overwrite=True)
            log(f"  [save] Model -> {save_file}")

        # Write subtracted result if requested
        if outfile is not None:
            if ndim == 2:
                hdul[0].data = corrected_planes[0]
            else:
                hdul[0].data = np.stack(corrected_planes, axis=0)
            hdul[0].header['HISTORY'] = f'autoflat: background subtracted (mode {config["mode"]})'
            hdul.writeto(outfile, overwrite=True)

    log(f"[file] Done {infile}")


# =========================
# Main
# =========================

def main():
    global DEBUG

    config = parse_args(sys.argv)
    DEBUG = config['debug']

    input_spec = config['input_spec']
    output_spec = config['output_spec']
    save_spec = config['save_spec']
    subtract = config['subtract']

    # Build file lists
    try:
        if subtract:
            io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
        else:
            # No subtraction, just input files
            in_files = batch_utils.expand_input_spec(input_spec)
            io_pairs = [(f, None) for f in in_files]
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    # Build save file list
    save_files = [None] * len(io_pairs)
    if save_spec is not None:
        try:
            # Use build_io_file_lists to expand save pattern
            save_pairs = batch_utils.build_io_file_lists(input_spec, save_spec)
            if len(save_pairs) != len(io_pairs):
                sys.stderr.write("Error: --save pattern produced different file count than input\n")
                sys.exit(1)
            save_files = [sp[1] for sp in save_pairs]
        except Exception as e:
            sys.stderr.write(f"Error expanding --save pattern: {e}\n")
            sys.exit(1)

    log(f"Found {len(io_pairs)} file(s). Mode={config['mode']}, poly={config['poly_order']}, debug={DEBUG}")

    total = len(io_pairs)
    for i, ((infile, outfile), save_file) in enumerate(zip(io_pairs, save_files), start=1):
        try:
            process_file(infile, outfile, save_file, config)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")
    sys.stderr.flush()


if __name__ == "__main__":
    main()
