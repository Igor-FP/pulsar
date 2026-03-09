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

Supports 2D (mono) and 3D (color, 3×H×W) FITS images.

Color modes (3D only):
  - Independent (default): MTF applied to each channel separately.
  - --preserve_color: MTF applied to per-pixel average luminance,
    channels scaled proportionally to preserve R/avg and B/avg ratios.
"""

import sys
import os
import numpy as np
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import shared utilities (two levels up from Personal/Mtf/)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils

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
        "\n"
        "Positional arguments:\n"
        "  input_spec       Input file(s): single file, wildcard (*.fit),\n"
        "                   numbered sequence (img0001.fit), or list (@list.txt)\n"
        "  output_spec      Output: single file, numbered pattern (out0001.fit),\n"
        "                   or directory\n"
        "  K                Midtones balance (0..1 exclusive)\n"
        "                     K < 0.5 -> brighten midtones\n"
        "                     K = 0.5 -> no change (identity)\n"
        "                     K > 0.5 -> darken midtones\n"
        "\n"
        "Options:\n"
        "  --preserve_color  For color (3-channel) images: apply MTF to\n"
        "                    per-pixel average, scale R/G/B proportionally\n"
        "                    (default: independent per-channel MTF)\n"
        "  --k2 K2           Second MTF parameter for dual-MTF blending.\n"
        "                    Result = mtf(x,K)*(1-blend) + mtf(x,K2)*blend\n"
        "  --blend B         Blend factor for dual MTF (default: 0.5).\n"
        "                    0.0 = only K, 1.0 = only K2\n"
        "  --clip [P]        After MTF, clip black/white by percentile P\n"
        "                    (default P=0.001). Ignores pixels <= 0.\n"
        "                    Remaps [black_P, white_P] to [0, 1].\n"
        "\n"
        "Examples:\n"
        "  mtf.py img.fit out.fit 0.2\n"
        "  mtf.py img.fit out.fit 0.15 --k2 0.45\n"
        "  mtf.py img.fit out.fit 0.15 --k2 0.45 --blend 0.3\n"
        "  mtf.py img0001.fit out0001.fit 0.25 --preserve_color\n"
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

    preserve_color = False
    k2 = None
    blend = 0.5
    clip = None  # None = disabled, float = percentile

    # Extract option flags
    i = 0
    positional = []
    while i < len(args):
        if args[i] == "--preserve_color":
            preserve_color = True
            i += 1
        elif args[i] == "--clip":
            clip = 0.001  # default percentile
            i += 1
            # Optional percentile value
            if i < len(args) and not args[i].startswith("--"):
                try:
                    clip = float(args[i])
                    i += 1
                except ValueError:
                    pass
        elif args[i] == "--k2":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --k2 requires a value.\n")
                sys.exit(1)
            k2 = _parse_k(args[i+1], "K2")
            i += 2
        elif args[i] == "--blend":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --blend requires a value.\n")
                sys.exit(1)
            try:
                blend = float(args[i+1])
            except ValueError:
                sys.stderr.write("Error: --blend value must be a number.\n")
                sys.exit(1)
            if blend < 0.0 or blend > 1.0:
                sys.stderr.write("Error: --blend must be in range [0, 1].\n")
                sys.exit(1)
            i += 2
        else:
            positional.append(args[i])
            i += 1

    if len(positional) != 3:
        usage()

    input_spec = positional[0]
    output_spec = positional[1]
    m = _parse_k(positional[2], "K")

    return input_spec, output_spec, m, preserve_color, k2, blend, clip


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

    # Safe division: den=0 only at extreme m values with x=0 or x=1
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(den != 0.0, num / den, np.where(x <= 0.0, 0.0, 1.0))

    return np.clip(result, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Normalization helpers
# ---------------------------------------------------------------------------

def normalize_to_01(data, orig_dtype):
    """Normalize data to [0, 1] range based on dtype."""
    work = data.astype(np.float64)
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        work = (work - info.min) / (info.max - info.min)
    else:
        # Float: clip to [0, 1]
        np.clip(work, 0.0, 1.0, out=work)
    return work


def denormalize_from_01(data, orig_dtype):
    """Convert [0, 1] data back to original dtype range."""
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        work = data * (info.max - info.min) + info.min
        work = np.clip(work, info.min, info.max)
        work = np.rint(work)
        return work.astype(orig_dtype)
    if np.issubdtype(orig_dtype, np.floating):
        arr = data.astype(orig_dtype)
        bad = ~np.isfinite(arr)
        if np.any(bad):
            arr[bad] = 0
        return arr
    return data.astype(np.float32)


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------

def apply_dual_mtf(x, m1, m2, blend):
    """Apply two MTFs and blend: result = mtf(x,m1)*(1-blend) + mtf(x,m2)*blend."""
    r1 = apply_mtf(x, m1)
    r2 = apply_mtf(x, m2)
    return r1 * (1.0 - blend) + r2 * blend


def apply_mtf_preserve_color(work_3d, m, k2, blend):
    """Apply MTF to 3-channel image preserving color ratios.

    MTF is applied to the per-pixel average (luminance).
    Each channel is then scaled by new_avg/avg to preserve R/avg, B/avg ratios.
    """
    # Per-pixel average across channels (shape: H, W)
    avg = np.mean(work_3d, axis=0)

    # Apply MTF (single or dual) to average
    if k2 is not None:
        new_avg = apply_dual_mtf(avg, m, k2, blend)
    else:
        new_avg = apply_mtf(avg, m)

    # Scale factor: new_avg / avg, protect division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        scale = np.where(avg > 0.0, new_avg / avg, 0.0)

    # Scale all channels by the same factor
    result = work_3d * scale[np.newaxis, :, :]
    return np.clip(result, 0.0, 1.0)


def clip_black_white(work, clip_pct, orig_dtype):
    """Clip black/white levels by percentile and remap with margin.

    Ignores pixels <= 0 (zero-padded regions from crop).
    Works on normalized [0, 1] data. Applied per-channel for 3D.

    Remaps [black_pct, white_pct] → [margin_lo, 1 - margin_hi]
    so the clipped percentile is not lost but pushed to a small offset.
    """
    # Compute margins in normalized [0, 1] space
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

        # Preserve zeros, remap [black, white] → [dst_lo, dst_hi]
        zero_mask = plane <= 0
        result = dst_lo + (plane - black) / (white - black) * (dst_hi - dst_lo)
        result = np.clip(result, 0.0, 1.0)
        result[zero_mask] = 0.0
        return result

    if work.ndim == 2:
        return _clip_plane(work)

    # 3D color: compute black/white from average across channels,
    # apply same remap to all channels to preserve color balance
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


def process_file(infile, outfile, m, preserve_color, k2, blend, clip):
    """Process a single FITS file."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")

        data = hdul[0].data
        header = hdul[0].header.copy()
        orig_dtype = data.dtype

    # Normalize to [0, 1]
    work = normalize_to_01(data, orig_dtype)

    # Apply MTF
    is_color = (work.ndim == 3 and work.shape[0] == 3)
    if is_color and preserve_color:
        work = apply_mtf_preserve_color(work, m, k2, blend)
    elif k2 is not None:
        work = apply_dual_mtf(work, m, k2, blend)
    else:
        work = apply_mtf(work, m)

    # Optional black/white clip
    if clip is not None:
        work = clip_black_white(work, clip, orig_dtype)

    # Convert back
    out_data = denormalize_from_01(work, orig_dtype)

    # Update header
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    mode_str = "preserve_color" if (is_color and preserve_color) else "independent"
    if k2 is not None:
        header["HISTORY"] = (
            f"MTF applied by mtf.py: K={m}, K2={k2}, blend={blend}, mode={mode_str}")
    else:
        header["HISTORY"] = f"MTF applied by mtf.py: m={m}, mode={mode_str}"
    if clip is not None:
        header["HISTORY"] = f"  clip={clip}%"

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fits.PrimaryHDU(out_data, header=header).writeto(outfile, overwrite=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    input_spec, output_spec, m, preserve_color, k2, blend, clip = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    parts = [f"K={m}"]
    if k2 is not None:
        parts.append(f"K2={k2}, blend={blend}")
    if preserve_color:
        parts.append("preserve_color")
    if clip is not None:
        parts.append(f"clip={clip}%")
    ncpu = os.cpu_count() or 1
    workers = min(max(1, ncpu - 1), total)
    print(f"Applying MTF ({', '.join(parts)}) to {total} file(s), threads={workers}")

    done = 0
    errors = 0

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(process_file, infile, outfile, m, preserve_color, k2, blend, clip):
                os.path.basename(infile)
            for infile, outfile in io_pairs
        }

        for future in as_completed(futures):
            fname = futures[future]
            done += 1
            try:
                future.result()
                sys.stderr.write(f"\r{done}/{total}")
                sys.stderr.flush()
            except Exception as e:
                errors += 1
                sys.stderr.write(f"\n  Error '{fname}': {e}\n")

    sys.stderr.write(f"\nDone. {done - errors} OK, {errors} errors.\n")


if __name__ == "__main__":
    main()
