#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fits2tiff - Convert FITS images to TIFF, JPEG, or PNG.

Output format determined by:
  1. --format option (priority)
  2. Output file extension (.tif/.tiff, .jpg/.jpeg, .png)

Supports 2D (mono) and 3D RGB (3×H×W) FITS images.
JPEG and PNG are always 8-bit. TIFF supports 8/16/32-bit.

Optional --stretch for auto screen transfer (background-aware nonlinear stretch).
"""

import sys
import os
import re
import warnings
import numpy as np
from astropy.io import fits

try:
    from PIL import Image
except ImportError:
    sys.stderr.write(
        "Error: Pillow is required but not installed.\n"
        "Install it with: pip install Pillow\n"
    )
    sys.exit(1)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils

# Suppress numpy warnings from NaN operations
warnings.filterwarnings('ignore', message='.*All-NaN slice.*', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*Degrees of freedom <= 0.*', category=RuntimeWarning)


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  fits2tiff.py [options] input_spec output_spec\n"
        "\n"
        "    input_spec   - single file, wildcard (*.fit), numbered, @list.txt\n"
        "    output_spec  - single file, numbered pattern, or wildcard (*.tif)\n"
        "\n"
        "Format (auto-detected from extension, or override):\n"
        "  --format F     Force output format: tiff, jpeg, png\n"
        "                 Default: by extension (.tif/.tiff, .jpg/.jpeg, .png)\n"
        "\n"
        "Bit depth (TIFF only, JPEG/PNG always 8-bit):\n"
        "  --bits 8       Linear stretch [min,max] -> [0,255]\n"
        "  --bits 16      Clamp to [0,65535], uint16\n"
        "  --bits 32      Float32, no scaling\n"
        "  (default)      Auto: uint8/int8->8, uint16/int16->16, else->32\n"
        "\n"
        "Auto stretch (screen transfer function):\n"
        "  --stretch [T]  Auto black/white + MTF to place median at T\n"
        "                 (default T=0.1). Forces 8-bit output.\n"
        "  --keepcolor [K]  Preserve color during stretch (default K=1.0).\n"
        "                 K=1: full preservation via HSL. K=0: per-channel.\n"
        "\n"
        "Downscale:\n"
        "  --bin N        Bin NxN pixels by averaging (2 or 4)\n"
        "\n"
        "JPEG quality:\n"
        "  --jpeg Q       JPEG compression quality 1-100 (default: 99)\n"
        "\n"
        "Other:\n"
        "  --flip         Flip vertically (FITS bottom-left to image top-left)\n"
        "\n"
        "Examples:\n"
        "  fits2tiff img.fit img.tif\n"
        "  fits2tiff img.fit img.jpg --stretch\n"
        "  fits2tiff img.fit img.png --stretch 0.15 --keepcolor\n"
        "  fits2tiff *.fit *.jpg --stretch\n"
        "  fits2tiff img.fit img.tif --bits 16\n"
        "  fits2tiff img.fit img.jpg --format jpeg\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    bits = None
    flip = False
    fmt = None
    stretch_target = None  # None = disabled
    keepcolor_k = None
    bin_factor = None      # None = no binning
    jpeg_quality = 99

    positional = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == '--flip':
            flip = True
            i += 1
        elif a == '--bits':
            if i + 1 >= len(args):
                sys.stderr.write("Error: --bits requires a value (8, 16, or 32).\n")
                sys.exit(1)
            if args[i+1] not in ('8', '16', '32'):
                sys.stderr.write("Error: --bits must be 8, 16, or 32.\n")
                sys.exit(1)
            bits = int(args[i+1])
            i += 2
        elif a == '--format':
            if i + 1 >= len(args):
                sys.stderr.write("Error: --format requires a value.\n")
                sys.exit(1)
            fmt = args[i+1].lower()
            if fmt not in ('tiff', 'jpeg', 'png'):
                sys.stderr.write("Error: --format must be tiff, jpeg, or png.\n")
                sys.exit(1)
            i += 2
        elif a == '--stretch':
            stretch_target = 0.1  # default
            i += 1
            if i < len(args) and not args[i].startswith('-'):
                try:
                    stretch_target = float(args[i])
                    i += 1
                except ValueError:
                    pass
        elif a == '--keepcolor':
            keepcolor_k = 1.0  # default
            i += 1
            if i < len(args) and not args[i].startswith('-'):
                try:
                    keepcolor_k = float(args[i])
                    i += 1
                except ValueError:
                    pass
        elif a == '--bin':
            if i + 1 >= len(args):
                sys.stderr.write("Error: --bin requires a value (2 or 4).\n")
                sys.exit(1)
            if args[i+1] not in ('2', '4'):
                sys.stderr.write("Error: --bin must be 2 or 4.\n")
                sys.exit(1)
            bin_factor = int(args[i+1])
            i += 2
        elif a == '--jpeg':
            if i + 1 >= len(args):
                sys.stderr.write("Error: --jpeg requires a quality value (1-100).\n")
                sys.exit(1)
            try:
                jpeg_quality = int(args[i+1])
            except ValueError:
                sys.stderr.write("Error: --jpeg quality must be an integer.\n")
                sys.exit(1)
            jpeg_quality = max(1, min(100, jpeg_quality))
            i += 2
        elif a.startswith('-'):
            sys.stderr.write(f"Error: unknown option '{a}'\n")
            usage()
        else:
            positional.append(a)
            i += 1

    if len(positional) != 2:
        usage()

    return {
        'input_spec': positional[0],
        'output_spec': positional[1],
        'bits': bits,
        'flip': flip,
        'format': fmt,
        'stretch_target': stretch_target,
        'keepcolor_k': keepcolor_k,
        'bin_factor': bin_factor,
        'jpeg_quality': jpeg_quality,
    }


# =========================================================================
# Format detection
# =========================================================================

_EXT_FORMAT = {
    '.tif': 'tiff', '.tiff': 'tiff',
    '.jpg': 'jpeg', '.jpeg': 'jpeg',
    '.png': 'png',
}

def detect_format(outfile, fmt_override):
    """Determine output format from override or file extension."""
    if fmt_override:
        return fmt_override
    ext = os.path.splitext(outfile)[1].lower()
    return _EXT_FORMAT.get(ext, 'tiff')


# =========================================================================
# Output file pattern (supports tif/tiff/jpg/jpeg/png extensions)
# =========================================================================

_OUT_SEQ_RE = re.compile(r"^(.*?)(\d+)(\.\w+)$", re.IGNORECASE)


def _apply_wildcard_output(inputs, output_spec):
    """Replace * in output_spec with each input file's stem."""
    out_dir = os.path.dirname(output_spec)
    out_pattern = os.path.basename(output_spec)
    pairs = []
    for inp in inputs:
        stem = os.path.splitext(os.path.basename(inp))[0]
        outname = out_pattern.replace("*", stem)
        outpath = os.path.join(out_dir, outname) if out_dir else outname
        pairs.append((inp, outpath))
    return pairs


def build_io_pairs(input_spec, output_spec):
    """Build (input, output) pairs."""
    inputs = batch_utils.expand_input_spec(input_spec)
    if not inputs:
        raise ValueError("No input files found.")

    if "*" in output_spec:
        return _apply_wildcard_output(inputs, output_spec)

    if len(inputs) == 1:
        return [(inputs[0], output_spec)]

    base = os.path.basename(output_spec)
    m = _OUT_SEQ_RE.match(base)
    if not m:
        raise ValueError(
            "Output pattern must contain a numeric field when multiple "
            "input files are provided (e.g. out0001.tif), or use "
            "wildcard (e.g. *.jpg) to preserve input names.")

    prefix, digits, ext = m.group(1), m.group(2), m.group(3)
    width = len(digits)
    start_index = int(digits)
    out_dir = os.path.dirname(os.path.abspath(output_spec)) or "."

    pairs = []
    for i, inp in enumerate(inputs):
        fname = f"{prefix}{str(start_index + i).zfill(width)}{ext}"
        pairs.append((inp, os.path.join(out_dir, fname)))
    return pairs


# =========================================================================
# Auto stretch (screen transfer function)
# =========================================================================

def _compute_autoblack(plane):
    """Black level from darkest 3% of background cells. Returns value in data units."""
    from background import _build_cell_grid
    grid = _build_cell_grid(plane, 64, clip_k=1.7)
    valid = grid[grid > 0]
    if valid.size == 0:
        return 0.0
    n_dark = max(1, int(valid.size * 0.03))
    threshold = np.partition(valid, n_dark)[n_dark - 1]
    dark_values = valid[valid <= threshold]
    med = float(np.median(dark_values))
    mad = float(np.median(np.abs(dark_values - med)))
    return max(0.0, med - 5.0 * mad)


def _compute_autowhite(plane):
    """White level from brightest 0.01% of pixels."""
    valid = plane[plane > 0].ravel()
    if valid.size == 0:
        return 1.0
    n_bright = max(1, int(valid.size * 0.0001))
    idx = valid.size - n_bright
    partitioned = np.partition(valid, idx)
    return float(np.median(partitioned[idx:]))


def _apply_mtf(x, m):
    """Midtone transfer function on [0,1] array."""
    if abs(m - 0.5) < 1e-10:
        return x.copy()
    num = (1.0 - m) * x
    den = m + x * (1.0 - 2.0 * m)
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(den != 0.0, num / den, np.where(x <= 0.0, 0.0, 1.0))
    return np.clip(result, 0.0, 1.0)


def _find_mtf_k(median_norm, target):
    """Find MTF parameter K that maps median_norm to target.

    Binary search: K < 0.5 brightens, K > 0.5 darkens.
    We want mtf(median_norm, K) = target.
    """
    if median_norm <= 0 or median_norm >= 1:
        return 0.5
    # Analytical solution from MTF formula:
    # target = (1-K)*med / (K + med*(1-2K))
    # Solve for K:
    # target*(K + med - 2*K*med) = (1-K)*med
    # target*K + target*med - 2*target*K*med = med - K*med
    # K*(target - 2*target*med + med) = med - target*med
    # K = med*(1 - target) / (target - 2*target*med + med)
    denom = target - 2.0 * target * median_norm + median_norm
    if abs(denom) < 1e-15:
        return 0.5
    k = median_norm * (1.0 - target) / denom
    return np.clip(k, 0.001, 0.999)


def _rgb_to_hsl_pixel(r, g, b):
    """Vectorized RGB->HSL for stretch color preservation."""
    cmax = np.maximum(np.maximum(r, g), b)
    cmin = np.minimum(np.minimum(r, g), b)
    delta = cmax - cmin
    l = (cmax + cmin) / 2.0

    s = np.zeros_like(l)
    mask = delta > 0
    low = mask & (l <= 0.5)
    high = mask & (l > 0.5)
    dl = cmax[low] + cmin[low]
    dh = 2.0 - cmax[high] - cmin[high]
    s[low] = np.where(dl > 0, delta[low] / dl, 0.0)
    s[high] = np.where(dh > 0, delta[high] / dh, 0.0)

    h = np.zeros_like(l)
    with np.errstate(divide='ignore', invalid='ignore'):
        inv_d = np.where(delta > 0, 1.0 / delta, 0.0)
    rm = mask & (cmax == r)
    h[rm] = ((g[rm] - b[rm]) * inv_d[rm]) % 6.0
    gm = mask & (cmax == g) & ~rm
    h[gm] = (b[gm] - r[gm]) * inv_d[gm] + 2.0
    bm = mask & ~rm & ~gm
    h[bm] = (r[bm] - g[bm]) * inv_d[bm] + 4.0
    h /= 6.0
    h = h % 1.0
    return h, s, l


def _hsl_to_rgb_pixel(h, s, l):
    """Vectorized HSL->RGB."""
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

    return np.clip(r + m, 0, 1), np.clip(g + m, 0, 1), np.clip(b + m, 0, 1)


def auto_stretch(data, target, keepcolor_k):
    """Apply auto screen transfer function.

    1. Compute black/white (linked for RGB)
    2. Remap to [0,1]
    3. Find MTF K to place median at target
    4. Apply MTF (with optional color preservation)
    5. Scale to [0, 255] uint8

    data: 2D (H,W) or 3D (3,H,W) float64 array.
    Returns: uint8 array (H,W) or (H,W,3).
    """
    is_rgb = data.ndim == 3 and data.shape[0] == 3

    if is_rgb:
        # Linked black/white
        blacks = [_compute_autoblack(data[ch]) for ch in range(3)]
        whites = [_compute_autowhite(data[ch]) for ch in range(3)]
        black = min(blacks)
        white = max(whites)
    else:
        black = _compute_autoblack(data)
        white = _compute_autowhite(data)

    bw_range = white - black
    if bw_range <= 0:
        bw_range = 1.0

    if is_rgb:
        # Normalize to [0,1]
        channels = []
        for ch in range(3):
            norm = np.clip((data[ch] - black) / bw_range, 0.0, 1.0)
            channels.append(norm)

        # Compute median from luminance
        lum = (channels[0] + 2.0 * channels[1] + channels[2]) / 4.0
        med_valid = lum[lum > 0]
        median_norm = float(np.median(med_valid)) if med_valid.size > 0 else 0.5

        # Find MTF K
        k = _find_mtf_k(median_norm, target)

        if keepcolor_k is not None and keepcolor_k > 0:
            # HSL color preservation
            h, s, l = _rgb_to_hsl_pixel(channels[0], channels[1], channels[2])
            l_new = _apply_mtf(l, k)

            if keepcolor_k >= 1.0:
                r_out, g_out, b_out = _hsl_to_rgb_pixel(h, s, l_new)
            else:
                # Blend: per-channel vs HSL
                r_hsl, g_hsl, b_hsl = _hsl_to_rgb_pixel(h, s, l_new)
                r_pc = _apply_mtf(channels[0], k)
                g_pc = _apply_mtf(channels[1], k)
                b_pc = _apply_mtf(channels[2], k)
                kk = keepcolor_k
                r_out = r_pc * (1 - kk) + r_hsl * kk
                g_out = g_pc * (1 - kk) + g_hsl * kk
                b_out = b_pc * (1 - kk) + b_hsl * kk
        else:
            # Per-channel MTF
            r_out = _apply_mtf(channels[0], k)
            g_out = _apply_mtf(channels[1], k)
            b_out = _apply_mtf(channels[2], k)

        # To uint8 H×W×3
        result = np.stack([
            np.rint(np.clip(r_out * 255, 0, 255)).astype(np.uint8),
            np.rint(np.clip(g_out * 255, 0, 255)).astype(np.uint8),
            np.rint(np.clip(b_out * 255, 0, 255)).astype(np.uint8),
        ], axis=-1)

    else:
        # Mono
        norm = np.clip((data - black) / bw_range, 0.0, 1.0)
        med_valid = norm[norm > 0]
        median_norm = float(np.median(med_valid)) if med_valid.size > 0 else 0.5
        k = _find_mtf_k(median_norm, target)
        stretched = _apply_mtf(norm, k)
        result = np.rint(np.clip(stretched * 255, 0, 255)).astype(np.uint8)

    return result


def _bin_2d(arr, factor):
    """Bin a 2D array by factor (mean pooling). Trims to exact multiple."""
    h, w = arr.shape
    nh = (h // factor) * factor
    nw = (w // factor) * factor
    trimmed = arr[:nh, :nw]
    return trimmed.reshape(nh // factor, factor, nw // factor, factor).mean(axis=(1, 3))


def _bin_data(data, factor):
    """Bin 2D or 3D (C,H,W) data by factor."""
    if data.ndim == 2:
        return _bin_2d(data, factor)
    return np.stack([_bin_2d(data[ch], factor) for ch in range(data.shape[0])], axis=0)


# =========================================================================
# Bit depth helpers
# =========================================================================

def auto_bits(dtype):
    """Choose output bit depth based on FITS data type."""
    if dtype in (np.uint8, np.int8):
        return 8
    if dtype in (np.uint16, np.int16):
        return 16
    return 32


def _scale_channel(work, out_bits):
    """Scale a single 2D channel to the target bit depth."""
    if out_bits == 8:
        vmin = work.min()
        vmax = work.max()
        if vmax > vmin:
            scaled = (work - vmin) / (vmax - vmin) * 255.0
        else:
            scaled = np.zeros_like(work)
        return np.rint(np.clip(scaled, 0, 255)).astype(np.uint8)
    elif out_bits == 16:
        return np.rint(np.clip(work, 0, 65535)).astype(np.uint16)
    else:
        vmin = work.min()
        vmax = work.max()
        if vmax > vmin:
            work = (work - vmin) / (vmax - vmin)
        else:
            work = np.zeros_like(work)
        return work.astype(np.float32)


# =========================================================================
# Custom TIFF writers (Pillow can't handle 16-bit or 32-bit RGB)
# =========================================================================

def _save_tiff_rgb16(result, outfile):
    """Write 16-bit RGB TIFF manually."""
    import struct

    h, w = result.shape[:2]
    pixels = np.ascontiguousarray(result).tobytes()
    pixel_bytes = len(pixels)

    tags = [
        (256, 4, 1, w),
        (257, 4, 1, h),
        (258, 3, 3, 0),
        (259, 3, 1, 1),
        (262, 3, 1, 2),
        (273, 4, 1, 0),
        (277, 3, 1, 3),
        (278, 4, 1, h),
        (279, 4, 1, pixel_bytes),
        (284, 3, 1, 1),
    ]
    n_tags = len(tags)

    ifd_offset = 8
    ifd_size = 2 + n_tags * 12 + 4
    bps_offset = ifd_offset + ifd_size
    data_offset = bps_offset + 6

    fixed_tags = []
    for tag, typ, cnt, val in tags:
        if tag == 258:
            val = bps_offset
        elif tag == 273:
            val = data_offset
        fixed_tags.append((tag, typ, cnt, val))

    buf = bytearray()
    buf += struct.pack('<2sHI', b'II', 42, ifd_offset)
    buf += struct.pack('<H', n_tags)
    for tag, typ, cnt, val in fixed_tags:
        buf += struct.pack('<HHII', tag, typ, cnt, val)
    buf += struct.pack('<I', 0)
    buf += struct.pack('<HHH', 16, 16, 16)
    buf += pixels

    with open(outfile, 'wb') as f:
        f.write(buf)


def _save_tiff_rgb32f(result, outfile):
    """Write 32-bit float RGB TIFF manually."""
    import struct

    h, w = result.shape[:2]
    pixels = np.ascontiguousarray(result.astype(np.float32)).tobytes()
    pixel_bytes = len(pixels)

    tags = [
        (256, 4, 1, w),
        (257, 4, 1, h),
        (258, 3, 3, 0),
        (259, 3, 1, 1),
        (262, 3, 1, 2),
        (273, 4, 1, 0),
        (277, 3, 1, 3),
        (278, 4, 1, h),
        (279, 4, 1, pixel_bytes),
        (284, 3, 1, 1),
        (339, 3, 1, 3),
    ]
    n_tags = len(tags)

    ifd_offset = 8
    ifd_size = 2 + n_tags * 12 + 4
    bps_offset = ifd_offset + ifd_size
    data_offset = bps_offset + 6

    fixed_tags = []
    for tag, typ, cnt, val in tags:
        if tag == 258:
            val = bps_offset
        elif tag == 273:
            val = data_offset
        fixed_tags.append((tag, typ, cnt, val))

    buf = bytearray()
    buf += struct.pack('<2sHI', b'II', 42, ifd_offset)
    buf += struct.pack('<H', n_tags)
    for tag, typ, cnt, val in fixed_tags:
        buf += struct.pack('<HHII', tag, typ, cnt, val)
    buf += struct.pack('<I', 0)
    buf += struct.pack('<HHH', 32, 32, 32)
    buf += pixels

    with open(outfile, 'wb') as f:
        f.write(buf)


# =========================================================================
# File conversion
# =========================================================================

def _save_image(img, outfile, out_fmt, jpeg_quality):
    """Save PIL Image in the requested format."""
    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    if out_fmt == 'jpeg':
        img.save(outfile, 'JPEG', quality=jpeg_quality)
    elif out_fmt == 'png':
        img.save(outfile, 'PNG')
    else:
        img.save(outfile, 'TIFF')


def convert_file(infile, outfile, config):
    """Convert a single FITS file to TIFF/JPEG/PNG."""
    bits = config['bits']
    flip = config['flip']
    stretch_target = config['stretch_target']
    keepcolor_k = config['keepcolor_k']
    bin_factor = config['bin_factor']
    jpeg_quality = config['jpeg_quality']

    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"File '{infile}' has no primary image data.")
        data = hdul[0].data
        orig_dtype = data.dtype

    is_rgb = data.ndim == 3 and data.shape[0] == 3
    if data.ndim != 2 and not is_rgb:
        raise ValueError(f"File '{infile}': unsupported shape {data.shape}.")

    out_fmt = detect_format(outfile, config['format'])

    # Flip if requested
    if flip:
        if is_rgb:
            data = np.stack([np.flipud(data[ch]) for ch in range(3)], axis=0)
        else:
            data = np.flipud(data)

    # Binning (downscale by averaging NxN blocks)
    if bin_factor is not None:
        data = _bin_data(data.astype(np.float64), bin_factor)

    # Auto stretch mode
    if stretch_target is not None:
        work = data.astype(np.float64)
        work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)
        result = auto_stretch(work, stretch_target, keepcolor_k)
        img = Image.fromarray(result)
        _save_image(img, outfile, out_fmt, jpeg_quality)
        return

    # Standard conversion (no stretch)
    out_bits = bits if bits is not None else auto_bits(orig_dtype)

    # Force 8-bit for JPEG/PNG
    if out_fmt in ('jpeg', 'png'):
        out_bits = 8

    if is_rgb:
        channels = []
        for ch in range(3):
            work = data[ch].astype(np.float64)
            work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)
            channels.append(work)

        if out_bits == 8:
            vmin = min(ch.min() for ch in channels)
            vmax = max(ch.max() for ch in channels)
            result_channels = []
            for ch in channels:
                if vmax > vmin:
                    scaled = (ch - vmin) / (vmax - vmin) * 255.0
                else:
                    scaled = np.zeros_like(ch)
                result_channels.append(np.rint(np.clip(scaled, 0, 255)).astype(np.uint8))
            result = np.stack(result_channels, axis=-1)
        elif out_bits == 16:
            result = np.stack(
                [np.rint(np.clip(ch, 0, 65535)).astype(np.uint16) for ch in channels],
                axis=-1)
        else:
            vmin = min(ch.min() for ch in channels)
            vmax = max(ch.max() for ch in channels)
            result_channels = []
            for ch in channels:
                if vmax > vmin:
                    ch = (ch - vmin) / (vmax - vmin)
                else:
                    ch = np.zeros_like(ch)
                result_channels.append(ch.astype(np.float32))
            result = np.stack(result_channels, axis=-1)

        if out_bits == 16 and out_fmt == 'tiff':
            _save_tiff_rgb16(result, outfile)
            return
        if out_bits == 32 and out_fmt == 'tiff':
            _save_tiff_rgb32f(result, outfile)
            return

        img = Image.fromarray(result)
    else:
        work = data.astype(np.float64)
        work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)
        result = _scale_channel(work, out_bits)
        img = Image.fromarray(result)

    _save_image(img, outfile, out_fmt, jpeg_quality)


# =========================================================================
# Main
# =========================================================================

def main():
    config = parse_args(sys.argv)

    try:
        io_pairs = build_io_pairs(config['input_spec'], config['output_spec'])
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)

    if total == 1:
        infile, outfile = io_pairs[0]
        try:
            convert_file(infile, outfile, config)
            print(f"  {os.path.basename(outfile)}")
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")
    else:
        from concurrent.futures import ThreadPoolExecutor, as_completed
        max_workers = max(1, (os.cpu_count() or 2) - 1)
        done = 0
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for infile, outfile in io_pairs:
                fut = executor.submit(convert_file, infile, outfile, config)
                futures[fut] = infile
            for fut in as_completed(futures):
                infile = futures[fut]
                done += 1
                try:
                    fut.result()
                    sys.stdout.write(f"\r{done}/{total}")
                    sys.stdout.flush()
                except Exception as e:
                    sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    print(f"\nDone.")


if __name__ == "__main__":
    main()
