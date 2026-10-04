#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
animate - build a video animation from a sequence of FITS frames.

Collects an input sequence, sorts it by acquisition time (DATE-OBS), optionally
crops, then stretches every frame identically (percentile black/white points +
a Midtones Transfer Function that puts each frame's median at a target) and writes
a video. Mono frames -> grayscale video; RGB frames are background-neutralized so
the colour does not flicker. Source files are only read, never modified.

Two typical jobs:
  * Publication - smooth animations of objects/comets from aligned frames
    (--quality youtube for upload, or lossless for archival/editing).
  * Blinking / transient review - spotting moving or changing objects. The
    lossless default matters here: codec artifacts must not masquerade as, or
    hide, real changes. Use a low --fps (2-4) and --boomerang for classic blink.

Video is written with PyAV (the 'av' package): ffmpeg libraries ship inside the
wheel, so no external ffmpeg install is needed. If 'av' is missing, animate offers
to pip-install it (interactive runs only).

The per-frame stretch mirrors mtf.py's median-target mode:
  clip  = rescale [percentile(black), percentile(white)] -> [0, 1]
  m     = midtones so that MTF(median(clip), m) = --median target
  out   = MTF(clip, m)                       (PixInsight MTF, clipped to [0, 1])
then quantized to 8-bit for the video.
"""

import sys
import os
import math

import numpy as np

# Import shared utilities (one level up from Animate/).
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEF_FPS = 30.0
DEF_MEDIAN = 0.2          # target median for the MTF stretch
DEF_BLACK = 0.1           # default black-point percentile (%) for --black
DEF_AUTOBLACK_K = 5.0     # default MAD factor N for --autoblack [N] (median - N*MAD)
DEF_WHITE = 99.9          # white-point percentile (%)
DEF_QUALITY = "lossless"  # lossless | youtube | uncompressed
DEF_SORT = "date"         # date | name

QUALITIES = ("lossless", "youtube", "uncompressed")
SORTS = ("date", "name")
LABEL_TZS = ("utc", "local")


# 8x8 pixel font: public-domain font8x8_basic (Daniel Hepper), ASCII 32..126.
# Each glyph is 8 bytes (rows top->bottom); within a row byte, bit c (LSB first)
# is column c from the left. Crisp, integer-scaled, no anti-aliasing, no Pillow.
_FONT8X8 = bytes.fromhex(
    "0000000000000000183c3c1818001800363600000000000036367f367f3636000c3e031e301f0c00006333180c6663001c361c6e3b336e000606030000000000180c0606060c1800060c1818180c060000663cff3c660000000c0c3f0c0c00"
    "0000000000000c0c060000003f0000000000000000000c0c006030180c060301003e63737b6f673e000c0e0c0c0c0c3f001e33301c06333f001e33301c30331e00383c36337f3078003f031f3030331e001c06031f33331e003f3330180c0c"
    "0c001e33331e33331e001e33333e30180e00000c0c00000c0c00000c0c00000c0c06180c0603060c180000003f00003f0000060c1830180c06001e3330180c000c003e637b7b7b031e000c1e33333f3333003f66663e66663f003c66030303"
    "663c001f36666666361f007f46161e16467f007f46161e16060f003c66030373667c003333333f333333001e0c0c0c0c0c1e007830303033331e006766361e366667000f06060646667f0063777f7f6b63630063676f7b736363001c366363"
    "63361c003f66663e06060f001e3333333b1e38003f66663e366667001e33070e38331e003f2d0c0c0c0c1e003333333333333f0033333333331e0c006363636b7f7763006363361c1c3663003333331e0c0c1e007f6331184c667f001e0606"
    "0606061e0003060c18306040001e18181818181e00081c36630000000000000000000000ff0c0c18000000000000001e303e336e000706063e66663b0000001e3303331e003830303e33336e0000001e333f031e001c36060f06060f000000"
    "6e33333e301f0706366e666667000c000e0c0c0c1e00300030303033331e070666361e3667000e0c0c0c0c0c1e000000337f7f6b630000001f333333330000001e3333331e0000003b66663e060f00006e33333e307800003b6e66060f0000"
    "003e031e301f00080c3e0c0c2c18000000333333336e0000003333331e0c000000636b7f7f3600000063361c36630000003333333e301f00003f190c263f00380c0c070c0c38001818180018181800070c0c380c0c07006e3b000000000000")
_GLYPH_W = 8
_GLYPH_H = 8


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "animate - build a video animation from a sequence of FITS frames.\n"
        "\n"
        "Usage:\n"
        "  animate.py input_spec output [options]\n"
        "\n"
        "  input_spec   directory, wildcard (*.fit), numbered (img0001.fit) or @list.txt\n"
        "  output       output video path (.mkv for lossless/uncompressed, .mp4 for youtube)\n"
        "\n"
        "Stretch (applied per frame, identically):\n"
        "  --median T       target median for the MTF, in (0,1). Default %.2g\n"
        "  --black P        black-point percentile in %% (over non-zero pixels). Default %g\n"
        "  --autoblack [N]  alternative black point: median - N*MAD of the non-zero pixels,\n"
        "                   signed (robust, like mtf --autoblack); N default 5. Mutually\n"
        "                   exclusive with --black.\n"
        "  --white P        white-point percentile in %%. Default %.4g\n"
        "\n"
        "Video:\n"
        "  --fps N          frames per second. Default %g (use 2-4 for blinking)\n"
        "  --quality Q      lossless (default, FFV1/.mkv) | youtube (H.264 crf16/.mp4)\n"
        "                   | uncompressed (rawvideo). Lossless is artifact-free, best\n"
        "                   for transient review; youtube is for upload.\n"
        "  --loop N         repeat the whole sequence N times. Default 1\n"
        "  --boomerang      append the reversed sequence (ping-pong, seamless loop)\n"
        "  --label [TZ]     timezone for the burned label: utc (default, suffix 'Z') or\n"
        "                   local (machine timezone). The label (DATE-OBS timestamp + the\n"
        "                   source filename) is burned into every output by default.\n"
        "  --nostamp        do not burn any label (no timestamp, no filename).\n"
        "  (output .ser)    an output path ending in .ser writes one uncompressed 16-bit\n"
        "                   SER file (mono/RGB) with a per-frame UTC timestamp trailer,\n"
        "                   playable in SER Player / PIPP / AutoStakkert. No PyAV needed.\n"
        "  --png            write a numbered PNG series (<output>_0001.png, ...) with the\n"
        "                   DATE-OBS timestamp burned in, instead of a video. Mono -> 8-bit\n"
        "                   grayscale; colour -> 8-bit RGBA (opaque). Needs Pillow.\n"
        "\n"
        "Ordering / cropping:\n"
        "  --sort MODE      date (default, by DATE-OBS) | name (natural filename order).\n"
        "                   Missing DATE-OBS -> warn and fall back to name order.\n"
        "  --center W H     crop a WxH region at the image centre (pixels). (With --width/\n"
        "                   --height below, --center X Y is instead the centre point.)\n"
        "  --width W --height H [--center X Y]   crop WxH about a centre point (default:\n"
        "                   image centre), like crop.py.\n"
        "  --top N --bottom N --left N --right N  crop by trimming margins (pixels).\n"
        "\n"
        "Other:\n"
        "  --threads N      worker threads for frame processing. Default CPU cores - 1\n"
        "  --probe FILE     diagnostic: print what a produced video contains (frame\n"
        "                   count, timestamps, brightness) and exit\n"
        "  -h, --help       show this help\n"
        "\n"
        "Source files are only read, never modified. RGB frames get per-frame background\n"
        "neutralization; mono frames become a grayscale video. Needs the 'av' package\n"
        "(offered for pip-install on first run if missing).\n"
        % (DEF_MEDIAN, DEF_BLACK, DEF_WHITE, DEF_FPS)
    )
    sys.exit(1)


def _need_value(args, i, flag):
    if i + 1 >= len(args):
        sys.stderr.write("Error: %s requires a value.\n" % flag)
        sys.exit(1)
    return args[i + 1]


def _as_float(flag, s):
    try:
        return float(s)
    except ValueError:
        sys.stderr.write("Error: %s requires a numeric value, got '%s'.\n" % (flag, s))
        sys.exit(1)


def _as_int(flag, s):
    try:
        return int(s)
    except ValueError:
        sys.stderr.write("Error: %s requires an integer, got '%s'.\n" % (flag, s))
        sys.exit(1)


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def parse_args(argv):
    args = argv[1:]
    cfg = {
        "fps": DEF_FPS, "median": DEF_MEDIAN, "black": None, "autoblack": None, "white": DEF_WHITE,
        "quality": DEF_QUALITY, "sort": DEF_SORT, "loop": 1, "boomerang": False,
        "label": "utc", "nostamp": False, "threads": None, "probe": False, "png": False,
        "width": None, "height": None, "center": None,
        "top": 0, "bottom": 0, "left": 0, "right": 0,
    }
    positional = []
    i = 0
    while i < len(args):
        a = args[i]
        if a in ("-h", "--help"):
            usage()
        elif a == "--median" or a == "-m":
            cfg["median"] = _as_float(a, _need_value(args, i, a)); i += 2
        elif a == "--black":
            cfg["black"] = _as_float(a, _need_value(args, i, a)); i += 2
        elif a == "--autoblack":
            if i + 1 < len(args) and _is_number(args[i + 1]) and not args[i + 1].startswith("-"):
                cfg["autoblack"] = float(args[i + 1]); i += 2
            else:
                cfg["autoblack"] = DEF_AUTOBLACK_K; i += 1
        elif a == "--white":
            cfg["white"] = _as_float(a, _need_value(args, i, a)); i += 2
        elif a == "--fps":
            cfg["fps"] = _as_float(a, _need_value(args, i, a)); i += 2
        elif a == "--quality":
            cfg["quality"] = _need_value(args, i, a).lower(); i += 2
        elif a == "--sort":
            cfg["sort"] = _need_value(args, i, a).lower(); i += 2
        elif a == "--loop":
            cfg["loop"] = _as_int(a, _need_value(args, i, a)); i += 2
        elif a == "--boomerang":
            cfg["boomerang"] = True; i += 1
        elif a == "--label":
            # Optional value: 'utc' or 'local'; a bare --label means utc.
            if i + 1 < len(args) and args[i + 1].lower() in LABEL_TZS:
                cfg["label"] = args[i + 1].lower(); i += 2
            else:
                cfg["label"] = "utc"; i += 1
        elif a == "--nostamp":
            cfg["nostamp"] = True; i += 1
        elif a == "--threads":
            cfg["threads"] = _as_int(a, _need_value(args, i, a)); i += 2
        elif a == "--probe":
            cfg["probe"] = True; i += 1
        elif a == "--png":
            cfg["png"] = True; i += 1
        elif a == "--width":
            cfg["width"] = _as_int(a, _need_value(args, i, a)); i += 2
        elif a == "--height":
            cfg["height"] = _as_int(a, _need_value(args, i, a)); i += 2
        elif a == "--center":
            if i + 2 >= len(args):
                sys.stderr.write("Error: --center requires two values (W H, or X Y with --width/--height).\n"); sys.exit(1)
            cfg["center"] = (_as_int("--center", args[i + 1]), _as_int("--center", args[i + 2])); i += 3
        elif a == "--top":
            cfg["top"] = _as_int(a, _need_value(args, i, a)); i += 2
        elif a == "--bottom":
            cfg["bottom"] = _as_int(a, _need_value(args, i, a)); i += 2
        elif a == "--left":
            cfg["left"] = _as_int(a, _need_value(args, i, a)); i += 2
        elif a == "--right":
            cfg["right"] = _as_int(a, _need_value(args, i, a)); i += 2
        elif a.startswith("-") and a != "-":
            sys.stderr.write("Error: unknown option: %s\n" % a)
            usage()
        else:
            positional.append(a); i += 1

    if cfg["probe"]:
        if len(positional) != 1:
            sys.stderr.write("Error: --probe takes exactly one file to inspect.\n"); sys.exit(1)
        return positional[0], None, cfg
    if len(positional) != 2:
        usage()
    # Validate
    if not (0.0 < cfg["median"] < 1.0):
        sys.stderr.write("Error: --median must be in (0, 1).\n"); sys.exit(1)
    if cfg["black"] is not None and cfg["autoblack"] is not None:
        sys.stderr.write("Error: --black and --autoblack are mutually exclusive.\n"); sys.exit(1)
    if cfg["black"] is not None and not (0.0 <= cfg["black"] < cfg["white"]):
        sys.stderr.write("Error: --black must be a percentile in [0, --white).\n"); sys.exit(1)
    if cfg["autoblack"] is not None and (not math.isfinite(cfg["autoblack"]) or cfg["autoblack"] < 0.0):
        sys.stderr.write("Error: --autoblack N (MAD factor) must be a finite number >= 0.\n"); sys.exit(1)
    if not (0.0 < cfg["white"] <= 100.0):
        sys.stderr.write("Error: --white must be a percentile in (0, 100].\n"); sys.exit(1)
    if not math.isfinite(cfg["fps"]) or cfg["fps"] <= 0.0:
        sys.stderr.write("Error: --fps must be a finite positive number.\n"); sys.exit(1)
    if cfg["quality"] not in QUALITIES:
        sys.stderr.write("Error: --quality must be one of %s.\n" % ", ".join(QUALITIES)); sys.exit(1)
    if cfg["sort"] not in SORTS:
        sys.stderr.write("Error: --sort must be one of %s.\n" % ", ".join(SORTS)); sys.exit(1)
    if cfg["loop"] < 1:
        sys.stderr.write("Error: --loop must be >= 1.\n"); sys.exit(1)
    if cfg["threads"] is not None and cfg["threads"] < 1:
        sys.stderr.write("Error: --threads must be >= 1.\n"); sys.exit(1)
    uses_wh = cfg["width"] is not None or cfg["height"] is not None or cfg["center"] is not None
    uses_margin = any(cfg[k] for k in ("top", "bottom", "left", "right"))
    if uses_wh and uses_margin:
        sys.stderr.write("Error: mix of --width/--height/--center and --top/--bottom/--left/--right.\n")
        sys.exit(1)
    if (cfg["width"] is None) != (cfg["height"] is None):
        sys.stderr.write("Error: --width and --height must be given together.\n"); sys.exit(1)
    if any(cfg[k] < 0 for k in ("top", "bottom", "left", "right")):
        sys.stderr.write("Error: crop margins must be non-negative.\n"); sys.exit(1)

    return positional[0], positional[1], cfg


# ---------------------------------------------------------------------------
# Dependency: PyAV (offer to install on first interactive run)
# ---------------------------------------------------------------------------

def ensure_av():
    """Import PyAV ('av'); if missing, offer to pip-install it (interactive only)."""
    try:
        import av
        return av
    except ImportError:
        pass
    cmd = "%s -m pip install av" % sys.executable
    sys.stderr.write("animate: the 'av' package (PyAV) is required for video output "
                     "but is not installed.\n")
    interactive = bool(getattr(sys, "stdin", None)) and sys.stdin.isatty()
    if not interactive:
        sys.stderr.write("  Install it with:\n    %s\n" % cmd)
        sys.exit(1)
    try:
        ans = input("Install it now with pip? [y/N] ").strip().lower()
    except EOFError:
        ans = ""
    if ans not in ("y", "yes"):
        sys.stderr.write("Aborted. Install manually:\n    %s\n" % cmd)
        sys.exit(1)
    import subprocess
    if subprocess.run([sys.executable, "-m", "pip", "install", "av"]).returncode != 0:
        sys.stderr.write("Error: install failed. Install manually:\n    %s\n" % cmd)
        sys.exit(1)
    try:
        import av
        return av
    except ImportError:
        sys.stderr.write("Error: 'av' still not importable after install.\n")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Collect + sort
# ---------------------------------------------------------------------------

def _natural_key(path):
    """Filename sort key that orders img2 before img10 (numeric-aware)."""
    name = os.path.basename(path)
    import re
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", name)]


def _read_dateobs(path):
    """Return a sortable key from DATE-OBS (plus DATE-OBS string), or (None, None)."""
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdul:
        hdr = hdul[0].header
        val = hdr.get("DATE-OBS") or hdr.get("DATE_OBS")
    if not val:
        return None, None
    return _parse_time(str(val)), str(val)


def _parse_time(s):
    """Parse an ISO-ish DATE-OBS into a datetime (UTC-naive), or None."""
    from datetime import datetime
    t = s.strip().rstrip("Zz")
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S",
                "%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M"):
        try:
            return datetime.strptime(t, fmt)
        except ValueError:
            continue
    try:
        return datetime.fromisoformat(t)
    except ValueError:
        return None


def collect_and_sort(input_spec, sort_mode):
    """Expand the input spec and return a list of (path, dateobs_str) in order.

    sort=date orders by DATE-OBS; any missing/unparseable DATE-OBS triggers a loud
    warning and a fall back to natural filename order (sort=name forces that)."""
    if os.path.isdir(input_spec):
        import glob
        files = sorted(glob.glob(os.path.join(input_spec, "*.fit")) +
                       glob.glob(os.path.join(input_spec, "*.fits")))
        files = [os.path.abspath(f) for f in files]
        if not files:
            raise FileNotFoundError("no FITS files in directory: %s" % input_spec)
    else:
        files = batch_utils.expand_input_spec(input_spec)

    stamps = {}
    if sort_mode == "date":
        missing = []
        for f in files:
            key, raw = _read_dateobs(f)
            stamps[f] = raw
            if key is None:
                missing.append(f)
        if missing:
            sys.stderr.write("Warning: %d of %d files have no usable DATE-OBS; "
                             "falling back to filename order.\n" % (len(missing), len(files)))
            sort_mode = "name"
    if sort_mode == "date":
        files.sort(key=lambda f: _parse_time(stamps[f]))
    else:
        files.sort(key=_natural_key)
        if not stamps:
            # still read DATE-OBS for optional labels, best-effort
            for f in files:
                _, stamps[f] = _read_dateobs(f)
    return [(f, stamps.get(f)) for f in files]


def build_order(n, loop, boomerang):
    """Index order into the frame list, applying --loop and --boomerang."""
    base = list(range(n))
    if boomerang and n > 2:
        base = base + list(range(n - 2, 0, -1))
    return base * max(1, loop)


# ---------------------------------------------------------------------------
# Crop
# ---------------------------------------------------------------------------

def compute_crop(h, w, cfg):
    """Return (y0, y1, x0, x1) for the crop, or None for no crop. Validated here
    so a bad crop fails loudly before any frame is encoded. With --width/--height,
    --center X Y is the centre point; without them, --center W H is the central
    crop size."""
    if cfg["width"] is not None:
        cw, ch = cfg["width"], cfg["height"]
        if cw <= 0 or ch <= 0:
            raise ValueError("--width/--height must be positive")
        if cfg["center"] is not None:
            cx, cy = cfg["center"]
        else:
            cx, cy = w // 2, h // 2
        x0 = cx - cw // 2
        y0 = cy - ch // 2
        x1, y1 = x0 + cw, y0 + ch
        if x0 < 0 or y0 < 0 or x1 > w or y1 > h:
            raise ValueError("crop %dx%d about (%d,%d) is outside the %dx%d frame"
                             % (cw, ch, cx, cy, w, h))
        return y0, y1, x0, x1
    if cfg["center"] is not None:
        # No explicit --width/--height: --center W H is the SIZE of a central crop.
        cw, ch = cfg["center"]
        if cw <= 0 or ch <= 0:
            raise ValueError("--center W H must be positive")
        if cw > w or ch > h:
            raise ValueError("--center %dx%d is larger than the %dx%d frame" % (cw, ch, w, h))
        x0 = (w - cw) // 2
        y0 = (h - ch) // 2
        return y0, y0 + ch, x0, x0 + cw
    if any(cfg[k] for k in ("top", "bottom", "left", "right")):
        y0, y1 = cfg["top"], h - cfg["bottom"]
        x0, x1 = cfg["left"], w - cfg["right"]
        if y1 <= y0 or x1 <= x0:
            raise ValueError("crop margins remove the whole frame")
        return y0, y1, x0, x1
    return None


# ---------------------------------------------------------------------------
# Stretch (per frame) -> 8-bit
# ---------------------------------------------------------------------------

def _black_point(values, k):
    """Robust signed black level: median - k*MAD of the values (k = MAD factor,
    'how many MADs below the background'). Signed - may be negative. Same family
    as mtf.py --autoblack (median - 5*MAD)."""
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med)))
    return med - k * mad


def _black_level(values, cfg):
    """Black point for the stretch over the (non-zero) values: --autoblack N gives
    median - N*MAD (robust, signed); otherwise the --black percentile (default)."""
    if cfg["autoblack"] is not None:
        return _black_point(values, cfg["autoblack"])
    pct = cfg["black"] if cfg["black"] is not None else DEF_BLACK
    return float(np.percentile(values, pct))


def solve_midtones(median, target):
    """PixInsight MTF midtones m so that MTF(median, m) = target (both in (0,1))."""
    d = min(max(float(median), 1e-6), 1.0 - 1e-6)
    t = min(max(float(target), 1e-6), 1.0 - 1e-6)
    denom = d + t - 2.0 * d * t
    return d * (1.0 - t) / denom if denom > 0.0 else 0.5


def apply_mtf(x, m):
    """PixInsight MTF: (1-m)*x / (m + x*(1-2m)); m=0.5 identity. x in [0,1]."""
    if abs(m - 0.5) < 1e-9:
        return np.clip(x, 0.0, 1.0)
    denom = m + x * (1.0 - 2.0 * m)
    out = (1.0 - m) * x / np.where(denom != 0.0, denom, 1.0)
    return np.clip(out, 0.0, 1.0)


def _balance_background(rgb):
    """Neutralize an RGB frame's background: shift each channel so its background
    level (30th percentile, background-dominated) matches the channels' mean. In
    place on a float64 (C,H,W) array."""
    bg = []
    for c in range(3):
        nz = rgb[c][rgb[c] != 0.0]                          # ignore the zero alignment border
        bg.append(float(np.percentile(nz if nz.size else rgb[c], 30.0)))
    target = sum(bg) / 3.0
    for c in range(3):
        rgb[c] += (target - bg[c])
    return rgb


def stretch_to_u8(data, cfg, crop, maxval=255):
    """Crop, stretch (percentile/MAD black -> MTF median->target) and quantize one
    frame. maxval 255 -> uint8, 65535 -> uint16. Returns (H,W) mono or (H,W,3) RGB
    with even dimensions (H.264/yuv420p needs them). Colour channels share one
    midtones value (from luminance) so hue is preserved."""
    arr = np.asarray(data, dtype=np.float64)
    out_dtype = np.uint8 if maxval <= 255 else np.uint16

    # Determine layout: mono (H,W) or RGB as (3,H,W) / (H,W,3).
    if arr.ndim == 2:
        channels, color = [arr], False
    elif arr.ndim == 3 and arr.shape[0] == 3:
        channels, color = [arr[0], arr[1], arr[2]], True
    elif arr.ndim == 3 and arr.shape[2] == 3:
        channels, color = [arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]], True
    else:
        raise ValueError("unsupported frame shape %s (need 2D or 3-channel RGB)" % (arr.shape,))

    if crop is not None:
        y0, y1, x0, x1 = crop
        channels = [c[y0:y1, x0:x1] for c in channels]

    if color:
        stack = np.stack(channels, axis=0)
        border = ~np.any(stack != 0.0, axis=0)             # all-channel-zero pixels = alignment border
        _balance_background(stack)
        lum = stack.mean(axis=0)
        valid = lum[~border] if (~border).any() else lum.ravel()
        lo = _black_level(valid, cfg)                      # --autoblack: median-N*MAD; else --black percentile (zeros ignored)
        hi = float(np.percentile(valid, cfg["white"]))
        span = (hi - lo) if hi > lo else 1.0
        norm = np.clip((stack - lo) / span, 0.0, 1.0)      # (3,H,W)
        m = solve_midtones(float(np.median(np.clip((valid - lo) / span, 0.0, 1.0))), cfg["median"])
        out = apply_mtf(norm, m)
        out = np.nan_to_num(out, nan=0.0, posinf=1.0, neginf=0.0)
        u8 = np.rint(out * maxval).astype(out_dtype)       # (3,H,W)
        u8 = np.ascontiguousarray(np.transpose(u8, (1, 2, 0)))  # -> (H,W,3)
    else:
        g = channels[0]
        src = g[g != 0.0]                                   # ignore the zero alignment border
        if src.size == 0:
            src = g.ravel()
        lo = _black_level(src, cfg)                         # --autoblack: median-N*MAD; else --black percentile
        hi = float(np.percentile(src, cfg["white"]))
        span = (hi - lo) if hi > lo else 1.0
        norm = np.clip((g - lo) / span, 0.0, 1.0)
        m = solve_midtones(float(np.median(np.clip((src - lo) / span, 0.0, 1.0))), cfg["median"])
        out = apply_mtf(norm, m)
        out = np.nan_to_num(out, nan=0.0, posinf=1.0, neginf=0.0)
        u8 = np.rint(out * maxval).astype(out_dtype)       # (H,W)

    # Enforce even dimensions (drop a trailing row/column if odd).
    h, w = u8.shape[0], u8.shape[1]
    if h % 2:
        u8 = u8[:h - 1]
    if w % 2:
        u8 = u8[:, :w - 1]
    return np.ascontiguousarray(u8)


# ---------------------------------------------------------------------------
# Label (pixel font)
# ---------------------------------------------------------------------------

def make_label_text(dateobs, tz, filename=None):
    """Build the burned-in label '<timestamp>  <source filename>'. The timestamp is
    DATE-OBS in the chosen tz (UTC with a 'Z', or machine-local); if DATE-OBS is
    missing the label is just the filename. The source filename (basename) lets a
    peculiar artifact be traced back to its exact frame."""
    ts = ""
    if dateobs:
        dt = _parse_time(str(dateobs))
        if dt is not None:
            if tz == "local":
                from datetime import timezone
                ts = dt.replace(tzinfo=timezone.utc).astimezone().strftime("%Y-%m-%d %H:%M:%S")
            else:
                ts = dt.strftime("%Y-%m-%d %H:%M:%S") + "Z"
    name = os.path.basename(filename) if filename else ""
    if ts and name:
        return ts + "  " + name
    return ts or name or None


def draw_label(u8, text, scale, white=255):
    """Draw 'text' with the 8x8 pixel font near the bottom-left corner, on a dark
    backdrop for contrast. Pixel-perfect (integer scale, no anti-aliasing). 'white'
    is the lit-pixel value (255 for 8-bit, 65535 for 16-bit). Text wider than the
    frame is truncated, not dropped."""
    if not text:
        return u8
    pad = scale
    adv = (_GLYPH_W + 1) * scale                        # per-character advance
    H, W = u8.shape[0], u8.shape[1]
    text_h = _GLYPH_H * scale
    y0 = H - pad - text_h
    if y0 < 0:                                          # frame too short: skip
        return u8
    maxchars = max(0, (W - 2 * pad) // adv)             # truncate to fit the width
    if maxchars == 0:
        return u8
    text = text[:maxchars]
    x0 = pad
    text_w = len(text) * adv
    by0, by1 = max(0, y0 - pad), min(H, y0 + text_h + pad)
    bx0, bx1 = max(0, x0 - pad), min(W, x0 + text_w + pad)
    if u8.ndim == 2:
        u8[by0:by1, bx0:bx1] = 0                        # dark backdrop
    else:
        u8[by0:by1, bx0:bx1, :] = 0
    cx = x0
    for ch in text:
        o = ord(ch)
        if 32 <= o <= 126:
            base = (o - 32) * 8
            for r in range(_GLYPH_H):
                row = _FONT8X8[base + r]
                for c in range(_GLYPH_W):
                    if (row >> c) & 1:                  # LSB = leftmost column
                        yy, xx = y0 + r * scale, cx + c * scale
                        if u8.ndim == 2:
                            u8[yy:yy + scale, xx:xx + scale] = white
                        else:
                            u8[yy:yy + scale, xx:xx + scale, :] = white
        cx += adv
    return u8


def label_scale(width):
    """Integer font scale chosen from the frame width (readable but small)."""
    return max(1, int(round(width / 1200.0)))


# ---------------------------------------------------------------------------
# Encoding (PyAV)
# ---------------------------------------------------------------------------

def run_png(entries, output, cfg):
    """--png output: write a numbered PNG series instead of a video - one PNG per
    frame in sorted order, with the DATE-OBS timestamp burned in. Mono -> 8-bit
    grayscale ('L'); colour -> 8-bit RGBA ('RGBA', opaque, alpha=255 => 32 bit/px).
    Uses Pillow (not opencv). Frames are processed in parallel; each is saved on its
    own, so any number of frames works (only disk limits the result)."""
    try:
        from PIL import Image
    except ImportError:
        sys.stderr.write("animate: --png needs Pillow. Install it with:\n"
                         "    %s -m pip install pillow\n" % sys.executable)
        sys.exit(1)
    from astropy.io import fits

    # Crop geometry + label scale from the first frame.
    with fits.open(entries[0][0], memmap=False) as hdul:
        first = np.asarray(hdul[0].data)
    if first.ndim == 3 and first.shape[0] == 3:
        fh, fw = first.shape[1], first.shape[2]
    elif first.ndim == 3 and first.shape[2] == 3:
        fh, fw = first.shape[0], first.shape[1]
    elif first.ndim == 2:
        fh, fw = first.shape
    else:
        sys.stderr.write("Error: unsupported frame shape %s.\n" % (first.shape,)); sys.exit(1)
    try:
        crop = compute_crop(fh, fw, cfg)
    except ValueError as e:
        sys.stderr.write("Error: %s\n" % e); sys.exit(1)
    lscale = label_scale(crop[3] - crop[2] if crop else fw)

    stem = os.path.splitext(output)[0]
    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    total = len(entries)
    digits = max(4, len(str(total)))          # leading zeros, sized to the frame count
    threads = cfg["threads"] if cfg["threads"] else max(1, (os.cpu_count() or 2) - 1)

    def work(item):
        idx, entry = item
        u8 = process_one(entry, cfg, crop, lscale)
        if u8.ndim == 2:
            img = Image.fromarray(u8, mode="L")                 # 8-bit grayscale
        else:
            alpha = np.full(u8.shape[:2], 255, dtype=np.uint8)  # opaque
            rgba = np.ascontiguousarray(np.dstack([u8, alpha]))
            img = Image.fromarray(rgba, mode="RGBA")            # 8-bit RGBA (32 bit/px)
        img.save("%s_%0*d.png" % (stem, digits, idx + 1))

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=threads) as ex:
        for k, _ in enumerate(ex.map(work, enumerate(entries)), start=1):
            _progress(k, total)
    sys.stderr.write("\nanimate: wrote %d PNG frame(s): %s_%0*d.png .. %s_%0*d.png\n"
                     % (total, stem, digits, 1, stem, digits, total))


def run_ser(entries, output, cfg):
    """SER output (output path ends in .ser): write all frames into one uncompressed
    16-bit .ser via lib/ser_writer, one frame at a time (any count, disk-limited).
    Each frame carries the burned timestamp+filename label; the SER trailer also
    records per-frame UTC time from DATE-OBS."""
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
    from ser_writer import SerWriter
    from astropy.io import fits

    with fits.open(entries[0][0], memmap=False) as hdul:
        f0 = np.asarray(hdul[0].data)
    if f0.ndim == 3 and f0.shape[0] == 3:
        fh, fw = f0.shape[1], f0.shape[2]
    elif f0.ndim == 3 and f0.shape[2] == 3:
        fh, fw = f0.shape[0], f0.shape[1]
    elif f0.ndim == 2:
        fh, fw = f0.shape
    else:
        sys.stderr.write("Error: unsupported frame shape %s.\n" % (f0.shape,)); sys.exit(1)
    try:
        crop = compute_crop(fh, fw, cfg)
    except ValueError as e:
        sys.stderr.write("Error: %s\n" % e); sys.exit(1)
    lscale = label_scale(crop[3] - crop[2] if crop else fw)

    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    def ts_of(entry):
        return _parse_time(str(entry[1])) if entry[1] else None

    first = process_one(entries[0], cfg, crop, lscale, maxval=65535)   # 16-bit
    color = first.ndim == 3
    vh, vw = first.shape[0], first.shape[1]
    total = len(entries)
    threads = cfg["threads"] if cfg["threads"] else max(1, (os.cpu_count() or 2) - 1)

    writer = SerWriter(output, vw, vh, color=color, depth=16)
    try:
        writer.add_frame(first, timestamp=ts_of(entries[0]))
        _progress(1, total)
        from concurrent.futures import ThreadPoolExecutor
        window = max(2, threads * 2)
        with ThreadPoolExecutor(max_workers=threads) as ex:
            inflight = {}
            nxt = 1
            while nxt < total and nxt <= window:
                inflight[nxt] = ex.submit(process_one, entries[nxt], cfg, crop, lscale, 65535); nxt += 1
            for pos in range(1, total):
                u = inflight.pop(pos).result()
                if u.shape != first.shape:
                    raise ValueError("frame %d shape %s != first frame %s "
                                     "(frames must share width, height and channel count; "
                                     "crop to a common size)" % (pos, u.shape, first.shape))
                writer.add_frame(u, timestamp=ts_of(entries[pos]))
                _progress(pos + 1, total)
                if nxt < total:
                    inflight[nxt] = ex.submit(process_one, entries[nxt], cfg, crop, lscale, 65535); nxt += 1
    finally:
        writer.close()
    sys.stderr.write("\nanimate: wrote %d frames -> %s (16-bit %s SER)\n"
                     % (total, output, "RGB" if color else "mono"))


def run_probe(av, path):
    """Diagnostic: report what a produced video actually contains - the packet
    count (frames written), their timestamps, and the first decoded frames'
    brightness (to see whether frames are distinct). Read-only."""
    try:
        probe = av.open(path)
    except Exception as exc:
        raise RuntimeError("cannot open '%s' as a media file (%s)" % (path, exc))
    with probe as c:
        if not c.streams.video:
            raise RuntimeError("'%s' has no video stream" % path)
        vs = c.streams.video[0]
        cc = vs.codec_context
        sys.stderr.write("probe: %s\n" % path)
        sys.stderr.write("  codec=%s pix_fmt=%s size=%dx%d avg_rate=%s time_base=%s\n"
                         % (cc.name, cc.pix_fmt, cc.width or 0, cc.height or 0,
                            vs.average_rate, vs.time_base))
    npkt = 0
    pkt_pts = []
    with av.open(path) as c:
        for p in c.demux(video=0):
            if p.size == 0:
                continue
            if npkt < 12:
                pkt_pts.append(p.pts)
            npkt += 1
    rows = []
    with av.open(path) as c:
        for i, fr in enumerate(c.decode(video=0)):
            rows.append((i, fr.pts, round(float(np.asarray(fr.to_ndarray()).mean()), 3)))
            if i >= 11:
                break
    sys.stderr.write("  packets (frames written): %d\n" % npkt)
    sys.stderr.write("  first packet pts: %s\n" % pkt_pts)
    sys.stderr.write("  first decoded frames (index, pts, mean-brightness):\n")
    for i, pts, m in rows:
        sys.stderr.write("    %3d  pts=%s  mean=%s\n" % (i, pts, m))


def open_encoder(av, output, quality, width, height, color, fps):
    """Open an output container + video stream for the chosen quality mode."""
    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    try:
        container = av.open(output, mode="w")
    except Exception as exc:
        raise RuntimeError(
            "cannot open '%s' for writing (%s). Use a known extension: .mp4 (youtube), "
            ".mkv (lossless/uncompressed), .ser (SER), or --png for a PNG series." % (output, exc))
    # PyAV's add_stream wants an int/Fraction rate, not a float (a float raises
    # "'float' object has no attribute 'numerator'").
    from fractions import Fraction
    rate = Fraction(fps).limit_denominator(1000000)
    if quality == "youtube":
        stream = container.add_stream("libx264", rate=rate)
        stream.pix_fmt = "yuv420p"
        stream.options = {"crf": "16", "preset": "slow"}
    elif quality == "uncompressed":
        stream = container.add_stream("rawvideo", rate=rate)
        stream.pix_fmt = "rgb24" if color else "gray"
    else:  # lossless
        stream = container.add_stream("ffv1", rate=rate)
        stream.pix_fmt = "gbrp" if color else "gray"
    stream.width = width
    stream.height = height
    return container, stream


def encode_frame(av, container, stream, u8):
    """Encode one 8-bit frame ((H,W) gray or (H,W,3) rgb), reformatting to the
    stream's pixel format as needed."""
    src_fmt = "rgb24" if u8.ndim == 3 else "gray"
    frame = av.VideoFrame.from_ndarray(u8, format=src_fmt)
    if stream.pix_fmt and stream.pix_fmt != src_fmt:
        frame = frame.reformat(format=stream.pix_fmt)
    for packet in stream.encode(frame):
        container.mux(packet)


# ---------------------------------------------------------------------------
# Processing + progress
# ---------------------------------------------------------------------------

def _progress(done, total):
    pct = int(round(100.0 * done / total)) if total else 100
    sys.stderr.write("\ranimate: %d/%d frames (%3d%%)" % (done, total, pct))
    sys.stderr.flush()


def process_one(entry, cfg, crop, label_scale_px, maxval=255):
    """Load one FITS, stretch and quantize (maxval 255 -> 8-bit, 65535 -> 16-bit),
    draw the label (timestamp + source filename). Returns the frame array."""
    from astropy.io import fits
    path, dateobs = entry
    with fits.open(path, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("no image data in %s" % path)
        data = hdul[0].data
        u = stretch_to_u8(data, cfg, crop, maxval)
    if not cfg["nostamp"]:
        draw_label(u, make_label_text(dateobs, cfg["label"], path), label_scale_px, white=maxval)
    return u


def main():
    input_spec, output, cfg = parse_args(sys.argv)

    if cfg["probe"]:
        run_probe(ensure_av(), input_spec)
        return

    try:
        entries = collect_and_sort(input_spec, cfg["sort"])
    except Exception as e:
        sys.stderr.write("Error: %s\n" % e); sys.exit(1)
    if not entries:
        sys.stderr.write("Error: no input frames.\n"); sys.exit(1)
    sys.stderr.write("animate: %d frames; sort=%s quality=%s fps=%g\n"
                     % (len(entries), cfg["sort"], cfg["quality"], cfg["fps"]))

    if cfg["png"]:
        run_png(entries, output, cfg)
        return
    if output.lower().endswith(".ser"):
        run_ser(entries, output, cfg)
        return

    # Video path. Give the output a container extension if it has none, so PyAV can
    # pick a muxer (lossless/uncompressed -> .mkv, youtube -> .mp4).
    if not os.path.splitext(output)[1]:
        output += ".mp4" if cfg["quality"] == "youtube" else ".mkv"
        sys.stderr.write("animate: no output extension given; writing %s\n" % output)

    av = ensure_av()

    # Crop geometry + reference size come from the FIRST frame; the encoder is then
    # opened and the first frame encoded, so a size mismatch fails before encoding.
    from astropy.io import fits
    with fits.open(entries[0][0], memmap=False) as hdul:
        first = np.asarray(hdul[0].data)
    if first.ndim == 3 and first.shape[0] == 3:
        fh, fw = first.shape[1], first.shape[2]
    elif first.ndim == 3 and first.shape[2] == 3:
        fh, fw = first.shape[0], first.shape[1]
    elif first.ndim == 2:
        fh, fw = first.shape
    else:
        sys.stderr.write("Error: unsupported frame shape %s.\n" % (first.shape,)); sys.exit(1)
    try:
        crop = compute_crop(fh, fw, cfg)
    except ValueError as e:
        sys.stderr.write("Error: %s\n" % e); sys.exit(1)

    lscale = label_scale(crop[3] - crop[2] if crop else fw)

    first_u8 = process_one(entries[0], cfg, crop, lscale)
    color = first_u8.ndim == 3
    vh, vw = first_u8.shape[0], first_u8.shape[1]

    order = build_order(len(entries), cfg["loop"], cfg["boomerang"])
    threads = cfg["threads"] if cfg["threads"] else max(1, (os.cpu_count() or 2) - 1)

    container, stream = open_encoder(av, output, cfg["quality"], vw, vh, color, cfg["fps"])
    try:
        encode_frame(av, container, stream, first_u8)

        from concurrent.futures import ThreadPoolExecutor
        window = max(2, threads * 2)
        total = len(order)

        def frame_for(pos):
            idx = order[pos]
            if idx == 0:
                return first_u8 if pos == 0 else process_one(entries[0], cfg, crop, lscale)
            return process_one(entries[idx], cfg, crop, lscale)

        done = 1
        _progress(done, total)
        with ThreadPoolExecutor(max_workers=threads) as ex:
            inflight = {}
            submit_pos = 1
            while submit_pos < total and submit_pos <= window:
                inflight[submit_pos] = ex.submit(frame_for, submit_pos); submit_pos += 1
            for pos in range(1, total):
                u8 = inflight.pop(pos).result()
                if u8.shape != first_u8.shape:
                    raise ValueError("frame %d shape %s != first frame %s "
                                     "(frames must share width, height and channel count; "
                                     "crop to a common size)" % (order[pos], u8.shape, first_u8.shape))
                encode_frame(av, container, stream, u8)
                done += 1
                _progress(done, total)
                if submit_pos < total:
                    inflight[submit_pos] = ex.submit(frame_for, submit_pos); submit_pos += 1

        for packet in stream.encode():     # flush
            container.mux(packet)
    finally:
        container.close()

    sys.stderr.write("\nanimate: wrote %s (%dx%d, %d frames @ %g fps)\n"
                     % (output, vw, vh, len(order), cfg["fps"]))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nanimate: interrupted\n"); sys.exit(1)
    except (RuntimeError, ValueError, OSError) as exc:
        sys.stderr.write("animate: %s\n" % exc); sys.exit(1)
