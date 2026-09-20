#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
makemask - Build processing masks from FITS images.

Optional pipeline, applied in this order:
  1. colour -> greyscale via g = (R + 2*G + B) / 4  (only if input is 3-channel)
  2. black/white clip + linear stretch (--black / --white); each endpoint is
     given either as a PERCENTILE (a trailing '%', e.g. 90%) or as an ABSOLUTE
     brightness (a bare 0..1 value, a fraction of the format's full scale). The
     two forms may be mixed and are both computed on the ORIGINAL frame, before
     the stretch, which maps [black_level, white_level] onto the full range.
  3. morphological grow/shrink by a circular aperture (--grow R): R>0 dilates
     (bright stars grow, max filter), R<0 erodes (min filter). Applied BEFORE
     inversion.
  4. inversion (--invert): negative image, done last.

Full scale is the dtype maximum for integer frames (e.g. 65535 for uint16) and
1.0 for float frames. So an absolute --black 0.01 means 0.01 on a float frame and
0.01*65535 on a uint16 frame; a percentile --white 90% means the 90th percentile
(the brightest 10% saturate to white). --black 0 --white 1 (the defaults) is a
true no-op; --black 0% --white 100% clips nothing but min-max normalizes to the
full range. Output is ALWAYS single-channel and the input dtype is preserved.
"""

import sys
import os
import numpy as np
from astropy.io import fits

# Import shared utilities
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "makemask - Build processing masks from FITS images.\n"
        "\n"
        "Usage:\n"
        "  makemask.py input_spec output_spec [options]\n"
        "\n"
        "  input_spec   - single file, wildcard (*.fit), numbered (img0001.fit),\n"
        "                 or @list.txt\n"
        "  output_spec  - single file, numbered pattern, or directory\n"
        "\n"
        "Options (all optional; applied in this order):\n"
        "  --black LEVEL - black point (maps to 0). Two forms:\n"
        "                    '90%'  -> a PERCENTILE (0..100) of the frame\n"
        "                    '0.01' -> an ABSOLUTE brightness, a fraction (0..1)\n"
        "                              of full scale (dtype max for int, 1.0 for\n"
        "                              float)\n"
        "  --white LEVEL - white point (maps to full scale), same two forms.\n"
        "                  Both endpoints are measured from 0 and computed on the\n"
        "                  original frame, BEFORE the stretch, which then maps\n"
        "                  [black_level, white_level] onto the full range. The two\n"
        "                  forms may be mixed (e.g. --black 1% --white 0.95).\n"
        "                  Giving either enables the stretch; the missing side\n"
        "                  falls back to its default. Defaults --black 0 --white 1\n"
        "                  are a true no-op; --black 0% --white 100% clips nothing\n"
        "                  but min-max normalizes to the full range. The resolved\n"
        "                  black level must be < white.\n"
        "  --grow R      - grow/shrink by a circular aperture of radius R px:\n"
        "                  R>0 dilates (stars grow, max filter), R<0 erodes\n"
        "                  (min filter). Applied before inversion.\n"
        "  --invert      - negative image (done last)\n"
        "\n"
        "Full scale: dtype maximum for integer frames (65535 for uint16, ...),\n"
        "1.0 for float frames. So absolute --black 0.01 is 0.01 on a float frame\n"
        "and 0.01*65535 on a uint16 frame. A percentile --white 90% clips the\n"
        "brightest 10% to white.\n"
        "\n"
        "Examples:\n"
        "  makemask.py src.fit mask.fit --black 0.01 --white 90%\n"
        "      black at absolute 0.01 (0.01*full scale), white at the 90th\n"
        "      percentile (brightest 10% saturate to white).\n"
        "  makemask.py src.fit mask.fit --black 1% --white 0.95\n"
        "      black at the 1st percentile, white at absolute 0.95 of full scale.\n"
        "  makemask.py stars.fit halo.fit --white 80% --grow 6 --invert\n"
        "      stretch, grow bright stars by 6 px, then invert.\n"
        "\n"
        "Colour input is reduced to grey with g = (R + 2*G + B) / 4 before\n"
        "processing; the output is always single-channel.\n"
    )
    sys.exit(1)


def _opt_value(args, idx, name):
    if idx + 1 >= len(args):
        sys.stderr.write(f"Error: {name} requires a value.\n")
        sys.exit(1)
    try:
        return float(args[idx + 1])
    except ValueError:
        sys.stderr.write(f"Error: {name} must be a number.\n")
        sys.exit(1)


def _parse_level(spec, name):
    """Parse a black/white level. A trailing '%' means a percentile (0..100);
    otherwise it is an absolute brightness as a fraction of full scale (0..1).
    Returns ('pct', p) or ('val', v)."""
    s = str(spec).strip()
    if s.endswith("%"):
        try:
            p = float(s[:-1])
        except ValueError:
            sys.stderr.write(f"Error: {name} percentile must be a number before '%'.\n")
            sys.exit(1)
        if not (0.0 <= p <= 100.0):
            sys.stderr.write(f"Error: {name} percentile must be in 0..100.\n")
            sys.exit(1)
        return ("pct", p)
    try:
        v = float(s)
    except ValueError:
        sys.stderr.write(f"Error: {name} must be a value in 0..1 or a percentile like '90%'.\n")
        sys.exit(1)
    if not (0.0 <= v <= 1.0):
        sys.stderr.write(f"Error: {name} brightness value must be in 0..1 "
                         f"(or use a percentile like '90%').\n")
        sys.exit(1)
    return ("val", v)


def _fmt_level(level):
    """Format a parsed level for HISTORY: ('pct', 90) -> '90%', ('val', 0.01) -> '0.01'."""
    mode, num = level
    return f"{num:g}%" if mode == "pct" else f"{num:g}"


def parse_args(argv):
    args = argv[1:]
    black = None
    white = None
    grow = 0.0
    invert = False

    positional = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--invert":
            invert = True
            i += 1
            continue
        if a == "--black":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --black requires a value.\n")
                sys.exit(1)
            black = _parse_level(args[i + 1], "--black")
            i += 2
            continue
        if a == "--white":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --white requires a value.\n")
                sys.exit(1)
            white = _parse_level(args[i + 1], "--white")
            i += 2
            continue
        if a == "--grow":
            grow = _opt_value(args, i, "--grow")
            i += 2
            continue
        if a.startswith("--"):
            sys.stderr.write(f"Error: unknown option: {a}\n")
            usage()
        positional.append(a)
        i += 1

    if len(positional) < 2:
        usage()

    if black is not None or white is not None:
        # Giving either level enables the stretch; the missing side falls back to
        # its no-op default (black at 0.0 of full scale, white at 1.0 of full
        # scale). Levels cross-check happens in build_mask once resolved, because
        # a percentile and an absolute value cannot be compared before the frame
        # is known.
        if black is None:
            black = ("val", 0.0)
        if white is None:
            white = ("val", 1.0)

    return positional[0], positional[1], black, white, grow, invert


# ---------------------------------------------------------------------------
# Processing
# ---------------------------------------------------------------------------

def to_gray(data):
    """Reduce to a 2D float64 grey image. 3-channel input -> (R+2G+B)/4."""
    if data.ndim == 2:
        return data.astype(np.float64), False
    if data.ndim == 3:
        if data.shape[0] == 3:            # (C, H, W) - PULSAR convention
            r, g, b = data[0], data[1], data[2]
        elif data.shape[2] == 3:          # (H, W, C)
            r, g, b = data[..., 0], data[..., 1], data[..., 2]
        else:
            raise ValueError(f"Unsupported 3D shape {data.shape}: need 3 channels.")
        r = r.astype(np.float64)
        gray = (r + 2.0 * g.astype(np.float64) + b.astype(np.float64)) / 4.0
        return gray, True
    raise ValueError(f"Unsupported FITS shape: {data.shape}")


def disk_footprint(radius):
    """Circular boolean footprint of the given integer radius."""
    y, x = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return (x * x + y * y) <= radius * radius


def _resolve_level(level, finite_vals, full):
    """Resolve a parsed level to an absolute pixel value on the input scale.
    Percentile ('pct', p) -> p-th percentile of the finite pixels; absolute
    ('val', v) -> v * full (v is a fraction of full scale)."""
    mode, num = level
    if mode == "pct":
        return float(np.percentile(finite_vals, num))
    return float(num) * full


def build_mask(gray, orig_dtype, black, white, grow, invert):
    """Apply the mask pipeline to a 2D float64 grey image; return (work, stretched).

    black/white are parsed levels (see _parse_level): each is either a percentile
    ('pct', p) or an absolute fraction of full scale ('val', v). They are resolved
    independently and may mix modes. Both are computed from the ORIGINAL grey
    frame, before any stretch; the stretch then maps [black_level, white_level]
    onto the full output range. Non-finite pixels (NaN/Inf, e.g. out-of-footprint
    regions) are excluded from the percentiles and pinned to the black point, so a
    single bad pixel cannot poison the stretch and zero the whole frame."""
    is_int = np.issubdtype(orig_dtype, np.integer)
    full = float(np.iinfo(orig_dtype).max) if is_int else 1.0

    finite_mask = np.isfinite(gray)
    if not finite_mask.any():
        raise RuntimeError("input frame has no finite pixels")
    all_finite = bool(finite_mask.all())
    finite_vals = gray if all_finite else gray[finite_mask]

    work = gray
    stretched = False
    if black is not None:                 # (white is also set by parse_args)
        vmin = _resolve_level(black, finite_vals, full)
        vmax = _resolve_level(white, finite_vals, full)
        if vmax <= vmin:
            raise RuntimeError(
                f"resolved black level ({vmin:g}) >= white level ({vmax:g}); "
                f"check --black/--white (a constant frame or crossed levels)")
        # non-finite pixels -> black point, so finite pixels stretch correctly
        src = gray if all_finite else np.where(finite_mask, gray, vmin)
        work = (src - vmin) / (vmax - vmin)
        np.clip(work, 0.0, 1.0, out=work)
        work = work * full
        stretched = True
    elif not all_finite:
        # no stretch: keep the pipeline finite for morphology / inversion
        work = np.where(finite_mask, gray, float(finite_vals.min()))

    radius = int(round(abs(grow)))
    if grow != 0.0 and radius >= 1:
        from scipy.ndimage import maximum_filter, minimum_filter
        footprint = disk_footprint(radius)
        if grow > 0.0:
            work = maximum_filter(work, footprint=footprint)
        else:
            work = minimum_filter(work, footprint=footprint)

    if invert:
        white_ref = full if stretched else float(work.max())  # work is finite here
        work = white_ref - work

    return work, stretched


def to_output(work, orig_dtype):
    """Sanitize NaN/Inf and cast to the original dtype (clamp+round for int)."""
    work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        return np.clip(np.rint(work), info.min, info.max).astype(orig_dtype)
    return work.astype(orig_dtype)


def process_file(infile, outfile, black, white, grow, invert):
    with fits.open(infile, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
    if data is None:
        raise RuntimeError("no image data")

    orig_dtype = data.dtype
    gray, was_colour = to_gray(data)
    work, stretched = build_mask(gray, orig_dtype, black, white, grow, invert)
    out = to_output(work, orig_dtype)

    # Output is single-channel: drop cube/scaling keys.
    for key in ("NAXIS3", "BSCALE", "BZERO"):
        if key in header:
            del header[key]
    header["NAXIS"] = 2

    parts = []
    if was_colour:
        parts.append("gray=(R+2G+B)/4")
    if stretched:
        parts.append(f"stretch black={_fmt_level(black)} white={_fmt_level(white)}")
    radius = int(round(abs(grow)))
    if grow != 0.0 and radius >= 1:
        parts.append(f"grow={'+' if grow > 0 else '-'}{radius}px")
    if invert:
        parts.append("invert")
    header["HISTORY"] = "makemask.py: " + (", ".join(parts) if parts else "passthrough (mono)")

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    fits.PrimaryHDU(out, header=header).writeto(outfile, overwrite=True)


def main():
    input_spec, output_spec, black, white, grow, invert = parse_args(sys.argv)

    io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
    total = len(io_pairs)
    if total == 0:
        sys.stderr.write("Error: no input files.\n")
        sys.exit(1)

    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            process_file(infile, outfile, black, white, grow, invert)
            sys.stderr.write(f"\rProcessed {i} / {total}")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
