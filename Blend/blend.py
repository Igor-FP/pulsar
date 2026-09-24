#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
blend - Combine two FITS images through an opacity mask.

    blend.py source output mask operand [options]

The mask is a greyscale image whose value is the OPACITY of the operand at
each pixel (0 = fully source, full scale = fully operand, half = arithmetic
mean). With the normalized mask m in [0, 1]:

    output = source * (1 - m) + operand * m

  white mask pixel (m = 1) -> output = operand
  black mask pixel (m = 0) -> output = source
  grey  mask pixel (m = 0.5) -> output = (source + operand) / 2

Mask normalization: the mask is scaled to [0, 1] by its OWN full scale - the
dtype maximum for integer masks (65535 for uint16, 255 for uint8, ...) and 1.0
for float masks (a float mask must therefore already be in [0, 1]). A colour
mask is reduced to grey with (R + 2*G + B) / 4.

Options:
  --mtf [K]   Apply the PixInsight midtone transfer function to the mask (in
              [0, 1]) before blending, and before --invert. K is the midtones
              balance, same meaning as mtf.py (0 < K < 1): K < 0.5 brightens the
              mask midtones (more operand shows through), K > 0.5 darkens them
              (more source). K defaults to 0.25 when omitted.
  --invert    Use the inverted mask (m -> 1 - m). Applied AFTER --mtf.

The operand may be a numeric constant (a "virtual file" filled with that value)
to blend toward a flat level. source and operand must share the same shape and
data scale; the output keeps the source dtype and header. 2D (mono) and
3-channel colour images are supported; a mono mask is broadcast across channels.
Standard batch_utils input/output specs apply (single file, wildcard, numbered
sequence, @list.txt); mask and operand may be a single file (broadcast to all
sources) or a sequence matching the source count.
"""

import sys
import os
import numpy as np
from astropy.io import fits

# Import shared utilities
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


DEFAULT_MTF_K = 0.25
DEFAULT_OPACITY = 0.5


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "blend - Combine two FITS images through an opacity mask and/or a flat opacity.\n"
        "\n"
        "Usage:\n"
        "  blend.py source output mask operand [--mtf [K]] [--invert] [-o N]\n"
        "  blend.py source output operand -o N [--mtf [K]] [--invert]      (no mask)\n"
        "\n"
        "Positional arguments (in this order):\n"
        "  source    - base image(s): single file, wildcard (*.fit), numbered\n"
        "              sequence (img0001.fit), or @list.txt\n"
        "  output    - single file, numbered pattern (out0001.fit), or directory\n"
        "  mask      - greyscale opacity image (FITS file, or a sequence matching the\n"
        "              source count). NOT a numeric constant. May be omitted only\n"
        "              when -o/--opacity is given (then a flat opacity is used).\n"
        "  operand   - image blended in where opacity is high: a FITS file, a\n"
        "              matching sequence, or a numeric constant (flat level).\n"
        "\n"
        "The operand's opacity at each pixel is m in [0, 1]:\n"
        "  output = source * (1 - m) + operand * m    (m=1 -> operand, m=0 -> source)\n"
        "\n"
        "m is built in THIS ORDER:\n"
        "  1. base: the mask, scaled to [0,1] by its own full scale (dtype maximum\n"
        "     for integer masks, 1.0 for float masks - a float mask must be in\n"
        "     [0,1]; a colour mask -> grey (R+2*G+B)/4). With no mask, the base is a\n"
        "     flat 1.0 everywhere.\n"
        "  2. --mtf [K]: reshape the base opacity's midtones (no effect on a flat\n"
        "     base, which has no midtones).\n"
        "  3. --invert: m -> 1 - m.\n"
        "  4. -o N / --opacity N: multiply the whole result by N, applied LAST -\n"
        "     'do everything as usual, but at N strength'.\n"
        "\n"
        "So -o only scales the FINAL opacity down:\n"
        "  with a mask:    m = (mask, after MTF and invert) * N   -> a weaker mask\n"
        "  without a mask: m = N on every pixel                   -> a flat N overlay\n"
        "\n"
        "Options:\n"
        "  -o N, --opacity N  scale the final opacity by N, applied LAST (after --mtf\n"
        "              and --invert). N is a fraction in [0,1], or a percent with a\n"
        "              '%' suffix (e.g. 0.3 or 30%); N defaults to 0.5 when the flag\n"
        "              is given without a value. With a mask this weakens it to N\n"
        "              strength; with the mask argument omitted it gives a flat N\n"
        "              opacity everywhere.\n"
        "  --mtf [K] - apply the MTF (midtone transfer function) to the base opacity\n"
        "              in [0,1], BEFORE --invert and BEFORE -o. K is the midtones\n"
        "              balance (0<K<1), same as mtf.py: K<0.5 raises opacity (more\n"
        "              operand), K>0.5 lowers it. K defaults to 0.25 when omitted.\n"
        "  --invert  - use the inverted opacity (m -> 1 - m); applied AFTER --mtf\n"
        "              and BEFORE -o.\n"
        "\n"
        "Examples:\n"
        "  blend base.fit out.fit mask.fit stars.fit\n"
        "      blend stars.fit over base.fit where mask.fit is bright (full mask).\n"
        "  blend base.fit out.fit mask.fit stars.fit -o 50%\n"
        "      apply the same mask at half strength (mask * 0.5).\n"
        "  blend base.fit out.fit stars.fit -o 30%\n"
        "      flat 30% overlay of stars.fit (no mask).\n"
        "  blend base.fit out.fit hi.fit -o 0.5\n"
        "      an even 50/50 average of the two images.\n"
        "  blend base.fit out.fit mask.fit hi.fit --mtf 0.2 --invert\n"
        "      reshape then invert the mask before blending.\n"
        "\n"
        "source and operand must share shape and data scale; the output keeps the\n"
        "source dtype and header. 2D and 3-channel colour images are supported\n"
        "(a mono mask or the flat opacity is broadcast across channels).\n"
    )
    sys.exit(1)


def _validate_mtf_k(value):
    """Validate an MTF midtones balance: a float strictly in (0, 1)."""
    if not (0.0 < value < 1.0):     # also rejects NaN (all comparisons False)
        sys.stderr.write("Error: --mtf K must be a finite value in range (0, 1) exclusive.\n")
        sys.exit(1)
    return float(value)


def _looks_numeric(token):
    """True if the token parses as an opacity value (bare number or N%)."""
    try:
        float(token[:-1] if token.endswith("%") else token)
        return True
    except ValueError:
        return False


def _parse_opacity(spec):
    """Flat opacity: a fraction in [0,1], or a percent with a '%' suffix."""
    s = str(spec).strip()
    if s.endswith("%"):
        try:
            value = float(s[:-1]) / 100.0
        except ValueError:
            sys.stderr.write("Error: --opacity percent must be a number before '%'.\n")
            sys.exit(1)
    else:
        try:
            value = float(s)
        except ValueError:
            sys.stderr.write("Error: --opacity must be a number in [0,1] or a percent like '30%'.\n")
            sys.exit(1)
    if not (0.0 <= value <= 1.0):      # also rejects NaN
        sys.stderr.write("Error: --opacity must resolve to [0,1] (0..100%).\n")
        sys.exit(1)
    return value


def parse_args(argv):
    args = argv[1:]
    invert = False
    mtf_k = None                 # None = MTF disabled
    opacity = None               # None = mask mode; float in [0,1] = flat opacity
    positional = []

    i = 0
    while i < len(args):
        a = args[i]
        if a == "--invert":
            invert = True
            i += 1
            continue
        if a == "--mtf":
            mtf_k = DEFAULT_MTF_K
            i += 1
            # optional K: consume the next token only if it parses as a number
            if i < len(args) and not args[i].startswith("--"):
                try:
                    val = float(args[i])
                except ValueError:
                    val = None
                if val is not None:
                    mtf_k = _validate_mtf_k(val)
                    i += 1
            continue
        if a in ("-o", "--opacity"):
            opacity = DEFAULT_OPACITY
            i += 1
            # optional N: consume the next token only if it looks like a value (N or N%)
            if i < len(args) and not args[i].startswith("-") and _looks_numeric(args[i]):
                opacity = _parse_opacity(args[i])
                i += 1
            continue
        if a.startswith("--"):
            sys.stderr.write(f"Error: unknown option: {a}\n")
            usage()
        positional.append(a)
        i += 1

    if opacity is not None:
        # -o scales the final opacity; the mask is OPTIONAL in this mode.
        #   4 positionals -> weaken a mask (mask * N)
        #   3 positionals -> flat N opacity (no mask)
        if len(positional) == 4:
            source, output, mask, operand = positional
        elif len(positional) == 3:
            source, output, operand = positional
            mask = None
        else:
            sys.stderr.write(
                "Error: with -o/--opacity need 4 positional arguments "
                "(source output mask operand) to weaken a mask, or 3 "
                "(source output operand) for a flat opacity; "
                f"got {len(positional)}.\n")
            usage()
        return {
            "source": source,
            "output": output,
            "mask": mask,
            "operand": operand,
            "invert": invert,
            "mtf_k": mtf_k,
            "opacity": opacity,
        }

    if len(positional) != 4:
        sys.stderr.write(
            "Error: need exactly 4 positional arguments: source output mask operand "
            f"(got {len(positional)}); or use -o/--opacity N for a flat opacity "
            "(source output operand).\n")
        if mtf_k is not None:
            sys.stderr.write(
                "Hint: if the operand is a numeric constant, give --mtf an explicit K "
                "(e.g. --mtf 0.3) so the constant is not read as K.\n")
        usage()

    return {
        "source": positional[0],
        "output": positional[1],
        "mask": positional[2],
        "operand": positional[3],
        "invert": invert,
        "mtf_k": mtf_k,
        "opacity": None,
    }


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def apply_mtf(x, m):
    """Midtone transfer function on a normalized [0,1] array (as in mtf.py).

    mtf(x, m) = (1-m)*x / (m + x*(1-2*m)); copied here to keep the tool
    self-contained (tools must not import from other tool modules)."""
    if abs(m - 0.5) < 1e-10:
        return x.copy()
    num = (1.0 - m) * x
    den = m + x * (1.0 - 2.0 * m)
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.where(den != 0.0, num / den, np.where(x <= 0.0, 0.0, 1.0))
    return np.clip(result, 0.0, 1.0)


def to_gray(data):
    """Reduce mask data to a 2D float64 grey image. 3-channel -> (R+2G+B)/4."""
    if data.ndim == 2:
        return data.astype(np.float64)
    if data.ndim == 3:
        if data.shape[0] == 3:              # (C, H, W) - PULSAR convention
            r, g, b = data[0], data[1], data[2]
        elif data.shape[2] == 3:            # (H, W, C)
            r, g, b = data[..., 0], data[..., 1], data[..., 2]
        else:
            raise ValueError(f"unsupported mask 3D shape {data.shape}: need 3 channels")
        return (r.astype(np.float64) + 2.0 * g.astype(np.float64)
                + b.astype(np.float64)) / 4.0
    raise ValueError(f"unsupported mask shape {data.shape}")


def load_mask(mask_path, spatial_shape):
    """Load a mask and return opacity m in [0,1] with shape spatial_shape (2D).

    Integer masks are scaled by the dtype maximum; float masks are used as-is
    (and must already lie in [0,1]). Non-finite pixels map to 0 (fully source)."""
    with fits.open(mask_path, memmap=False) as hdul:
        mdata = hdul[0].data
    if mdata is None:
        raise ValueError(f"mask '{mask_path}' has no image data")
    mdtype = mdata.dtype

    gray = to_gray(mdata)
    if gray.shape != spatial_shape:
        raise ValueError(
            f"mask shape {gray.shape} does not match image spatial shape "
            f"{spatial_shape}")

    finite = np.isfinite(gray)
    if not finite.any():
        raise ValueError(f"mask '{mask_path}' has no finite pixels")

    if np.issubdtype(mdtype, np.integer):
        full = float(np.iinfo(mdtype).max)
        m = np.nan_to_num(gray, nan=0.0, posinf=0.0, neginf=0.0) / full
    else:
        vals = gray[finite]
        fmin, fmax = float(vals.min()), float(vals.max())
        if fmax > 1.5 or fmin < -0.5:
            raise ValueError(
                f"float mask '{os.path.basename(mask_path)}' looks un-normalized "
                f"(range [{fmin:g}, {fmax:g}]); a float opacity mask must be in "
                f"[0,1] - rescale it or use an integer mask")
        m = np.nan_to_num(gray, nan=0.0, posinf=0.0, neginf=0.0)

    np.clip(m, 0.0, 1.0, out=m)
    return m


def blend_arrays(source, operand, m):
    """Blend: source*(1-m) + operand*m in float64. A 2D mask m is broadcast
    across the channel axis of a 3-channel source/operand."""
    src = source.astype(np.float64)
    op = operand.astype(np.float64)
    if src.ndim == 3:
        if src.shape[0] == 3:               # (C, H, W)
            mm = m[np.newaxis, :, :]
        else:                               # (H, W, C)
            mm = m[:, :, np.newaxis]
    else:
        mm = m
    return src * (1.0 - mm) + op * mm


def to_output(work, orig_dtype):
    """Sanitize NaN/Inf and cast to the original dtype (clamp+round for int)."""
    work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        return np.clip(np.rint(work), info.min, info.max).astype(orig_dtype)
    return work.astype(orig_dtype)


def _spatial_shape(data):
    """(H, W) spatial shape of a 2D or 3-channel image; raises otherwise."""
    if data.ndim == 2:
        return data.shape
    if data.ndim == 3:
        if data.shape[0] == 3:
            return data.shape[1:]
        if data.shape[2] == 3:
            return data.shape[:2]
    raise ValueError(f"unsupported source shape {data.shape}")


def process_file(sfile, ofile, mask_spec, operand_spec, index, invert, mtf_k, opacity):
    with fits.open(sfile, memmap=False) as hdul:
        sdata = hdul[0].data
        header = hdul[0].header.copy()
    if sdata is None:
        raise ValueError("source has no image data")

    orig_dtype = sdata.dtype
    spatial = _spatial_shape(sdata)

    # Base opacity m in [0,1]: the mask (scaled), or a flat 1.0 when no mask is
    # given. Pipeline order: MTF -> invert -> scale by the flat opacity N (LAST).
    if mask_spec is not None:
        mask_path = batch_utils.get_operand_for_file(mask_spec, index)
        m = load_mask(mask_path, spatial)
        m_desc = f"mask={os.path.basename(mask_path)}"
    else:
        m = np.full(spatial, 1.0, dtype=np.float64)
        m_desc = "mask=none"
    if mtf_k is not None:
        m = apply_mtf(m, mtf_k)
    if invert:
        m = 1.0 - m
    if opacity is not None:
        m = m * float(opacity)          # weaken the final opacity to N strength (last)

    # Operand (constant or file), matched to the full source shape.
    operand_raw = batch_utils.get_operand_for_file(operand_spec, index)
    operand = batch_utils.resolve_operand_value(operand_raw, sdata.shape, orig_dtype)

    out = to_output(blend_arrays(sdata, operand, m), orig_dtype)

    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    if isinstance(operand_raw, (int, float)):
        op_desc = f"{float(operand_raw):g}"
    else:
        op_desc = os.path.basename(operand_raw)
    parts = [m_desc]
    if mtf_k is not None:
        parts.append(f"mtf K={mtf_k:g}")
    if invert:
        parts.append("inverted")
    if opacity is not None:
        parts.append(f"opacity x{opacity:g}")
    parts.append(f"operand={op_desc}")
    header["HISTORY"] = ("blend.py: out = source*(1-m) + operand*m; "
                         + ", ".join(parts))

    out_dir = os.path.dirname(ofile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    fits.PrimaryHDU(out, header=header).writeto(ofile, overwrite=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    cfg = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(cfg["source"], cfg["output"])
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)
    if not io_pairs:
        sys.stderr.write("Error: no source files to process.\n")
        sys.exit(1)

    total = len(io_pairs)

    mask_spec = None
    if cfg["mask"] is not None:
        try:
            mask_spec = batch_utils.build_operand_spec(cfg["mask"], total)
        except Exception as e:
            sys.stderr.write(f"Error (mask): {e}\n")
            sys.exit(1)
        if isinstance(mask_spec, float):
            sys.stderr.write(
                "Error: mask must be a FITS file (or a matching sequence), "
                "not a numeric constant.\n")
            sys.exit(1)

    try:
        operand_spec = batch_utils.build_operand_spec(cfg["operand"], total)
    except Exception as e:
        sys.stderr.write(f"Error (operand): {e}\n")
        sys.exit(1)

    for i, (sfile, ofile) in enumerate(io_pairs, start=1):
        try:
            process_file(sfile, ofile, mask_spec, operand_spec, i - 1,
                         cfg["invert"], cfg["mtf_k"], cfg["opacity"])
            sys.stderr.write(f"\rProcessed {i} / {total}")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{sfile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
