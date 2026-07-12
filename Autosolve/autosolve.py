#!/usr/bin/env python3
# autosolve.py — WSL-aware astrometric solver + optional rectification + subpixel fine align + optional WCS refit from .corr
# Depends on: astropy, numpy, scipy (fine-align), reproject (interp/exact), astrometry.net (solve-field, wcs-resample via WSL)
# Uses ../lib/batch_utils.py for input/pattern handling.
# !!! Code comments in ENGLISH ONLY please !!!

import argparse
import os
import re
import sys
import glob
import shlex
import subprocess
import urllib.request
import urllib.parse
from typing import List, Tuple, Optional

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

# ---- optional reproject imports (lazy) ----
_have_reproject = True
try:
    from reproject import reproject_interp, reproject_exact
except Exception:
    _have_reproject = False

# ---- scipy for spline shift (fine align) ----
_have_scipy = True
try:
    from scipy.ndimage import shift as ndi_shift, median_filter as ndi_median_filter
except Exception:
    _have_scipy = False

# ---- optional Pillow (JPEG input/output) ----
_have_pil = True
try:
    from PIL import Image
except Exception:
    _have_pil = False

# ---- import batch_utils from ../lib ----
HERE = os.path.abspath(os.path.dirname(__file__))
LIBDIR = os.path.abspath(os.path.join(HERE, "..", "lib"))
if LIBDIR not in sys.path:
    sys.path.insert(0, LIBDIR)
try:
    import batch_utils as bu
except Exception as e:
    print(f"Failed to import batch_utils from {LIBDIR}: {e}", file=sys.stderr)
    sys.exit(1)

# FITS-header <-> JPEG codec (shared with fits2tiff; single source of truth for
# the on-disk format, so a JPEG solved here reads back identically elsewhere).
try:
    from image_fits_header import (fits_header_to_json, json_to_fits_header,
                                   embed_header_in_jpeg, read_fits_header_from_jpeg)
    _have_jpeg_codec = True
except Exception:
    _have_jpeg_codec = False


# ---------------------- RGB support ----------------------
def _normalize_rgb_data(data):
    """Normalize FITS data to standard layout. Returns (data, is_rgb).

    - 2D (ny, nx) -> unchanged, is_rgb=False
    - 3D (3, ny, nx) -> unchanged, is_rgb=True
    - 3D (ny, nx, 3) -> transposed to (3, ny, nx), is_rgb=True
    """
    if data.ndim == 2:
        return data, False
    if data.ndim == 3:
        if data.shape[0] == 3:
            return data, True
        if data.shape[2] == 3:
            return np.transpose(data, (2, 0, 1)), True
    return data, False


# ---------------------- JPEG input/output (solve round-trip) ----------------------
# autosolve can take JPEGs (e.g. produced by fits2tiff) and return JPEGs with the
# solved WCS packed into a COM marker, via the shared lib/image_fits_header codec.

# WCS/SIP keyword prefixes stripped from a recovered header before re-solving, so
# the output JPEG carries only the fresh WCS produced by solve-field.
_WCS_STRIP_PREFIXES = (
    "WCSAXES", "CRVAL", "CRPIX", "CDELT", "CD", "PC", "CTYPE", "CUNIT",
    "CROTA", "PV", "PS", "A_", "B_", "AP_", "BP_", "LONPOLE", "LATPOLE",
    "RADESYS", "EQUINOX",
)


def _is_jpeg(path: str) -> bool:
    return os.path.splitext(path)[1].lower() in (".jpg", ".jpeg")


def _safe_remove(path):
    if path and os.path.exists(path):
        try:
            os.remove(path)
        except OSError:
            pass


def _clean_recovered_header(hdr):
    """Strip WCS/SIP from a JPEG-recovered header so the solved output carries a
    fresh WCS only. Preserve the field center as an RA/DEC hint when the source
    had only a prior WCS (CRVAL) and no explicit RA/DEC."""
    if hdr is None:
        return None
    if ("RA" not in hdr) and ("CRVAL1" in hdr):
        try:
            hdr["RA"] = float(hdr["CRVAL1"])
        except Exception:
            pass
    if ("DEC" not in hdr) and ("CRVAL2" in hdr):
        try:
            hdr["DEC"] = float(hdr["CRVAL2"])
        except Exception:
            pass
    for key in list(hdr.keys()):
        if not key:
            continue
        ku = str(key).upper()
        if any(ku.startswith(p) for p in _WCS_STRIP_PREFIXES):
            try:
                del hdr[key]
            except Exception:
                pass
    return hdr


def jpeg_to_temp_fits(jpeg_path: str, temp_fits_path: str) -> None:
    """Decode a JPEG into a temporary FITS for solving. Recovers any embedded
    FITS header (structural + WCS stripped) so acquisition hints (FOCALLEN,
    RA/DEC) guide the solve. Mono 'L' -> (ny,nx); colour -> (3,ny,nx)."""
    if not _have_pil:
        raise RuntimeError("Pillow is required for JPEG input (pip install Pillow)")
    with Image.open(jpeg_path) as img:
        if img.mode == "L":
            arr = np.array(img, dtype=np.uint8)
        else:
            arr = np.transpose(np.array(img.convert("RGB"), dtype=np.uint8), (2, 0, 1))

    hdr = None
    if _have_jpeg_codec:
        try:
            hdr_dict = read_fits_header_from_jpeg(jpeg_path)
            if hdr_dict:
                hdr = _clean_recovered_header(json_to_fits_header(hdr_dict))
        except Exception:
            hdr = None

    fits.PrimaryHDU(data=arr, header=hdr).writeto(temp_fits_path, overwrite=True)


def fits_to_jpeg(fits_path: str, jpeg_path: str, quality: int) -> None:
    """Write a solved FITS (mono or RGB) to JPEG, embedding its header (WCS
    included) in a COM marker. Pixel values are treated as 8-bit display data
    (the input JPEG was 8-bit); non-uint8 arrays are clipped to [0,255]."""
    if not _have_pil:
        raise RuntimeError("Pillow is required for JPEG output (pip install Pillow)")
    with fits.open(fits_path, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
    if data is None:
        raise ValueError(f"'{fits_path}' has no image data to write as JPEG.")

    data, is_rgb = _normalize_rgb_data(data)
    work = np.asarray(data)
    if work.dtype != np.uint8:
        work = np.nan_to_num(work.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)
        work = np.clip(np.rint(work), 0, 255).astype(np.uint8)

    if is_rgb:
        img = Image.fromarray(np.transpose(work, (1, 2, 0)))   # (3,ny,nx)->(ny,nx,3)
    else:
        img = Image.fromarray(work)

    out_dir = os.path.dirname(jpeg_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    img.save(jpeg_path, "JPEG", quality=int(quality))
    if _have_jpeg_codec:
        embed_header_in_jpeg(jpeg_path, fits_header_to_json(header))


# ---------------------- helpers: OS / paths ----------------------
def is_windows() -> bool:
    return os.name == "nt"


def win_to_wsl_path(path: str) -> str:
    """Convert Windows path (S:\\dir\\file) to WSL path (/mnt/s/dir/file)."""
    p = path.replace("\\", "/")
    if p.startswith("/mnt/") or p.startswith("/"):
        return p
    m = re.match(r"^([A-Za-z]):/(.*)$", p)
    if m:
        drive = m.group(1).lower()
        rest = m.group(2)
        return f"/mnt/{drive}/{rest}"
    return p


# ---------------------- helpers: console ----------------------
def print_if_verbose(msg: str, verbose: bool):
    if verbose:
        print(msg, file=sys.stderr)


def progress_bar(i: int, n: int, width: int = 40):
    """Render a single-line progress bar."""
    done = int(width * i / max(1, n))
    bar = "#" * done + "-" * (width - done)
    print(f"[{bar}] {i}/{n}", end="\r", file=sys.stderr)


# ---------------------- helpers: WCS hints ----------------------
def sexagesimal_to_deg(sex, is_ra=True):
    """Parse a sexagesimal angle to degrees. RA is read in hours (H M S -> deg
    via x15); Dec in degrees. Accepts space / colon / 'h m s' / comma separators.
    Returns None on failure. Sign is taken from the string (so '-00 30 00' works)."""
    s = str(sex).strip().lower()
    neg = s.startswith("-")
    s = s.replace("h", " ").replace("m", " ").replace("s", " ")
    s = s.replace(":", " ").replace(",", " ")
    parts = s.split()
    if len(parts) < 3:
        return None
    try:
        a, b, c = abs(float(parts[0])), float(parts[1]), float(parts[2])
    except Exception:
        return None
    val = a + b / 60.0 + c / 3600.0
    if is_ra:
        return val * 15.0
    return -val if neg else val


def _parse_angle_cli(val, is_ra):
    """Parse an --ra/--dec CLI value: a bare decimal is DEGREES; otherwise a
    sexagesimal string (RA in HOURS, Dec in degrees). Raises ValueError."""
    s = str(val).strip()
    try:
        deg = float(s)           # plain decimal degrees
    except ValueError:
        deg = sexagesimal_to_deg(s, is_ra=is_ra)
        if deg is None:
            fmt = "HH MM SS" if is_ra else "[+-]DD MM SS"
            raise ValueError(f"cannot parse {'RA' if is_ra else 'Dec'} value '{val}' "
                             f"(use decimal degrees or sexagesimal '{fmt}')")
    if not np.isfinite(deg):
        raise ValueError(f"{'RA' if is_ra else 'Dec'} value '{val}' is not a finite number")
    return deg


_SESAME_URL = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/SNV?"


def resolve_name_sesame(name, timeout=20):
    """Resolve an object name to (ra_deg, dec_deg) via CDS Sesame (SIMBAD / NED /
    VizieR). Raises RuntimeError on failure.

    Note: resolution follows the service's disambiguation -- e.g. 'Abell 35' is
    the galaxy cluster ACO 35; the planetary nebula is 'PN A66 35'."""
    url = _SESAME_URL + urllib.parse.quote(str(name).strip())
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", "replace")
    except Exception as e:
        raise RuntimeError(f"Sesame lookup failed for '{name}': {e}")
    for line in text.splitlines():
        if line.startswith("%J "):
            parts = line.split()
            try:
                return float(parts[1]), float(parts[2])
            except (IndexError, ValueError):
                break
    raise RuntimeError(f"Sesame could not resolve name '{name}' "
                       f"(try a more specific identifier, e.g. 'PN A66 35')")


def parse_ra_dec_from_header(hdr) -> Tuple[Optional[float], Optional[float]]:
    """Try CRVAL, then RA/DEC in degrees, then OBJCTRA/OBJCTDEC sexagesimal."""
    if "CRVAL1" in hdr and "CRVAL2" in hdr:
        try:
            return float(hdr["CRVAL1"]), float(hdr["CRVAL2"])
        except Exception:
            pass

    ra = hdr.get("RA")
    dec = hdr.get("DEC")
    if ra is not None and dec is not None:
        try:
            return float(ra), float(dec)
        except Exception:
            pass

    ra_s = hdr.get("OBJCTRA")
    dec_s = hdr.get("OBJCTDEC")
    if ra_s and dec_s:
        ra = sexagesimal_to_deg(ra_s, True)
        dec = sexagesimal_to_deg(dec_s, False)
        if ra is not None and dec is not None:
            return ra, dec

    return None, None


def estimate_scale_from_header(hdr):
    """Estimate arcsec/pixel from FOCALLEN (mm) and X/YPIXSZ (microns)."""
    fl = hdr.get("FOCALLEN") or hdr.get("FOCLEN") or hdr.get("FOCUSLEN")
    xpix = hdr.get("XPIXSZ") or hdr.get("XPIXSIZE") or hdr.get("PIXSIZE1")
    ypix = hdr.get("YPIXSZ") or hdr.get("YPIXSIZE") or hdr.get("PIXSIZE2")
    if not fl or (not xpix and not ypix):
        return None, None, None
    try:
        fl = float(fl)
        if xpix:
            xpix = float(xpix)
        if ypix:
            ypix = float(ypix)
    except Exception:
        return None, None, None
    pix = (xpix + ypix) * 0.5 if (xpix and ypix) else (xpix or ypix)
    scale = 206.264806 * (pix / fl)
    return scale, scale * 0.8, scale * 1.2


# Search-cone radius (deg) used when the field centre is known but the pixel
# scale is not (no FOCALLEN/pixel size, so no FOV to size the cone from). Wide
# enough that slightly imprecise header coordinates still fall inside it, yet far
# tighter than a blind all-sky search.
_NOSCALE_RADIUS_DEG = 12.0


def estimate_radius_from_header(hdr, scale_arcsec):
    """Search radius (deg) from image size and scale. Falls back to a wide cone
    (_NOSCALE_RADIUS_DEG) when the scale or image size is unknown."""
    nx = hdr.get("NAXIS1")
    ny = hdr.get("NAXIS2")
    if not nx or not ny or not scale_arcsec:
        return _NOSCALE_RADIUS_DEG
    try:
        nx = int(nx)
        ny = int(ny)
    except Exception:
        return _NOSCALE_RADIUS_DEG
    fov_x = nx * scale_arcsec / 3600.0
    fov_y = ny * scale_arcsec / 3600.0
    r = max(fov_x, fov_y) * 0.75
    return max(r, 1.0)


# ---------------------- scale / FOV helpers ----------------------
def _write_scale_fov(wcs_fit_path, verbose):
    """Compute and write pixel scale and FOV into the solved WCS header."""
    try:
        with fits.open(wcs_fit_path, mode="update") as hdul:
            hdr = hdul[0].header
            data = hdul[0].data
            ny, nx = data.shape[-2:]
            w = WCS(hdr, naxis=2)

            # Extract CD matrix
            if w.wcs.cd is not None and w.wcs.cd.size == 4:
                cd = np.array(w.wcs.cd, dtype=float)
            else:
                pc = np.array(w.wcs.get_pc(), dtype=float)
                cdelt = np.array(w.wcs.cdelt, dtype=float)
                cd = pc * cdelt.reshape(2, 1)

            # Pixel scale (arcsec/pixel) from CD matrix
            sx = np.hypot(cd[0, 0], cd[1, 0]) * 3600.0  # deg -> arcsec
            sy = np.hypot(cd[0, 1], cd[1, 1]) * 3600.0
            scale = (sx + sy) / 2.0

            # FOV in arcminutes
            fov_x = nx * sx / 60.0
            fov_y = ny * sy / 60.0

            hdr['SCALE'] = (round(scale, 4), 'Pixel scale [arcsec/pixel]')
            hdr['FOV1'] = (round(fov_x, 2), 'Field of view X [arcmin]')
            hdr['FOV2'] = (round(fov_y, 2), 'Field of view Y [arcmin]')

            hdul.flush()
            print_if_verbose(
                f"Scale: {scale:.4f} arcsec/pix, "
                f"FOV: {fov_x:.1f}' x {fov_y:.1f}'", verbose)
    except Exception as e:
        print_if_verbose(f"Warning: could not write scale/FOV: {e}", verbose)


# ---------------------- astrometry calls ----------------------
def run_solve_field(input_fits: str,
                    output_wcs: str,
                    tweak_order: int,
                    use_wsl: bool,
                    verbose: bool,
                    ra_override: Optional[float] = None,
                    dec_override: Optional[float] = None,
                    radius_override: Optional[float] = None,
                    fovx: Optional[float] = None,
                    fovy: Optional[float] = None,
                    downsample: int = 1) -> None:
    """Call solve-field; write directly to output_wcs via --new-fits.

    For RGB input: extracts green channel to a temp mono file for solving,
    then restores full RGB data into the _wcs.fit output.

    Optional hints override header-derived ones: ra/dec_override (deg) set the
    search centre; fovx/fovy (deg) give the field of view, from which the pixel
    scale is computed; radius_override (deg) sets the search radius; downsample
    is passed to solve-field's source extraction.
    """
    with fits.open(input_fits, memmap=False) as hdul:
        hdr = hdul[0].header
        orig_data = hdul[0].data
        orig_dtype = orig_data.dtype

    orig_data, is_rgb = _normalize_rgb_data(orig_data)
    ny, nx = orig_data.shape[-2:]

    # For RGB: extract green channel to temp mono file for solve-field
    solve_input = input_fits
    green_tmp = None
    if is_rgb:
        green_tmp = os.path.splitext(input_fits)[0] + "_green_tmp.fit"
        green_data = orig_data[1]  # G channel (index 1)
        green_hdr = hdr.copy()
        # Ensure 2D header for solve-field
        for key in ('NAXIS3',):
            if key in green_hdr:
                del green_hdr[key]
        green_hdr['NAXIS'] = 2
        fits.PrimaryHDU(data=green_data, header=green_hdr).writeto(
            green_tmp, overwrite=True)
        solve_input = green_tmp
        print_if_verbose("RGB input: extracting G channel for solving", verbose)

    # Position hint: explicit CLI override, else derive from the header.
    if ra_override is not None and dec_override is not None:
        ra_hint, dec_hint = ra_override, dec_override
    else:
        ra_hint, dec_hint = parse_ra_dec_from_header(hdr)

    # Scale hint: from a user-supplied field of view (deg), else from the header.
    if fovx is not None or fovy is not None:
        if fovx is not None and nx:
            scale = float(fovx) * 3600.0 / float(nx)
        elif fovy is not None and ny:
            scale = float(fovy) * 3600.0 / float(ny)
        else:
            scale = None
        if scale and scale > 0:
            # Generous +-25% window: a FOV given by hand is only approximate.
            scale_low, scale_high = scale * 0.75, scale * 1.25
            print_if_verbose(f"FOV hint -> pixel scale {scale:.4f} arcsec/pix", verbose)
        else:
            scale = scale_low = scale_high = None
    else:
        scale, scale_low, scale_high = estimate_scale_from_header(hdr)

    # Search radius: explicit override, else sized from the (possibly FOV-derived) scale.
    radius = radius_override if radius_override is not None else estimate_radius_from_header(hdr, scale)

    if use_wsl:
        input_cmd = win_to_wsl_path(solve_input)
        output_cmd = win_to_wsl_path(output_wcs)
    else:
        input_cmd = solve_input
        output_cmd = output_wcs

    solve_args = [
        "solve-field",
        "--no-plots",
        "--overwrite",
        "--new-fits", output_cmd,
        "--tweak-order", str(tweak_order),
        "--downsample", str(int(downsample)),
    ]
    if scale_low and scale_high:
        solve_args += ["--scale-low", f"{scale_low}", "--scale-high", f"{scale_high}", "--scale-units", "arcsecperpix"]
    if ra_hint is not None and dec_hint is not None:
        # A position hint does not require a known scale: estimate_radius_from_header
        # falls back to a wide cone (_NOSCALE_RADIUS_DEG) when scale is absent, so a
        # centre from --ra/--dec/--name (or the header) still constrains the search.
        solve_args += ["--ra", f"{ra_hint}", "--dec", f"{dec_hint}", "--radius", f"{radius}"]
    solve_args.append(input_cmd)

    if use_wsl:
        solve_cmd_str = " ".join(shlex.quote(a) for a in solve_args)
        cmd = ["wsl", "bash", "-lc", solve_cmd_str]
    else:
        cmd = solve_args

    print_if_verbose("Running: " + " ".join(cmd), verbose)
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if verbose:
        print(res.stdout, file=sys.stderr)

    if res.returncode != 0:
        raise RuntimeError("solve-field failed")
    if not os.path.exists(output_wcs):
        raise RuntimeError(
            "solve-field did not find a solution (no WCS written). Narrow the search "
            "with --fovx/--fovy (scale) and --name or --ra/--dec (position), and/or use "
            "--downsample 2 on large frames.")

    # For RGB: replace green-only data in _wcs.fit with full RGB,
    # keeping the solved WCS header
    if is_rgb:
        with fits.open(output_wcs, mode="update") as hdul:
            hdul[0].data = orig_data.astype(orig_dtype)
            hdul.flush()
        print_if_verbose("RGB: restored full 3-channel data into WCS output", verbose)

    # Write human-readable scale and FOV into header
    _write_scale_fov(output_wcs, verbose)

    # Clean up temp green file
    if green_tmp and os.path.exists(green_tmp):
        try:
            os.remove(green_tmp)
        except OSError:
            pass


def cleanup_astrometry_side_products(input_fits: str, verbose: bool):
    """Remove side products emitted by astrometry.net for this input base.

    For RGB input, solve-field runs on the extracted green temp, so its side
    products carry the '<base>_green_tmp' base -- clean that base too.
    """
    try:
        base = os.path.splitext(input_fits)[0]
        exts = ["-indx.xyls", ".corr", ".wcs", ".rdls", ".match", ".solved", ".axy"]
        extras = [b + e for b in (base, base + "_green_tmp") for e in exts]
        extras.append(f"{base}_green_tmp.fit")
        for f in extras:
            if os.path.exists(f):
                os.remove(f)
    except Exception as e:
        print_if_verbose(f"Warning: cleanup failed: {e}", verbose)


def read_cd_center_scale(wcs_fit_path):
    """Read CD matrix, image center (RA,Dec), and shape from a solved WCS FITS."""
    with fits.open(wcs_fit_path) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data
    ny, nx = data.shape[-2:]
    w = WCS(hdr, naxis=2)
    if w.wcs.cd is not None and w.wcs.cd.size == 4:
        cd = np.array(w.wcs.cd, dtype=float)
    else:
        pc = np.array(w.wcs.get_pc(), dtype=float)
        cdelt = np.array(w.wcs.cdelt, dtype=float)
        cd = pc * cdelt.reshape(2, 1)
    ra_c, dec_c = w.wcs_pix2world([[nx / 2.0, ny / 2.0]], 0)[0]
    return cd, (ra_c, dec_c), (ny, nx)


def read_full_wcs_header(wcs_fit_path):
    """Return a COPY of the full primary header (including SIP, CD, CRPIX...)."""
    with fits.open(wcs_fit_path) as hdul:
        hdr = hdul[0].header.copy()
        data = hdul[0].data
    ny, nx = data.shape[-2:]
    return hdr, (ny, nx)


def build_rectified_wcs_header(wcs_fits: str,
                               center_ra_dec,
                               pixscale_arcsec,
                               base_cd: Optional[np.ndarray],
                               base_shape: Optional[Tuple[int, int]],
                               base_hdr_full: Optional[fits.Header],
                               no_rotate: bool = False):
    """
    Create target WCS header.

    If no_rotate: use the current file's own CD matrix (preserving rotation)
    but create a pure TAN projection (strip SIP distortion only).

    Otherwise:
    - If base_hdr_full is given: use it AS-IS (full WCS of the first frame).
    - Else if base_cd is given: use TAN with that CD (same scale+angle signs).
    - Else: fallback to diagonal TAN with signs/scale estimated from current file.
    """
    # --no-rotate: strip SIP distortion but keep original rotation
    if no_rotate:
        with fits.open(wcs_fits) as hdul:
            hdr = hdul[0].header
            data = hdul[0].data
        ny, nx = data.shape[-2:]
        w = WCS(hdr, naxis=2)

        # Extract CD matrix WITH rotation from solved WCS
        if w.wcs.cd is not None and w.wcs.cd.size == 4:
            cd = np.array(w.wcs.cd, dtype=float)
        else:
            pc = np.array(w.wcs.get_pc(), dtype=float)
            cdelt = np.array(w.wcs.cdelt, dtype=float)
            cd = pc * cdelt.reshape(2, 1)

        ra_c, dec_c = w.wcs_pix2world([[nx / 2.0, ny / 2.0]], 0)[0]

        # Pure TAN with original CD (rotation preserved, SIP stripped)
        w_tan = WCS(naxis=2)
        w_tan.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w_tan.wcs.crval = [ra_c, dec_c]
        w_tan.wcs.crpix = [nx / 2.0, ny / 2.0]
        w_tan.wcs.cd = cd

        hdr_tan = w_tan.to_header()
        hdr_tan["NAXIS"] = 2
        hdr_tan["NAXIS1"] = nx
        hdr_tan["NAXIS2"] = ny
        return hdr_tan, (ny, nx)

    if base_hdr_full is not None:
        hdr_tan = base_hdr_full.copy()
        nx = hdr_tan.get("NAXIS1")
        ny = hdr_tan.get("NAXIS2")
        if nx is None or ny is None:
            if base_shape is not None:
                ny, nx = base_shape
            else:
                with fits.open(wcs_fits) as hdul:
                    ny, nx = hdul[0].data.shape[-2:]
            hdr_tan["NAXIS"] = 2
            hdr_tan["NAXIS1"] = int(nx)
            hdr_tan["NAXIS2"] = int(ny)
        return hdr_tan, (int(ny), int(nx))

    with fits.open(wcs_fits) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data
    ny, nx = data.shape[-2:]
    w = WCS(hdr, naxis=2)

    if center_ra_dec is not None:
        ra_c, dec_c = center_ra_dec
    else:
        ra_c, dec_c = w.wcs_pix2world([[nx / 2.0, ny / 2.0]], 0)[0]

    if base_cd is not None:
        cd_use = np.array(base_cd, dtype=float)
        if base_shape is not None:
            ny, nx = base_shape[0], base_shape[1]
    else:
        if w.wcs.cd is not None and w.wcs.cd.size == 4:
            cd_orig = np.array(w.wcs.cd, dtype=float)
        else:
            pc = np.array(w.wcs.get_pc(), dtype=float)
            cdelt = np.array(w.wcs.cdelt, dtype=float)
            cd_orig = pc * cdelt.reshape(2, 1)
        sx = np.hypot(cd_orig[0, 0], cd_orig[1, 0])
        sy = np.hypot(cd_orig[0, 1], cd_orig[1, 1])
        scale_deg = (float(pixscale_arcsec) / 3600.0) if pixscale_arcsec else (sx + sy) / 2.0
        sign_ra = 1.0 if cd_orig[0, 0] > 0 else -1.0
        sign_dec = 1.0 if cd_orig[1, 1] > 0 else -1.0
        cd_use = np.array([[sign_ra * scale_deg, 0.0],
                           [0.0,                 sign_dec * scale_deg]])

    w_tan = WCS(naxis=2)
    w_tan.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w_tan.wcs.crval = [ra_c, dec_c]
    w_tan.wcs.crpix = [nx / 2.0, ny / 2.0]
    w_tan.wcs.cd = cd_use

    hdr_tan = w_tan.to_header()
    hdr_tan["NAXIS"] = 2
    hdr_tan["NAXIS1"] = nx
    hdr_tan["NAXIS2"] = ny
    return hdr_tan, (ny, nx)


def merge_headers_preserve_wcs(rect_hdr, src_hdr):
    """Copy non-WCS metadata from src_hdr into rect_hdr, not overwriting WCS/structural keys."""
    block_prefixes = ("CRVAL", "CRPIX", "CD", "PC", "CDELT", "CTYPE", "CROTA",
                      "PV", "A_", "B_", "AP_", "BP_", "LONPOLE", "LATPOLE", "RADESYS")
    block_exact = {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND"}
    new_hdr = rect_hdr.copy()
    for card in src_hdr.cards:
        key = card.keyword
        if key in block_exact:
            continue
        if any(key.startswith(pfx) for pfx in block_prefixes):
            continue
        if key in new_hdr:
            continue
        new_hdr.append(card)
    return new_hdr


def _sanitize_2d(data_2d):
    """Replace NaN/Inf with fill value in a 2D float32 array. Returns sanitized copy."""
    data_2d = np.asarray(data_2d, dtype=np.float32)
    finite_mask = np.isfinite(data_2d)
    if not finite_mask.any():
        return None  # signal: all NaN
    if not finite_mask.all():
        fill = float(np.median(data_2d[finite_mask]))
        data_2d[~finite_mask] = fill
    return data_2d


def _reproject_2d(data_2d, wcs_in, wcs_out, ny, nx, rect_method, verbose):
    """Reproject a single 2D plane using interp or exact method."""
    try:
        if wcs_in.to_header(relax=True) == wcs_out.to_header(relax=True):
            return data_2d.copy()
    except Exception:
        pass

    try:
        if rect_method == "interp":
            result, _ = reproject_interp(
                (data_2d, wcs_in), wcs_out,
                shape_out=(ny, nx), return_footprint=True)
        else:
            result, _ = reproject_exact(
                (data_2d, wcs_in), wcs_out,
                shape_out=(ny, nx), return_footprint=True)
    except Exception as e:
        print_if_verbose(
            f"Warning: reproject_{rect_method} failed ({e}); "
            f"falling back to reproject_interp.", verbose)
        result, _ = reproject_interp(
            (data_2d, wcs_in), wcs_out,
            shape_out=(ny, nx), return_footprint=True)
    return result


def _wcsresample_2d(input_2d_fits, rect_wcs_path, output_path, use_wsl, verbose):
    """Run wcs-resample on a single 2D FITS. Returns True on success."""
    if use_wsl:
        inp = win_to_wsl_path(input_2d_fits)
        wcs_p = win_to_wsl_path(rect_wcs_path)
        out = win_to_wsl_path(output_path)
        cmd = ["wsl", "wcs-resample", inp, wcs_p, out]
    else:
        cmd = ["wcs-resample", input_2d_fits, rect_wcs_path, output_path]

    print_if_verbose("Running: " + " ".join(cmd), verbose)
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if verbose:
        print(res.stdout, file=sys.stderr)
    return res.returncode == 0 and os.path.exists(output_path)


def run_wcs_resample(wcs_fits: str,
                     output_rect: str,
                     center_ra_dec,
                     pixscale_arcsec,
                     use_wsl: bool,
                     verbose: bool,
                     base_cd: Optional[np.ndarray],
                     base_shape: Optional[Tuple[int, int]],
                     base_hdr_full: Optional[fits.Header],
                     rect_method: str,
                     no_rotate: bool = False) -> None:
    """Resample into target grid using chosen method, then merge metadata.

    Supports both 2D (mono) and 3D (RGB) input in wcs_fits.
    For RGB: reprojects each channel separately, then stacks.
    If no_rotate: strips SIP distortion but preserves field rotation.
    """
    rect_wcs_hdr, (ny, nx) = build_rectified_wcs_header(
        wcs_fits, center_ra_dec=center_ra_dec, pixscale_arcsec=pixscale_arcsec,
        base_cd=base_cd, base_shape=base_shape, base_hdr_full=base_hdr_full,
        no_rotate=no_rotate
    )

    # Read input data and detect RGB
    with fits.open(wcs_fits, memmap=False) as src_hdul:
        data_raw = src_hdul[0].data
        src_wcs_hdr = src_hdul[0].header

    data_raw, is_rgb = _normalize_rgb_data(data_raw)

    if rect_method == "wcsresample":
        # Create target WCS template (always 2D)
        rect_wcs_path = os.path.splitext(output_rect)[0] + "_rectwcs_tmp.fit"
        fits.PrimaryHDU(data=np.zeros((ny, nx), dtype=np.float32),
                        header=rect_wcs_hdr).writeto(rect_wcs_path, overwrite=True)

        if is_rgb:
            # Per-channel wcs-resample
            ch_results = []
            base_tmp = os.path.splitext(output_rect)[0]
            for ch_idx, ch_name in enumerate(["R", "G", "B"]):
                ch_tmp_in = f"{base_tmp}_ch{ch_name}_tmp.fit"
                ch_tmp_out = f"{base_tmp}_ch{ch_name}_rect_tmp.fit"
                # Write mono channel with WCS header
                ch_hdr = src_wcs_hdr.copy()
                if 'NAXIS3' in ch_hdr:
                    del ch_hdr['NAXIS3']
                ch_hdr['NAXIS'] = 2
                fits.PrimaryHDU(data=data_raw[ch_idx].astype(np.float32),
                                header=ch_hdr).writeto(ch_tmp_in, overwrite=True)
                if not _wcsresample_2d(ch_tmp_in, rect_wcs_path, ch_tmp_out, use_wsl, verbose):
                    raise RuntimeError(f"wcs-resample failed for {ch_name} channel")
                with fits.open(ch_tmp_out, memmap=False) as h:
                    ch_results.append(h[0].data.astype(np.float32))
                # Clean up temps
                for f in (ch_tmp_in, ch_tmp_out):
                    try:
                        os.remove(f)
                    except OSError:
                        pass

            data_out = np.stack(ch_results, axis=0)
        else:
            # Mono: single wcs-resample call (original behavior)
            if use_wsl:
                inp = win_to_wsl_path(wcs_fits)
                wcs_p = win_to_wsl_path(rect_wcs_path)
                out = win_to_wsl_path(output_rect)
                cmd = ["wsl", "wcs-resample", inp, wcs_p, out]
            else:
                cmd = ["wcs-resample", wcs_fits, rect_wcs_path, output_rect]

            print_if_verbose("Running: " + " ".join(cmd), verbose)
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if verbose:
                print(res.stdout, file=sys.stderr)

            if res.returncode != 0:
                raise RuntimeError("wcs-resample failed")
            if not os.path.exists(output_rect):
                raise FileNotFoundError(f"Expected rectified FITS not found: {output_rect}")
            data_out = None  # already written by wcs-resample

        try:
            os.remove(rect_wcs_path)
        except OSError:
            pass

        # For RGB: write stacked result
        if is_rgb and data_out is not None:
            fits.PrimaryHDU(data=data_out,
                            header=rect_wcs_hdr).writeto(output_rect, overwrite=True)

    else:
        # interp / exact methods
        if not _have_reproject:
            raise RuntimeError("reproject is not installed. Install with: pip install reproject")

        wcs_in = WCS(src_wcs_hdr, naxis=2)
        wcs_out = WCS(rect_wcs_hdr, naxis=2)

        if is_rgb:
            # Per-channel reprojection
            ch_results = []
            for ch_idx in range(3):
                ch_in = data_raw[ch_idx].astype(np.float32)
                ch_out = _reproject_2d(ch_in, wcs_in, wcs_out, ny, nx, rect_method, verbose)
                ch_out = _sanitize_2d(ch_out)
                if ch_out is None:
                    print_if_verbose(
                        f"Warning: reproject produced only NaNs for channel {ch_idx}; "
                        f"copying original.", verbose)
                    ch_out = ch_in.copy()
                ch_results.append(ch_out)
            data_out = np.stack(ch_results, axis=0)
        else:
            # Mono reprojection (original behavior)
            data_in = data_raw.astype(np.float32)
            data_out = _reproject_2d(data_in, wcs_in, wcs_out, ny, nx, rect_method, verbose)
            data_out = _sanitize_2d(data_out)
            if data_out is None:
                print_if_verbose(
                    f"Warning: reproject produced only NaNs for {os.path.basename(wcs_fits)}; "
                    f"copying original data.", verbose)
                data_out = data_in.copy()

        fits.PrimaryHDU(data=data_out,
                        header=rect_wcs_hdr).writeto(output_rect, overwrite=True)

    # Merge non-WCS metadata from source into output
    with fits.open(wcs_fits) as src, fits.open(output_rect, mode="update") as dst:
        new_hdr = merge_headers_preserve_wcs(dst[0].header, src[0].header)
        dst[0].header = new_hdr
        dst.flush()


# ---------------------- WCS refit using .corr ----------------------
def _load_corr_matches(corr_path: str):
    """Return (xy: Nx2 pixels, sky: SkyCoord) from a .corr file (FITS bin table or ASCII)."""
    try:
        with fits.open(corr_path) as hdul:
            hdu = None
            for h in hdul:
                if isinstance(h, fits.BinTableHDU):
                    hdu = h
                    break
            if hdu is not None:
                cols = {c.name.lower(): c.name for c in hdu.columns}

                def pick(*names):
                    for nm in names:
                        if nm in cols:
                            return cols[nm]
                    return None

                cx = pick("field_x", "x", "x_image")
                cy = pick("field_y", "y", "y_image")
                cra = pick("index_ra", "ra", "ref_ra", "ra_icrs")
                cdec = pick("index_dec", "dec", "ref_dec", "dec_icrs")
                if cx and cy and cra and cdec:
                    x = np.asarray(hdu.data[cx], dtype=float)
                    y = np.asarray(hdu.data[cy], dtype=float)
                    ra = np.asarray(hdu.data[cra], dtype=float)
                    dec = np.asarray(hdu.data[cdec], dtype=float)
                    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(ra) & np.isfinite(dec)
                    x, y, ra, dec = x[m], y[m], ra[m], dec[m]
                    if x.size >= 6:
                        xy = np.vstack([x, y]).T
                        sky = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
                        return xy, sky
    except Exception:
        pass

    try:
        arr = np.loadtxt(corr_path, comments="#")
        if arr.ndim == 1:
            arr = arr[None, :]
        if arr.shape[1] >= 4:
            x = arr[:, 0].astype(float)
            y = arr[:, 1].astype(float)
            ra = arr[:, -2].astype(float)
            dec = arr[:, -1].astype(float)
            m = np.isfinite(x) & np.isfinite(y) & np.isfinite(ra) & np.isfinite(dec)
            x, y, ra, dec = x[m], y[m], ra[m], dec[m]
            if x.size >= 6:
                xy = np.vstack([x, y]).T
                sky = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
                return xy, sky
    except Exception:
        pass

    return None, None


def refit_wcs_from_corr(input_fits: str,
                        wcs_fits: str,
                        sip_degree: int,
                        verbose: bool) -> bool:
    """
    Refit WCS using matched pairs from <input>.corr.
    Custom least-squares TAN (gnomonic) fit (affine) without fit_wcs_from_points.
    Returns True if refit applied, otherwise False.
    """
    corr = os.path.splitext(input_fits)[0] + ".corr"
    if not os.path.exists(corr):
        print_if_verbose(f"Refit skipped: not found {corr}", verbose)
        return False

    xy, sky = _load_corr_matches(corr)
    if xy is None or sky is None or xy.shape[0] < 6:
        print_if_verbose(f"Refit skipped: cannot parse enough matches in {corr}", verbose)
        return False

    try:
        sky_icrs = sky.icrs
        ra = sky_icrs.ra.to(u.radian).value
        dec = sky_icrs.dec.to(u.radian).value

        ra0 = np.median(ra)
        dec0 = np.median(dec)

        dra = ra - ra0
        dra = (dra + np.pi) % (2.0 * np.pi) - np.pi

        sin_dec = np.sin(dec)
        cos_dec = np.cos(dec)
        sin_dec0 = np.sin(dec0)
        cos_dec0 = np.cos(dec0)

        cosc = sin_dec0 * sin_dec + cos_dec0 * cos_dec * np.cos(dra)
        cosc[cosc == 0] = 1e-12

        xi = cos_dec * np.sin(dra) / cosc
        eta = (cos_dec0 * sin_dec - sin_dec0 * cos_dec * np.cos(dra)) / cosc

        xi_deg = np.degrees(xi)
        eta_deg = np.degrees(eta)

        x = xy[:, 0]
        y = xy[:, 1]

        A = np.stack([x, y, np.ones_like(x)], axis=1)
        coeff_xi, _, _, _ = np.linalg.lstsq(A, xi_deg, rcond=None)
        coeff_eta, _, _, _ = np.linalg.lstsq(A, eta_deg, rcond=None)

        M = np.vstack([coeff_xi[:2], coeff_eta[:2]])
        c = np.array([coeff_xi[2], coeff_eta[2]])

        try:
            Minv = np.linalg.inv(M)
        except np.linalg.LinAlgError:
            print_if_verbose("Refit failed: singular transform matrix", verbose)
            return False

        crpix_vec = -Minv @ c
        crpix = (float(crpix_vec[0]), float(crpix_vec[1]))

        with fits.open(wcs_fits, mode="update") as hdul:
            data = hdul[0].data
            orig_hdr = hdul[0].header.copy()
            ny, nx = data.shape[-2:]

            w = WCS(naxis=2)
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            w.wcs.crval = [np.degrees(ra0), np.degrees(dec0)]
            w.wcs.crpix = [crpix[0], crpix[1]]
            w.wcs.cd = M

            new_hdr = w.to_header(relax=True)
            new_hdr["NAXIS"] = 2
            new_hdr["NAXIS1"] = nx
            new_hdr["NAXIS2"] = ny

            merged = merge_headers_preserve_wcs(new_hdr, orig_hdr)
            hdul[0].header = merged
            hdul.flush()

        print_if_verbose(
            f"Refit applied from {os.path.basename(corr)} "
            f"(TAN, affine, matches={xy.shape[0]})", verbose
        )
        return True

    except Exception as e:
        print_if_verbose(f"Refit failed: {e}", verbose)
        return False


# ---------------------- fine alignment (subpixel, no ringing) ----------------------
def _hann2d(h, w):
    """2D Hann window to reduce edge effects in FFT-based registration."""
    wy = np.hanning(h)
    wx = np.hanning(w)
    return np.outer(wy, wx)


def _despike_inplace(a: np.ndarray, thresh: float = 6.0):
    """Replace strong outliers (hot pixels / CR hits) by local 3x3 median."""
    if not _have_scipy:
        return
    med = np.median(a)
    mad = np.median(np.abs(a - med)) + 1e-6
    z = (a - med) / (1.4826 * mad)
    mask = np.abs(z) > thresh
    if mask.any():
        mf = ndi_median_filter(a, size=3)
        a[mask] = mf[mask]


def _phase_corr_shift(ref, img):
    """
    Estimate (dy, dx) shift that aligns 'img' to 'ref' using phase correlation.
    Returns subpixel shift in pixels (positive dy=down, dx=right).
    """
    ref = np.asarray(ref, dtype=np.float32)
    img = np.asarray(img, dtype=np.float32)
    h, w = ref.shape
    win = _hann2d(h, w)
    R = np.fft.rfftn(ref * win, axes=(0, 1))
    I = np.fft.rfftn(img * win, axes=(0, 1))

    CPS = R * np.conj(I)
    denom = np.abs(CPS)
    denom = np.where(denom == 0, 1.0, denom)
    CPS = np.nan_to_num(CPS / denom, nan=0.0, posinf=0.0, neginf=0.0)
    cc = np.fft.irfftn(CPS, s=ref.shape, axes=(0, 1))

    maxpos = np.unravel_index(np.argmax(cc), cc.shape)
    peak_y, peak_x = maxpos

    def quad_refine(values):
        a = (values[2] + values[0] - 2 * values[1]) / 2.0
        b = (values[2] - values[0]) / 2.0
        if a == 0:
            return 0.0
        return -b / (2 * a)

    def wrap_to_signed(p, n):
        return p if p <= n // 2 else p - n

    dy_i = wrap_to_signed(peak_y, h)
    dx_i = wrap_to_signed(peak_x, w)

    y_prev = (peak_y - 1) % h
    y_next = (peak_y + 1) % h
    x_prev = (peak_x - 1) % w
    x_next = (peak_x + 1) % w

    sub_y = quad_refine([cc[y_prev, peak_x], cc[peak_y, peak_x], cc[y_next, peak_x]])
    sub_x = quad_refine([cc[peak_y, x_prev], cc[peak_y, peak_x], cc[peak_y, x_next]])

    return dy_i + sub_y, dx_i + sub_x


def _spline_shift_2d(img, dy, dx):
    """Subpixel shift using cubic spline (order=3), no ringing, no wrap-around."""
    if not _have_scipy:
        raise RuntimeError("SciPy is required for --fine-align (pip install scipy)")
    return ndi_shift(img, shift=(dy, dx), order=3, mode="nearest", prefilter=True).astype(np.float32, copy=False)


def fine_align_inplace(rect_paths: List[str], verbose: bool):
    """
    Load first rect image as reference; align all subsequent rect images to it.
    Shift data by subpixel (cubic spline). Keep headers unchanged.

    For RGB: computes shift from green channel only, applies to all 3 channels.
    """
    if len(rect_paths) <= 1:
        return

    with fits.open(rect_paths[0], memmap=False) as ref_hdul:
        ref_data = ref_hdul[0].data.astype(np.float32)

    ref_data, ref_is_rgb = _normalize_rgb_data(ref_data)

    # Use green channel for correlation (or mono as-is)
    if ref_is_rgb:
        ref_plane = ref_data[1].copy()  # G channel
    else:
        ref_plane = ref_data.copy()
    ref_plane = ref_plane - np.median(ref_plane)
    _despike_inplace(ref_plane)

    for p in rect_paths[1:]:
        try:
            with fits.open(p, mode="update") as hdul:
                raw = hdul[0].data.astype(np.float32)
                raw, is_rgb = _normalize_rgb_data(raw)
                orig_dtype = hdul[0].data.dtype

                # Compute shift from green channel (or mono)
                if is_rgb:
                    corr_plane = raw[1].copy()
                else:
                    corr_plane = raw.copy()
                corr_norm = corr_plane - np.median(corr_plane)
                _despike_inplace(corr_norm)
                dy, dx = _phase_corr_shift(ref_plane, corr_norm)

                # Apply shift to all channels
                if is_rgb:
                    channels = []
                    for ch in range(3):
                        ch_data = raw[ch]
                        ch_med = np.median(ch_data)
                        ch_norm = ch_data - ch_med
                        ch_shifted = _spline_shift_2d(ch_norm, dy, dx)
                        ch_shifted += ch_med
                        channels.append(ch_shifted)
                    shifted = np.stack(channels, axis=0)
                else:
                    med = np.median(raw)
                    raw_norm = raw - med
                    shifted = _spline_shift_2d(raw_norm, dy, dx)
                    shifted += med

                if np.issubdtype(orig_dtype, np.integer):
                    info = np.iinfo(orig_dtype)
                    out = np.clip(np.round(shifted), info.min, info.max).astype(orig_dtype)
                else:
                    out = shifted.astype(orig_dtype)

                hdul[0].data = out
                hdul.flush()
            print_if_verbose(f"Fine-align {os.path.basename(p)} : dy={dy:.3f} dx={dx:.3f}", verbose)
        except Exception as e:
            print_if_verbose(f"Fine-align failed for {p}: {e}", verbose)


# ---------------------- input expansion (AstroBatch style) ----------------------
def expand_io_pairs(input_spec: str,
                    output_spec: Optional[str]) -> List[Tuple[str, str]]:
    """
    Build list of (infile, base_out) pairs.
    If output_spec is provided, use batch_utils.build_io_file_lists() and take its outfile as base.
    If output_spec is None:
        - expand input_spec to a list of files (wildcards or numbered sequence or single file),
        - use base_out == infile path (outputs go near inputs with suffixes).
    """
    if output_spec:
        return bu.build_io_file_lists(input_spec, output_spec)

    pairs = []
    if bu.has_wildcards(input_spec):
        matched = sorted(f for f in glob.glob(input_spec) if os.path.isfile(f))
        if not matched:
            raise FileNotFoundError(f"No input files match pattern '{input_spec}'.")
        for f in matched:
            pairs.append((os.path.abspath(f), os.path.abspath(f)))
        return pairs

    seq = bu.find_numbered_sequence(input_spec)
    if seq:
        for (name, _, _) in seq:
            pairs.append((os.path.abspath(name), os.path.abspath(name)))
        return pairs

    if not os.path.exists(input_spec):
        raise FileNotFoundError(f"Input file '{input_spec}' not found.")
    f = os.path.abspath(input_spec)
    pairs.append((f, f))
    return pairs


# ---------------------- utilities ----------------------
def get_center_from_wcs(wcs_fit_path: str) -> Tuple[float, float]:
    """Return (RA, Dec) in degrees at image center from solved WCS FITS."""
    with fits.open(wcs_fit_path) as hdul:
        data = hdul[0].data
        hdr = hdul[0].header
    ny, nx = data.shape[-2:]
    w = WCS(hdr, naxis=2)
    ra_c, dec_c = w.wcs_pix2world([[nx / 2.0, ny / 2.0]], 0)[0]
    return float(ra_c), float(dec_c)


# ---------------------- main batch driver ----------------------
def main():
    parser = argparse.ArgumentParser(
        description="Astrometric solve (WSL-aware) + optional WCS refit + optional rectification "
                    "(+subpixel fine align). Batch-friendly. Accepts FITS or JPEG input; a JPEG "
                    "input yields a JPEG output with the solved WCS packed into its header."
    )
    parser.add_argument("input", help="Input spec: file | numbered pattern | wildcard mask (FITS or JPEG)")
    parser.add_argument("output", nargs="?", default=None,
                        help="Optional output. If given, the solved result takes this EXACT name; "
                             "if omitted, outputs go near the inputs with a _wcs (and _rect) suffix.")

    parser.add_argument("-r", "--rectify", action="store_true",
                        help="Produce <base>_rect.fit (reprojection). If omitted, only <base>_wcs.fit is created.")
    parser.add_argument("-i", "--individual", action="store_true",
                        help="Per-file rectification centers (ignore 'first file as base' logic).")
    parser.add_argument("--rect-center-ra", type=float, help="Rectified center RA (deg).")
    parser.add_argument("--rect-center-dec", type=float, help="Rectified center Dec (deg).")
    parser.add_argument("--rect-pixscale", type=float, help="Rectified pixel scale (arcsec/pix).")
    parser.add_argument("--rect-method", choices=["wcsresample", "interp", "exact"], default="wcsresample",
                        help="Resampling engine: wcsresample (fast), interp (bilinear), exact (highest fidelity).")

    parser.add_argument("--no-rotate", action="store_true",
                        help="Rectify without derotation: strip SIP distortion but preserve "
                             "field rotation. Image dimensions unchanged, no black borders. "
                             "Useful for mosaics, surveys, and tile pyramids.")
    parser.add_argument("--fine-align", action="store_true",
                        help="After rectification, subpixel-align all _rect.fit to the first (phase correlation + cubic spline shift).")

    parser.add_argument("--refit-wcs", action="store_true",
                        help="After solve-field, refit WCS from <input>.corr via least-squares and rewrite <base>_wcs.fit.")
    parser.add_argument("--refit-sip-degree", type=int, default=2,
                        help="SIP polynomial degree for WCS refit (1=affine). Default: 2 (currently only affine used).")

    parser.add_argument("--tweak-order", type=int, default=3, help="SIP polynomial order for solve-field.")

    # Solve hints (useful for header-less inputs such as processed JPEGs).
    parser.add_argument("--fovx", type=float, metavar="DEG",
                        help="Field-of-view WIDTH in degrees; the pixel scale is computed from it "
                             "(alternative to reading FOCALLEN from the header). Use --fovx OR --fovy.")
    parser.add_argument("--fovy", type=float, metavar="DEG",
                        help="Field-of-view HEIGHT in degrees (alternative to --fovx).")
    parser.add_argument("--ra", metavar="RA",
                        help="Search-centre RA. Decimal = DEGREES (e.g. 193.39); sexagesimal = HOURS "
                             "(e.g. \"12 53 32\" = 193.39 deg). Requires --dec.")
    parser.add_argument("--dec", metavar="DEC",
                        help="Search-centre Dec in DEGREES: decimal (e.g. -22.87) or sexagesimal "
                             "(e.g. \"-22 52 23\"). Requires --ra.")
    parser.add_argument("--name", metavar="ID",
                        help="Resolve the search centre from an object name via CDS Sesame/SIMBAD "
                             "(e.g. --name \"PN A66 35\"). Explicit --ra/--dec take precedence. Note: "
                             "\"Abell 35\" is the galaxy cluster; the planetary nebula is \"PN A66 35\".")
    parser.add_argument("--radius", type=float, metavar="DEG",
                        help="Search radius in degrees around the centre hint (default: from FOV, else wide).")
    parser.add_argument("--downsample", type=int, default=1, metavar="N",
                        help="solve-field source-extraction downsampling (default 1). 2-4 speeds up large frames.")

    parser.add_argument("--jpeg", type=int, nargs="?", const=99, default=99, metavar="Q",
                        help="JPEG output quality 1-100 (default 99). Applies when the input is a "
                             "JPEG: output is <base>_wcs.jpg (and <base>_rect.jpg with -r) with the "
                             "solved WCS packed into a COM marker (shared fits2tiff format).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output (show solver logs).")

    args = parser.parse_args()

    try:
        io_pairs = expand_io_pairs(args.input, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Resolve a one-time centre hint from the CLI (applies to every frame).
    cli_ra = cli_dec = None
    if (args.ra is not None) or (args.dec is not None):
        if (args.ra is None) or (args.dec is None):
            print("Error: --ra and --dec must be given together.", file=sys.stderr)
            sys.exit(1)
        try:
            cli_ra = _parse_angle_cli(args.ra, is_ra=True)
            cli_dec = _parse_angle_cli(args.dec, is_ra=False)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.name:
        try:
            cli_ra, cli_dec = resolve_name_sesame(args.name)
        except RuntimeError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        print(f"Resolved '{args.name}' -> RA {cli_ra:.5f}, Dec {cli_dec:+.5f} (deg)", file=sys.stderr)

    for fov_name, fov_val in (("--fovx", args.fovx), ("--fovy", args.fovy)):
        if fov_val is not None and not (np.isfinite(fov_val) and fov_val > 0):
            print(f"Error: {fov_name} must be a positive number of degrees.", file=sys.stderr)
            sys.exit(1)
    if args.fovx is not None and args.fovy is not None:
        print("Note: both --fovx and --fovy given; using --fovx.", file=sys.stderr)

    use_wsl = is_windows()
    n = len(io_pairs)
    ok = 0
    fail = 0

    use_shared_center = False
    shared_center = None
    base_cd = None
    base_shape = None
    base_hdr_full = None
    rect_paths_for_fine = []
    jpeg_rect_tasks = []   # (rect_fits, out_jpg, quality): converted after fine-align

    if args.rectify:
        if args.no_rotate:
            # --no-rotate: each file keeps its own orientation, no shared grid
            use_shared_center = False
        elif (args.rect_center_ra is not None) and (args.rect_center_dec is not None):
            shared_center = (args.rect_center_ra, args.rect_center_dec)
            use_shared_center = True
        elif not args.individual:
            use_shared_center = True

    for idx, (infile, outbase) in enumerate(io_pairs, start=1):
        progress_bar(idx - 1, n)
        temp_solve_fits = None
        output_wcs = None
        solve_input = None
        try:
            abs_out = os.path.abspath(outbase)
            base_out_noext, _ = os.path.splitext(abs_out)
            infile_is_jpeg = _is_jpeg(infile)

            # The _wcs/_rect suffix is only appended when the output name is
            # auto-derived. With an explicit output, its exact name goes to the
            # primary product (rectified if -r, else the solved image); the other
            # / intermediate files keep derived names.
            explicit_out = args.output is not None
            if explicit_out and (not infile_is_jpeg) and (not args.rectify):
                output_wcs = abs_out                      # FITS deliverable = exact name
            else:
                output_wcs = base_out_noext + "_wcs.fit"
            if explicit_out and (not infile_is_jpeg) and args.rectify:
                output_rect = abs_out                     # rectified FITS = exact name
            else:
                output_rect = base_out_noext + "_rect.fit"

            # JPEG input: decode to a temp FITS (carrying recovered hints) and
            # solve that; the output is converted back to JPEG below.
            solve_input = infile
            if infile_is_jpeg:
                temp_solve_fits = base_out_noext + "_solvetmp.fit"
                jpeg_to_temp_fits(infile, temp_solve_fits)
                solve_input = temp_solve_fits

            run_solve_field(
                input_fits=solve_input,
                output_wcs=output_wcs,
                tweak_order=args.tweak_order,
                use_wsl=use_wsl,
                verbose=args.verbose,
                ra_override=cli_ra,
                dec_override=cli_dec,
                radius_override=args.radius,
                fovx=args.fovx,
                fovy=args.fovy,
                downsample=args.downsample,
            )

            if args.refit_wcs:
                applied = refit_wcs_from_corr(
                    input_fits=solve_input,
                    wcs_fits=output_wcs,
                    sip_degree=args.refit_sip_degree,
                    verbose=args.verbose
                )
                if applied and args.rectify and (not args.individual) and (not args.no_rotate) and base_hdr_full is None and \
                   (args.rect_center_ra is None and args.rect_center_dec is None and args.rect_pixscale is None):
                    try:
                        base_hdr_full, base_shape = read_full_wcs_header(output_wcs)
                    except Exception as e:
                        print_if_verbose(f"Warning: cannot read base full WCS from {output_wcs}: {e}", args.verbose)

            if args.rectify and use_shared_center and shared_center is None and \
               not ((args.rect_center_ra is not None) and (args.rect_center_dec is not None)):
                try:
                    shared_center = get_center_from_wcs(output_wcs)
                except Exception as e:
                    print_if_verbose(f"Warning: could not read center from {output_wcs}: {e}", args.verbose)

            if args.rectify and (not args.individual) and (not args.no_rotate) and base_hdr_full is None and \
               (args.rect_center_ra is None and args.rect_center_dec is None and args.rect_pixscale is None):
                try:
                    base_hdr_full, base_shape = read_full_wcs_header(output_wcs)
                except Exception as e:
                    print_if_verbose(f"Warning: cannot read base full WCS from {output_wcs}: {e}", args.verbose)

            if args.rectify:
                center = shared_center if use_shared_center else None
                run_wcs_resample(
                    wcs_fits=output_wcs,
                    output_rect=output_rect,
                    center_ra_dec=center,
                    pixscale_arcsec=args.rect_pixscale,
                    use_wsl=use_wsl,
                    verbose=args.verbose,
                    base_cd=base_cd,
                    base_shape=base_shape,
                    base_hdr_full=base_hdr_full,
                    rect_method=args.rect_method,
                    no_rotate=args.no_rotate
                )
                rect_paths_for_fine.append(output_rect)

            cleanup_astrometry_side_products(solve_input, verbose=args.verbose)

            # JPEG input -> JPEG output(s) with embedded WCS. The solve product
            # (_wcs) is converted now; the rectified product is deferred until
            # after fine-align (which edits _rect.fit in place).
            if infile_is_jpeg:
                q = max(1, min(100, args.jpeg))
                if args.rectify:
                    # Primary = rectified jpg (exact name if explicit); the solved
                    # jpg is a secondary byproduct with a derived name.
                    wcs_jpg = base_out_noext + "_wcs.jpg"
                    rect_jpg = abs_out if explicit_out else base_out_noext + "_rect.jpg"
                    jpeg_rect_tasks.append((output_rect, rect_jpg, q))
                else:
                    wcs_jpg = abs_out if explicit_out else base_out_noext + "_wcs.jpg"
                fits_to_jpeg(output_wcs, wcs_jpg, q)
                _safe_remove(output_wcs)
                _safe_remove(temp_solve_fits)

            ok += 1
            progress_bar(idx, n)
            print(f"  OK: {os.path.basename(infile)}" + (' ' * 30), file=sys.stderr)

        except Exception as e:
            fail += 1
            # Don't leave astrometry side products behind when a later step
            # (e.g. rectify) fails; for JPEG input also drop the temp + orphan _wcs.fit.
            if solve_input:
                cleanup_astrometry_side_products(solve_input, verbose=args.verbose)
            if temp_solve_fits is not None:
                _safe_remove(output_wcs)
                _safe_remove(temp_solve_fits)
            progress_bar(idx, n)
            print(f" FAIL: {os.path.basename(infile)}  [{e}]" + (' ' * 30), file=sys.stderr)

    if args.rectify and args.fine_align and len(rect_paths_for_fine) > 1:
        if not _have_scipy:
            print("Warning: SciPy not installed, skipping --fine-align. Install with: pip install scipy", file=sys.stderr)
        else:
            print_if_verbose("Starting fine alignment (phase correlation + cubic spline shift)...", args.verbose)
            fine_align_inplace(rect_paths_for_fine, verbose=args.verbose)

    # JPEG rectified outputs: convert after fine-align has finished editing _rect.fit.
    for (rect_fits, out_jpg, q) in jpeg_rect_tasks:
        try:
            fits_to_jpeg(rect_fits, out_jpg, q)
            _safe_remove(rect_fits)
        except Exception as e:
            print(f"Warning: JPEG rectified output failed for {os.path.basename(rect_fits)}: {e}", file=sys.stderr)

    progress_bar(n, n)
    print("", file=sys.stderr)
    print(f"Done. OK: {ok}, FAIL: {fail}", file=sys.stderr)


if __name__ == "__main__":
    main()
