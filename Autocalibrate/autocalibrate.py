#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import glob
import traceback
from datetime import datetime, timezone, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from scipy.ndimage import median_filter, gaussian_filter
from astropy.io import fits

# Local paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

LIB_DIR = os.path.join(SCRIPT_DIR, "..", "lib")
if os.path.isdir(LIB_DIR) and LIB_DIR not in sys.path:
    sys.path.insert(0, LIB_DIR)

COSME_DIR = os.path.join(SCRIPT_DIR, "..", "Cosme")
if os.path.isdir(COSME_DIR) and COSME_DIR not in sys.path:
    sys.path.insert(0, COSME_DIR)

import batch_utils  # wildcard helpers etc.

# Cosmetic correction (optional)
try:
    from cosme import read_cosme_list, apply_cosmetic_correction
    HAVE_COSME = True
except Exception:
    HAVE_COSME = False

FITS_EXTS = (".fit", ".fits", ".fts")

# =============================================================================
# BESTFLAT ALGORITHM CONSTANTS
# =============================================================================
#
# Best flat selection algorithm (--bestflat mode):
#
# 1. PREPARE LIGHT PREVIEWS (per session):
#    - Load each light frame, subtract dark
#    - Downscale by BESTFLAT_DOWNSCALE_FACTOR (8x8 averaging)
#    - Apply median filter (size=BESTFLAT_LIGHT_MEDIAN_SIZE) to remove stars
#    - Combine all session lights via pixel-wise median
#
# 2. PREPARE FLAT PREVIEWS (cached):
#    - Downscale by BESTFLAT_DOWNSCALE_FACTOR (no median filter - preserve dust!)
#    - Normalize to median 5000 for consistent comparison
#
# 3. EVALUATE EACH CANDIDATE FLAT:
#    - Divide: light_preview / flat_preview
#    - Apply high-pass filter (removes gradients, keeps dust/artifacts)
#      Sigma = BESTFLAT_HIGHPASS_SIGMA_FACTOR * mean(width, height)
#    - Apply gaussian blur (size=BESTFLAT_SIGMA_BLUR_SIZE) to suppress noise
#    - Compute sigma (std deviation) of the result
#
# 4. SELECT BEST:
#    - Flat with lowest sigma = best dust/vignetting match
#
# =============================================================================

BESTFLAT_DOWNSCALE_FACTOR = 8       # Downscale factor for previews (8 = 8x8 averaging)
BESTFLAT_LIGHT_MEDIAN_SIZE = 16     # Median filter size for lights (removes stars)
BESTFLAT_HIGHPASS_SIGMA_FACTOR = 0.3  # High-pass gaussian sigma as fraction of image size
BESTFLAT_SIGMA_BLUR_SIZE = 5        # Gaussian blur before sigma calculation (removes HF noise)
BESTFLAT_DEBUG_SAVE_BLURRED = True  # Save debug files with blur applied (same as sigma input)
BESTFLAT_DEBUG_CALIBRATION_LINE = True  # Add calibration line (0 | 32767) to fix viewer scaling


def print_usage():
    print(
        "Usage:\n"
        "  autocalibrate.py [options] rawfiles.fit out_path dark_path flat_path\n\n"
        "Options:\n"
        "  --bestflat   Auto-select best flat per observing session by minimizing\n"
        "               stddev after flat division. Sessions are grouped by filter\n"
        "               and local noon-to-noon periods.\n"
        "  --debug      Save downscaled test images to ./debug/ folder for inspection.\n"
        "               Creates subfolders per session with divided previews.\n"
        "  --flat-future-days N\n"
        "               Max days in future to accept a flat (default: 2).\n"
        "               Allows using flats taken shortly after observations\n"
        "               (e.g., after clear sky session before new calibration).\n"
        "               Example: --flat-future-days 3\n"
        "  --flatlog FILE\n"
        "               CSV file with flat renewal timestamps.\n"
        "               Format: DATETIME_UTC,CAMERA_ID[,COMMENT]\n"
        "               Example: 2024-05-18T14:30:00,2600MM\n"
        "               When specified, prioritizes flats from current interval:\n"
        "               - First searches current interval (between log entries)\n"
        "               - If multiple flats in interval, picks nearest within +N days\n"
        "               - If none within +N days but interval has flats, uses latest\n"
        "               - If interval empty, falls back to global search (+N days)\n"
        "               Camera ID is matched as substring against INSTRUME header.\n\n"
        "rawfiles.fit : raw input spec (file, list file, or wildcard mask)\n"
        "out_path     : base output file name pattern (e.g. out.fit or path/out)\n"
        "               Result names:\n"
        "                 <base>_exp<EXPTIME>_<FILTER>_<N>.ext\n"
        "               where N is counted separately for each (EXPTIME,FILTER).\n"
        "               If extension is missing, .fit is used.\n"
        "dark_path    : directory tree with master darks and cosmetic *.lst\n"
        "flat_path    : directory tree with master flats\n\n"
        "Calibration formula:\n"
        "  RESULT = ((RAW - DARK) * MULT) / FLAT\n"
        "where MULT = median(FLAT), computed automatically.\n"
        "Cosmetic correction applied after dark subtraction,\n"
        "then flat division, computed in float64,\n"
        "then clamped back to source data type range.\n"
    )


def parse_args(argv):
    bestflat = False
    debug = False
    flat_future_days = 2.0
    flatlog = None
    args = argv[1:]

    # Parse options
    while args and args[0].startswith("--"):
        if args[0] == "--bestflat":
            bestflat = True
            args = args[1:]
        elif args[0] == "--debug":
            debug = True
            args = args[1:]
        elif args[0] == "--flat-future-days":
            if len(args) < 2:
                print("ERROR: --flat-future-days requires a value")
                sys.exit(1)
            try:
                flat_future_days = float(args[1])
            except ValueError:
                print(f"ERROR: --flat-future-days must be a number, got: {args[1]}")
                sys.exit(1)
            args = args[2:]
        elif args[0] == "--flatlog":
            if len(args) < 2:
                print("ERROR: --flatlog requires a file path")
                sys.exit(1)
            flatlog = args[1]
            if not os.path.isfile(flatlog):
                print(f"ERROR: Flat log file not found: {flatlog}")
                sys.exit(1)
            args = args[2:]
        else:
            print(f"ERROR: Unknown option: {args[0]}")
            print_usage()
            sys.exit(1)

    if len(args) < 4:
        print_usage()
        sys.exit(1)

    raw_spec = args[0]
    out_spec = args[1]
    dark_path = args[2]
    flat_path = args[3]

    if not os.path.isdir(dark_path):
        print(f"ERROR: dark_path is not a directory: {dark_path}")
        sys.exit(1)
    if not os.path.isdir(flat_path):
        print(f"ERROR: flat_path is not a directory: {flat_path}")
        sys.exit(1)

    return raw_spec, out_spec, dark_path, flat_path, bestflat, debug, flat_future_days, flatlog


def read_list_file(path):
    files = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith(";"):
                continue
            files.append(line)
    return files


def has_wildcards(s):
    """Check if string contains glob wildcard characters."""
    return any(c in s for c in '*?[')


def resolve_raw_files(raw_spec):
    # Wildcards
    if has_wildcards(raw_spec):
        files = sorted(glob.glob(raw_spec))
        if not files:
            raise FileNotFoundError(f"No files match pattern: {raw_spec}")
        return files

    # Direct file or list file
    if os.path.isfile(raw_spec):
        ext = os.path.splitext(raw_spec)[1].lower()
        if ext in FITS_EXTS:
            return [raw_spec]

        # Treat as list file, paths relative to list location
        entries = read_list_file(raw_spec)
        if not entries:
            raise FileNotFoundError(f"No entries in list file: {raw_spec}")
        base_dir = os.path.dirname(raw_spec)
        files = []
        for p in entries:
            full = p if os.path.isabs(p) else os.path.join(base_dir, p)
            if not os.path.exists(full):
                raise FileNotFoundError(f"File from list not found: {full}")
            files.append(full)
        return files

    raise FileNotFoundError(f"Input spec not recognized or file not found: {raw_spec}")


def find_cosme_for_dark(dark_path):
    """
    Find cosmetic file for a dark.

    Expects cosme file in same folder with name pattern:
    dark<X>.fit -> cosme<X>.lst

    Returns cosmetic coords list or None.
    """
    if not HAVE_COSME:
        return None

    dark_dir = os.path.dirname(dark_path)
    dark_name = os.path.basename(dark_path)
    dark_stem = os.path.splitext(dark_name)[0]

    # Replace 'dark' with 'cosme' in filename
    if dark_stem.lower().startswith("dark"):
        cosme_stem = "cosme" + dark_stem[4:]
    else:
        cosme_stem = dark_stem.replace("dark", "cosme")

    cosme_path = os.path.join(dark_dir, cosme_stem + ".lst")

    if os.path.isfile(cosme_path):
        try:
            coords = read_cosme_list(cosme_path)
            return coords
        except Exception as e:
            sys.stderr.write(f"WARNING: Failed to read cosme '{cosme_path}': {e}\n")
            return None

    return None


def load_darks(dark_path):
    """Recursively load all dark frames as float64 with EXPTIME, JD, and cosmetic coords."""
    darks = []

    # Collect all .lst files for fallback (load lazily)
    all_lst_files = glob.glob(os.path.join(dark_path, "**", "*.lst"), recursive=True)
    fallback_lst = max(all_lst_files, key=os.path.getmtime) if all_lst_files else None
    fallback_cosme = None  # Will be loaded on first use
    fallback_used = False

    for root, _, files in os.walk(dark_path):
        for name in files:
            if os.path.splitext(name)[1].lower() not in FITS_EXTS:
                continue
            path = os.path.join(root, name)

            try:
                with fits.open(path, memmap=False) as hdul:
                    if hdul[0].data is None:
                        continue
                    hdr = hdul[0].header
                    data = hdul[0].data.astype(np.float64, copy=True)
            except Exception as e:
                sys.stderr.write(f"WARNING: Failed to read dark '{path}': {e}, skipped.\n")
                continue

            exptime = hdr.get("EXPTIME", hdr.get("EXPOSURE"))
            if exptime is None:
                sys.stderr.write(f"WARNING: Dark '{path}' has no EXPTIME, skipped.\n")
                continue

            jd = hdr.get("JD")
            if jd is None:
                sys.stderr.write(f"WARNING: Dark '{path}' has no JD, skipped.\n")
                continue

            # Find corresponding cosmetic file
            cosme_coords = find_cosme_for_dark(path)
            if cosme_coords is None and fallback_lst and HAVE_COSME:
                # Load fallback on first use
                if fallback_cosme is None and not fallback_used:
                    try:
                        fallback_cosme = read_cosme_list(fallback_lst)
                    except Exception as e:
                        sys.stderr.write(f"WARNING: Failed to read fallback cosme '{fallback_lst}': {e}\n")
                    fallback_used = True

                if fallback_cosme is not None:
                    sys.stderr.write(f"WARNING: No cosme for '{name}', using fallback: {os.path.basename(fallback_lst)}\n")
                    cosme_coords = fallback_cosme

            darks.append({
                "path": path,
                "data": data,
                "exptime": float(exptime),
                "jd": float(jd),
                "cosme_coords": cosme_coords,
            })

    if not darks:
        print("ERROR: No usable darks found in dark_path tree.")
        sys.exit(1)

    return darks


def load_flats(flat_path):
    """Recursively load all flats as float64 with FILTER, JD, and median (mult)."""
    flats = []
    for root, _, files in os.walk(flat_path):
        for name in files:
            if os.path.splitext(name)[1].lower() not in FITS_EXTS:
                continue
            path = os.path.join(root, name)

            try:
                with fits.open(path, memmap=False) as hdul:
                    if hdul[0].data is None:
                        continue
                    hdr = hdul[0].header
                    data = hdul[0].data.astype(np.float64, copy=True)
            except Exception as e:
                sys.stderr.write(f"WARNING: Failed to read flat '{path}': {e}, skipped.\n")
                continue

            flt = hdr.get("FILTERS", hdr.get("FILTER"))
            if flt is None:
                sys.stderr.write(f"WARNING: Flat '{path}' has no FILTER/FILTERS, skipped.\n")
                continue
            flt = str(flt).strip()
            if not flt:
                sys.stderr.write(f"WARNING: Flat '{path}' has empty FILTER/FILTERS, skipped.\n")
                continue

            jd = hdr.get("JD")
            if jd is None:
                sys.stderr.write(f"WARNING: Flat '{path}' has no JD, skipped.\n")
                continue

            # Compute multiplier as fast median of the flat
            mult = batch_utils.fast_median(data)
            if mult <= 0:
                sys.stderr.write(f"WARNING: Flat '{path}' has zero or negative median, skipped.\n")
                continue

            flats.append({
                "path": path,
                "data": data,
                "filter": flt,
                "jd": float(jd),
                "mult": mult,
            })

    if not flats:
        print("ERROR: No usable flats found in flat_path tree.")
        sys.exit(1)

    return flats




def clamp_to_dtype(data, dtype):
    """Clamp float64 to original dtype range."""
    if np.issubdtype(dtype, np.integer):
        info = np.iinfo(dtype)
        return np.clip(data, info.min, info.max).astype(dtype)
    elif np.issubdtype(dtype, np.floating):
        return data.astype(dtype)
    else:
        return data.astype(np.float64)


def jd_to_local_datetime(jd):
    """Convert Julian Date to local datetime."""
    # JD 2440587.5 = Unix epoch (1970-01-01 00:00:00 UTC)
    unix_seconds = (jd - 2440587.5) * 86400.0
    utc_dt = datetime.fromtimestamp(unix_seconds, tz=timezone.utc)
    # Convert to local time
    local_dt = utc_dt.astimezone()
    return local_dt


def jd_to_session_key(jd):
    """
    Convert JD to observing session key.

    Sessions are defined as noon-to-noon in local time (24h format):
    - From 12:01 today to 12:00 tomorrow = today's session
    - 12:00 exactly is the END of the previous session
    - 12:01 is the START of the new session

    Returns session date as string 'YYYY-MM-DD'.
    """
    local_dt = jd_to_local_datetime(jd)

    # Calculate total minutes from midnight
    total_minutes = local_dt.hour * 60 + local_dt.minute

    # 12:00 (720 min) or earlier → previous day's session
    # 12:01 (721 min) or later → current day's session
    if total_minutes <= 720:
        session_date = local_dt.date() - timedelta(days=1)
    else:
        session_date = local_dt.date()

    return session_date.isoformat()


def datetime_to_jd(dt):
    """Convert datetime to Julian Date."""
    # Ensure UTC
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)

    unix_seconds = dt.timestamp()
    # JD 2440587.5 = Unix epoch (1970-01-01 00:00:00 UTC)
    return unix_seconds / 86400.0 + 2440587.5


def read_maintenance_log(log_path):
    """
    Read maintenance log CSV file.

    Format: DATETIME_UTC,CAMERA_ID[,COMMENT]
    Lines starting with # are comments.

    Returns list of dicts: [{"jd": float, "camera": str}, ...]
    sorted by JD ascending.
    """
    entries = []

    with open(log_path, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split(",")
            if len(parts) < 2:
                sys.stderr.write(
                    f"WARNING: Maintenance log line {line_num}: invalid format, skipping.\n"
                )
                continue

            datetime_str = parts[0].strip()
            camera = parts[1].strip()

            try:
                # Parse ISO datetime, handle Z suffix
                dt = datetime.fromisoformat(datetime_str.replace("Z", "+00:00"))
                # If no timezone, assume UTC
                if dt.tzinfo is None:
                    dt = dt.replace(tzinfo=timezone.utc)
                jd = datetime_to_jd(dt)
            except ValueError as e:
                sys.stderr.write(
                    f"WARNING: Maintenance log line {line_num}: invalid datetime '{datetime_str}', skipping.\n"
                )
                continue

            entries.append({
                "jd": jd,
                "camera": camera,
                "datetime_str": datetime_str,
            })

    # Sort by JD ascending
    entries.sort(key=lambda e: e["jd"])
    return entries


def get_camera_from_header(header):
    """Extract camera identifier from FITS header (INSTRUME field)."""
    instrume = header.get("INSTRUME", "")
    return str(instrume).strip()


def find_maintenance_interval(jd, camera, maintenance_entries):
    """
    Find the maintenance interval for a given JD and camera.

    Returns tuple (start_jd, end_jd) where:
    - start_jd: JD of maintenance event at start of interval (or None if before first)
    - end_jd: JD of maintenance event at end of interval (or None if after last)

    Camera matching is substring-based: camera ID from log is searched
    as substring in INSTRUME value.
    """
    # Filter entries for this camera (substring match)
    camera_entries = [
        e for e in maintenance_entries
        if e["camera"] in camera or camera in e["camera"]
    ]

    if not camera_entries:
        # No maintenance entries for this camera - entire history is one interval
        return (None, None)

    # Find which interval the JD falls into
    start_jd = None
    end_jd = None

    for entry in camera_entries:
        if entry["jd"] <= jd:
            start_jd = entry["jd"]
        elif end_jd is None:
            end_jd = entry["jd"]
            break

    return (start_jd, end_jd)


def is_flat_in_interval(flat_jd, interval):
    """
    Check if a flat JD falls within the maintenance interval.

    interval is (start_jd, end_jd) where either can be None.
    """
    start_jd, end_jd = interval

    if start_jd is not None and flat_jd < start_jd:
        return False
    if end_jd is not None and flat_jd >= end_jd:
        return False

    return True


def downscale_8x8(data):
    """
    Downscale image by 8x8 averaging.

    Reduces image size by factor of 8 in each dimension.
    Truncates to fit exact 8x8 blocks.
    """
    h, w = data.shape
    new_h = h // 8
    new_w = w // 8

    # Truncate to exact multiple of 8
    truncated = data[:new_h * 8, :new_w * 8]

    # Reshape to (new_h, 8, new_w, 8) and average
    reshaped = truncated.reshape(new_h, 8, new_w, 8)
    return reshaped.mean(axis=(1, 3))


def apply_highpass_filter(data):
    """
    Apply high-pass filter by subtracting gaussian-blurred copy.

    Blur sigma = BESTFLAT_HIGHPASS_SIGMA_FACTOR * ((width + height) / 2)
    This removes large-scale gradients while preserving local variations
    (dust spots, flat field artifacts).

    Returns:
        (highpassed, blur_mean) - highpassed data and mean of the blurred image
                                  (for debug file offset restoration)
    """
    h, w = data.shape
    sigma = BESTFLAT_HIGHPASS_SIGMA_FACTOR * ((w + h) / 2)
    blurred = gaussian_filter(data, sigma=sigma)
    blur_mean = np.mean(blurred)
    return data - blurred, blur_mean


def normalize_to_median(data, target_median=5000):
    """
    Normalize image to target median (ngain-style).

    Scales so that median becomes target_median.
    Clamps to uint16 range [0, 32767].
    """
    # Shift to positive range first (high-pass can produce negatives)
    min_val = np.min(data)
    shifted = data - min_val

    # Scale to target median
    shifted_median = np.median(shifted)
    if shifted_median > 0:
        scaled = shifted * (target_median / shifted_median)
    else:
        # shifted_median == 0 means median(data) == min(data)
        # Can't scale meaningfully, just shift to put median at target
        # (this will make the output have median = target_median when min_val == median)
        scaled = shifted + target_median

    # Clamp to 0-32767
    return np.clip(scaled, 0, 32767)


def choose_dark(raw_exptime, raw_jd, darks):
    """
    Select DARK:
      1) Same EXPTIME (with small tolerance), choose latest by JD.
      2) Otherwise: closest EXPTIME, then latest JD as tie-breaker.
    """

    def eq(a, b):
        if a == b:
            return True
        ref = max(1.0, abs(a), abs(b))
        return abs(a - b) <= 1e-6 * ref

    same = [d for d in darks if eq(d["exptime"], raw_exptime)]
    if same:
        return max(same, key=lambda d: d["jd"])

    return min(
        darks,
        key=lambda d: (abs(d["exptime"] - raw_exptime), -d["jd"])
    )


def choose_flat(raw_filter, raw_jd, flats, flat_future_days=2.0, maintenance_interval=None):
    """
    Select FLAT:
      - Use flats with FILTER == raw_filter, JD_flat <= JD_raw + flat_future_days.
      - If maintenance_interval is provided:
        1. First try current interval
        2. If not found, search in earlier intervals (before interval start)
      - If none found, try FILTER == 'L' with same logic.
      - If still none, raise error.

    Args:
        raw_filter: filter name to match
        raw_jd: Julian Date of the observation
        flats: list of flat dicts
        flat_future_days: max days in future to accept a flat
        maintenance_interval: tuple (start_jd, end_jd) from maintenance log, or None
    """
    if not flats:
        raise RuntimeError("No flats available.")

    raw_filter = (raw_filter or "").strip()
    limit_jd = raw_jd + flat_future_days

    def select_in_interval(candidates, interval):
        """Select best (latest) flat within given interval.

        If interval is set:
          - First try with +N days limit within interval
          - If nothing found but interval has flats, take latest from interval (priority)
          - If interval is empty, return None (caller will try without interval)
        If interval is None:
          - Apply +N days limit only
        """
        if interval is not None:
            start_jd, end_jd = interval
            if start_jd is not None or end_jd is not None:
                # Get all flats in interval
                in_interval = [f for f in candidates if is_flat_in_interval(f["jd"], interval)]
                if not in_interval:
                    return None  # Nothing in interval, caller will try without

                # Try with +N days limit first
                with_limit = [f for f in in_interval if f["jd"] <= limit_jd]
                if with_limit:
                    return max(with_limit, key=lambda f: f["jd"])

                # Interval has flats but none within +N days - take latest anyway (interval priority)
                return max(in_interval, key=lambda f: f["jd"])
            else:
                # Interval is (None, None) - treat as no interval
                valid = [f for f in candidates if f["jd"] <= limit_jd]
        else:
            # No interval: apply future days limit only
            valid = [f for f in candidates if f["jd"] <= limit_jd]

        if not valid:
            return None
        return max(valid, key=lambda f: f["jd"])

    def search_with_fallback(candidates):
        """Search current interval, then without interval if not found."""
        # Try current interval first
        best = select_in_interval(candidates, maintenance_interval)
        if best is not None:
            return best, False  # found in current interval

        # If interval was set but nothing found, try without interval (global search)
        if maintenance_interval is not None:
            best = select_in_interval(candidates, None)
            if best is not None:
                return best, True  # found outside interval

        return None, False

    # Try exact filter match
    exact = [f for f in flats if f["filter"] == raw_filter]
    best, from_outside = search_with_fallback(exact)
    if best is not None:
        if from_outside:
            flat_date = jd_to_local_datetime(best["jd"]).strftime("%Y-%m-%d")
            sys.stderr.write(
                f"WARNING: No flat for filter '{raw_filter}' in current interval, "
                f"using flat from {flat_date}.\n"
            )
        return best, False

    # Try L fallback
    l_cands = [f for f in flats if f["filter"] == "L"]
    best, from_outside = search_with_fallback(l_cands)
    if best is not None:
        msg = f"WARNING: No flat for filter '{raw_filter}', using L-flat '{os.path.basename(best['path'])}'"
        if from_outside:
            flat_date = jd_to_local_datetime(best["jd"]).strftime("%Y-%m-%d")
            msg += f" from outside interval ({flat_date})"
        sys.stderr.write(msg + ".\n")
        return best, True

    # Build descriptive error message
    interval_info = ""
    if maintenance_interval is not None:
        start_jd, end_jd = maintenance_interval
        if start_jd is not None:
            start_dt = jd_to_local_datetime(start_jd).strftime("%Y-%m-%d %H:%M")
            interval_info += f" after {start_dt}"
        if end_jd is not None:
            end_dt = jd_to_local_datetime(end_jd).strftime("%Y-%m-%d %H:%M")
            interval_info += f" before {end_dt}"

    raise RuntimeError(
        f"No suitable flat for filter '{raw_filter}'{interval_info} (and no L-flat)."
    )


def select_best_flat_for_session(light_previews, candidate_flats, flat_previews_cache, session_key, filter_name, debug_dir=None, maintenance_interval=None):
    """
    Select optimal flat for a session by minimizing stddev after flat division.

    Args:
        light_previews: list of downscaled dark-subtracted light arrays
        candidate_flats: list of flat dicts matching the filter
        flat_previews_cache: dict to cache downscaled flats {path: array}
        session_key: session identifier string for logging
        filter_name: filter name for logging
        debug_dir: if set, save debug images to this directory
        maintenance_interval: tuple (start_jd, end_jd) from maintenance log, or None

    Returns:
        (best_flat_dict, used_L_fallback)
    """
    if not light_previews:
        raise RuntimeError(f"No light previews for session {session_key}")

    if not candidate_flats:
        raise RuntimeError(f"No candidate flats for filter '{filter_name}' in session {session_key}")

    # Filter by maintenance interval if provided
    if maintenance_interval is not None:
        candidate_flats = [
            f for f in candidate_flats
            if is_flat_in_interval(f["jd"], maintenance_interval)
        ]
        if not candidate_flats:
            interval_info = ""
            start_jd, end_jd = maintenance_interval
            if start_jd is not None:
                start_dt = jd_to_local_datetime(start_jd).strftime("%Y-%m-%d %H:%M")
                interval_info += f" after {start_dt}"
            if end_jd is not None:
                end_dt = jd_to_local_datetime(end_jd).strftime("%Y-%m-%d %H:%M")
                interval_info += f" before {end_dt}"
            raise RuntimeError(
                f"No candidate flats for filter '{filter_name}' in session {session_key} "
                f"within maintenance interval{interval_info}"
            )

    # Median combine light previews
    if len(light_previews) == 1:
        light_median = light_previews[0]
    else:
        stack = np.array(light_previews)
        light_median = np.median(stack, axis=0)

    # Create debug subfolder for this session if needed
    session_debug_dir = None
    if debug_dir:
        safe_session = session_key.replace("-", "_")
        safe_filter = filter_name.replace("/", "_").replace("\\", "_")
        session_debug_dir = os.path.join(debug_dir, f"{safe_session}_{safe_filter}")
        os.makedirs(session_debug_dir, exist_ok=True)

        # Save light median preview
        light_debug_path = os.path.join(session_debug_dir, "light_median.fit")
        fits.writeto(light_debug_path, light_median.astype(np.float32), overwrite=True)

    # Evaluate each flat
    results = []
    for flat in candidate_flats:
        flat_path = flat["path"]

        # Get or compute downscaled flat (NO median filter - preserve dust!)
        if flat_path not in flat_previews_cache:
            # Skip if flat has invalid multiplier
            if flat["mult"] <= 0:
                sys.stderr.write(f"WARNING: Flat '{os.path.basename(flat_path)}' has invalid mult, skipping.\n")
                continue

            # Downscale only (no median filter - we need to see dust spots)
            downscaled = downscale_8x8(flat["data"])
            # Normalize to median 5000 for consistent comparison
            downscaled = downscaled * (5000.0 / flat["mult"])
            flat_previews_cache[flat_path] = downscaled

            # Save downscaled flat to debug root if debug enabled
            if debug_dir:
                flat_basename = os.path.splitext(os.path.basename(flat_path))[0]
                flat_date_str = jd_to_local_datetime(flat["jd"]).strftime("%Y%m%d")
                flat_debug_name = f"flat_{flat_basename}_{flat_date_str}.fit"
                flat_debug_path = os.path.join(debug_dir, flat_debug_name)
                fits.writeto(flat_debug_path, flat_previews_cache[flat_path].astype(np.float32), overwrite=True)

        flat_preview = flat_previews_cache[flat_path]

        # Ensure same shape
        if flat_preview.shape != light_median.shape:
            sys.stderr.write(
                f"WARNING: Flat '{os.path.basename(flat_path)}' preview shape mismatch, skipping.\n"
            )
            continue

        # Divide and compute sigma
        # flat_preview is normalized to median 5000, so use 5000 as multiplier
        safe_flat = flat_preview.copy()
        safe_flat[safe_flat <= 0] = np.nan
        divided = light_median * 5000.0 / safe_flat
        divided = np.nan_to_num(divided, nan=0.0, posinf=0.0, neginf=0.0)

        # Apply high-pass filter before computing sigma
        # (removes large-scale gradients while preserving dust/vignetting artifacts)
        highpassed, blur_mean = apply_highpass_filter(divided)

        # Apply gaussian blur to suppress high-frequency noise before sigma calculation
        blurred_for_sigma = gaussian_filter(highpassed, sigma=BESTFLAT_SIGMA_BLUR_SIZE)

        # Standard deviation (after high-pass and blur)
        sigma = np.std(blurred_for_sigma)

        # Get date from flat JD for display
        flat_date = jd_to_local_datetime(flat["jd"]).strftime("%Y-%m-%d %H:%M")
        flat_date_str = jd_to_local_datetime(flat["jd"]).strftime("%Y%m%d")

        # Store data for debug output (will be saved after sorting to mark best)
        if session_debug_dir:
            debug_data = blurred_for_sigma if BESTFLAT_DEBUG_SAVE_BLURRED else highpassed
        else:
            debug_data = None
        results.append({
            "flat": flat,
            "sigma": sigma,
            "date": flat_date,
            "date_str": flat_date_str,
            "debug_data": debug_data,
            "blur_mean": blur_mean,
        })

    if not results:
        raise RuntimeError(f"No valid flats for session {session_key} filter '{filter_name}'")

    # Sort by sigma ascending
    results.sort(key=lambda r: r["sigma"])

    # Save debug files with _best marker for the best one
    if session_debug_dir:
        for i, r in enumerate(results):
            flat_basename = os.path.splitext(os.path.basename(r["flat"]["path"]))[0]
            best_marker = "_best" if i == 0 else ""

            debug_data = r["debug_data"]
            if debug_data is None:
                sys.stderr.write(f"DEBUG: debug_data is None for flat {flat_basename}\n")
                continue

            # Shift by blur_mean to restore offset, scale down, clamp negatives to 0
            shifted = (debug_data + r["blur_mean"]) * 0.5
            shifted = np.clip(shifted, 0, None)

            # Add calibration line to fix viewer scaling (first row: left half=0, right half=32767)
            if BESTFLAT_DEBUG_CALIBRATION_LINE:
                h, w = shifted.shape
                shifted[0, :w // 2] = 0
                shifted[0, w // 2:] = 32767

            debug_name = f"hp_{flat_basename}_{r['date_str']}_sigma{r['sigma']:.1f}{best_marker}.fit"
            debug_path = os.path.join(session_debug_dir, debug_name)
            fits.writeto(debug_path, shifted.astype(np.float32), overwrite=True)

    # Print results table
    print(f"\n--- Best flat selection for session {session_key}, filter '{filter_name}' ---")
    print(f"{'Flat file':<40} {'Date':<18} {'Sigma':>12}")
    print("-" * 72)
    for i, r in enumerate(results):
        marker = " <-- BEST" if i == 0 else ""
        print(f"{os.path.basename(r['flat']['path']):<40} {r['date']:<18} {r['sigma']:>12.2f}{marker}")
    print()

    best = results[0]["flat"]
    print(f"Selected: {best['path']}")
    print()

    return best, False


def load_raw_header(path):
    """Load only header and dtype of RAW."""
    with fits.open(path, memmap=False) as hdul:
        if hdul[0].data is None:
            raise RuntimeError(f"No image data in file: {path}")
        return hdul[0].header, hdul[0].data.dtype


def sanitize_filter_for_name(f):
    """Make filter name safe for filenames. Default to 'L' if empty."""
    if not f:
        f = "L"
    f = str(f).strip()
    if not f:
        f = "L"
    f = f.replace(" ", "")
    for ch in '\\/:*?"<>|':
        f = f.replace(ch, "")
    return f or "L"


def build_jobs(raw_files, out_spec, darks, flats, bestflat=False, debug_dir=None, flat_future_days=2.0, maintenance_entries=None):
    """
    Build job list.

    out_spec -> dir/base.ext

    For each RAW file:
      - Read JD (required).
      - Read EXPTIME/EXPOSURE (required).
      - Read FILTERS/FILTER:
          if missing/empty -> treat as 'L' (and use as such in names and selection).
      - Select DARK and FLAT as per rules.
      - Output filename:
          <base>_exp<EXPTIME>_<FILTER>_<N>.ext
        where:
          EXPTIME: integer if integral; else dot replaced by 'p' (e.g. 0p5).
          FILTER: sanitized string.
          N: independent counter for each (EXPTIME,FILTER) pair.

    If bestflat=True:
      - Group files by (filter, session) where session is noon-to-noon local time
      - For each group, compute optimal flat by minimizing stddev

    If maintenance_entries is provided:
      - Determine maintenance interval for each observation based on camera
      - Use strict flat selection (exact filter, no L fallback)
      - Search current interval first, then earlier intervals
    """
    total = len(raw_files)

    out_dir = os.path.dirname(out_spec) or "."
    base = os.path.basename(out_spec)
    stem, ext = os.path.splitext(base)
    if ext.lower() not in FITS_EXTS:
        ext = ".fit"
    if not stem:
        stem = "out"

    os.makedirs(out_dir, exist_ok=True)

    # First pass: collect metadata for all files
    file_info = []
    for raw_path in raw_files:
        hdr, src_dtype = load_raw_header(raw_path)

        jd = hdr.get("JD")
        if jd is None:
            raise RuntimeError(f"Raw '{raw_path}' has no JD.")
        jd = float(jd)

        exptime = hdr.get("EXPTIME", hdr.get("EXPOSURE"))
        if exptime is None:
            raise RuntimeError(f"Raw '{raw_path}' has no EXPTIME/EXPOSURE.")
        exptime = float(exptime)

        flt = hdr.get("FILTERS", hdr.get("FILTER"))
        flt = (flt or "").strip()
        if not flt:
            flt = "L"

        dark = choose_dark(exptime, jd, darks)
        session_key = jd_to_session_key(jd)

        # Get camera and maintenance interval if maintenance log is provided
        camera = None
        maintenance_interval = None
        if maintenance_entries is not None:
            camera = get_camera_from_header(hdr)
            if camera:
                maintenance_interval = find_maintenance_interval(jd, camera, maintenance_entries)
            else:
                sys.stderr.write(
                    f"WARNING: No INSTRUME in '{raw_path}', cannot use maintenance log.\n"
                )

        file_info.append({
            "raw_path": raw_path,
            "src_dtype": src_dtype,
            "jd": jd,
            "exptime": exptime,
            "filter": flt,
            "dark": dark,
            "session": session_key,
            "camera": camera,
            "maintenance_interval": maintenance_interval,
        })

    # Determine flat selection method
    if bestflat:
        # Group by (filter, session)
        groups = {}
        for i, info in enumerate(file_info):
            key = (info["filter"], info["session"])
            if key not in groups:
                groups[key] = []
            groups[key].append(i)

        # Cache for downscaled flats
        flat_previews_cache = {}

        # For each group, select best flat
        session_flat_map = {}  # (filter, session) -> (flat, used_L)

        for (filter_name, session_key), indices in groups.items():
            print(f"Processing session {session_key}, filter '{filter_name}' ({len(indices)} files)...")

            # Get maintenance interval from first file in group
            # (all files in same session should have same interval)
            first_info = file_info[indices[0]]
            group_maintenance_interval = first_info.get("maintenance_interval")

            # Find candidate flats for this filter
            candidate_flats = [f for f in flats if f["filter"] == filter_name]
            used_L_fallback = False

            if not candidate_flats:
                # Try L fallback
                candidate_flats = [f for f in flats if f["filter"] == "L"]
                if candidate_flats:
                    sys.stderr.write(
                        f"WARNING: No flats for filter '{filter_name}', trying L-flats.\n"
                    )
                    used_L_fallback = True
                else:
                    raise RuntimeError(
                        f"No flats for filter '{filter_name}' and no L-flats available."
                    )

            # Load and process light previews for this group
            light_previews = []
            for idx in indices:
                info = file_info[idx]
                raw_path = info["raw_path"]
                dark = info["dark"]

                # Load raw data
                with fits.open(raw_path, memmap=False) as hdul:
                    raw = hdul[0].data.astype(np.float64, copy=True)

                # Subtract dark
                if raw.shape != dark["data"].shape:
                    raise RuntimeError(
                        f"Shape mismatch RAW {raw.shape} vs DARK {dark['data'].shape} "
                        f"for '{raw_path}'"
                    )
                work = raw - dark["data"]

                # Apply cosmetic correction (from dark's cosme_coords)
                work = apply_cosmetic_if_needed(work, dark.get("cosme_coords"))

                # Downscale and median filter to remove stars
                preview = downscale_8x8(work)
                preview = median_filter(preview, size=BESTFLAT_LIGHT_MEDIAN_SIZE)
                light_previews.append(preview)

            # Select best flat
            try:
                best_flat, _ = select_best_flat_for_session(
                    light_previews, candidate_flats, flat_previews_cache,
                    session_key, filter_name, debug_dir, group_maintenance_interval
                )
            except ZeroDivisionError as e:
                sys.stderr.write(f"\n=== ZeroDivisionError in select_best_flat_for_session ===\n")
                sys.stderr.write(f"Session: {session_key}, Filter: {filter_name}\n")
                sys.stderr.write(f"Candidate flats ({len(candidate_flats)}):\n")
                for f in candidate_flats:
                    sys.stderr.write(f"  {f['path']}: mult={f['mult']}\n")
                sys.stderr.write(f"\nFull traceback:\n")
                traceback.print_exc(file=sys.stderr)
                raise
            session_flat_map[(filter_name, session_key)] = (best_flat, used_L_fallback)

    # Build jobs
    per_combo_counts = {}
    jobs = []

    for idx, info in enumerate(file_info, start=1):
        raw_path = info["raw_path"]
        flt = info["filter"]
        jd = info["jd"]
        exptime = info["exptime"]
        dark = info["dark"]
        src_dtype = info["src_dtype"]
        session_key = info["session"]

        if bestflat:
            flat, used_L = session_flat_map[(flt, session_key)]
        else:
            maintenance_interval = info.get("maintenance_interval")
            flat, used_L = choose_flat(flt, jd, flats, flat_future_days, maintenance_interval)

        safe_filter = sanitize_filter_for_name(flt)

        # Format exposure for filename
        if abs(exptime - int(exptime)) < 1e-6:
            exp_str = str(int(exptime))
        else:
            exp_str = str(exptime).replace(".", "p")

        combo_key = (exp_str, safe_filter)
        count = per_combo_counts.get(combo_key, 0) + 1
        per_combo_counts[combo_key] = count

        out_name = f"{stem}_exp{exp_str}_{safe_filter}_{count}{ext}"
        out_file = os.path.join(out_dir, out_name)

        jobs.append({
            "index": idx,
            "total": total,
            "raw": raw_path,
            "out": out_file,
            "src_dtype": src_dtype,
            "raw_jd": jd,
            "raw_exptime": exptime,
            "raw_filter": flt,
            "dark": dark,
            "flat": flat,
            "mult": flat["mult"],
            "used_L_fallback": used_L,
        })

    return jobs


def apply_cosmetic_if_needed(data, cosme_coords):
    if cosme_coords is None or not HAVE_COSME:
        return data
    return apply_cosmetic_correction(data, cosme_coords)


def process_one(job):
    """
    Worker: calibrate single file.
    Runs in threads. No printing here.

    Calibration order: dark subtraction -> cosmetic -> flat division
    """
    raw_path = job["raw"]
    out_path = job["out"]
    dark = job["dark"]
    flat = job["flat"]
    mult = job["mult"]
    src_dtype = job["src_dtype"]
    cosme_coords = dark.get("cosme_coords")

    with fits.open(raw_path, memmap=False) as hdul:
        if hdul[0].data is None:
            raise RuntimeError(f"No image data in file: {raw_path}")
        header = hdul[0].header.copy()
        raw = hdul[0].data.astype(np.float64, copy=True)

    dark_data = dark["data"]
    flat_data = flat["data"]

    if raw.shape != dark_data.shape:
        raise RuntimeError(
            f"Shape mismatch RAW {raw.shape} vs DARK {dark_data.shape} for '{raw_path}'"
        )
    if raw.shape != flat_data.shape:
        raise RuntimeError(
            f"Shape mismatch RAW {raw.shape} vs FLAT {flat_data.shape} for '{raw_path}'"
        )

    # Step 1: Dark subtraction
    work = raw - dark_data

    # Step 2: Cosmetic correction (before flat division)
    work = apply_cosmetic_if_needed(work, cosme_coords)

    # Step 3: Flat division - RESULT = work * MULT / FLAT
    safe_flat = flat_data.copy()
    safe_flat[safe_flat <= 0] = np.nan
    work = work * mult / safe_flat

    # Replace NaN/inf by finite values
    work = np.nan_to_num(work, copy=False)

    out_data = clamp_to_dtype(work, src_dtype)

    # Preserve original header, add calibration info
    header["HISTORY"] = "Calibrated by autocalibrate.py"
    header["HISTORY"] = f"  DARK = {os.path.basename(dark['path'])}"
    header["HISTORY"] = f"  FLAT = {os.path.basename(flat['path'])}"
    header["HISTORY"] = f"  MULT = {mult}"
    if job.get("used_L_fallback"):
        header["HISTORY"] = "  FLAT L fallback used"
    header["CALIB"] = ("AUTO", "Calibrated by autocalibrate.py")

    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fits.writeto(out_path, out_data, header, overwrite=True)

    return (
        job["index"],
        job["total"],
        raw_path,
        out_path,
        dark["path"],
        flat["path"],
    )


def main(argv=None):
    if argv is None:
        argv = sys.argv

    raw_spec, out_spec, dark_path, flat_path, bestflat, debug, flat_future_days, flatlog = parse_args(argv)

    if bestflat:
        print("Mode: --bestflat (auto-select optimal flat per session)")

    # Setup debug directory if requested
    debug_dir = None
    if debug:
        debug_dir = os.path.join(os.getcwd(), "debug")
        os.makedirs(debug_dir, exist_ok=True)
        print(f"Debug mode: saving previews to {debug_dir}")

    # Load maintenance log if provided
    maintenance_entries = None
    if flatlog:
        print(f"Using maintenance log: {flatlog}")
        maintenance_entries = read_maintenance_log(flatlog)
        if maintenance_entries:
            print(f"  Loaded {len(maintenance_entries)} maintenance entries")
            for entry in maintenance_entries:
                print(f"    {entry['datetime_str']} - {entry['camera']}")
        else:
            print("  WARNING: Maintenance log is empty")

    # Resolve RAW files
    try:
        raw_files = resolve_raw_files(raw_spec)
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    if not raw_files:
        print("ERROR: No raw files resolved.")
        sys.exit(1)

    print(f"Found {len(raw_files)} raw file(s).")

    # Load masters once (darks include their cosmetic coords)
    darks = load_darks(dark_path)
    flats = load_flats(flat_path)

    # Build jobs (per-file DARK/FLAT selection, output names)
    try:
        jobs = build_jobs(raw_files, out_spec, darks, flats, bestflat, debug_dir, flat_future_days, maintenance_entries)
    except Exception as e:
        print(f"ERROR while preparing jobs: {e}")
        sys.exit(1)

    total = len(jobs)

    # Thread pool size: up to CPU-1, not more than number of jobs
    try:
        cpu_count = os.cpu_count() or 1
    except Exception:
        cpu_count = 1
    max_workers = max(1, cpu_count - 1)
    max_workers = min(max_workers, total)

    print(f"Using {max_workers} worker thread(s).")

    # Multithreaded processing; logging only from main thread
    futures = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for job in jobs:
            fut = executor.submit(process_one, job)
            futures[fut] = job

        done = 0
        for fut in as_completed(futures):
            job = futures[fut]
            try:
                idx, total_count, raw_path, out_path_file, dark_used, flat_used = fut.result()
            except Exception as e:
                sys.stderr.write("\n")
                print(f"ERROR processing '{job['raw']}': {e}")
                sys.exit(1)

            done += 1
            print(
                f"[{done}/{total_count}] "
                f"{os.path.basename(raw_path)} -> {os.path.basename(out_path_file)} "
                f"| DARK={os.path.basename(dark_used)} "
                f"| FLAT={os.path.basename(flat_used)}"
            )

    print("Auto-calibration completed.")


if __name__ == "__main__":
    main()
