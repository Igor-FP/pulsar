#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
absession - Generate AstroBin acquisition session CSV from FITS files.

Scans a directory (recursively) for FITS files, groups them by date, filter,
and exposure time, then outputs a CSV compatible with AstroBin's CSV import.

Two modes of operation:
  - FITS header reading (default): reads FILTER, EXPTIME, DATE-OBS from headers
  - Filename parsing (--parsename): extracts filter and exposure from filename
    following the convention: Name_Type_FILTER_..._EXPOSURE_..._Seq.fit

Session date is determined by DATE-OBS header (default) or file modification
time (--parsename / fallback). Images taken before noon are counted as the
previous evening's session (astronomers image from dusk to dawn).

Requires: astropy (pip install astropy) for default mode.
Filename mode (--parsename) has no external dependencies.

Part of PULSAR (Portable Utility Library for Shell Astrophotography Routines).
"""

import sys
import os
import re
import datetime
from collections import defaultdict

# ---------------------------------------------------------------
# AstroBin filter IDs
# ---------------------------------------------------------------
# *** IMPORTANT: EDIT THESE VALUES! ***
# The IDs below are EXAMPLES. Each AstroBin filter has a unique numeric ID.
# To find your filter IDs:
#   1. Search for your filter at:
#      https://app.astrobin.com/equipment/explorer/filter?page=1
#   2. Open the filter page — the ID is the number in the URL, e.g.:
#      https://app.astrobin.com/equipment/explorer/filter/4388/filter-name-...
#      → ID = 4388
#   3. Replace the values below with YOUR filter IDs
# ---------------------------------------------------------------

ASTROBIN_FILTER_IDS = {
    'L':    5652,   # CHANGE to your Luminance filter ID
    'R':    5656,   # CHANGE to your Red filter ID
    'G':    5646,   # CHANGE to your Green filter ID
    'B':    5642,   # CHANGE to your Blue filter ID
    'Ha':   4388,   # CHANGE to your H-alpha filter ID
    'SII':  4396,   # CHANGE to your S-II filter ID
    'OIII': 4392,   # CHANGE to your O-III filter ID
}

# Common aliases for filter names found in FITS headers
FILTER_ALIASES = {
    # Luminance
    'L':          'L',
    'LUM':        'L',
    'LUMINANCE':  'L',
    'CLEAR':      'L',
    'CLR':        'L',
    # Red
    'R':          'R',
    'RED':        'R',
    # Green
    'G':          'G',
    'GREEN':      'G',
    # Blue
    'B':          'B',
    'BLUE':       'B',
    # Narrowband
    'HA':         'Ha',
    'H-ALPHA':    'Ha',
    'HALPHA':     'Ha',
    'H_ALPHA':    'Ha',
    'SII':        'SII',
    'S-II':       'SII',
    'S2':         'SII',
    'SULPHUR':    'SII',
    'SULFUR':     'SII',
    'OIII':       'OIII',
    'O-III':      'OIII',
    'O3':         'OIII',
    'OXYGEN':     'OIII',
}

# Channel display order
CHANNEL_ORDER = ['L', 'R', 'G', 'B', 'Ha', 'SII', 'OIII']


# ---------------------------------------------------------------
# CLI
# ---------------------------------------------------------------

def usage():
    sys.stderr.write(
        "Usage:\n"
        "  absession.py [options] [directory]\n"
        "\n"
        "  Scans directory for FITS files (.fit, .fits) and generates\n"
        "  an AstroBin-compatible CSV acquisition session list.\n"
        "\n"
        "  Default mode reads FITS headers (requires astropy).\n"
        "\n"
        "Options:\n"
        "  --parsename     Extract filter/exposure from filenames instead\n"
        "                  of FITS headers (no dependencies required)\n"
        "  --out FILE      Write CSV to file instead of stdout\n"
        "  --flat          Do not recurse into subdirectories\n"
        "\n"
        "  directory       Path to scan (default: current directory)\n"
        "\n"
        "Filename convention (--parsename mode):\n"
        "  Name_Type_FILTER_SeqNum_Binning_EXPOSURE_Temp.fit\n"
        "  Example: NGC253_L_R_54205_Bin1x1_120s_-10C.fit\n"
        "\n"
        "  Only FILTER and EXPOSURE are used; other fields are ignored.\n"
        "  FILTER:   R, G, B, L, Ha, SII, OIII (default: L)\n"
        "  EXPOSURE: integer seconds, with optional 's' suffix (e.g. 600, 120s)\n"
        "\n"
        "  Mosaic variant (extra field after Type):\n"
        "  Name_Type_Mosaic_FILTER_SeqNum_Binning_EXPOSURE_Temp.fit\n"
        "\n"
        "  Configure this naming pattern in APT, N.I.N.A., SGPro,\n"
        "  or your imaging capture software.\n"
        "\n"
        "Output CSV format (AstroBin import):\n"
        "  date,filter,number,duration\n"
        "  2024-03-15,5656,12,600\n"
        "\n"
        "IMPORTANT: Filter IDs in this script are EXAMPLES.\n"
        "  You MUST edit ASTROBIN_FILTER_IDS in absession.py to match\n"
        "  YOUR filters on AstroBin.\n"
        "  Search: https://app.astrobin.com/equipment/explorer/filter?page=1\n"
        "  The ID is the number in the filter page URL, e.g.:\n"
        "  .../filter/4388/filter-name-...  ->  ID = 4388\n"
        "\n"
        "Default IDs (examples):\n"
        "  L=5652, R=5656, G=5646, B=5642, Ha=4388, SII=4396, OIII=4392\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]
    parsename = False
    out_file = None
    flat = False
    scan_dir = "."

    while args:
        if args[0] == '--parsename':
            parsename = True
            args.pop(0)
        elif args[0] == '--flat':
            flat = True
            args.pop(0)
        elif args[0] == '--out' and len(args) > 1:
            out_file = args[1]
            args = args[2:]
        elif args[0] in ('-h', '--help', '-?'):
            usage()
        elif args[0].startswith('-'):
            sys.stderr.write(f"Unknown option: {args[0]}\n")
            usage()
        else:
            scan_dir = args[0]
            args.pop(0)

    return scan_dir, parsename, out_file, flat


# ---------------------------------------------------------------
# Session date from file modification time
# ---------------------------------------------------------------

def get_session_date(filepath):
    """
    Get session date from file modification time.
    If the file was modified before noon, count it as the previous day's
    session (overnight imaging session).
    """
    mtime = os.path.getmtime(filepath)
    dt = datetime.datetime.fromtimestamp(mtime)
    if dt.hour < 12:
        dt -= datetime.timedelta(hours=12)
    return dt.date()


# ---------------------------------------------------------------
# Filename-based parsing (--parsename mode)
# ---------------------------------------------------------------

def _split_filename(filename):
    """Split filename by underscores, dots, and quotes (like C++ version)."""
    parts = []
    current = []
    for ch in filename:
        if ch in ('_', '"', '.'):
            if current:
                parts.append(''.join(current))
                current = []
        else:
            current.append(ch)
    if current:
        parts.append(''.join(current))
    return parts


def parse_filename(filepath):
    """
    Parse FITS filename to extract filter and exposure time.
    Returns (filter_name, exposure_seconds) or None if unparseable.

    Expected filename format (7 fields + extension, underscore-separated):

      NGC253_L_R_54205_Bin1x1_120s_-10C.fit
      |      | | |     |      |    |    |
      |      | | |     |      |    |    +-- extension (.fit/.fits)
      |      | | |     |      |    +------- temperature (ignored)
      |      | | |     |      +------------ exposure: 120 seconds
      |      | | |     +------------------- binning (ignored)
      |      | | +------------------------- sequence number (ignored)
      |      | +--------------------------- FILTER: R, G, B, L, Ha, SII, OIII
      |      +----------------------------- exposure type (ignored)
      +------------------------------------ object name

    Mosaic variant (8 fields + extension): panel number after object name.
    In capture software, set object name to e.g. "M8_02" — the panel number
    becomes an extra field, shifting all subsequent positions by 1:

      M8_02_L_G_11788_Bin1x1_120s_-10C.fit
      |  |  | | |     |      |    |
      |  |  | | |     |      |    +-- temperature (ignored)
      |  |  | | |     |      +------ exposure: 120 seconds
      |  |  | | |     +------------- binning (ignored)
      |  |  | | +------------------- sequence number (ignored)
      |  |  | +--------------------- FILTER: G
      |  |  +----------------------- exposure type (ignored)
      |  +-------------------------- mosaic panel number (ignored)
      +----------------------------- object name

    Detection: 8 parts (excl. extension) = mosaic, 7 parts = normal.

    This naming pattern can be configured in APT, N.I.N.A., Sequence
    Generator Pro, or other imaging capture software.

    Only FILTER and EXPOSURE fields are used; all other fields are ignored.
    """
    basename = os.path.basename(filepath)
    parts = _split_filename(basename)

    # Must end with 'fit' or 'fits'
    if not parts or parts[-1].lower() not in ('fit', 'fits'):
        return None

    # Need 8 or 9 parts (including extension)
    if len(parts) not in (8, 9):
        return None

    # Mosaic shift: filter position moves by 1
    shift = 1 if len(parts) == 9 else 0

    # Filter at position [2 + shift]
    filter_pos = 2 + shift
    filter_str = parts[filter_pos] if filter_pos < len(parts) else ''

    # Map filter string to canonical name
    if filter_str in ('R', 'G', 'B', 'Ha', 'SII', 'OIII'):
        filter_name = filter_str
    elif filter_str.upper() in FILTER_ALIASES:
        filter_name = FILTER_ALIASES[filter_str.upper()]
    else:
        filter_name = 'L'  # Default to Luminance

    # Exposure at position [5 + shift]
    # Value may have 's' suffix (e.g. "120s" or "120")
    exp_pos = 5 + shift
    try:
        exp_str = parts[exp_pos].rstrip('sS')
        exposure = int(exp_str)
    except (ValueError, IndexError):
        return None

    if exposure <= 0:
        return None

    return filter_name, exposure


# ---------------------------------------------------------------
# FITS header-based parsing (default mode)
# ---------------------------------------------------------------

def _normalize_filter(filter_raw):
    """Map raw filter string to canonical name."""
    filter_upper = filter_raw.upper().strip()
    if filter_upper in FILTER_ALIASES:
        return FILTER_ALIASES[filter_upper]
    if filter_raw in ASTROBIN_FILTER_IDS:
        return filter_raw
    return 'L'


def parse_fits_header(filepath):
    """
    Read FITS header to extract filter, exposure, observation date, and object.
    Returns (filter_name, exposure_seconds, obs_date_or_None, object_or_None)
    or None.  Requires astropy.
    """
    try:
        from astropy.io import fits
    except ImportError:
        sys.stderr.write(
            "Error: astropy is required for FITS header reading (default mode).\n"
            "Install with: pip install astropy\n"
            "Or use --parsename to extract info from filenames instead.\n"
        )
        sys.exit(1)

    try:
        with fits.open(filepath, memmap=False) as hdul:
            header = hdul[0].header
    except Exception:
        return None

    # Exposure time
    exptime = header.get('EXPTIME') or header.get('EXPOSURE')
    if exptime is None:
        return None
    try:
        exposure = int(round(float(exptime)))
    except (ValueError, TypeError):
        return None
    if exposure <= 0:
        return None

    # Filter
    filter_raw = str(header.get('FILTER', 'L')).strip()
    filter_name = _normalize_filter(filter_raw)

    # Observation date (optional — fallback to file mtime in caller)
    obs_date = None
    date_obs = header.get('DATE-OBS')
    if date_obs:
        try:
            dt_str = str(date_obs).strip()
            # Parse ISO format: YYYY-MM-DDTHH:MM:SS or YYYY-MM-DD
            if 'T' in dt_str:
                dt = datetime.datetime.fromisoformat(dt_str.split('.')[0])
            else:
                dt = datetime.datetime.strptime(dt_str[:10], '%Y-%m-%d')
                dt = dt.replace(hour=23)  # Assume evening if no time
            # Apply noon rollback for overnight sessions
            if dt.hour < 12:
                dt -= datetime.timedelta(hours=12)
            obs_date = dt.date()
        except (ValueError, TypeError):
            pass

    # Object name (optional)
    obj_name = header.get('OBJECT')
    if obj_name:
        obj_name = str(obj_name).strip()

    return filter_name, exposure, obs_date, obj_name


# ---------------------------------------------------------------
# Directory scanning
# ---------------------------------------------------------------

_FITS_RE = re.compile(r'\.fits?$', re.IGNORECASE)


def scan_files(scan_dir, recursive=True):
    """Find all FITS files in directory."""
    result = []
    if recursive:
        for root, dirs, files in os.walk(scan_dir):
            for fname in files:
                if _FITS_RE.search(fname):
                    result.append(os.path.join(root, fname))
    else:
        try:
            for entry in os.scandir(scan_dir):
                if entry.is_file() and _FITS_RE.search(entry.name):
                    result.append(entry.path)
        except OSError as e:
            sys.stderr.write(f"Error scanning '{scan_dir}': {e}\n")
    return result


# ---------------------------------------------------------------
# Grouping and output
# ---------------------------------------------------------------

def format_time(total_seconds):
    """Format seconds as 'Xh Ym'."""
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    return f"{hours}h {minutes}m"


def generate_report(groups, out, object_name=None):
    """
    Generate AstroBin CSV and summary statistics.

    groups: dict of (date, filter_name, exposure) -> count
    out: file-like object for CSV output
    object_name: object name string (from FITS OBJECT or filename), or None
    """
    if not groups:
        sys.stderr.write("No valid FITS files found.\n")
        return

    # Sort by (date, filter order, exposure)
    def sort_key(item):
        (date, filt, exp), count = item
        filt_idx = CHANNEL_ORDER.index(filt) if filt in CHANNEL_ORDER else 99
        return (date, filt_idx, exp)

    sorted_groups = sorted(groups.items(), key=sort_key)

    # CSV header
    out.write("date,filter,number,duration\n")

    # Statistics accumulators
    total_files = 0
    total_seconds = 0
    # Aggregate by (filter, exposure) across all dates
    agg = defaultdict(int)  # (filter, exposure) -> total count

    for (date, filt, exp), count in sorted_groups:
        filt_id = ASTROBIN_FILTER_IDS.get(filt, 5652)  # Default L
        out.write(f"{date.isoformat()},{filt_id},{count},{exp}\n")

        total_files += count
        total_seconds += count * exp
        agg[(filt, exp)] += count

    # Summary to stderr (so CSV on stdout stays clean)
    sys.stderr.write("\n")
    if object_name:
        sys.stderr.write(f"Object: {object_name}\n")
    sys.stderr.write(f"Total files: {total_files}\n")
    sys.stderr.write(f"Exposure: {total_seconds} ( {format_time(total_seconds)} )\n")
    sys.stderr.write("\n")

    for filt in CHANNEL_ORDER:
        # Collect all exposures for this filter, sorted by exposure time
        filt_entries = sorted(
            [(exp, cnt) for (f, exp), cnt in agg.items() if f == filt]
        )
        for exp, count in filt_entries:
            filt_total = count * exp
            if exp % 60 == 0:
                exp_str = f"{exp // 60}"
            else:
                exp_str = f"{exp / 60:.1f}"
            sys.stderr.write(
                f"{filt}: {count} x {exp_str} min"
                f" ( {format_time(filt_total)} )\n"
            )


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

def main():
    scan_dir, parsename, out_file, flat = parse_args(sys.argv)

    if not os.path.isdir(scan_dir):
        sys.stderr.write(f"Error: '{scan_dir}' is not a directory.\n")
        sys.exit(1)

    # Find files
    files = scan_files(scan_dir, recursive=not flat)
    if not files:
        sys.stderr.write(f"No FITS files found in '{scan_dir}'.\n")
        sys.exit(0)

    sys.stderr.write(f"Scanning {len(files)} FITS files...\n")

    # Parse and group
    groups = defaultdict(int)  # (date, filter, exposure) -> count
    object_names = set()
    skipped = 0

    for filepath in files:
        try:
            if parsename:
                result = parse_filename(filepath)
                if result is None:
                    skipped += 1
                    continue
                filter_name, exposure = result
                obs_date = get_session_date(filepath)
            else:
                result = parse_fits_header(filepath)
                if result is None:
                    skipped += 1
                    continue
                filter_name, exposure, obs_date, obj_name = result
                if obs_date is None:
                    obs_date = get_session_date(filepath)
                if obj_name:
                    object_names.add(obj_name)

            groups[(obs_date, filter_name, exposure)] += 1

        except Exception:
            skipped += 1

    if skipped > 0:
        sys.stderr.write(f"Skipped {skipped} files (unrecognized format).\n")

    # Determine object name (use most common if multiple)
    object_name = None
    if object_names:
        if len(object_names) == 1:
            object_name = next(iter(object_names))
        else:
            object_name = ', '.join(sorted(object_names))

    # Output
    if out_file:
        with open(out_file, 'w', newline='') as f:
            generate_report(groups, f, object_name)
        sys.stderr.write(f"CSV written to: {out_file}\n")
    else:
        generate_report(groups, sys.stdout, object_name)


if __name__ == "__main__":
    main()
