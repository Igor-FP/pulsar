#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import shutil
from typing import List, Tuple, Optional

import numpy as np
from astropy.io import fits
from astropy.time import Time

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils  # noqa: E402
from batch_utils import expand_input_spec as bu_expand_input_spec


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  sortfits.py input_spec output_pattern [options]\n"
        "\n"
        "input_spec:\n"
        "  firstNNNN.fits          numbered sequence starting at NNNN\n"
        "  @list.txt | list.txt    text file with paths (one per line)\n"
        "  single.fits             single file (treated as 1 item)\n"
        "  *.fit, img???.fits      wildcard masks (via batch_utils)\n"
        "\n"
        "output_pattern:\n"
        "  Without --auto:  must contain a numeric field (e.g. out0001.fit)\n"
        "  With --auto:     output directory (created if needed)\n"
        "  With --sessions: base name used; files named as\n"
        "                   <basename>_Sssss_Fffff.fit\n"
        "\n"
        "Modes:\n"
        "  (default)                numbered rename: out0001.fit, out0002.fit, ...\n"
        "  -s, --sessions           session/frame numbering: base_S0001_F0001.fit\n"
        "  --auto                   auto-naming from FITS headers:\n"
        "                           {OBJECT}_exp{TIME}_{FILTER}[_S{SSSS}]_{NNNN}.fit\n"
        "\n"
        "Options:\n"
        "  --auto                 auto-name from FITS headers (OBJECT, EXPTIME, FILTER)\n"
        "  --name NAME            override object name (default: OBJECT header or 'obj')\n"
        "  --group-num            number within each group instead of globally\n"
        "  -s, --sessions         enable session splitting (works with --auto too)\n"
        "  --gap-hours H          gap in hours to split sessions (default: 1.0)\n"
        "\n"
        "Examples:\n"
        "  sortfits *.fit out0001.fit\n"
        "  sortfits *.fit out.fit --sessions --gap-hours 2\n"
        "  sortfits *.fit sorted/ --auto\n"
        "  sortfits *.fit sorted/ --auto --sessions --name NGC1097\n"
        "  sortfits *.fit sorted/ --auto --group-num\n"
        "\n"
    )
    sys.exit(1)


# ---------------- Small helpers ---------------- #


def read_list_file(path: str) -> List[str]:
    out = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            out.append(line)
    return out


def expand_input_spec(input_spec: str) -> List[str]:
    """
    Expand input specification:
      - @list.txt / list.txt / list.lst
      - numbered sequence via batch_utils.find_numbered_sequence
      - single file
    """
    if input_spec.startswith("@"):
        return read_list_file(input_spec[1:])

    if input_spec.lower().endswith((".txt", ".lst")) and os.path.isfile(input_spec):
        return read_list_file(input_spec)

    if os.path.isfile(input_spec):
        # Try numbered sequence like prefixNNNN.fit(s)
        base = os.path.basename(input_spec)
        m = re.match(r"^(.*?)(\d+)(\.(?:fit|fits))$", base, flags=re.IGNORECASE)
        if m:
            seq = batch_utils.find_numbered_sequence(os.path.abspath(input_spec))
            if seq:
                return [name for (name, _, _) in seq]
        # Fallback: single file
        return [os.path.abspath(input_spec)]

    raise FileNotFoundError(f"Input spec not recognized or file not found: {input_spec}")


# ---------------- Time extraction ---------------- #


def _get_float(header, key) -> Optional[float]:
    try:
        v = header.get(key)
        if v is None:
            return None
        return float(v)
    except Exception:
        return None


def _get_str(header, key) -> Optional[str]:
    try:
        v = header.get(key)
        if v is None:
            return None
        s = str(v).strip()
        return s if s else None
    except Exception:
        return None


def _parse_date_obs(header) -> Optional[float]:
    """
    Try to parse DATE-OBS using astropy, return MJD (float) or None.
    """
    date_obs = header.get("DATE-OBS") or header.get("DATE")
    if not date_obs:
        return None

    try:
        # Astropy Time can handle many formats
        # We convert to MJD in days.
        t = Time(date_obs, format="isot", scale="utc")
        return float(t.mjd)
    except Exception:
        pass

    try:
        t = Time(date_obs, scale="utc")
        return float(t.mjd)
    except Exception:
        return None


def extract_mjd(path: str) -> Optional[float]:
    """
    Extract MJD for exposure start or mid-exposure if possible.

    Strategy:
      1) EXPMJD (if present)
      2) MJD-OBS
      3) DATE-OBS + (EXPTIME / 2)/86400  (mid-exposure)
      4) DATE-OBS as MJD
    """
    try:
        with fits.open(path, memmap=True) as hdul:
            header = hdul[0].header
    except Exception:
        return None

    # 1) EXPMJD - direct MJD
    expmjd = _get_float(header, "EXPMJD")
    if expmjd is not None:
        return expmjd

    # 2) MJD-OBS
    mjd_obs = _get_float(header, "MJD-OBS")
    if mjd_obs is not None:
        return mjd_obs

    # 3) DATE-OBS + half exposure
    date_mjd = _parse_date_obs(header)
    if date_mjd is not None:
        exptime = _get_float(header, "EXPTIME") or _get_float(header, "EXPOSURE")
        if exptime is not None:
            # Mid-exposure: start + exptime/2
            return date_mjd + (exptime / 2.0) / 86400.0
        return date_mjd

    return None


def extract_metadata(path: str) -> Tuple[Optional[float], Optional[str], Optional[float]]:
    """
    Extract (mjd, filter_name, exptime) from FITS header.
    """
    try:
        with fits.open(path, memmap=True) as hdul:
            header = hdul[0].header
    except Exception:
        return None, None, None

    # MJD
    mjd = None
    expmjd = _get_float(header, "EXPMJD")
    if expmjd is not None:
        mjd = expmjd
    else:
        mjd_obs = _get_float(header, "MJD-OBS")
        if mjd_obs is not None:
            mjd = mjd_obs
        else:
            date_mjd = _parse_date_obs(header)
            if date_mjd is not None:
                exptime = _get_float(header, "EXPTIME") or _get_float(header, "EXPOSURE")
                if exptime is not None:
                    mjd = date_mjd + (exptime / 2.0) / 86400.0
                else:
                    mjd = date_mjd

    # Filter
    filt = _get_str(header, "FILTER") or _get_str(header, "FILTER1")

    # Exptime
    exptime = _get_float(header, "EXPTIME") or _get_float(header, "EXPOSURE")

    # Object name
    obj = _get_str(header, "OBJECT")

    return mjd, filt, exptime, obj


# ---------------- Sorting logic ---------------- #


def sort_by_jd(files: List[str]) -> List[Tuple[str, float]]:
    """
    Return list of (path, mjd) sorted by mjd.
    Files with no time info are placed at the end, in original order.
    """
    items: List[Tuple[str, Optional[float]]] = []
    for path in files:
        mjd = extract_mjd(path)
        items.append((path, mjd))

    # Keep original order for those with no MJD
    with_mjd = [(p, mjd) for (p, mjd) in items if mjd is not None]
    without_mjd = [(p, mjd) for (p, mjd) in items if mjd is None]

    with_mjd.sort(key=lambda x: x[1])

    # For those without MJD, just keep original order
    return with_mjd + without_mjd


def detect_sessions(sorted_items: List[Tuple[str, float]], gap_hours: float) -> List[List[Tuple[str, float]]]:
    """
    Split into sessions by gap in time.
    Only items with a valid MJD are considered for gaps; items without MJD
    are kept in the same session as the previous one (if any), or in the first.
    """
    if not sorted_items:
        return []

    sessions: List[List[Tuple[str, float]]] = []
    current: List[Tuple[str, float]] = []
    last_mjd: Optional[float] = None

    gap_days = gap_hours / 24.0

    for (path, mjd) in sorted_items:
        if mjd is None:
            # No time info: keep in current session
            current.append((path, mjd))
            continue

        if last_mjd is None:
            # First with MJD
            current.append((path, mjd))
            last_mjd = mjd
            continue

        # Check gap
        if mjd - last_mjd > gap_days:
            # New session
            sessions.append(current)
            current = [(path, mjd)]
        else:
            current.append((path, mjd))

        last_mjd = mjd

    if current:
        sessions.append(current)

    return sessions


# ---------------- Output naming ---------------- #


def validate_output_pattern(pattern: str) -> Tuple[str, int, int]:
    """
    Validate output pattern for non-session mode:
    It must contain a sequence of digits (e.g. out0001.fit).
    Return (pattern, index_of_digits_start, num_digits).
    """
    base = os.path.basename(pattern)
    m = re.search(r"(\d+)", base)
    if not m:
        raise ValueError(
            "Output pattern must contain a numeric field, e.g. out0001.fit"
        )
    num_digits = len(m.group(1))
    start_idx = m.start(1)
    return pattern, start_idx, num_digits


def make_output_name(pattern: str, index: int) -> str:
    """
    Substitute index into numeric field of pattern (non-session mode).
    """
    base = os.path.basename(pattern)
    m = re.search(r"(\d+)", base)
    assert m is not None
    num_digits = len(m.group(1))
    start_idx = m.start(1)
    end_idx = m.end(1)

    number_str = f"{index:0{num_digits}d}"
    new_base = base[:start_idx] + number_str + base[end_idx:]
    return new_base


def make_session_name(base_pattern: str, session_idx: int, frame_idx: int) -> str:
    """
    Session/frame mode: generate <basename>_Sssss_Fffff.fit
    with 4-digit session and 4-digit frame numbers.
    """
    base = os.path.basename(base_pattern)
    root, ext = os.path.splitext(base)
    if not ext:
        ext = ".fit"

    s_str = f"{session_idx:04d}"
    f_str = f"{frame_idx:04d}"
    return f"{root}_S{s_str}_F{f_str}{ext}"


def format_exptime(exptime: Optional[float]) -> str:
    """Format exposure time for auto-naming: exp120, exp500ms, exp0."""
    if exptime is None:
        return "exp0"
    if exptime < 1.0:
        ms = int(round(exptime * 1000))
        return f"exp{ms}ms"
    return f"exp{int(round(exptime))}"


def sanitize_name(name: str) -> str:
    """Remove characters unsafe for filenames."""
    # Keep alphanumeric, dash, underscore, dot, plus
    s = re.sub(r'[^\w\-+.]', '_', name)
    # Collapse multiple underscores
    s = re.sub(r'_+', '_', s)
    return s.strip('_')


def make_auto_name(obj_name: str, exptime_str: str, filt: str,
                   session_idx: Optional[int], num: int, num_width: int,
                   ext: str) -> str:
    """
    Build auto filename: {obj}_{exptime}_{filter}[_S{session}]_{num}.{ext}
    """
    parts = [obj_name, exptime_str, filt]
    if session_idx is not None:
        parts.append(f"S{session_idx:04d}")
    parts.append(f"{num:0{num_width}d}")
    return "_".join(parts) + ext


def copy_or_link(src: str, dst: str) -> None:
    """
    Copy file (using shutil.copy2). Could be replaced by hardlink or symlink
    in the future if desired.
    """
    os.makedirs(os.path.dirname(dst) if os.path.dirname(dst) else ".", exist_ok=True)
    shutil.copy2(src, dst)


# ---------------- Num width helper ---------------- #

def num_width(count: int) -> int:
    """Determine zero-padding width for count items (minimum 4)."""
    if count <= 0:
        return 4
    w = len(str(count))
    return max(w, 4)


# ---------------- Main script ---------------- #


def main():
    # Parse args
    args = sys.argv[1:]

    sessions_mode = False
    auto_mode = False
    group_num = False
    gap_hours = 1.0
    override_name = None

    positional = []
    i = 0
    while i < len(args):
        a = args[i]
        if a in ("-s", "--sessions"):
            sessions_mode = True
            i += 1
        elif a == "--auto":
            auto_mode = True
            i += 1
        elif a == "--group-num":
            group_num = True
            i += 1
        elif a == "--name":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --name requires a value.\n")
                sys.exit(1)
            override_name = args[i + 1]
            i += 2
        elif a == "--gap-hours":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --gap-hours requires a value.\n")
                sys.exit(1)
            try:
                gap_hours = float(args[i + 1])
            except ValueError:
                sys.stderr.write("Error: --gap-hours must be a number.\n")
                sys.exit(1)
            i += 2
        elif a.startswith("--") or a == "-s":
            sys.stderr.write(f"Unknown option: {a}\n")
            usage()
        else:
            positional.append(a)
            i += 1

    if auto_mode:
        if len(positional) < 2:
            sys.stderr.write("Error: --auto requires input_spec and output_dir.\n")
            usage()
        input_spec = positional[0]
        output_dir = positional[1]
    elif sessions_mode:
        if len(positional) < 2:
            sys.stderr.write("Error: --sessions requires input_spec and output_pattern.\n")
            usage()
        input_spec = positional[0]
        output_pattern = positional[1]
    else:
        if len(positional) < 2:
            usage()
        input_spec = positional[0]
        output_pattern = positional[1]

    # Validate incompatible options
    if not auto_mode and (group_num or override_name is not None):
        sys.stderr.write("Error: --group-num and --name require --auto mode.\n")
        sys.exit(1)

    # Expand input files
    try:
        files = bu_expand_input_spec(input_spec)
    except Exception as e:
        sys.stderr.write(f"Error expanding input: {e}\n")
        sys.exit(1)

    if not files:
        sys.stderr.write("No input files.\n")
        sys.exit(1)

    print(f"Sorting {len(files)} file(s)...")

    # ---- AUTO MODE ----
    if auto_mode:
        _run_auto_mode(files, output_dir, override_name, sessions_mode,
                       gap_hours, group_num)
        return

    # ---- LEGACY MODES ----

    # Sort by JD
    try:
        items = sort_by_jd(files)
    except Exception as e:
        sys.stderr.write(f"Error reading times from FITS files: {e}\n")
        sys.exit(1)

    if not sessions_mode:
        # Non-session mode: just rename/copy in sorted order
        try:
            _, _, _ = validate_output_pattern(output_pattern)
        except Exception as e:
            sys.stderr.write(f"Invalid output pattern: {e}\n")
            sys.exit(1)

        outdir = os.path.dirname(output_pattern) or "."
        index = 1
        for (path, mjd) in items:
            out_name = make_output_name(output_pattern, index)
            dst = os.path.join(outdir, out_name)
            print(f"{index:04d}  {path}  ->  {dst}  (MJD={mjd})")
            copy_or_link(path, dst)
            index += 1
    else:
        # Session mode: split into sessions, number sessions and frames
        sorted_with_mjd = [(p, mjd) for (p, mjd) in items if mjd is not None]
        if not sorted_with_mjd:
            sys.stderr.write(
                "Warning: no valid MJD found, all files have no time information.\n"
                "Session mode will treat all files as one session.\n"
            )
            sessions = [items]
        else:
            sessions = detect_sessions(items, gap_hours)

        outdir = os.path.dirname(output_pattern) or "."
        base_pattern = output_pattern

        session_idx = 1
        for session in sessions:
            frame_idx = 1
            print(f"Session {session_idx:04d}: {len(session)} files")
            for (path, mjd) in session:
                out_name = make_session_name(base_pattern, session_idx, frame_idx)
                dst = os.path.join(outdir, out_name)
                print(
                    f"  S{session_idx:04d} F{frame_idx:04d}: {path} -> {dst}  (MJD={mjd})"
                )
                copy_or_link(path, dst)
                frame_idx += 1
            session_idx += 1

    print("Done.")


# ---------------- Auto mode ---------------- #


def _run_auto_mode(files, output_dir, override_name, sessions_mode,
                   gap_hours, group_num):
    """
    Auto-naming mode: read FITS headers, build names from metadata.
    Format: {OBJECT}_exp{TIME}_{FILTER}[_S{SSSS}]_{NNNN}.fit
    """
    # Read metadata from all files
    file_meta = []  # [(path, mjd, filter, exptime, object), ...]
    for idx, fpath in enumerate(files, 1):
        mjd, filt, exptime, obj = extract_metadata(fpath)
        file_meta.append((fpath, mjd, filt, exptime, obj))
        sys.stderr.write(f"\r  Reading headers: {idx}/{len(files)}")
        sys.stderr.flush()
    sys.stderr.write("\n")

    # Determine object name
    if override_name:
        obj_name = sanitize_name(override_name)
    else:
        # Use most common OBJECT value from headers
        obj_counts = {}
        for _, _, _, _, obj in file_meta:
            if obj:
                key = obj.strip()
                obj_counts[key] = obj_counts.get(key, 0) + 1
        if obj_counts:
            obj_name = sanitize_name(max(obj_counts, key=obj_counts.get))
        else:
            obj_name = "obj"

    # Sort by MJD
    file_meta.sort(key=lambda x: (x[1] is None, x[1] or 0))

    # Build group key for each file: (exptime_str, filter)
    def group_key(meta):
        _, _, filt, exptime, _ = meta
        return (format_exptime(exptime), filt or "nofilter")

    # Determine extension from first input file
    _, ext = os.path.splitext(files[0])
    if not ext:
        ext = ".fit"

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    if not sessions_mode:
        # No sessions: just sort by JD and number
        if group_num:
            # Number within each (exptime, filter) group
            groups = {}
            for meta in file_meta:
                gk = group_key(meta)
                if gk not in groups:
                    groups[gk] = []
                groups[gk].append(meta)

            total = 0
            for gk in sorted(groups.keys()):
                members = groups[gk]
                nw = num_width(len(members))
                exptime_str, filt = gk
                for i, meta in enumerate(members, 1):
                    fpath = meta[0]
                    out_name = make_auto_name(obj_name, exptime_str, filt,
                                              None, i, nw, ext)
                    dst = os.path.join(output_dir, out_name)
                    print(f"  {fpath}  ->  {out_name}")
                    copy_or_link(fpath, dst)
                    total += 1
        else:
            # Global numbering across all files
            nw = num_width(len(file_meta))
            for i, meta in enumerate(file_meta, 1):
                fpath, mjd, filt, exptime, _ = meta
                exptime_str = format_exptime(exptime)
                filt_str = filt or "nofilter"
                out_name = make_auto_name(obj_name, exptime_str, filt_str,
                                          None, i, nw, ext)
                dst = os.path.join(output_dir, out_name)
                print(f"  {fpath}  ->  {out_name}")
                copy_or_link(fpath, dst)

    else:
        # With sessions: split by time gaps, then name with session number
        sorted_items = [(m[0], m[1]) for m in file_meta]
        sessions = detect_sessions(sorted_items, gap_hours)

        # Re-index file_meta by path for quick lookup
        meta_by_path = {m[0]: m for m in file_meta}

        if group_num:
            # Number within each (exptime, filter, session) group
            total = 0
            for s_idx, session in enumerate(sessions, 1):
                groups = {}
                for (path, mjd) in session:
                    meta = meta_by_path[path]
                    gk = group_key(meta)
                    if gk not in groups:
                        groups[gk] = []
                    groups[gk].append(meta)

                print(f"Session {s_idx:04d}: {len(session)} files")
                for gk in sorted(groups.keys()):
                    members = groups[gk]
                    nw = num_width(len(members))
                    exptime_str, filt = gk
                    for i, meta in enumerate(members, 1):
                        fpath = meta[0]
                        out_name = make_auto_name(obj_name, exptime_str, filt,
                                                  s_idx, i, nw, ext)
                        dst = os.path.join(output_dir, out_name)
                        print(f"    {fpath}  ->  {out_name}")
                        copy_or_link(fpath, dst)
                        total += 1
        else:
            # Global numbering across all files (but with session tag)
            nw = num_width(len(file_meta))
            global_idx = 1
            for s_idx, session in enumerate(sessions, 1):
                print(f"Session {s_idx:04d}: {len(session)} files")
                for (path, mjd) in session:
                    meta = meta_by_path[path]
                    _, _, filt, exptime, _ = meta
                    exptime_str = format_exptime(exptime)
                    filt_str = filt or "nofilter"
                    out_name = make_auto_name(obj_name, exptime_str, filt_str,
                                              s_idx, global_idx, nw, ext)
                    dst = os.path.join(output_dir, out_name)
                    print(f"    {fpath}  ->  {out_name}")
                    copy_or_link(fpath, dst)
                    global_idx += 1

    # Summary
    groups_found = set()
    for meta in file_meta:
        gk = group_key(meta)
        groups_found.add(gk)

    print(f"\nDone. {len(file_meta)} files -> {output_dir}/")
    print(f"  Object: {obj_name}")
    print(f"  Groups: {len(groups_found)}")
    for exptime_str, filt in sorted(groups_found):
        count = sum(1 for m in file_meta if group_key(m) == (exptime_str, filt))
        print(f"    {exptime_str}_{filt}: {count} files")


if __name__ == "__main__":
    main()
