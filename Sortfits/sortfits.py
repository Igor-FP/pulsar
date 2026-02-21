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
        "  sortfits.py input_spec output_pattern [-s|--sessions] [--gap-hours H]\n"
        "\n"
        "input_spec:\n"
        "  firstNNNN.fits          numbered sequence starting at NNNN\n"
        "  @list.txt | list.txt    text file with paths (one per line)\n"
        "  single.fits             single file (treated as 1 item)\n"
        "  *.fit, img???.fits      wildcard masks (via batch_utils)\n"
        "\n"
        "output_pattern:\n"
        "  Without --sessions: must contain a numeric field (e.g. out0001.fit)\n"
        "  With --sessions:     base name used; files will be named as\n"
        "                       <basename>_Sssss_Fffff.fit\n"
        "                       ssss - session number (0001..9999)\n"
        "                       ffff - frame number in session (0001..9999)\n"
        "\n"
        "Options:\n"
        "  -s, --sessions         enable session/frame numbering mode\n"
        "  --gap-hours H          gap in hours to split sessions (default: 1.0)\n"
        "\n"
    )


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
            if not current:
                current.append((path, mjd))
            else:
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


def copy_or_link(src: str, dst: str) -> None:
    """
    Copy file (using shutil.copy2). Could be replaced by hardlink or symlink
    in the future if desired.
    """
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)


# ---------------- Main script ---------------- #


def main():
    if len(sys.argv) < 3:
        usage()

    # Parse args (simple manual parsing to keep small and robust)
    args = sys.argv[1:]
    input_spec = args[0]
    output_pattern = args[1]

    sessions_mode = False
    gap_hours = 1.0

    i = 2
    while i < len(args):
        a = args[i]
        if a in ("-s", "--sessions"):
            sessions_mode = True
            i += 1
            continue
        if a == "--gap-hours":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --gap-hours requires a value.\n")
                sys.exit(1)
            try:
                gap_hours = float(args[i + 1])
            except ValueError:
                sys.stderr.write("Error: --gap-hours must be a number.\n")
                sys.exit(1)
            i += 2
            continue
        sys.stderr.write(f"Unknown argument: {a}\n")
        usage()

    try:
        files = bu_expand_input_spec(input_spec)
    except Exception as e:
        sys.stderr.write(f"Error expanding input: {e}\n")
        sys.exit(1)

    if not files:
        sys.stderr.write("No input files.\n")
        sys.exit(1)

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


if __name__ == "__main__":
    main()
