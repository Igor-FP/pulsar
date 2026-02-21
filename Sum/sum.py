#!/usr/bin/env python3
# sum.py

import sys
import os
import re
import glob
from typing import List, Optional

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
import warnings
from astropy.io.fits.verify import VerifyWarning

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils

# Suppress FITS "card too long" warnings
warnings.filterwarnings("ignore", category=VerifyWarning)


def parse_time(s: str) -> Optional[Time]:
    if not s:
        return None
    s = s.strip()
    try:
        return Time(s, format="isot", scale="utc")
    except Exception:
        try:
            return Time(s, scale="utc")
        except Exception:
            return None


def get_start_time(hdr) -> Optional[Time]:
    if "DATE-OBS" in hdr:
        t = parse_time(hdr["DATE-OBS"])
        if t is not None:
            return t
    d = hdr.get("DATE-OBS")
    tm = hdr.get("TIME-OBS")
    if d and tm:
        t = parse_time(f"{d}T{tm}")
        if t is not None:
            return t
    cam = hdr.get("HIERARCH CAMERA-DATE-OBS")
    if cam:
        t = parse_time(cam)
        if t is not None:
            return t
    return None


def get_exptime(hdr) -> Optional[float]:
    for k in ("EXPTIME", "EXPOSURE"):
        if k in hdr:
            try:
                return float(hdr[k])
            except Exception:
                pass
    return None


def bitpix_to_dtype(bitpix: int):
    if bitpix == 8:
        return np.uint8
    if bitpix == 16:
        return np.int16
    if bitpix == 32:
        return np.int32
    if bitpix == 64:
        return np.int64
    if bitpix == -32:
        return np.float32
    if bitpix == -64:
        return np.float64
    raise ValueError(f"Unsupported BITPIX={bitpix}")


def print_progress(idx: int, total: int, fname: str, bar_len: int = 30):
    frac = idx / total if total > 0 else 1.0
    fill = int(bar_len * frac)
    bar = "#" * fill + "-" * (bar_len - fill)
    msg = f"\r[{bar}] {idx}/{total} {fname}"
    sys.stdout.write(msg + " " * max(0, 80 - len(msg)))
    sys.stdout.flush()


def collect_input_files(first_input: str) -> List[str]:
    """
    Collect input files for summation.

    New behavior:
      - If first_input contains '*' or '?', treat as wildcard mask
        (glob all .fit/.fits files, sorted alphabetically).

    Old behavior (preserved):
      - If name contains digits, use numeric sequence via batch_utils.find_numbered_sequence().
      - If file exists, use single file.
    """
    # Wildcard mask
    if "*" in first_input or "?" in first_input:
        matched = sorted(
            f for f in glob.glob(first_input)
            if os.path.isfile(f) and f.lower().endswith((".fit", ".fits"))
        )
        if not matched:
            raise FileNotFoundError(f"No input files match pattern '{first_input}'")
        return matched

    # Numeric sequence (old logic)
    if re.search(r"\d+", first_input):
        seq = batch_utils.find_numbered_sequence(first_input)
        if seq:
            return [name for (name, _, _) in seq]
        if os.path.exists(first_input):
            return [first_input]
        raise FileNotFoundError(
            f"No input files found matching sequence starting from {first_input}"
        )

    # Fallback: single file
    if os.path.exists(first_input):
        return [first_input]

    raise FileNotFoundError(f"Input file not found: {first_input}")


def main():
    if len(sys.argv) != 3:
        print("Usage: sum.py <first_input.fit | mask*.fit> <output.fit>", file=sys.stderr)
        sys.exit(1)

    first_input = sys.argv[1]
    output_name = sys.argv[2]
    if not re.search(r"\.[Ff][Ii][Tt][Ss]?$", output_name):
        output_name += ".fit"

    try:
        files = collect_input_files(first_input)
    except Exception as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)

    if not files:
        print(f"No input files found matching {first_input}", file=sys.stderr)
        sys.exit(1)

    total_files = len(files)

    with fits.open(files[0], memmap=False) as hdul:
        if hdul[0].data is None:
            print("First FITS file has no primary image data.", file=sys.stderr)
            sys.exit(1)
        data0 = hdul[0].data
        base_header = hdul[0].header.copy()

    bitpix = int(base_header.get("BITPIX"))
    base_dtype = bitpix_to_dtype(bitpix)
    is_int = bitpix > 0

    sum_data = np.zeros_like(data0, dtype=np.float64)

    total_exptime = 0.0
    have_exptime = False
    start_times: List[Time] = []
    end_times: List[Time] = []

    for idx, fname in enumerate(files, start=1):
        with fits.open(fname, memmap=False) as hdul:
            data = hdul[0].data
            hdr = hdul[0].header

            if data is None:
                print(f"\nWarning: {fname} has no primary image data, skipping.", file=sys.stderr)
                continue

            if data.shape != sum_data.shape:
                print(f"\nError: shape mismatch in {fname}: {data.shape} != {sum_data.shape}", file=sys.stderr)
                sys.exit(1)

            bpx = int(hdr.get("BITPIX"))
            if bpx != bitpix:
                print(f"\nError: BITPIX mismatch in {fname}: {bpx} != {bitpix}", file=sys.stderr)
                sys.exit(1)

            sum_data += data.astype(np.float64, copy=False)

            exptime = get_exptime(hdr)
            if exptime is not None:
                total_exptime += exptime
                have_exptime = True

            st = get_start_time(hdr)
            if st is not None:
                start_times.append(st)

            et = None
            if st is not None and exptime is not None:
                try:
                    et = st + TimeDelta(exptime, format="sec")
                except Exception:
                    et = None
            if et is None:
                for k in ("DATE-END", "END-OBS", "DATE_OBS_END"):
                    if k in hdr:
                        t = parse_time(hdr[k])
                        if t is not None:
                            et = t
                            break
            if et is not None:
                end_times.append(et)

        print_progress(idx, total_files, os.path.basename(fname))

    print()
    n = total_files
    out_header = base_header.copy()

    if have_exptime:
        for k in ("EXPTIME", "EXPOSURE"):
            if k in out_header:
                out_header[k] = total_exptime
    if "SNAPSHOT" in out_header:
        out_header["SNAPSHOT"] = n

    earliest_start = min(start_times) if start_times else None
    latest_end = max(end_times) if end_times else None
    if earliest_start is not None:
        iso = earliest_start.isot
        out_header["DATE-OBS"] = iso
        if "TIME-OBS" in out_header:
            out_header["TIME-OBS"] = iso.split("T", 1)[1]
        if "HIERARCH CAMERA-DATE-OBS" in out_header:
            out_header["HIERARCH CAMERA-DATE-OBS"] = iso
    if earliest_start and latest_end:
        mid = earliest_start + 0.5 * (latest_end - earliest_start)
        mid_iso = mid.isot
        end_iso = latest_end.isot
        for k in ("DATE-MID", "MID-OBS", "MIDTIME", "DATE-AVG"):
            if k in out_header:
                out_header[k] = mid_iso
        for k in ("DATE-END", "END-OBS", "DATE_OBS_END"):
            if k in out_header:
                out_header[k] = end_iso

    if is_int:
        dtype_out = base_dtype
        info = np.iinfo(dtype_out)
        data_scaled = sum_data.copy()
        if info.min < 0:
            data_scaled[data_scaled < 0.0] = 0.0
        cur_max = float(data_scaled.max()) if data_scaled.size > 0 else 0.0
        if cur_max <= 0.0:
            out_data = np.zeros_like(sum_data, dtype=dtype_out)
        else:
            scale = info.max / cur_max
            data_scaled *= scale
            low = 0 if info.min < 0 else info.min
            data_scaled = np.clip(data_scaled, low, info.max)
            out_data = data_scaled.astype(dtype_out)
        out_header["BITPIX"] = bitpix
        for k in ("BZERO", "BSCALE"):
            out_header.pop(k, None)
    else:
        avg_data = sum_data / float(n) if n > 0 else sum_data
        if bitpix == -32:
            out_data = avg_data.astype(np.float32)
        elif bitpix == -64:
            out_data = avg_data.astype(np.float64)
        else:
            out_data = avg_data
        out_header["BITPIX"] = bitpix
        for k in ("BZERO", "BSCALE"):
            out_header.pop(k, None)

    fits.HDUList([fits.PrimaryHDU(data=out_data, header=out_header)]).writeto(
        output_name, overwrite=True, output_verify="silentfix"
    )

    print(f"Processed {n} files into {output_name}")


if __name__ == "__main__":
    main()
