#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
raw2fits - Convert Camera RAW files to FITS with raw Bayer data. Currently supports Canon CR2/CR3.

Default mode: raw CFA Bayer mosaic (2D uint16, no processing).
With --debaer: linear demosaicing to 3-channel RGB FITS (uint16).
"""

import sys
import os
import re
import glob
import io
import struct
import contextlib
import numpy as np
from astropy.io import fits

try:
    import rawpy
except ImportError:
    sys.stderr.write(
        "Error: rawpy is required but not installed.\n"
        "Install it with: pip install rawpy\n"
    )
    sys.exit(1)

try:
    import exifread
except ImportError:
    sys.stderr.write(
        "Error: exifread is required but not installed.\n"
        "Install it with: pip install exifread\n"
    )
    sys.exit(1)

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


# ---------------------------------------------------------------
# Constants
# ---------------------------------------------------------------

_COLOR_CHAR = {0: 'R', 1: 'G', 2: 'B'}

_RAW_SEQ_RE = re.compile(r"^(.*?)(\d+)(\.cr[23])$", re.IGNORECASE)
_FITS_SEQ_RE = re.compile(r"^(.*?)(\d+)(\.(?:fit|fits))$", re.IGNORECASE)


# ---------------------------------------------------------------
# CLI
# ---------------------------------------------------------------

def usage():
    sys.stderr.write(
        "Usage:\n"
        "  raw2fits.py [--debaer] [--debug] input_spec output_spec\n"
        "\n"
        "    input_spec   - single RAW file, wildcard (*.cr2), numbered pattern\n"
        "                   (img0001.cr2), or @list.txt\n"
        "    output_spec  - single file (e.g. image.fit) OR numbered pattern\n"
        "                   (e.g. out0001.fit) OR wildcard (e.g. *.fit)\n"
        "                   Wildcard: * is replaced with input filename stem\n"
        "\n"
        "    --debaer     - demosaic (debayer) to linear RGB instead of raw CFA\n"
        "    --debug      - write detailed log to <input_name>.log for each file\n"
        "\n"
        "  Default: raw Bayer CFA mosaic (2D monochrome uint16, original sensor values)\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]
    debaer = False
    debug = False

    if "--debaer" in args:
        debaer = True
        args.remove("--debaer")

    if "--debug" in args:
        debug = True
        args.remove("--debug")

    if len(args) != 2:
        usage()

    return args[0], args[1], debaer, debug


# ---------------------------------------------------------------
# RAW input expansion (batch_utils only handles .fit/.fits sequences)
# ---------------------------------------------------------------

def _find_raw_sequence(first_path):
    """Discover contiguous numbered RAW sequence starting from first_path."""
    base = os.path.basename(first_path)
    m = _RAW_SEQ_RE.match(base)
    if not m:
        return []

    prefix, digits, ext = m.group(1), m.group(2), m.group(3)
    width = len(digits)
    start_index = int(digits)
    dir_ = os.path.dirname(os.path.abspath(first_path)) or "."

    seq = []
    idx = start_index
    while True:
        name = f"{prefix}{str(idx).zfill(width)}{ext}"
        full = os.path.join(dir_, name)
        if not os.path.isfile(full):
            break
        seq.append(full)
        idx += 1

    return seq


def expand_raw_input_spec(spec):
    """Expand input spec supporting RAW file numbered sequences."""
    spec = spec.strip()

    # Wildcard
    if "*" in spec or "?" in spec:
        matches = glob.glob(spec)
        files = [os.path.abspath(m) for m in matches if os.path.isfile(m)]
        files.sort()
        if not files:
            raise FileNotFoundError(f"No files match wildcard pattern: {spec}")
        return files

    # @list.txt
    if spec.startswith("@"):
        return batch_utils._expand_list_file(spec[1:])

    # list.txt / list.lst
    if spec.lower().endswith((".txt", ".lst")) and os.path.isfile(spec):
        return batch_utils._expand_list_file(spec)

    # Single file or numbered sequence
    if os.path.isfile(spec):
        seq = _find_raw_sequence(spec)
        if seq:
            return seq
        return [os.path.abspath(spec)]

    raise FileNotFoundError(f"Input not recognized or file not found: {spec}")


# ---------------------------------------------------------------
# Output pairing
# ---------------------------------------------------------------

def _apply_wildcard_output(inputs, output_spec):
    """Replace * in output_spec with each input file's stem."""
    out_dir = os.path.dirname(output_spec)
    out_pattern = os.path.basename(output_spec)

    pairs = []
    for inp in inputs:
        stem = os.path.splitext(os.path.basename(inp))[0]
        outname = out_pattern.replace("*", stem)
        if out_dir:
            outpath = os.path.join(out_dir, outname)
        else:
            outpath = outname
        pairs.append((inp, outpath))
    return pairs


def build_io_pairs(input_spec, output_spec):
    """Build (input, output) pairs for RAW->FITS conversion."""
    inputs = expand_raw_input_spec(input_spec)
    if not inputs:
        raise ValueError("No input files found.")

    # Wildcard output: *.fit, prefix_*.fits, etc.
    if "*" in output_spec:
        return _apply_wildcard_output(inputs, output_spec)

    # Single input -> single output
    if len(inputs) == 1:
        return [(inputs[0], output_spec)]

    # Multiple inputs => try numbered FITS pattern
    base = os.path.basename(output_spec)
    m = _FITS_SEQ_RE.match(base)
    if m:
        prefix, digits, ext = m.group(1), m.group(2), m.group(3)
        width = len(digits)
        start_index = int(digits)
        out_dir = os.path.dirname(os.path.abspath(output_spec)) or "."

        pairs = []
        idx = start_index
        for inp in inputs:
            fname = f"{prefix}{str(idx).zfill(width)}{ext}"
            pairs.append((inp, os.path.join(out_dir, fname)))
            idx += 1
        return pairs

    raise ValueError(
        "Output pattern must contain a numeric field when multiple "
        "input files are provided (e.g. out0001.fit), or use "
        "wildcard (e.g. *.fit) to preserve input names."
    )


# ---------------------------------------------------------------
# Minimal TIFF IFD parser (fallback when exifread can't follow sub-IFDs)
# ---------------------------------------------------------------

# EXIF tag IDs
_TAG_MAKE = 0x010F
_TAG_MODEL = 0x0110
_TAG_DATETIME = 0x0132
_TAG_EXIF_IFD_PTR = 0x8769
_TAG_EXPOSURE_TIME = 0x829A
_TAG_FNUMBER = 0x829D
_TAG_ISO = 0x8827
_TAG_DATETIME_ORIG = 0x9003
_TAG_FOCAL_35MM = 0xA405

# Canon MakerNote tag IDs
_TAG_CANON_SHOT_INFO = 0x0004  # ShotInfo array, index 12 = CameraTemperature + 128


def _parse_tiff_tags(tiff_data, log=None):
    """
    Minimal TIFF IFD parser. Reads IFD0 and follows ExifIFD pointer.
    Returns dict with exifread-compatible key names and plain Python values
    (str, int, float — not exifread IfdTag objects).
    """
    def _log(msg):
        if log:
            log.write(msg + "\n")

    if len(tiff_data) < 8:
        return {}

    # Byte order
    bom = tiff_data[:2]
    if bom == b'II':
        endian = '<'
    elif bom == b'MM':
        endian = '>'
    else:
        _log(f"[tiff_parse] Unknown byte order: {bom.hex()}")
        return {}

    magic = struct.unpack(endian + 'H', tiff_data[2:4])[0]
    if magic != 42:
        _log(f"[tiff_parse] Bad TIFF magic: {magic}")
        return {}

    ifd0_offset = struct.unpack(endian + 'I', tiff_data[4:8])[0]
    _log(f"[tiff_parse] Byte order={'LE' if endian == '<' else 'BE'},"
         f" IFD0 at offset {ifd0_offset}")

    def _read_ifd_entries(offset):
        """Read raw IFD entries at given offset."""
        if offset + 2 > len(tiff_data):
            return {}
        count = struct.unpack(endian + 'H', tiff_data[offset:offset + 2])[0]
        entries = {}
        for i in range(count):
            eo = offset + 2 + i * 12
            if eo + 12 > len(tiff_data):
                break
            tag_id = struct.unpack(endian + 'H', tiff_data[eo:eo + 2])[0]
            type_id = struct.unpack(endian + 'H', tiff_data[eo + 2:eo + 4])[0]
            cnt = struct.unpack(endian + 'I', tiff_data[eo + 4:eo + 8])[0]
            val_raw = tiff_data[eo + 8:eo + 12]
            entries[tag_id] = (type_id, cnt, val_raw)
        return entries

    # TIFF type sizes: 1=BYTE 2=ASCII 3=SHORT 4=LONG 5=RATIONAL 7=UNDEF 10=SRATIONAL
    _type_sz = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 7: 1, 10: 8}

    def _read_value(type_id, cnt, val_raw):
        """Decode IFD entry value."""
        sz = _type_sz.get(type_id, 0)
        if sz == 0:
            return None
        total = sz * cnt
        if total <= 4:
            data = val_raw[:total]
        else:
            data_off = struct.unpack(endian + 'I', val_raw)[0]
            if data_off + total > len(tiff_data):
                return None
            data = tiff_data[data_off:data_off + total]

        if type_id == 2:  # ASCII
            return data.decode('ascii', errors='replace').rstrip('\x00')
        elif type_id == 3:  # SHORT
            vals = [struct.unpack(endian + 'H', data[i*2:i*2+2])[0]
                    for i in range(cnt)]
            return vals[0] if cnt == 1 else vals
        elif type_id == 4:  # LONG
            vals = [struct.unpack(endian + 'I', data[i*4:i*4+4])[0]
                    for i in range(cnt)]
            return vals[0] if cnt == 1 else vals
        elif type_id in (5, 10):  # RATIONAL / SRATIONAL
            fmt = endian + ('II' if type_id == 5 else 'ii')
            vals = []
            for i in range(cnt):
                num, den = struct.unpack(fmt, data[i*8:i*8+8])
                vals.append(num / den if den != 0 else 0.0)
            return vals[0] if cnt == 1 else vals
        return None

    tags = {}

    # Map of interesting tags (tag_id -> output key name)
    _IFD0_MAP = {
        _TAG_MAKE: 'Image Make',
        _TAG_MODEL: 'Image Model',
        _TAG_DATETIME: 'Image DateTime',
    }
    _EXIF_MAP = {
        _TAG_EXPOSURE_TIME: 'EXIF ExposureTime',
        _TAG_FNUMBER: 'EXIF FNumber',
        _TAG_ISO: 'EXIF ISOSpeedRatings',
        _TAG_DATETIME_ORIG: 'EXIF DateTimeOriginal',
        _TAG_FOCAL_35MM: 'EXIF FocalLengthIn35mmFilm',
    }
    _MAKERNOTE_MAP = {
        _TAG_CANON_SHOT_INFO: 'Canon ShotInfo',
    }

    def _extract_from_ifd(ifd_entries, tag_map, ifd_name, dump_all=False):
        """Extract known tags from IFD entries."""
        for tag_id, (type_id, cnt, val_raw) in ifd_entries.items():
            if tag_id in tag_map:
                val = _read_value(type_id, cnt, val_raw)
                if val is not None:
                    key = tag_map[tag_id]
                    tags[key] = val
                    _log(f"[tiff_parse]   {ifd_name}: {key} = {val}")
            elif dump_all:
                val = _read_value(type_id, cnt, val_raw)
                # Truncate long arrays for readability
                preview = repr(val)
                if len(preview) > 200:
                    preview = preview[:200] + "..."
                _log(f"[tiff_parse]   {ifd_name}: 0x{tag_id:04X}"
                     f" type={type_id} cnt={cnt} = {preview}")
            # Also follow ExifIFD pointer from any IFD
            if tag_id == _TAG_EXIF_IFD_PTR:
                exif_off = struct.unpack(endian + 'I', val_raw)[0]
                _log(f"[tiff_parse]   {ifd_name}: ExifIFD pointer -> offset {exif_off}")
                exif_ifd = _read_ifd_entries(exif_off)
                _log(f"[tiff_parse]   ExifIFD: {len(exif_ifd)} entries,"
                     f" tags: {[f'0x{t:04X}' for t in sorted(exif_ifd.keys())]}")
                _extract_from_ifd(exif_ifd, _EXIF_MAP, 'ExifIFD')

    def _get_next_ifd(offset, entry_count):
        """Read NextIFD pointer after IFD entries."""
        next_ptr_off = offset + 2 + entry_count * 12
        if next_ptr_off + 4 > len(tiff_data):
            return 0
        return struct.unpack(endian + 'I', tiff_data[next_ptr_off:next_ptr_off + 4])[0]

    # Walk IFD chain: IFD0 -> IFD1 -> IFD2 -> ...
    ifd_offset = ifd0_offset
    ifd_num = 0
    visited = set()
    while ifd_offset != 0 and ifd_offset not in visited and ifd_num < 10:
        visited.add(ifd_offset)
        ifd_entries = _read_ifd_entries(ifd_offset)
        ifd_name = f"IFD{ifd_num}"
        _log(f"[tiff_parse] {ifd_name} at offset {ifd_offset}: {len(ifd_entries)} entries,"
             f" tags: {[f'0x{t:04X}' for t in sorted(ifd_entries.keys())]}")

        # Extract known tags + follow ExifIFD pointer
        # Dump all values if this looks like MakerNote (has low tag IDs)
        is_makernote = any(t < 0x0100 for t in ifd_entries.keys())
        combined_map = {**_IFD0_MAP, **_EXIF_MAP, **_MAKERNOTE_MAP}
        _extract_from_ifd(ifd_entries, combined_map, ifd_name,
                          dump_all=is_makernote)

        # Follow NextIFD
        next_ifd = _get_next_ifd(ifd_offset, len(ifd_entries))
        _log(f"[tiff_parse] {ifd_name} NextIFD -> {next_ifd}")
        ifd_offset = next_ifd
        ifd_num += 1

    _log(f"[tiff_parse] Total parsed tags: {len(tags)}")
    return tags


# ---------------------------------------------------------------
# CR3 EXIF extraction (ISOBMFF container parsing)
# ---------------------------------------------------------------

# Canon CMT UUID boxes inside moov (differ only in last byte):
#   CMT1 (6a48): IFD0 — Make, Model, DateTime, Artist, Copyright
#   CMT2 (6a49): EXIF IFD — ExposureTime, ISOSpeedRatings, DateTimeOriginal
#   CMT3 (6a4a): MakerNote — Canon-specific (CameraTemperature, LensModel, etc.)
#   CMT4 (6a4b): GPS IFD
_CMT_UUID_PREFIX = bytes.fromhex('85c0b687820f11e08111f4ce462b6a')
# Last byte: 0x48=CMT1, 0x49=CMT2, 0x4a=CMT3, 0x4b=CMT4


def _iter_isobmff_boxes(data, start, end):
    """Iterate over ISOBMFF boxes in data[start:end]."""
    pos = start
    while pos < end - 8:
        size = struct.unpack('>I', data[pos:pos + 4])[0]
        btype = data[pos + 4:pos + 8]
        hdr = 8
        if size == 1 and pos + 16 <= end:
            size = struct.unpack('>Q', data[pos + 8:pos + 16])[0]
            hdr = 16
        elif size == 0:
            size = end - pos
        if size < hdr:
            break
        yield btype, pos + hdr, min(pos + size, end)
        pos += size


def _extract_cr3_exif(filepath, log=None):
    """
    Extract EXIF tags from CR3 (ISOBMFF) by finding Canon CMT UUID boxes
    inside moov. Each CMT box contains a self-contained TIFF:
      CMT1 (6a48): IFD0 — Make, Model, DateTime
      CMT2 (6a49): EXIF IFD — ExposureTime, ISO, DateTimeOriginal
      CMT3 (6a4a): MakerNote — Canon-specific (temperature, lens, etc.)
      CMT4 (6a4b): GPS IFD
    Tags from all CMT boxes are merged into one dict.
    Falls back to brute-force TIFF signature search.
    """
    _CMT_NAMES = {0x48: 'CMT1', 0x49: 'CMT2', 0x4a: 'CMT3', 0x4b: 'CMT4'}

    def _log(msg):
        if log:
            log.write(msg + "\n")

    with open(filepath, 'rb') as f:
        head = f.read(512 * 1024)

    _log(f"[cr3_exif] Read {len(head)} bytes from file header")
    _log(f"[cr3_exif] First 16 bytes: {head[:16].hex(' ')}")

    all_tags = {}
    tiff_blocks = []  # list of (name, tiff_bytes)

    # Strategy 1: parse ISOBMFF boxes -> moov -> all CMT uuid boxes
    try:
        for btype, dstart, dend in _iter_isobmff_boxes(head, 0, len(head)):
            _log(f"[cr3_exif] Top box: type={btype!r} offset={dstart} size={dend - dstart}")
            if btype == b'moov':
                _log(f"[cr3_exif] Found moov box, scanning children...")
                for btype2, dstart2, dend2 in _iter_isobmff_boxes(head, dstart, dend):
                    box_size = dend2 - dstart2
                    _log(f"[cr3_exif]   moov child: type={btype2!r} offset={dstart2} size={box_size}")
                    if btype2 == b'uuid' and box_size >= 16:
                        uuid_val = head[dstart2:dstart2 + 16]
                        uuid_prefix = uuid_val[:15]
                        uuid_last = uuid_val[15]
                        cmt_name = _CMT_NAMES.get(uuid_last)
                        _log(f"[cr3_exif]   UUID: {uuid_val.hex()}"
                             f" -> {cmt_name or 'unknown'}")
                        if uuid_prefix == _CMT_UUID_PREFIX and cmt_name:
                            payload = head[dstart2 + 16:dend2]
                            _log(f"[cr3_exif]   {cmt_name} payload: {len(payload)} bytes,"
                                 f" first 16: {payload[:16].hex(' ')}")
                            # Log sub-boxes inside payload (before TIFF)
                            try:
                                for stype, ss, se in _iter_isobmff_boxes(payload, 0, len(payload)):
                                    _log(f"[cr3_exif]   {cmt_name} sub-box:"
                                         f" type={stype!r} size={se - ss}")
                            except Exception:
                                pass
                            # Find ALL TIFF headers inside payload
                            # (Canon wraps TIFF in CNCV/CCTP sub-boxes)
                            for sig in (b'\x49\x49\x2a\x00', b'\x4d\x4d\x00\x2a'):
                                search_start = 0
                                while True:
                                    tiff_pos = payload.find(sig, search_start)
                                    if tiff_pos < 0:
                                        break
                                    # Find end: next TIFF sig or end of payload
                                    next_pos = payload.find(sig, tiff_pos + 4)
                                    if next_pos < 0:
                                        tiff_bytes = payload[tiff_pos:]
                                    else:
                                        tiff_bytes = payload[tiff_pos:next_pos]
                                    _log(f"[cr3_exif]   {cmt_name}: TIFF at offset"
                                         f" {tiff_pos}, {len(tiff_bytes)} bytes")
                                    tiff_blocks.append((f"{cmt_name}@{tiff_pos}", tiff_bytes))
                                    search_start = tiff_pos + 4
                            if not any(cmt_name in name for name, _ in tiff_blocks):
                                _log(f"[cr3_exif]   {cmt_name}: no TIFF signature found")
                break
    except Exception as e:
        _log(f"[cr3_exif] ISOBMFF parse error: {e}")

    # Parse each TIFF block: try own parser for all blocks
    # (exifread often fails to follow sub-IFDs in CR3 embedded TIFFs)
    for block_name, tiff_bytes in tiff_blocks:
        _log(f"[cr3_exif] Parsing {block_name} with own TIFF parser...")
        own_tags = _parse_tiff_tags(tiff_bytes, log=log)
        if own_tags:
            for k, v in own_tags.items():
                if k not in all_tags:
                    all_tags[k] = v

    # Strategy 2: brute-force if no tags found from ISOBMFF
    if not all_tags:
        _log("[cr3_exif] No tags from ISOBMFF, trying brute-force TIFF search...")
        for sig in (b'\x49\x49\x2a\x00', b'\x4d\x4d\x00\x2a'):
            pos = head.find(sig)
            if pos >= 0:
                _log(f"[cr3_exif] Found TIFF signature {sig.hex()} at offset {pos}")
                own_tags = _parse_tiff_tags(head[pos:], log=log)
                if own_tags:
                    all_tags.update(own_tags)
                break
        if not all_tags:
            _log("[cr3_exif] No TIFF/EXIF found anywhere")

    _log(f"[cr3_exif] Total merged tags: {len(all_tags)}")

    return all_tags


# ---------------------------------------------------------------
# EXIF -> FITS header
# ---------------------------------------------------------------

def _tag_to_float(tag):
    """Convert EXIF tag to float. Handles exifread IfdTag and plain values."""
    # Plain float/int from our own parser
    if isinstance(tag, (int, float)):
        return float(tag)
    # exifread IfdTag with .values list of Ratio objects
    if hasattr(tag, 'values'):
        vals = tag.values
        if not vals:
            return None
        r = vals[0]
        if hasattr(r, 'num') and hasattr(r, 'den'):
            if r.den == 0:
                return None
            return float(r.num) / float(r.den)
        return float(r)
    # Fallback: try converting string
    try:
        return float(str(tag))
    except (ValueError, TypeError):
        return None


def _parse_datetime(dt_str):
    """Convert EXIF datetime 'YYYY:MM:DD HH:MM:SS' to FITS 'YYYY-MM-DDTHH:MM:SS'."""
    dt_str = str(dt_str).strip()
    if len(dt_str) >= 10 and dt_str[4] == ':' and dt_str[7] == ':':
        date_part = dt_str[:10].replace(':', '-')
        time_part = dt_str[10:].strip()
        if time_part:
            return f"{date_part}T{time_part}"
        return date_part
    return dt_str


def _datetime_to_jd(dt_str):
    """Convert EXIF datetime string to Julian Date (approximate)."""
    import datetime
    try:
        dt_str = str(dt_str).strip()
        clean = dt_str.replace('T', ' ')
        if clean[4] == '-':
            dt = datetime.datetime.strptime(clean, "%Y-%m-%d %H:%M:%S")
        else:
            dt = datetime.datetime.strptime(clean, "%Y:%m:%d %H:%M:%S")
        epoch = datetime.datetime(1970, 1, 1)
        unix_seconds = (dt - epoch).total_seconds()
        return unix_seconds / 86400.0 + 2440587.5
    except Exception:
        return None


def build_fits_header(exif_tags, raw_pattern, black_levels, debaer):
    """Build FITS header from EXIF tags and RAW metadata."""
    header = fits.Header()

    # INSTRUME = Make + Model
    make = str(exif_tags.get('Image Make', '')).strip()
    model = str(exif_tags.get('Image Model', '')).strip()
    if make and model:
        if model.lower().startswith(make.lower()):
            instrume = model
        else:
            instrume = f"{make} {model}"
    elif model:
        instrume = model
    elif make:
        instrume = make
    else:
        instrume = 'Unknown'
    header['INSTRUME'] = (instrume, 'Camera make and model')

    # EXPTIME / EXPOSURE
    exp_tag = exif_tags.get('EXIF ExposureTime')
    if exp_tag is not None:
        exptime = _tag_to_float(exp_tag)
        if exptime is not None:
            header['EXPTIME'] = (exptime, '[s] Exposure time')
            header['EXPOSURE'] = (exptime, '[s] Exposure time')

    # GAIN (ISO)
    iso_tag = exif_tags.get('EXIF ISOSpeedRatings')
    if iso_tag:
        try:
            iso_val = int(str(iso_tag))
            header['GAIN'] = (iso_val, 'ISO speed')
        except (ValueError, TypeError):
            pass

    # DATE-OBS and JD (prefer DateTimeOriginal, fallback to DateTime)
    dt_tag = exif_tags.get('EXIF DateTimeOriginal') or exif_tags.get('Image DateTime')
    if dt_tag:
        dt_str = str(dt_tag)
        date_obs = _parse_datetime(dt_str)
        header['DATE-OBS'] = (date_obs, 'Date/time of observation')
        jd = _datetime_to_jd(dt_str)
        if jd is not None:
            header['JD'] = (jd, 'Julian Date of observation')

    # OBJECT (from ImageDescription if present)
    desc_tag = exif_tags.get('Image ImageDescription')
    if desc_tag:
        desc = str(desc_tag).strip()
        if desc:
            header['OBJECT'] = (desc, 'Object description')

    # Fixed metadata
    header['XPIXSZ'] = (5.36, '[um] Pixel size X (Canon APS-C default)')
    header['YPIXSZ'] = (5.36, '[um] Pixel size Y (Canon APS-C default)')
    header['FOCALLEN'] = (1000, '[mm] Focal length (default)')
    # CCD-TEMP from Canon ShotInfo[12] - 128 (camera body temperature)
    shot_info = exif_tags.get('Canon ShotInfo')
    if isinstance(shot_info, list) and len(shot_info) > 12:
        header['CCD-TEMP'] = (float(shot_info[12] - 128), '[C] Camera body temperature')
    else:
        header['CCD-TEMP'] = (0.0, '[C] Sensor temperature (unknown)')
    header['SWCREATE'] = ('raw2fits (AstroBatch)', 'Software that created this file')
    header['XBINNING'] = (1, 'Binning factor X')
    header['YBINNING'] = (1, 'Binning factor Y')

    if debaer:
        # Demosaiced RGB mode
        header['FILTER'] = ('RGB', 'Demosaiced linear RGB')
    else:
        # Raw CFA mode
        header['FILTER'] = ('CFA', 'Color filter array (raw Bayer)')
        # BAYERPAT from raw_pattern
        try:
            pattern_2x2 = raw_pattern[:2, :2]
            bayerpat = ''.join(_COLOR_CHAR[c] for c in pattern_2x2.flat)
            header['BAYERPAT'] = (bayerpat, 'Bayer color pattern')
        except (KeyError, IndexError):
            header['BAYERPAT'] = ('RGGB', 'Bayer color pattern (assumed)')

    # Black levels as info
    header['BLKLVL_R'] = (int(black_levels[0]), 'Black level Red channel')
    header['BLKLVL_G'] = (int(black_levels[1]), 'Black level Green1 channel')
    header['BLKLVL_B'] = (int(black_levels[2]), 'Black level Blue channel')
    header['BLKLVLG2'] = (int(black_levels[3]), 'Black level Green2 channel')

    return header


# ---------------------------------------------------------------
# Dark/Light auto-detection
# ---------------------------------------------------------------

def detect_image_type(bayer, black_levels):
    """Detect whether image is Dark Frame or Light Frame."""
    min_black = min(int(bl) for bl in black_levels)
    usable_range = int(bayer.max()) - min_black
    if usable_range <= 0:
        return "Dark Frame"

    threshold = 0.30 * usable_range + min_black
    bright_fraction = np.count_nonzero(bayer > threshold) / bayer.size
    return "Dark Frame" if bright_fraction < 0.0002 else "Light Frame"


# ---------------------------------------------------------------
# Conversion
# ---------------------------------------------------------------

def convert_file(infile, outfile, debaer, debug=False):
    """Convert a single RAW file to FITS."""
    log = None
    if debug:
        log_path = os.path.splitext(infile)[0] + ".log"
        log = open(log_path, 'w', encoding='utf-8')
        log.write(f"=== raw2fits debug log ===\n")
        log.write(f"Input:  {infile}\n")
        log.write(f"Output: {outfile}\n")
        log.write(f"Debaer: {debaer}\n\n")

    try:
        _convert_file_inner(infile, outfile, debaer, log)
    finally:
        if log:
            log.write("\n=== done ===\n")
            log.close()


def _convert_file_inner(infile, outfile, debaer, log):
    """Inner conversion with optional debug log."""
    def _log(msg):
        if log:
            log.write(msg + "\n")

    # --- rawpy ---
    _log("[rawpy] Opening file...")
    raw = rawpy.imread(infile)

    raw_pattern = raw.raw_pattern.copy()
    black_levels = list(raw.black_level_per_channel)

    _log(f"[rawpy] raw_pattern = {raw_pattern.tolist()}")
    _log(f"[rawpy] black_level_per_channel = {black_levels}")
    _log(f"[rawpy] color_desc = {raw.color_desc}")
    _log(f"[rawpy] num_colors = {raw.num_colors}")
    _log(f"[rawpy] raw_image shape = {raw.raw_image.shape}, dtype = {raw.raw_image.dtype}")
    _log(f"[rawpy] raw_image_visible shape = {raw.raw_image_visible.shape}")
    _log(f"[rawpy] camera_whitebalance = {raw.camera_whitebalance}")

    if debaer:
        _log("[rawpy] Postprocessing (linear demosaic)...")
        rgb = raw.postprocess(
            demosaic_algorithm=rawpy.DemosaicAlgorithm.AHD,
            output_bps=16,
            no_auto_bright=True,
            gamma=(1, 1),
            use_camera_wb=True,
            output_color=rawpy.ColorSpace.raw,
        )
        raw.close()
        _log(f"[rawpy] postprocess result: shape={rgb.shape}, dtype={rgb.dtype}")
        data = np.transpose(rgb, (2, 0, 1)).astype(np.uint16)
        detect_data = data[1]
    else:
        bayer = raw.raw_image_visible.copy()
        raw.close()
        _log(f"[rawpy] raw_image_visible: shape={bayer.shape}, dtype={bayer.dtype}")
        _log(f"[rawpy] min={bayer.min()}, max={bayer.max()}, median={np.median(bayer):.0f}")
        data = bayer.astype(np.uint16)
        detect_data = data

    _log(f"[data] Output array: shape={data.shape}, dtype={data.dtype}")

    # --- EXIF ---
    _log("\n[exif] Trying exifread on raw file...")
    stderr_capture = io.StringIO()
    with open(infile, 'rb') as f:
        with contextlib.redirect_stderr(stderr_capture):
            exif_tags = exifread.process_file(f, details=False)

    stderr_msg = stderr_capture.getvalue().strip()
    if stderr_msg:
        _log(f"[exif] exifread stderr: {stderr_msg}")
    _log(f"[exif] exifread returned {len(exif_tags)} tags")
    if exif_tags:
        for k, v in sorted(exif_tags.items()):
            _log(f"[exif]   {k} = {v}")

    # CR3 fallback
    if not exif_tags:
        _log("\n[exif] No tags from exifread, trying CR3 ISOBMFF extraction...")
        exif_tags = _extract_cr3_exif(infile, log=log)
        _log(f"[exif] CR3 extraction returned {len(exif_tags)} tags")

    # --- Header ---
    header = build_fits_header(exif_tags, raw_pattern, black_levels, debaer)

    image_type = detect_image_type(detect_data, black_levels)
    header['IMAGETYP'] = (image_type, 'Frame type (auto-detected)')
    _log(f"\n[detect] IMAGETYP = {image_type}")

    mode_str = "demosaiced RGB" if debaer else "raw CFA Bayer"
    header['HISTORY'] = f'Converted from {os.path.basename(infile)} ({mode_str}) by raw2fits'

    _log(f"\n[header] Final FITS header:")
    for key in header.keys():
        if key and key != 'HISTORY':
            _log(f"[header]   {key:10s} = {header[key]!r}")

    # --- Write ---
    hdu = fits.PrimaryHDU(data, header=header)
    hdu.writeto(outfile, overwrite=True)
    _log(f"\n[write] Written to {outfile}")


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

def main():
    input_spec, output_spec, debaer, debug = parse_args(sys.argv)

    try:
        io_pairs = build_io_pairs(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            convert_file(infile, outfile, debaer, debug)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
