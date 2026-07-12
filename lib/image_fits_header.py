#!/usr/bin/env python3
# image_fits_header.py — shared codec for embedding a FITS header inside 8-bit
# image containers as JSON. Single source of truth so every tool round-trips
# the same bytes.
#
# Format: a JSON object stored in a JPEG COM marker (0xFFFE) placed right after
# SOI (FF D8). COMMENT cards -> "_COMMENT": [...]; HISTORY cards -> "_HISTORY".
# (PNG stores the same JSON string in a "fits_header" text chunk.)
#
# Producers: fits2tiff.py (FITS -> JPEG/PNG). Consumers: autosolve.py (JPEG
# solve round-trip), and anything reading a header back out of a JPEG.
# !!! Code comments in ENGLISH ONLY please !!!

import os
import json
import struct

from astropy.io import fits

# FITS structural keywords: derived by astropy from the array itself. They must
# never be reconstructed from a serialized header — they would clash with the
# real image dtype/shape.
_STRUCTURAL_KEYS = frozenset([
    "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3",
    "BZERO", "BSCALE", "EXTEND", "PCOUNT", "GCOUNT",
])


def fits_header_to_json(header):
    """Serialize a FITS header to a JSON string (JSON-safe scalars only)."""
    result = {}
    comments = []
    history = []

    for card in header.cards:
        key = card.keyword
        val = card.value
        if key == 'COMMENT':
            comments.append(str(val))
        elif key == 'HISTORY':
            history.append(str(val))
        elif key == '' or key is None:
            continue
        else:
            # Convert to JSON-safe type
            if isinstance(val, (int, float, bool, str)):
                result[key] = val
            else:
                result[key] = str(val)

    if comments:
        result['_COMMENT'] = comments
    if history:
        result['_HISTORY'] = history

    return json.dumps(result, ensure_ascii=False)


def json_to_fits_header(data, strip_structural=True):
    """Rebuild a FITS header from a dict (or JSON string) produced by
    fits_header_to_json(). Inverse of fits_header_to_json().

    strip_structural: drop SIMPLE/BITPIX/NAXIS*/BZERO/BSCALE/... so the header
    can be safely attached to an array of different dtype/shape. Cards astropy
    rejects (invalid keyword/value) are skipped rather than raising.
    """
    if isinstance(data, (bytes, bytearray)):
        data = data.decode('utf-8')
    if isinstance(data, str):
        data = json.loads(data)

    hdr = fits.Header()
    if not isinstance(data, dict):
        return hdr

    for key, val in data.items():
        if key in ('_COMMENT', '_HISTORY'):
            continue
        if strip_structural and key in _STRUCTURAL_KEYS:
            continue
        try:
            hdr[key] = val
        except Exception:
            # Invalid FITS keyword or unrepresentable value -> skip.
            continue

    for c in (data.get('_COMMENT') or []):
        try:
            hdr['COMMENT'] = str(c)
        except Exception:
            pass
    for h in (data.get('_HISTORY') or []):
        try:
            hdr['HISTORY'] = str(h)
        except Exception:
            pass

    return hdr


def embed_header_in_jpeg(jpeg_path, header_json_str):
    """Inject a FITS-header JSON string as a JPEG COM marker after SOI."""
    header_bytes = header_json_str.encode('utf-8')

    # A single COM marker is limited to 65533 payload bytes (length field is
    # 16-bit and counts itself). Truncate if it somehow exceeds that.
    if len(header_bytes) > 65533:
        header_bytes = header_bytes[:65533]

    com_len = len(header_bytes) + 2  # +2 for the length field itself
    com_marker = b'\xff\xfe' + struct.pack('>H', com_len) + header_bytes

    with open(jpeg_path, 'rb') as f:
        jpeg_data = f.read()

    # Insert the COM marker immediately after SOI (FF D8).
    with open(jpeg_path, 'wb') as f:
        f.write(jpeg_data[:2])   # SOI
        f.write(com_marker)      # COM with header
        f.write(jpeg_data[2:])   # rest of JPEG


def read_fits_header_from_jpeg(jpeg_path):
    """Extract the FITS header dict from a JPEG COM marker. Returns dict or None."""
    with open(jpeg_path, 'rb') as f:
        data = f.read(min(os.path.getsize(jpeg_path), 131072))  # first 128 KB

    pos = 2  # skip SOI (FF D8)
    while pos < len(data) - 3:
        if data[pos] != 0xFF:
            break
        marker = data[pos + 1]
        if marker == 0xD9:  # EOI
            break
        if marker == 0xFE:  # COM
            length = struct.unpack('>H', data[pos + 2:pos + 4])[0]
            comment = data[pos + 4:pos + 2 + length]
            try:
                return json.loads(comment.decode('utf-8'))
            except (json.JSONDecodeError, UnicodeDecodeError):
                pass
        # Skip to next marker
        if marker in (0xD8, 0x00):
            pos += 2
        else:
            length = struct.unpack('>H', data[pos + 2:pos + 4])[0]
            pos += 2 + length

    return None
