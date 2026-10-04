#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ser_writer - minimal writer for the SER astro/planetary video format (Lucam
Recorder v3): a 178-byte header, raw frames concatenated, and an optional
per-frame UTC timestamp trailer.

Pure Python (struct + the array's own bytes) - no external dependencies. Frames
are written one at a time as they are produced, so any number of frames works
within the memory of a single frame; the output size is limited only by disk.

Format notes:
  * The header is little-endian. The 'LittleEndian' field follows the de-facto
    convention used by SER software (0 = little-endian pixel data, 1 = big-endian)
    - the opposite of the original spec. We write 16-bit data little-endian and set
    the field to 0.
  * Mono: ColorID = 0, one plane. RGB: ColorID = 100, three interleaved planes
    (R,G,B per pixel).
  * 8-bit frames are uint8; 16-bit frames are little-endian uint16.
  * The optional trailer holds FrameCount 64-bit UTC timestamps as .NET "ticks"
    (100-ns intervals since 0001-01-01 00:00:00 UTC). Written only if at least one
    timestamp was supplied (missing ones are written as 0).

Reference: SER format v3 (Grischa Hahn), as implemented by PIPP / SER Player / Siril.
"""

import struct
from datetime import datetime, timezone

import numpy as np

COLOR_MONO = 0
COLOR_RGB = 100
COLOR_BGR = 101

_HEADER_SIZE = 178
_FILE_ID = b"LUCAM-RECORDER"          # 14 bytes
_FRAMECOUNT_OFFSET = 38
_DATETIME_OFFSET = 162
_TICKS_PER_SECOND = 10_000_000        # .NET ticks: 100-ns intervals
_NET_EPOCH = datetime(1, 1, 1, tzinfo=timezone.utc)


def datetime_to_ticks(dt):
    """Convert a datetime (naive is treated as UTC) to .NET UTC ticks (int64)."""
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    delta = dt.astimezone(timezone.utc) - _NET_EPOCH
    return int(round(delta.total_seconds() * _TICKS_PER_SECOND))


def _pad_ascii(s, n):
    """Encode a string to exactly n bytes of ASCII (truncated / zero-padded)."""
    b = ("" if s is None else str(s)).encode("ascii", "replace")[:n]
    return b + b"\x00" * (n - len(b))


class SerWriter:
    """Streaming writer for a single SER file.

    with SerWriter(path, width, height, color=False, depth=16) as w:
        for frame in frames:
            w.add_frame(frame, timestamp=dt)     # dt: UTC datetime, or None

    add_frame() takes (H, W) mono or (H, W, 3) RGB arrays, uint8 for depth=8 or
    uint16 for depth=16; the shape must match width/height/color. FrameCount and
    the timestamp trailer are finalized in close()."""

    def __init__(self, path, width, height, color=False, depth=16,
                 observer="", instrument="", telescope="PULSAR animate"):
        if depth not in (8, 16):
            raise ValueError("SER depth must be 8 or 16, got %r" % depth)
        self.path = path
        self.width = int(width)
        self.height = int(height)
        self.color = bool(color)
        self.depth = int(depth)
        self._dtype = np.uint8 if depth == 8 else np.dtype("<u2")
        self._observer = observer
        self._instrument = instrument
        self._telescope = telescope
        self._count = 0
        self._timestamps = []
        self._any_ts = False
        self._f = open(path, "wb")
        self._f.write(self._header(0))            # placeholder FrameCount/DateTime

    def _header(self, frame_count, first_ts=None):
        color_id = COLOR_RGB if self.color else COLOR_MONO
        ticks = datetime_to_ticks(first_ts) if first_ts is not None else 0
        h = bytearray()
        h += _FILE_ID                             # FileID (14)
        h += struct.pack("<i", 0)                 # LuID
        h += struct.pack("<i", color_id)          # ColorID
        h += struct.pack("<i", 0)                 # LittleEndian (0 = LE data, de-facto)
        h += struct.pack("<i", self.width)        # ImageWidth
        h += struct.pack("<i", self.height)       # ImageHeight
        h += struct.pack("<i", self.depth)        # PixelDepthPerPlane
        h += struct.pack("<i", frame_count)       # FrameCount
        h += _pad_ascii(self._observer, 40)
        h += _pad_ascii(self._instrument, 40)
        h += _pad_ascii(self._telescope, 40)
        h += struct.pack("<q", ticks)             # DateTime (local; we store UTC ticks)
        h += struct.pack("<q", ticks)             # DateTimeUTC
        if len(h) != _HEADER_SIZE:
            raise RuntimeError("SER header is %d bytes, expected %d" % (len(h), _HEADER_SIZE))
        return bytes(h)

    def add_frame(self, arr, timestamp=None):
        """Append one frame. arr: (H,W) mono or (H,W,3) RGB, matching the writer's
        size/colour/depth. timestamp: a UTC datetime for the trailer, or None."""
        if self._f is None:
            raise RuntimeError("SerWriter is closed")
        a = np.asarray(arr)
        expected = (self.height, self.width, 3) if self.color else (self.height, self.width)
        if a.shape != expected:
            raise ValueError("frame shape %s != expected %s" % (a.shape, expected))
        a = np.ascontiguousarray(a.astype(self._dtype, copy=False))   # LE uint16 / uint8
        self._f.write(a.tobytes())
        self._timestamps.append(timestamp)
        if timestamp is not None:
            self._any_ts = True
        self._count += 1

    def close(self):
        """Finalize: append the timestamp trailer (if any), patch FrameCount and the
        header DateTime, and close the file."""
        if self._f is None:
            return
        if self._any_ts:
            self._f.write(b"".join(
                struct.pack("<q", datetime_to_ticks(t) if t is not None else 0)
                for t in self._timestamps))
        self._f.seek(_FRAMECOUNT_OFFSET)
        self._f.write(struct.pack("<i", self._count))
        if self._any_ts:
            first = next((t for t in self._timestamps if t is not None), None)
            if first is not None:
                ticks = struct.pack("<q", datetime_to_ticks(first))
                self._f.seek(_DATETIME_OFFSET)
                self._f.write(ticks + ticks)      # DateTime + DateTimeUTC
        self._f.close()
        self._f = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
