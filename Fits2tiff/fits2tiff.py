#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import numpy as np
from astropy.io import fits

try:
    from PIL import Image
except ImportError:
    sys.stderr.write(
        "Error: Pillow is required but not installed.\n"
        "Install it with: pip install Pillow\n"
    )
    sys.exit(1)

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  fits2tiff.py [--bits 8|16|32] [--flip] input_spec output_spec\n"
        "\n"
        "    input_spec   - single file OR numbered pattern (e.g. light0001.fit)\n"
        "                   OR wildcard mask (e.g. *.fit) OR @list.txt\n"
        "    output_spec  - single file (e.g. image.tif) OR numbered pattern\n"
        "                   (e.g. out0001.tif) OR wildcard (e.g. *.tif)\n"
        "                   Wildcard: * is replaced with input filename stem\n"
        "\n"
        "    --bits 8     - linear stretch [min,max] -> [0,255], uint8 output\n"
        "    --bits 16    - clamp to [0,65535], uint16 output\n"
        "    --bits 32    - float32 output, no scaling\n"
        "    (default)    - auto: uint8/int8->8, uint16/int16->16, else->32\n"
        "\n"
        "    --flip       - flip image vertically (FITS bottom-left to TIFF top-left)\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]
    bits = None
    flip = False

    # Extract --flip flag
    if "--flip" in args:
        flip = True
        args.remove("--flip")

    # Extract --bits option
    if "--bits" in args:
        idx = args.index("--bits")
        if idx + 1 >= len(args):
            sys.stderr.write("Error: --bits requires a value (8, 16, or 32).\n")
            sys.exit(1)
        bits_str = args[idx + 1]
        if bits_str not in ("8", "16", "32"):
            sys.stderr.write("Error: --bits must be 8, 16, or 32.\n")
            sys.exit(1)
        bits = int(bits_str)
        args = args[:idx] + args[idx + 2:]

    if len(args) != 2:
        usage()

    input_spec = args[0]
    output_spec = args[1]

    return input_spec, output_spec, bits, flip


# Numbered pattern for TIFF output: prefix + digits + .tif(f)
_TIFF_SEQ_RE = re.compile(r"^(.*?)(\d+)(\.tiff?)$", re.IGNORECASE)


def _apply_wildcard_output(inputs, output_spec):
    """Replace * in output_spec with each input file's stem (name without ext)."""
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


def build_tiff_io_pairs(input_spec, output_spec):
    """Build (input, output) pairs supporting .tif/.tiff output patterns."""
    inputs = batch_utils.expand_input_spec(input_spec)
    if not inputs:
        raise ValueError("No input files found.")

    # Wildcard output: *.tif, prefix_*.tiff, etc.
    if "*" in output_spec:
        return _apply_wildcard_output(inputs, output_spec)

    # Single input -> single output
    if len(inputs) == 1:
        return [(inputs[0], output_spec)]

    # Multiple inputs => output must be numbered pattern
    base = os.path.basename(output_spec)
    m = _TIFF_SEQ_RE.match(base)
    if not m:
        raise ValueError(
            "Output pattern must contain a numeric field when multiple "
            "input files are provided (e.g. out0001.tif), or use "
            "wildcard (e.g. *.tif) to preserve input names."
        )

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


def auto_bits(dtype):
    """Choose output bit depth based on FITS data type."""
    if dtype in (np.uint8, np.int8):
        return 8
    if dtype in (np.uint16, np.int16):
        return 16
    return 32


def _scale_channel(work, out_bits):
    """Scale a single 2D channel to the target bit depth. Returns ndarray."""
    if out_bits == 8:
        vmin = work.min()
        vmax = work.max()
        if vmax > vmin:
            scaled = (work - vmin) / (vmax - vmin) * 255.0
        else:
            scaled = np.zeros_like(work)
        return np.rint(np.clip(scaled, 0, 255)).astype(np.uint8)

    elif out_bits == 16:
        return np.rint(np.clip(work, 0, 65535)).astype(np.uint16)

    else:  # 32-bit float, normalize to [0, 1] for Photoshop compatibility
        vmin = work.min()
        vmax = work.max()
        if vmax > vmin:
            work = (work - vmin) / (vmax - vmin)
        else:
            work = np.zeros_like(work)
        return work.astype(np.float32)


def convert_file(infile, outfile, bits, flip):
    """Convert a single FITS file to TIFF. Supports 2D (mono) and 3D RGB (3×H×W)."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"File '{infile}' has no primary image data.")
        data = hdul[0].data
        orig_dtype = data.dtype

    is_rgb = data.ndim == 3 and data.shape[0] == 3
    if data.ndim != 2 and not is_rgb:
        raise ValueError(
            f"File '{infile}': unsupported shape {data.shape}. "
            f"Expected 2D (H×W) or 3D RGB (3×H×W)."
        )

    # Determine output bits
    out_bits = bits if bits is not None else auto_bits(orig_dtype)

    if is_rgb:
        # 3×H×W -> process each channel, then stack to H×W×3
        channels = []
        for ch in range(3):
            work = data[ch].astype(np.float64)
            work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)
            if flip:
                work = np.flipud(work)
            channels.append(work)

        if out_bits == 8:
            # Use global min/max across all channels for consistent stretch
            vmin = min(ch.min() for ch in channels)
            vmax = max(ch.max() for ch in channels)
            result_channels = []
            for ch in channels:
                if vmax > vmin:
                    scaled = (ch - vmin) / (vmax - vmin) * 255.0
                else:
                    scaled = np.zeros_like(ch)
                result_channels.append(np.rint(np.clip(scaled, 0, 255)).astype(np.uint8))
            result = np.stack(result_channels, axis=-1)  # H×W×3
        elif out_bits == 16:
            result = np.stack(
                [np.rint(np.clip(ch, 0, 65535)).astype(np.uint16) for ch in channels],
                axis=-1
            )
        else:
            vmin = min(ch.min() for ch in channels)
            vmax = max(ch.max() for ch in channels)
            result_channels = []
            for ch in channels:
                if vmax > vmin:
                    ch = (ch - vmin) / (vmax - vmin)
                else:
                    ch = np.zeros_like(ch)
                result_channels.append(ch.astype(np.float32))
            result = np.stack(result_channels, axis=-1)

        if out_bits == 16:
            _save_tiff_rgb16(result, outfile)
            return
        img = Image.fromarray(result)
    else:
        # 2D mono
        work = data.astype(np.float64)
        work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)
        if flip:
            work = np.flipud(work)
        result = _scale_channel(work, out_bits)
        img = Image.fromarray(result)

    img.save(outfile)


def _save_tiff_rgb16(result, outfile):
    """Write 16-bit RGB TIFF manually (Pillow doesn't support RGB uint16)."""
    import struct

    h, w = result.shape[:2]
    pixels = np.ascontiguousarray(result).tobytes()
    pixel_bytes = len(pixels)

    # IFD entries: tag(H), type(H), count(I), value/offset(I)
    # Types: 3=SHORT, 4=LONG
    tags = [
        (256, 4, 1, w),                # ImageWidth
        (257, 4, 1, h),                # ImageLength
        (258, 3, 3, 0),                # BitsPerSample (3 values -> offset, filled below)
        (259, 3, 1, 1),                # Compression = None
        (262, 3, 1, 2),                # PhotometricInterpretation = RGB
        (273, 4, 1, 0),                # StripOffsets (filled below)
        (277, 3, 1, 3),                # SamplesPerPixel = 3
        (278, 4, 1, h),                # RowsPerStrip = all rows
        (279, 4, 1, pixel_bytes),      # StripByteCounts
        (284, 3, 1, 1),                # PlanarConfiguration = chunky
    ]
    n_tags = len(tags)

    # Layout:
    #   0..7   : TIFF header (8 bytes)
    #   8      : IFD: count(2) + entries(n*12) + next_ifd(4)
    #   after IFD: BitsPerSample values (3 x uint16 = 6 bytes)
    #   then: pixel data
    ifd_offset = 8
    ifd_size = 2 + n_tags * 12 + 4
    bps_offset = ifd_offset + ifd_size         # BitsPerSample array
    data_offset = bps_offset + 6               # pixel data

    # Fix offsets in tags
    fixed_tags = []
    for tag, typ, cnt, val in tags:
        if tag == 258:    # BitsPerSample -> offset to 3 shorts
            val = bps_offset
        elif tag == 273:  # StripOffsets -> pixel data start
            val = data_offset
        fixed_tags.append((tag, typ, cnt, val))

    buf = bytearray()
    # TIFF header: little-endian, magic 42, IFD offset
    buf += struct.pack('<2sHI', b'II', 42, ifd_offset)
    # IFD entry count
    buf += struct.pack('<H', n_tags)
    for tag, typ, cnt, val in fixed_tags:
        buf += struct.pack('<HHII', tag, typ, cnt, val)
    # Next IFD = 0
    buf += struct.pack('<I', 0)
    # BitsPerSample values
    buf += struct.pack('<HHH', 16, 16, 16)
    # Pixel data
    buf += pixels

    with open(outfile, 'wb') as f:
        f.write(buf)


def main():
    input_spec, output_spec, bits, flip = parse_args(sys.argv)

    try:
        io_pairs = build_tiff_io_pairs(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)

    if total == 1:
        # Single file — no threading overhead
        infile, outfile = io_pairs[0]
        try:
            convert_file(infile, outfile, bits, flip)
            sys.stderr.write(f"\rProcessed 1 / 1 files")
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")
    else:
        from concurrent.futures import ThreadPoolExecutor, as_completed
        max_workers = max(1, os.cpu_count() - 1)
        done = 0
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for infile, outfile in io_pairs:
                fut = executor.submit(convert_file, infile, outfile, bits, flip)
                futures[fut] = infile
            for fut in as_completed(futures):
                infile = futures[fut]
                done += 1
                try:
                    fut.result()
                    sys.stderr.write(f"\rProcessed {done} / {total} files")
                    sys.stderr.flush()
                except Exception as e:
                    sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
