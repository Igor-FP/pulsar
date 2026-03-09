#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
debayer - Demosaic Bayer-pattern FITS images to RGB.

Input:  2D monochrome FITS with Bayer CFA pattern (H x W).
Output: 3-channel RGB FITS (3 x H x W).

Algorithms:
  bilinear  - Pure numpy bilinear interpolation (default, fast)
  vng       - Variable Number of Gradients via OpenCV (requires cv2)

Bayer pattern is read from BAYERPAT header keyword, or set via --pattern.
"""

import sys
import os
import numpy as np
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import shared utilities (two levels up from Personal/Debayer/)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "debayer - Demosaic Bayer-pattern FITS to RGB.\n"
        "\n"
        "Usage:\n"
        "  debayer.py input_spec output_spec [options]\n"
        "\n"
        "Positional arguments:\n"
        "  input_spec       Input 2D Bayer FITS: single file, wildcard (*.fit),\n"
        "                   numbered sequence (img0001.fit), or list (@list.txt)\n"
        "  output_spec      Output: single file, numbered pattern (out0001.fit),\n"
        "                   or directory\n"
        "\n"
        "Options:\n"
        "  --pattern PAT    Bayer pattern: RGGB, BGGR, GRBG, GBRG\n"
        "                   (default: from BAYERPAT header, or RGGB)\n"
        "  --method METHOD  Demosaicing algorithm (default: bilinear)\n"
        "                     bilinear  - fast, pure numpy\n"
        "                     vng       - high quality, requires OpenCV (cv2)\n"
        "\n"
        "Examples:\n"
        "  debayer.py img0001.fit out0001.fit\n"
        "  debayer.py *.fit rgb/ --method vng\n"
        "  debayer.py img.fit out.fit --pattern BGGR\n"
    )
    sys.exit(1)


VALID_PATTERNS = {'RGGB', 'BGGR', 'GRBG', 'GBRG'}


def parse_args(argv):
    args = argv[1:]

    pattern = None  # auto from header
    method = 'bilinear'

    i = 0
    positional = []
    while i < len(args):
        if args[i] == "--pattern":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --pattern requires a value.\n")
                sys.exit(1)
            pattern = args[i+1].upper()
            if pattern not in VALID_PATTERNS:
                sys.stderr.write(
                    f"Error: invalid pattern '{pattern}'. "
                    f"Must be one of: {', '.join(sorted(VALID_PATTERNS))}\n")
                sys.exit(1)
            i += 2
        elif args[i] == "--method":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --method requires a value.\n")
                sys.exit(1)
            method = args[i+1].lower()
            if method not in ('bilinear', 'vng'):
                sys.stderr.write(
                    "Error: --method must be 'bilinear' or 'vng'.\n")
                sys.exit(1)
            i += 2
        else:
            positional.append(args[i])
            i += 1

    if len(positional) < 2:
        usage()

    return positional[0], positional[1], pattern, method


# ---------------------------------------------------------------------------
# Bayer pattern helpers
# ---------------------------------------------------------------------------

def get_pattern(header, cli_pattern):
    """Get Bayer pattern from CLI override or FITS header."""
    if cli_pattern:
        return cli_pattern
    pat = header.get('BAYERPAT', '').strip().upper()
    if pat in VALID_PATTERNS:
        return pat
    return 'RGGB'


def parse_pattern(pattern):
    """Parse pattern string to color positions.

    Returns dict: color -> list of (row_offset, col_offset) within 2x2 block.
    Pattern string maps: [0]=(0,0), [1]=(0,1), [2]=(1,0), [3]=(1,1).
    """
    positions = {}
    for i, c in enumerate(pattern):
        positions.setdefault(c, []).append((i // 2, i % 2))
    return positions


# ---------------------------------------------------------------------------
# Bilinear demosaicing (pure numpy)
# ---------------------------------------------------------------------------

def debayer_bilinear(raw_2d, pattern):
    """Bilinear interpolation demosaicing.

    For each pixel site, missing color channels are interpolated from
    nearest neighbors of that color:
      - Same-color diagonal: average of 4 diagonal neighbors
      - Same-color cardinal: average of 4 cardinal neighbors
      - Same-color horizontal: average of left + right
      - Same-color vertical: average of top + bottom
    """
    h, w = raw_2d.shape
    h -= h % 2
    w -= w % 2
    img = raw_2d[:h, :w].astype(np.float64)

    pos = parse_pattern(pattern)
    r_pos = pos['R'][0]
    b_pos = pos['B'][0]
    g_pos = pos['G']  # two positions

    # Padded image for safe neighbor access
    p = np.pad(img, 1, mode='reflect')

    # Precompute neighbor averages (all shape h x w)
    cardinal = (p[0:h, 1:w+1] + p[2:h+2, 1:w+1] +
                p[1:h+1, 0:w] + p[1:h+1, 2:w+2]) * 0.25
    diagonal = (p[0:h, 0:w] + p[0:h, 2:w+2] +
                p[2:h+2, 0:w] + p[2:h+2, 2:w+2]) * 0.25
    horiz = (p[1:h+1, 0:w] + p[1:h+1, 2:w+2]) * 0.5
    vert = (p[0:h, 1:w+1] + p[2:h+2, 1:w+1]) * 0.5

    rgb = np.zeros((3, h, w), dtype=np.float64)

    # --- Red channel (index 0) ---
    # R at R sites: original value
    rgb[0, r_pos[0]::2, r_pos[1]::2] = img[r_pos[0]::2, r_pos[1]::2]
    # R at B sites: diagonal average (4 surrounding R pixels)
    rgb[0, b_pos[0]::2, b_pos[1]::2] = diagonal[b_pos[0]::2, b_pos[1]::2]
    # R at G sites: horiz or vert depending on G position relative to R
    for gr, gc in g_pos:
        if gr == r_pos[0]:
            # G in same row as R -> R neighbors are left and right
            rgb[0, gr::2, gc::2] = horiz[gr::2, gc::2]
        else:
            # G in same column as R -> R neighbors are top and bottom
            rgb[0, gr::2, gc::2] = vert[gr::2, gc::2]

    # --- Blue channel (index 2) ---
    # B at B sites: original value
    rgb[2, b_pos[0]::2, b_pos[1]::2] = img[b_pos[0]::2, b_pos[1]::2]
    # B at R sites: diagonal average (4 surrounding B pixels)
    rgb[2, r_pos[0]::2, r_pos[1]::2] = diagonal[r_pos[0]::2, r_pos[1]::2]
    # B at G sites
    for gr, gc in g_pos:
        if gr == b_pos[0]:
            # G in same row as B -> B neighbors are left and right
            rgb[2, gr::2, gc::2] = horiz[gr::2, gc::2]
        else:
            # G in same column as B -> B neighbors are top and bottom
            rgb[2, gr::2, gc::2] = vert[gr::2, gc::2]

    # --- Green channel (index 1) ---
    # G at G sites: original value
    for gr, gc in g_pos:
        rgb[1, gr::2, gc::2] = img[gr::2, gc::2]
    # G at R sites: cardinal average (4 surrounding G pixels)
    rgb[1, r_pos[0]::2, r_pos[1]::2] = cardinal[r_pos[0]::2, r_pos[1]::2]
    # G at B sites: cardinal average (4 surrounding G pixels)
    rgb[1, b_pos[0]::2, b_pos[1]::2] = cardinal[b_pos[0]::2, b_pos[1]::2]

    return rgb


# ---------------------------------------------------------------------------
# OpenCV VNG demosaicing
# ---------------------------------------------------------------------------

# OpenCV Bayer code mapping: pattern -> cv2 constant name
_CV2_BAYER_MAP = {
    'RGGB': 'COLOR_BayerRG2RGB',
    'BGGR': 'COLOR_BayerBG2RGB',
    'GRBG': 'COLOR_BayerGR2RGB',
    'GBRG': 'COLOR_BayerGB2RGB',
}

_CV2_BAYER_VNG_MAP = {
    'RGGB': 'COLOR_BayerRG2RGB_VNG',
    'BGGR': 'COLOR_BayerBG2RGB_VNG',
    'GRBG': 'COLOR_BayerGR2RGB_VNG',
    'GBRG': 'COLOR_BayerGB2RGB_VNG',
}


def debayer_vng(raw_2d, pattern):
    """VNG demosaicing via OpenCV."""
    try:
        import cv2
    except ImportError:
        sys.stderr.write(
            "Error: --method vng requires OpenCV. Install: pip install opencv-python\n")
        sys.exit(1)

    code_name = _CV2_BAYER_VNG_MAP.get(pattern)
    if code_name is None:
        sys.stderr.write(f"Error: unsupported pattern '{pattern}' for OpenCV.\n")
        sys.exit(1)

    code = getattr(cv2, code_name)

    h, w = raw_2d.shape
    h -= h % 2
    w -= w % 2
    img = raw_2d[:h, :w]

    # OpenCV demosaicing works with uint8 or uint16
    if img.dtype == np.uint16 or img.dtype == np.uint8:
        work = img
    else:
        # Convert to uint16 for processing
        work = img.astype(np.float64)
        work = np.clip(work, 0, 65535)
        work = np.rint(work).astype(np.uint16)

    # cv2.cvtColor returns (H, W, 3) in BGR order
    bgr = cv2.cvtColor(work, code)

    # Convert BGR -> RGB and to (3, H, W)
    rgb = np.transpose(bgr[:, :, ::-1], (2, 0, 1)).astype(np.float64)
    return rgb


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

METHODS = {
    'bilinear': debayer_bilinear,
    'vng': debayer_vng,
}


# ---------------------------------------------------------------------------
# dtype conversion
# ---------------------------------------------------------------------------

def convert_to_original_dtype(data, orig_dtype):
    """Convert float64 work array back to original FITS dtype."""
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        arr = np.clip(data, info.min, info.max)
        arr = np.rint(arr)
        return arr.astype(orig_dtype)
    if np.issubdtype(orig_dtype, np.floating):
        arr = data.astype(orig_dtype)
        bad = ~np.isfinite(arr)
        if np.any(bad):
            arr[bad] = 0
        return arr
    return data.astype(np.float32)


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------

def process_file(infile, outfile, cli_pattern, method):
    """Process a single FITS file."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")
        if hdul[0].data.ndim != 2:
            raise ValueError(
                f"Expected 2D Bayer image, got ndim={hdul[0].data.ndim}")

        data = hdul[0].data
        header = hdul[0].header.copy()
        orig_dtype = data.dtype

    pattern = get_pattern(header, cli_pattern)

    # Demosaic
    demosaic_fn = METHODS[method]
    rgb = demosaic_fn(data, pattern)

    # Convert to original dtype
    out_data = convert_to_original_dtype(rgb, orig_dtype)

    # Update header
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    # Update FILTER to indicate RGB
    header['FILTER'] = ('RGB', 'Demosaiced linear RGB')
    header['HISTORY'] = f"Demosaiced by debayer.py: pattern={pattern}, method={method}"

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fits.PrimaryHDU(out_data, header=header).writeto(outfile, overwrite=True)

    return pattern


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    input_spec, output_spec, cli_pattern, method = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    pat_str = cli_pattern if cli_pattern else "auto"

    ncpu = os.cpu_count() or 1
    workers = min(max(1, ncpu - 1), total)
    print(f"Debayer: {total} file(s), pattern={pat_str}, method={method}, "
          f"threads={workers}")

    done = 0
    errors = 0

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(process_file, infile, outfile, cli_pattern, method):
                os.path.basename(infile)
            for infile, outfile in io_pairs
        }

        for future in as_completed(futures):
            fname = futures[future]
            done += 1
            try:
                future.result()
                sys.stderr.write(f"\r{done}/{total}")
                sys.stderr.flush()
            except Exception as e:
                errors += 1
                sys.stderr.write(f"\n  Error '{fname}': {e}\n")

    sys.stderr.write(f"\nDone. {done - errors} OK, {errors} errors.\n")


if __name__ == "__main__":
    main()
