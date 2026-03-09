#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
crop - Crop FITS images.

Two modes:
  1. --width WW --height HH  [--center XX YY | --autocenter [PERCENTILE] | --positions FILE]
     Crop to fixed size around a center point.
  2. --top N --bottom N --left N --right N
     Trim N pixels from each side (margins).

If the crop box extends beyond image boundaries, missing pixels are filled
with zeros (zero-padding).

Supports 2D (mono) and 3D (color, N×H×W) FITS images.
"""

import sys
import os
import csv
import numpy as np
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import shared utilities (two levels up from Personal/Crop/)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AUTOCENTER_MEDIAN_RADIUS = 4  # median filter kernel = 2*R+1 = 17


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "crop - Crop FITS images.\n"
        "\n"
        "Usage:\n"
        "  crop.py input_spec output_spec [options]\n"
        "\n"
        "Positional arguments:\n"
        "  input_spec       Input file(s): single file, wildcard (*.fit),\n"
        "                   numbered sequence (img0001.fit), or list (@list.txt)\n"
        "  output_spec      Output: single file, numbered pattern (out0001.fit),\n"
        "                   or directory\n"
        "\n"
        "Mode 1 — crop to size:\n"
        "  --width WW         Output width in pixels (required with --height)\n"
        "  --height HH        Output height in pixels (required with --width)\n"
        "  --center XX YY     Center point (pixel coords, default: image center)\n"
        "  --autocenter [P]   Auto-detect center by brightness per file (default P=50)\n"
        "  --positions FILE   Per-file centers from CSV\n"
        "\n"
        "Positions CSV format:\n"
        "  filename,cx,cy\n"
        "  009A7220.fit,1960,1091\n"
        "  009A7221.fit,1960,1091\n"
        "  Extra columns (e.g. radius) are ignored.\n"
        "\n"
        "Mode 2 — trim margins:\n"
        "  --top N            Pixels to trim from top (default: 0)\n"
        "  --bottom N         Pixels to trim from bottom (default: 0)\n"
        "  --left N           Pixels to trim from left (default: 0)\n"
        "  --right N          Pixels to trim from right (default: 0)\n"
        "\n"
        "Modes are mutually exclusive. --center, --autocenter, and --positions\n"
        "are mutually exclusive (all require --width/--height).\n"
        "\n"
        "If crop box extends beyond image boundaries, missing pixels\n"
        "are filled with zeros.\n"
        "\n"
        "Examples:\n"
        "  crop.py img.fit out.fit --width 1000 --height 1000\n"
        "  crop.py img.fit out.fit --width 1000 --height 1000 --center 3000 2000\n"
        "  crop.py img.fit out.fit --width 1000 --height 1000 --autocenter\n"
        "  crop.py img.fit out.fit --width 1000 --height 1000 --autocenter 30\n"
        "  crop.py img0001.fit out0001.fit --width 1000 --height 1000 --positions pos.csv\n"
        "  crop.py img0001.fit out0001.fit --top 100 --bottom 100 --left 200 --right 200\n"
    )
    sys.exit(1)


def _parse_int(s, name):
    """Parse and validate a positive integer parameter."""
    try:
        v = int(s)
    except ValueError:
        sys.stderr.write(f"Error: {name} must be an integer.\n")
        sys.exit(1)
    if v < 0:
        sys.stderr.write(f"Error: {name} must be >= 0.\n")
        sys.exit(1)
    return v


def parse_args(argv):
    args = argv[1:]

    width = None
    height = None
    center_x = None
    center_y = None
    autocenter = False
    percentile = 50.0
    positions_file = None
    margin_top = 0
    margin_bottom = 0
    margin_left = 0
    margin_right = 0
    has_margins = False

    i = 0
    positional = []
    while i < len(args):
        if args[i] == "--positions":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --positions requires a file path.\n")
                sys.exit(1)
            positions_file = args[i + 1]
            i += 2
        elif args[i] == "--width":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --width requires a value.\n")
                sys.exit(1)
            width = _parse_int(args[i+1], "--width")
            if width <= 0:
                sys.stderr.write("Error: --width must be > 0.\n")
                sys.exit(1)
            i += 2
        elif args[i] == "--height":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --height requires a value.\n")
                sys.exit(1)
            height = _parse_int(args[i+1], "--height")
            if height <= 0:
                sys.stderr.write("Error: --height must be > 0.\n")
                sys.exit(1)
            i += 2
        elif args[i] == "--center":
            if i + 2 >= len(args):
                sys.stderr.write("Error: --center requires XX YY.\n")
                sys.exit(1)
            center_x = _parse_int(args[i+1], "--center XX")
            center_y = _parse_int(args[i+2], "--center YY")
            i += 3
        elif args[i] == "--autocenter":
            autocenter = True
            i += 1
            # Optional percentile argument (next arg if it's a number)
            if i < len(args) and not args[i].startswith("--"):
                try:
                    percentile = float(args[i])
                    i += 1
                except ValueError:
                    pass  # not a number, leave as positional
        elif args[i] == "--top":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --top requires a value.\n")
                sys.exit(1)
            margin_top = _parse_int(args[i+1], "--top")
            has_margins = True
            i += 2
        elif args[i] == "--bottom":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --bottom requires a value.\n")
                sys.exit(1)
            margin_bottom = _parse_int(args[i+1], "--bottom")
            has_margins = True
            i += 2
        elif args[i] == "--left":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --left requires a value.\n")
                sys.exit(1)
            margin_left = _parse_int(args[i+1], "--left")
            has_margins = True
            i += 2
        elif args[i] == "--right":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --right requires a value.\n")
                sys.exit(1)
            margin_right = _parse_int(args[i+1], "--right")
            has_margins = True
            i += 2
        else:
            positional.append(args[i])
            i += 1

    if len(positional) < 2:
        usage()

    has_size = (width is not None or height is not None)

    # Validate mode exclusivity
    if has_size and has_margins:
        sys.stderr.write("Error: --width/--height and --top/--bottom/--left/--right "
                         "are mutually exclusive.\n")
        sys.exit(1)

    if not has_size and not has_margins:
        usage()

    # Size mode validation
    if has_size:
        if width is None or height is None:
            sys.stderr.write("Error: both --width and --height are required.\n")
            sys.exit(1)
        # --center, --autocenter, --positions are mutually exclusive
        exclusive_count = sum([
            center_x is not None,
            autocenter,
            positions_file is not None,
        ])
        if exclusive_count > 1:
            sys.stderr.write("Error: --center, --autocenter, and --positions "
                             "are mutually exclusive.\n")
            sys.exit(1)

    # Margins mode: --center/--autocenter/--positions not allowed
    if has_margins:
        if center_x is not None or autocenter or positions_file is not None:
            sys.stderr.write("Error: --center/--autocenter/--positions "
                             "require --width/--height.\n")
            sys.exit(1)

    return (positional[0], positional[1],
            width, height, center_x, center_y,
            autocenter, percentile, positions_file,
            margin_top, margin_bottom, margin_left, margin_right)


# ---------------------------------------------------------------------------
# Autocenter
# ---------------------------------------------------------------------------

def find_autocenter(data, percentile):
    """Find brightness center of an image.

    1. Average channels for color images.
    2. Median-filter to suppress noise.
    3. Build max-projection rows/columns.
    4. Threshold at percentile% between min and max.
    5. Find first/last pixel above threshold -> center.
    """
    from scipy.ndimage import median_filter

    # Get 2D plane
    if data.ndim == 3:
        plane = np.mean(data.astype(np.float64), axis=0)
    else:
        plane = data.astype(np.float64)

    # Median filter
    kernel = 2 * AUTOCENTER_MEDIAN_RADIUS + 1
    filtered = median_filter(plane, size=kernel)

    # Center X: max-projection along rows -> 1D array of length W
    max_row = np.max(filtered, axis=0)
    val_max = max_row.max()
    val_min = max_row.min()
    threshold_x = val_min + (val_max - val_min) * percentile / 100.0
    above_x = np.where(max_row >= threshold_x)[0]
    if len(above_x) > 0:
        cx = int((above_x[0] + above_x[-1]) // 2)
    else:
        cx = len(max_row) // 2

    # Center Y: max-projection along columns -> 1D array of length H
    max_col = np.max(filtered, axis=1)
    val_max = max_col.max()
    val_min = max_col.min()
    threshold_y = val_min + (val_max - val_min) * percentile / 100.0
    above_y = np.where(max_col >= threshold_y)[0]
    if len(above_y) > 0:
        cy = int((above_y[0] + above_y[-1]) // 2)
    else:
        cy = len(max_col) // 2

    return cx, cy


# ---------------------------------------------------------------------------
# Crop logic
# ---------------------------------------------------------------------------

def compute_crop_box_size(img_h, img_w, crop_w, crop_h, cx, cy):
    """Compute crop box for size mode. Returns (x0, y0, x1, y1)."""
    x0 = cx - crop_w // 2
    y0 = cy - crop_h // 2
    x1 = x0 + crop_w
    y1 = y0 + crop_h
    return x0, y0, x1, y1


def compute_crop_box_margins(img_h, img_w, top, bottom, left, right):
    """Compute crop box for margin mode. Returns (x0, y0, x1, y1)."""
    x0 = left
    y0 = top
    x1 = img_w - right
    y1 = img_h - bottom
    if x1 <= x0 or y1 <= y0:
        raise ValueError(
            f"Margins too large: result would be {x1 - x0}x{y1 - y0} pixels.")
    return x0, y0, x1, y1


def crop_data(data, x0, y0, x1, y1):
    """Crop data array with zero-padding for out-of-bounds regions."""
    if data.ndim == 2:
        img_h, img_w = data.shape
    else:
        img_h, img_w = data.shape[-2], data.shape[-1]

    crop_w = x1 - x0
    crop_h = y1 - y0

    # Intersection with image
    src_x0 = max(0, x0)
    src_y0 = max(0, y0)
    src_x1 = min(img_w, x1)
    src_y1 = min(img_h, y1)

    # Position in output
    dst_x0 = src_x0 - x0
    dst_y0 = src_y0 - y0
    dst_x1 = dst_x0 + (src_x1 - src_x0)
    dst_y1 = dst_y0 + (src_y1 - src_y0)

    # Build output shape
    if data.ndim == 2:
        out = np.zeros((crop_h, crop_w), dtype=data.dtype)
        out[dst_y0:dst_y1, dst_x0:dst_x1] = data[src_y0:src_y1, src_x0:src_x1]
    else:
        out_shape = data.shape[:-2] + (crop_h, crop_w)
        out = np.zeros(out_shape, dtype=data.dtype)
        out[..., dst_y0:dst_y1, dst_x0:dst_x1] = data[..., src_y0:src_y1, src_x0:src_x1]

    return out


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------

def _write_cropped(data, header, outfile, x0, y0, x1, y1, history_str):
    """Crop data and write output FITS."""
    out_data = crop_data(data, x0, y0, x1, y1)

    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    header["HISTORY"] = history_str

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fits.PrimaryHDU(out_data, header=header).writeto(outfile, overwrite=True)


def process_file(infile, outfile, crop_box, history_str):
    """Process a single FITS file with a precomputed crop box."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")
        data = hdul[0].data
        header = hdul[0].header.copy()

    _write_cropped(data, header, outfile, *crop_box, history_str)


def process_file_autocenter(infile, outfile, width, height, percentile):
    """Process a single FITS file with per-file autocenter.

    Returns (cx, cy) for progress reporting.
    """
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")
        data = hdul[0].data
        header = hdul[0].header.copy()

    cx, cy = find_autocenter(data, percentile)
    crop_box = compute_crop_box_size(
        data.shape[-2], data.shape[-1], width, height, cx, cy)
    x0, y0, x1, y1 = crop_box

    history_str = (f"Cropped by crop.py: [{x0}:{x1}, {y0}:{y1}] "
                   f"{x1 - x0}x{y1 - y0} autocenter=({cx},{cy})")

    _write_cropped(data, header, outfile, x0, y0, x1, y1, history_str)
    return cx, cy


# ---------------------------------------------------------------------------
# Positions CSV
# ---------------------------------------------------------------------------

def load_positions_csv(csv_path):
    """Load positions from CSV. Returns dict {basename: (cx, cy)}.

    Expected CSV columns: filename, cx, cy
    Optional column: radius (ignored by crop, used by align.py)
    """
    positions = {}
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(
            (row for row in f if not row.startswith("#")),
        )
        for row in reader:
            positions[row["filename"]] = (
                int(row["cx"]),
                int(row["cy"]),
            )
    if not positions:
        raise ValueError(f"No positions found in {csv_path}")
    return positions


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    (input_spec, output_spec,
     width, height, center_x, center_y,
     autocenter, percentile, positions_file,
     margin_top, margin_bottom, margin_left, margin_right) = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    is_size_mode = (width is not None)

    ncpu = os.cpu_count() or 1
    workers = min(max(1, ncpu - 1), total)
    done = 0
    errors = 0

    if is_size_mode and positions_file:
        # --- Per-file positions from CSV ---
        try:
            pos_dict = load_positions_csv(positions_file)
        except Exception as e:
            sys.stderr.write(f"Error loading positions: {e}\n")
            sys.exit(1)

        # Validate: all input files must be in CSV
        missing = []
        for infile, _ in io_pairs:
            bname = os.path.basename(infile)
            if bname not in pos_dict:
                missing.append(bname)
        if missing:
            sys.stderr.write(f"Error: {len(missing)} file(s) not found in positions CSV:\n")
            for m in missing[:10]:
                sys.stderr.write(f"  {m}\n")
            if len(missing) > 10:
                sys.stderr.write(f"  ... and {len(missing) - 10} more\n")
            sys.exit(1)

        print(f"Crop to {width}x{height}, positions from {positions_file}, "
              f"{total} file(s), threads={workers}")

        # Pre-compute crop boxes and history strings for each file
        file_tasks = []
        for infile, outfile in io_pairs:
            bname = os.path.basename(infile)
            cx, cy = pos_dict[bname]
            crop_box = compute_crop_box_size(0, 0, width, height, cx, cy)
            x0, y0, x1, y1 = crop_box
            history_str = (f"Cropped by crop.py: [{x0}:{x1}, {y0}:{y1}] "
                           f"{x1 - x0}x{y1 - y0} positions=({cx},{cy})")
            file_tasks.append((infile, outfile, crop_box, history_str))

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_file, infile, outfile, crop_box, history_str):
                    os.path.basename(infile)
                for infile, outfile, crop_box, history_str in file_tasks
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

    elif is_size_mode and autocenter:
        # --- Per-file autocenter: each file computes its own center ---
        print(f"Crop to {width}x{height}, autocenter per file "
              f"(percentile={percentile}), {total} file(s), threads={workers}")

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_file_autocenter,
                            infile, outfile, width, height, percentile):
                    os.path.basename(infile)
                for infile, outfile in io_pairs
            }

            for future in as_completed(futures):
                fname = futures[future]
                done += 1
                try:
                    cx, cy = future.result()
                    sys.stderr.write(f"\r{done}/{total}  {fname} center=({cx},{cy})")
                    sys.stderr.flush()
                except Exception as e:
                    errors += 1
                    sys.stderr.write(f"\n  Error '{fname}': {e}\n")

    else:
        # --- Fixed crop box (computed once from first file) ---
        if is_size_mode:
            # --center or geometric center
            with fits.open(io_pairs[0][0], memmap=False) as hdul:
                if hdul[0].data is None:
                    sys.stderr.write("Error: no primary image data.\n")
                    sys.exit(1)
                img_h = hdul[0].data.shape[-2]
                img_w = hdul[0].data.shape[-1]

            if center_x is not None:
                cx, cy = center_x, center_y
            else:
                cx = img_w // 2
                cy = img_h // 2

            crop_box = compute_crop_box_size(img_h, img_w, width, height, cx, cy)
            x0, y0, x1, y1 = crop_box
            history_str = (f"Cropped by crop.py: [{x0}:{x1}, {y0}:{y1}] "
                           f"{x1 - x0}x{y1 - y0} center=({cx},{cy})")
            print(f"Crop to {width}x{height}, center=({cx},{cy}), "
                  f"box=[{x0}:{x1}, {y0}:{y1}], {total} file(s)")
        else:
            # Margins mode
            with fits.open(io_pairs[0][0], memmap=False) as hdul:
                if hdul[0].data is None:
                    sys.stderr.write("Error: no primary image data.\n")
                    sys.exit(1)
                img_h = hdul[0].data.shape[-2]
                img_w = hdul[0].data.shape[-1]

            try:
                crop_box = compute_crop_box_margins(
                    img_h, img_w,
                    margin_top, margin_bottom, margin_left, margin_right)
            except ValueError as e:
                sys.stderr.write(f"Error: {e}\n")
                sys.exit(1)

            x0, y0, x1, y1 = crop_box
            out_w = x1 - x0
            out_h = y1 - y0
            history_str = (f"Cropped by crop.py: margins T={margin_top} B={margin_bottom} "
                           f"L={margin_left} R={margin_right} -> {out_w}x{out_h}")
            print(f"Crop margins T={margin_top} B={margin_bottom} "
                  f"L={margin_left} R={margin_right} -> {out_w}x{out_h}, "
                  f"{total} file(s)")

        # crop_box and history_str computed BEFORE threads start
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_file, infile, outfile, crop_box, history_str):
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
