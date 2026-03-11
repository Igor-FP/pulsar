#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rgbbalance - RGB color balance and brightness normalization for color FITS.

Neutralizes background color, optionally applies white balance (auto or manual),
and normalizes brightness across a sequence of frames.

Input FITS must be 3-channel RGB: NAXIS=3, NAXIS3=3, shape (3, H, W).

Algorithm:
  1. Compute per-channel quantile medians (bottom kmin%, top kmax%)
  2. Shift each channel so dark background is neutral (common black level)
  3. Apply per-channel color balance coefficients (--auto or --rgb)
  4. Normalize brightness across frames (reference = first or specified file)

Modes:
  --auto [file]   Auto white balance with reference from file (default: first input)
  --autoeach      Auto white balance computed independently per file
  --rgb R G B     Manual per-channel coefficients
  (none)          Only background neutralization + brightness normalization
"""

import sys
import os
import numpy as np
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils
import star_utils


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "rgbbalance - RGB color balance and brightness normalization.\n"
        "\n"
        "Usage:\n"
        "  rgbbalance.py input_spec output_spec [options]\n"
        "\n"
        "Positional arguments:\n"
        "  input_spec       Input file(s): single file, wildcard (*.fit),\n"
        "                   numbered sequence (img0001.fit), or list (@list.txt)\n"
        "  output_spec      Output: single file, numbered pattern (out0001.fit),\n"
        "                   or directory\n"
        "\n"
        "Options:\n"
        "  --auto [file]    Auto white balance: scale R & B to match G range.\n"
        "                   Reference from file (default: first input file).\n"
        "  --autoeach       Auto white balance computed independently per file.\n"
        "  --autostar [file] Star-based white balance: detect stars on G channel,\n"
        "                   measure flux on R/G/B via aperture photometry,\n"
        "                   compute K from total star flux ratios.\n"
        "                   Reference from file (default: first input file).\n"
        "  --rgb R G B      Manual per-channel scaling coefficients.\n"
        "  --kmin N         Percent of darkest pixels for black level  (default: 5)\n"
        "  --kmax N         Percent of brightest pixels for white level (default: 5)\n"
        "  --snr N          SNR threshold for --autostar detection (default: 38)\n"
        "\n"
        "  --rgb, --auto, --autoeach, --autostar are mutually exclusive.\n"
        "  Without any: only background neutralization + brightness normalization.\n"
        "\n"
        "Examples:\n"
        "  rgbbalance.py img.fit out.fit --auto\n"
        "  rgbbalance.py img0001.fit out0001.fit --auto ref.fit\n"
        "  rgbbalance.py img0001.fit out0001.fit --autoeach\n"
        "  rgbbalance.py img0001.fit out0001.fit --autostar\n"
        "  rgbbalance.py img0001.fit out0001.fit --rgb 1.1 2.0 0.5\n"
        "  rgbbalance.py img0001.fit out0001.fit --auto --kmin 10 --kmax 3\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    rgb_k = None
    auto = False
    autoeach = False
    autostar = False
    auto_ref_file = None
    kmin = 5.0
    kmax = 5.0
    snr = 38.0

    # Parse option flags
    i = 0
    positional = []
    while i < len(args):
        if args[i] == "--rgb":
            if i + 3 >= len(args):
                sys.stderr.write("Error: --rgb requires 3 values (R G B).\n")
                sys.exit(1)
            try:
                rgb_k = (float(args[i+1]), float(args[i+2]), float(args[i+3]))
            except ValueError:
                sys.stderr.write("Error: --rgb values must be numbers.\n")
                sys.exit(1)
            i += 4
        elif args[i] == "--auto":
            auto = True
            i += 1
            # Optional reference filename (next arg if not a flag)
            if i < len(args) and not args[i].startswith("--"):
                candidate = args[i]
                if "." in os.path.basename(candidate):
                    auto_ref_file = candidate
                    i += 1
        elif args[i] == "--autoeach":
            autoeach = True
            i += 1
        elif args[i] == "--autostar":
            autostar = True
            i += 1
            # Optional reference filename
            if i < len(args) and not args[i].startswith("--"):
                candidate = args[i]
                if "." in os.path.basename(candidate):
                    auto_ref_file = candidate
                    i += 1
        elif args[i] == "--kmin":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --kmin requires a value.\n")
                sys.exit(1)
            try:
                kmin = float(args[i+1])
            except ValueError:
                sys.stderr.write("Error: --kmin value must be a number.\n")
                sys.exit(1)
            i += 2
        elif args[i] == "--kmax":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --kmax requires a value.\n")
                sys.exit(1)
            try:
                kmax = float(args[i+1])
            except ValueError:
                sys.stderr.write("Error: --kmax value must be a number.\n")
                sys.exit(1)
            i += 2
        elif args[i] == "--snr":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --snr requires a value.\n")
                sys.exit(1)
            try:
                snr = float(args[i+1])
            except ValueError:
                sys.stderr.write("Error: --snr value must be a number.\n")
                sys.exit(1)
            i += 2
        else:
            positional.append(args[i])
            i += 1

    # Validate mutual exclusivity
    mode_count = sum([rgb_k is not None, auto, autoeach, autostar])
    if mode_count > 1:
        sys.stderr.write("Error: --rgb, --auto, --autoeach, --autostar are mutually exclusive.\n")
        sys.exit(1)

    if len(positional) < 2:
        usage()

    input_spec = positional[0]
    output_spec = positional[1]

    return (input_spec, output_spec, rgb_k, auto, autoeach, autostar,
            auto_ref_file, kmin, kmax, snr)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def compute_channel_stats(channel_2d, kmin, kmax):
    """Compute quantile medians for a single 2D channel.

    Returns (min_median, max_median):
      min_median = median of the darkest kmin% of pixels
      max_median = median of the brightest kmax% of pixels
    """
    flat = channel_2d.ravel().astype(np.float64)

    # Exclude zero-padded pixels (from crop) for bottom percentile
    nonzero = flat[flat > 0]
    if len(nonzero) == 0:
        return 0.0, 0.0
    n = len(nonzero)

    # Bottom kmin%
    n_lo = max(1, int(n * kmin / 100.0))
    partitioned = np.partition(nonzero, n_lo)
    min_median = float(np.median(partitioned[:n_lo]))

    # Top kmax% (also on nonzero data for consistency)
    n_hi = max(1, int(n * kmax / 100.0))
    idx = n - n_hi
    partitioned = np.partition(nonzero, idx)
    max_median = float(np.median(partitioned[idx:]))

    return min_median, max_median


def compute_frame_stats(data_3d, kmin, kmax):
    """Compute (min_median, max_median) for each of the 3 RGB channels."""
    stats = []
    for ch in range(3):
        stats.append(compute_channel_stats(data_3d[ch], kmin, kmax))
    return stats  # [(min_R, max_R), (min_G, max_G), (min_B, max_B)]


# ---------------------------------------------------------------------------
# Reference calibration
# ---------------------------------------------------------------------------

def compute_reference(stats, rgb_k, auto):
    """Compute reference values from frame statistics.

    Returns (black, range_ref, K) where K = (K_R, K_G, K_B).
    """
    min_meds = [s[0] for s in stats]
    max_meds = [s[1] for s in stats]
    ranges = [mx - mn for mn, mx in stats]

    # Common black level = average of per-channel dark medians
    black = sum(min_meds) / 3.0

    # Reference dynamic range = average of per-channel ranges
    range_ref = sum(ranges) / 3.0

    # Color balance coefficients
    if rgb_k is not None:
        K = rgb_k
    elif auto:
        range_g = ranges[1]  # green is reference
        K_r = (range_g / ranges[0]) if ranges[0] > 0 else 1.0
        K_g = 1.0
        K_b = (range_g / ranges[2]) if ranges[2] > 0 else 1.0
        K = (K_r, K_g, K_b)
    else:
        K = (1.0, 1.0, 1.0)

    return black, range_ref, K


# ---------------------------------------------------------------------------
# Star-based color balance
# ---------------------------------------------------------------------------

def compute_star_balance(data_3d, snr=38.0):
    """Compute RGB balance coefficients from star photometry.

    Detects stars on the G channel, cross-matches with R and B detections
    (tolerance 1.5 px), filters out overexposed stars (peak > 50% range),
    and computes K from total flux ratios.

    Parameters
    ----------
    data_3d : 3D array (3, H, W) — RGB channels
    snr : SNR threshold for star detection

    Returns
    -------
    K : (K_R, K_G, K_B) coefficients
    n_stars : number of stars used
    """
    ch_r, ch_g, ch_b = data_3d[0], data_3d[1], data_3d[2]

    # Detect stars on G channel
    cat_g = star_utils.detect_stars(ch_g, snr=snr)
    if cat_g.n < 3:
        sys.stderr.write(f"  WARNING: only {cat_g.n} stars on G channel, "
                         "falling back to K=1.0\n")
        return (1.0, 1.0, 1.0), cat_g.n

    # Determine full range for overexposure filter (50% of dtype range)
    if np.issubdtype(data_3d.dtype, np.integer):
        max_range = float(np.iinfo(data_3d.dtype).max)
    else:
        max_range = float(np.max(data_3d))
    peak_limit = 0.5 * max_range

    # Filter: exclude overexposed stars (peak in absolute ADU)
    # cat_g.peak is above background, estimate absolute peak
    g_data = star_utils._ensure_sep_data(ch_g)
    g_bkg, _ = star_utils._estimate_background(g_data)
    # Approximate absolute peak = peak_above_bg + local background
    abs_peak = cat_g.peak.copy()
    for i in range(cat_g.n):
        iy = min(int(cat_g.y[i]), g_bkg.shape[0] - 1)
        ix = min(int(cat_g.x[i]), g_bkg.shape[1] - 1)
        abs_peak[i] += g_bkg[iy, ix]

    bright_ok = abs_peak < peak_limit

    if np.sum(bright_ok) < 3:
        sys.stderr.write(f"  WARNING: only {np.sum(bright_ok)} non-saturated stars, "
                         "using all stars\n")
        bright_ok = np.ones(cat_g.n, dtype=bool)

    # Positions and FWHM of usable stars
    x = cat_g.x[bright_ok]
    y = cat_g.y[bright_ok]
    fwhm = cat_g.fwhm[bright_ok]
    median_fwhm = float(np.median(fwhm))

    # Aperture radii: +2 px margin for chromatic shift tolerance
    r_aperture = 1.25 * median_fwhm + 2.0
    r_inner = max(2.0 * median_fwhm + 2.0, r_aperture + 2.0)
    r_outer = max(3.0 * median_fwhm + 2.0, r_inner + 3.0)

    # Ensure minimum sizes
    r_aperture = max(r_aperture, 4.0)
    r_inner = max(r_inner, r_aperture + 2.0)
    r_outer = max(r_outer, r_inner + 3.0)

    # Cross-match: detect on R and B, keep only stars near G positions
    cat_r = star_utils.detect_stars(ch_r, snr=snr)
    cat_b = star_utils.detect_stars(ch_b, snr=snr)

    tolerance = 1.5  # pixels

    # Find G stars that have matches in both R and B
    matched = np.zeros(len(x), dtype=bool)
    for i in range(len(x)):
        # Check R channel
        if cat_r.n > 0:
            dist_r = np.sqrt((cat_r.x - x[i])**2 + (cat_r.y - y[i])**2)
            has_r = np.min(dist_r) <= tolerance
        else:
            has_r = False
        # Check B channel
        if cat_b.n > 0:
            dist_b = np.sqrt((cat_b.x - x[i])**2 + (cat_b.y - y[i])**2)
            has_b = np.min(dist_b) <= tolerance
        else:
            has_b = False
        matched[i] = has_r and has_b

    n_matched = int(np.sum(matched))
    if n_matched < 3:
        sys.stderr.write(f"  WARNING: only {n_matched} cross-matched stars, "
                         "skipping match filter\n")
        matched = np.ones(len(x), dtype=bool)
        n_matched = len(x)

    x_use = x[matched]
    y_use = y[matched]

    # Aperture photometry on all 3 channels at same positions
    flux_r = star_utils.measure_flux_at(ch_r, x_use, y_use,
                                        r_aperture, r_inner, r_outer)
    flux_g = star_utils.measure_flux_at(ch_g, x_use, y_use,
                                        r_aperture, r_inner, r_outer)
    flux_b = star_utils.measure_flux_at(ch_b, x_use, y_use,
                                        r_aperture, r_inner, r_outer)

    # Use only stars with positive flux in all channels
    valid = (flux_r > 0) & (flux_g > 0) & (flux_b > 0)
    n_valid = int(np.sum(valid))
    if n_valid < 1:
        sys.stderr.write("  WARNING: no stars with positive flux in all channels\n")
        return (1.0, 1.0, 1.0), 0

    sum_r = np.sum(flux_r[valid])
    sum_g = np.sum(flux_g[valid])
    sum_b = np.sum(flux_b[valid])

    K_r = sum_g / sum_r
    K_b = sum_g / sum_b
    K = (K_r, 1.0, K_b)

    return K, n_valid


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
    # Fallback
    arr = data.astype(np.float32)
    bad = ~np.isfinite(arr)
    if np.any(bad):
        arr[bad] = 0
    return arr


# ---------------------------------------------------------------------------
# Frame processing
# ---------------------------------------------------------------------------

def _apply_balance(work, kmin, kmax, black, range_ref, K):
    """Apply balance to float64 3D array in-place. Returns brightness factor."""
    stats = compute_frame_stats(work, kmin, kmax)
    min_meds = [s[0] for s in stats]
    ranges = [s[1] - s[0] for s in stats]
    range_current = sum(ranges) / 3.0

    # Brightness normalization factor
    if range_ref > 0 and range_current > 0:
        brightness = range_ref / range_current
    else:
        brightness = 1.0

    # Process each channel
    for ch in range(3):
        work[ch] -= (min_meds[ch] - black)
        work[ch] = (work[ch] - black) * brightness * K[ch] + black

    return brightness


def process_file(infile, outfile, kmin, kmax, black, range_ref, K):
    """Process a single RGB FITS file with precomputed reference.

    Returns brightness factor.
    """
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")
        if hdul[0].data.ndim != 3 or hdul[0].data.shape[0] != 3:
            raise ValueError(
                f"Expected 3-channel RGB (3, H, W), got shape {hdul[0].data.shape}")

        orig_dtype = hdul[0].data.dtype
        header = hdul[0].header.copy()
        work = hdul[0].data.astype(np.float64)

    # Preserve zero-padded regions (from crop): mask pixels <= 0
    zero_mask = work <= 0

    brightness = _apply_balance(work, kmin, kmax, black, range_ref, K)

    # Restore zeros in padded regions
    work[zero_mask] = 0.0

    # Sanitize NaN/Inf
    work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)

    out_data = convert_to_original_dtype(work, orig_dtype)

    # Update header
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    header["HISTORY"] = "Color balanced by rgbbalance.py"
    header["HISTORY"] = f"  black={black:.2f}, brightness={brightness:.6f}"
    header["HISTORY"] = f"  K_R={K[0]:.6f}, K_G={K[1]:.6f}, K_B={K[2]:.6f}"

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fits.PrimaryHDU(out_data, header=header).writeto(outfile, overwrite=True)

    return brightness


def process_file_autoeach(infile, outfile, kmin, kmax):
    """Process a single RGB FITS file with independent auto balance.

    Computes its own reference from itself. Returns (K, brightness).
    """
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")
        if hdul[0].data.ndim != 3 or hdul[0].data.shape[0] != 3:
            raise ValueError(
                f"Expected 3-channel RGB (3, H, W), got shape {hdul[0].data.shape}")

        orig_dtype = hdul[0].data.dtype
        header = hdul[0].header.copy()
        work = hdul[0].data.astype(np.float64)

    # Preserve zero-padded regions (from crop)
    zero_mask = work <= 0

    # Compute reference from this file itself
    stats = compute_frame_stats(work, kmin, kmax)
    black, range_ref, K = compute_reference(stats, None, True)

    brightness = _apply_balance(work, kmin, kmax, black, range_ref, K)

    # Restore zeros in padded regions
    work[zero_mask] = 0.0

    # Sanitize NaN/Inf
    work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)

    out_data = convert_to_original_dtype(work, orig_dtype)

    # Update header
    for key in ("BSCALE", "BZERO"):
        if key in header:
            del header[key]

    header["HISTORY"] = "Color balanced by rgbbalance.py (autoeach)"
    header["HISTORY"] = f"  black={black:.2f}, brightness={brightness:.6f}"
    header["HISTORY"] = f"  K_R={K[0]:.6f}, K_G={K[1]:.6f}, K_B={K[2]:.6f}"

    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fits.PrimaryHDU(out_data, header=header).writeto(outfile, overwrite=True)

    return K, brightness


# ---------------------------------------------------------------------------
# Reference computation
# ---------------------------------------------------------------------------

def compute_ref_from_file(ref_file, kmin, kmax, rgb_k, auto):
    """Read a FITS file and compute reference values (black, range_ref, K).

    Does NOT write anything — only reads and computes stats.
    """
    with fits.open(ref_file, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError("No primary image data.")
        if hdul[0].data.ndim != 3 or hdul[0].data.shape[0] != 3:
            raise ValueError(
                f"Expected 3-channel RGB (3, H, W), got shape {hdul[0].data.shape}")
        data = hdul[0].data.astype(np.float64)

    stats = compute_frame_stats(data, kmin, kmax)
    black, range_ref, K = compute_reference(stats, rgb_k, auto)

    return black, range_ref, K, stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    (input_spec, output_spec,
     rgb_k, auto, autoeach, autostar, auto_ref_file,
     kmin, kmax, snr) = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    total = len(io_pairs)
    ncpu = os.cpu_count() or 1
    workers = min(max(1, ncpu - 1), total)

    if autostar:
        # --- Star-based white balance ---
        # Compute K from reference file, then apply as fixed coefficients
        ref_file = auto_ref_file if auto_ref_file else io_pairs[0][0]

        sys.stderr.write(f"Detecting stars on {os.path.basename(ref_file)} (snr={snr})...\n")
        with fits.open(ref_file, memmap=False) as hdul:
            ref_data = hdul[0].data.astype(np.float64)

        K_star, n_stars = compute_star_balance(ref_data, snr=snr)

        print(f"RGB balance: {total} file(s), mode=autostar, "
              f"ref={os.path.basename(ref_file)}")
        print(f"  Stars used: {n_stars}, "
              f"K_R={K_star[0]:.6f}, K_G={K_star[1]:.6f}, K_B={K_star[2]:.6f}")

        # Compute background reference for neutralization + brightness norm
        try:
            black, range_ref, _, ref_stats = compute_ref_from_file(
                ref_file, kmin, kmax, None, False)
        except Exception as e:
            sys.stderr.write(f"Error reading reference: {e}\n")
            sys.exit(1)

        # Use star-derived K
        K = K_star

        done = 0
        errors = 0
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_file, infile, outfile,
                            kmin, kmax, black, range_ref, K):
                    os.path.basename(infile)
                for infile, outfile in io_pairs
            }
            for future in as_completed(futures):
                fname = futures[future]
                done += 1
                try:
                    brightness = future.result()
                    sys.stderr.write(f"\r{done}/{total}  bright={brightness:.4f}")
                    sys.stderr.flush()
                except Exception as e:
                    errors += 1
                    sys.stderr.write(f"\n  Error '{fname}': {e}\n")

    elif autoeach:
        # --- Per-file independent auto balance ---
        print(f"RGB balance: {total} file(s), mode=autoeach, "
              f"kmin={kmin}%, kmax={kmax}%, threads={workers}")

        done = 0
        errors = 0

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_file_autoeach, infile, outfile, kmin, kmax):
                    os.path.basename(infile)
                for infile, outfile in io_pairs
            }

            for future in as_completed(futures):
                fname = futures[future]
                done += 1
                try:
                    K, brightness = future.result()
                    sys.stderr.write(
                        f"\r{done}/{total}  {fname}"
                        f"  K_R={K[0]:.3f} K_B={K[2]:.3f} bright={brightness:.4f}")
                    sys.stderr.flush()
                except Exception as e:
                    errors += 1
                    sys.stderr.write(f"\n  Error '{fname}': {e}\n")

    else:
        # --- Reference-based balance (auto/rgb/neutral) ---
        mode = "auto" if auto else ("rgb" if rgb_k is not None else "neutral")

        # Determine reference file
        if auto_ref_file:
            ref_file = auto_ref_file
        else:
            ref_file = io_pairs[0][0]

        # Compute reference BEFORE launching threads
        try:
            black, range_ref, K, ref_stats = compute_ref_from_file(
                ref_file, kmin, kmax, rgb_k, auto)
        except Exception as e:
            sys.stderr.write(f"Error reading reference file '{ref_file}': {e}\n")
            sys.exit(1)

        ranges_ref = [s[1] - s[0] for s in ref_stats]
        ref_label = os.path.basename(ref_file)
        print(f"RGB balance: {total} file(s), mode={mode}, ref={ref_label}, "
              f"threads={workers}")
        print(f"  black={black:.2f}, range: R={ranges_ref[0]:.2f} "
              f"G={ranges_ref[1]:.2f} B={ranges_ref[2]:.2f} (avg={range_ref:.2f})")
        if auto or rgb_k:
            print(f"  K_R={K[0]:.6f}, K_G={K[1]:.6f}, K_B={K[2]:.6f}")

        # Process all files with thread pool
        done = 0
        errors = 0

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_file, infile, outfile,
                            kmin, kmax, black, range_ref, K):
                    os.path.basename(infile)
                for infile, outfile in io_pairs
            }

            for future in as_completed(futures):
                fname = futures[future]
                done += 1
                try:
                    brightness = future.result()
                    sys.stderr.write(f"\r{done}/{total}  bright={brightness:.4f}")
                    sys.stderr.flush()
                except Exception as e:
                    errors += 1
                    sys.stderr.write(f"\n  Error '{fname}': {e}\n")

    sys.stderr.write(f"\nDone. {done - errors} OK, {errors} errors.\n")


if __name__ == "__main__":
    main()
