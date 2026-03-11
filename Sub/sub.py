#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sub - Subtract FITS images or constants, with optional continuum subtraction.

Continuum subtraction (--continuum) is used to isolate narrowband emission
(e.g. H-alpha) by subtracting a scaled broadband image (e.g. R channel).
Stars are detected in both images, cross-matched, and the broadband image
is scaled so that star flux matches — stars subtract to zero, leaving only
emission line signal.
"""

import sys
import os
from astropy.io import fits
import numpy as np

# Add path to shared utilities (batch_utils.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


def usage():
    sys.stderr.write(
        "sub - Subtract FITS images or constants.\n"
        "\n"
        "Usage:\n"
        "  sub.py input.fit output.fit operand [offset] [options]\n"
        "\n"
        "    input.fit   - single file OR numbered pattern (e.g. light0001.fit)\n"
        "    output.fit  - single file OR numbered pattern for results\n"
        "    operand     - number OR FITS file OR numbered FITS pattern\n"
        "                  (value to subtract from input)\n"
        "    offset      - optional numeric value added to the result\n"
        "\n"
        "Options:\n"
        "  --continuum [snr]  Continuum subtraction mode. Detects stars in both\n"
        "                     input and operand, cross-matches them, and scales\n"
        "                     the operand so star flux matches before subtraction.\n"
        "                     Result = input - K * operand (+ offset)\n"
        "                     SNR threshold for star detection (default: 38).\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    continuum = False
    continuum_snr = 38.0

    # Extract options
    i = 0
    positional = []
    while i < len(args):
        if args[i] in ("-h", "--help"):
            usage()
        elif args[i] == "--continuum":
            continuum = True
            i += 1
            # Optional SNR value
            if i < len(args) and not args[i].startswith("--"):
                try:
                    continuum_snr = float(args[i])
                    i += 1
                except ValueError:
                    pass  # not a number, treat as positional
        elif args[i].startswith("--"):
            sys.stderr.write(f"Error: unknown option '{args[i]}'.\n\n")
            usage()
        else:
            positional.append(args[i])
            i += 1

    if len(positional) < 3:
        usage()

    input_pattern = positional[0]
    output_pattern = positional[1]
    operand_str = positional[2]

    offset = 0.0
    if len(positional) >= 4:
        try:
            offset = float(positional[3])
        except ValueError:
            sys.stderr.write("Error: offset must be a number.\n")
            sys.exit(1)

    return input_pattern, output_pattern, operand_str, offset, continuum, continuum_snr


def compute_continuum_scale(data_ha, data_broad, snr=38.0):
    """Compute continuum scaling factor from star flux ratios.

    Detects stars in both images, cross-matches by position (tolerance 1.5 px),
    measures aperture photometry, and returns K = sum(flux_Ha) / sum(flux_broad)
    so that: Ha - K * Broad removes stars to zero.

    Returns (K, n_stars_used).
    """
    import star_utils

    # Detect stars on both images
    cat_ha = star_utils.detect_stars(data_ha, snr=snr)
    cat_br = star_utils.detect_stars(data_broad, snr=snr)

    if cat_ha.n < 3:
        sys.stderr.write(f"  WARNING: only {cat_ha.n} stars in input, "
                         "using K=1.0\n")
        return 1.0, cat_ha.n
    if cat_br.n < 3:
        sys.stderr.write(f"  WARNING: only {cat_br.n} stars in operand, "
                         "using K=1.0\n")
        return 1.0, cat_br.n

    # Overexposure filter: 50% of range
    if np.issubdtype(data_ha.dtype, np.integer):
        peak_limit = 0.5 * float(np.iinfo(data_ha.dtype).max)
    else:
        peak_limit = 0.5 * float(np.max(data_ha))

    # Use Ha catalog as reference, match against broadband
    tolerance = 1.5  # pixels

    # Cross-match: for each Ha star, find nearest broadband star
    matched_ha = []
    matched_br = []
    for i in range(cat_ha.n):
        dist = np.sqrt((cat_br.x - cat_ha.x[i])**2 +
                       (cat_br.y - cat_ha.y[i])**2)
        j = np.argmin(dist)
        if dist[j] <= tolerance:
            matched_ha.append(i)
            matched_br.append(j)

    n_matched = len(matched_ha)
    if n_matched < 3:
        sys.stderr.write(f"  WARNING: only {n_matched} cross-matched stars, "
                         "using K=1.0\n")
        return 1.0, n_matched

    matched_ha = np.array(matched_ha)
    matched_br = np.array(matched_br)

    # Use median FWHM from Ha for aperture sizing
    median_fwhm = float(np.median(cat_ha.fwhm))

    # Aperture radii: +2 px margin for alignment tolerance
    r_aperture = 1.25 * median_fwhm + 2.0
    r_inner = max(2.0 * median_fwhm + 2.0, r_aperture + 2.0)
    r_outer = max(3.0 * median_fwhm + 2.0, r_inner + 3.0)
    r_aperture = max(r_aperture, 4.0)
    r_inner = max(r_inner, r_aperture + 2.0)
    r_outer = max(r_outer, r_inner + 3.0)

    # Use Ha positions for photometry on both images (consistent measurement)
    x_use = cat_ha.x[matched_ha]
    y_use = cat_ha.y[matched_ha]

    flux_ha = star_utils.measure_flux_at(data_ha, x_use, y_use,
                                         r_aperture, r_inner, r_outer)
    flux_br = star_utils.measure_flux_at(data_broad, x_use, y_use,
                                         r_aperture, r_inner, r_outer)

    # Filter: positive flux in both, not overexposed on either image
    valid = (flux_ha > 0) & (flux_br > 0)

    # Peak limit for both images
    if np.issubdtype(data_broad.dtype, np.integer):
        peak_limit_br = 0.5 * float(np.iinfo(data_broad.dtype).max)
    else:
        peak_limit_br = 0.5 * float(np.max(data_broad))

    # Exclude overexposed on either Ha or broadband
    ha_data64 = star_utils._ensure_sep_data(data_ha)
    ha_bkg, _ = star_utils._estimate_background(ha_data64)
    br_data64 = star_utils._ensure_sep_data(data_broad)
    br_bkg, _ = star_utils._estimate_background(br_data64)

    for idx in range(len(matched_ha)):
        if not valid[idx]:
            continue
        iy = min(int(y_use[idx]), ha_bkg.shape[0] - 1)
        ix = min(int(x_use[idx]), ha_bkg.shape[1] - 1)
        # Check Ha
        abs_peak_ha = cat_ha.peak[matched_ha[idx]] + ha_bkg[iy, ix]
        if abs_peak_ha >= peak_limit:
            valid[idx] = False
            continue
        # Check broadband
        iy_b = min(iy, br_bkg.shape[0] - 1)
        ix_b = min(ix, br_bkg.shape[1] - 1)
        abs_peak_br = cat_br.peak[matched_br[idx]] + br_bkg[iy_b, ix_b]
        if abs_peak_br >= peak_limit_br:
            valid[idx] = False

    n_valid = int(np.sum(valid))
    if n_valid < 1:
        sys.stderr.write("  WARNING: no valid matched stars, using K=1.0\n")
        return 1.0, 0

    sum_ha = np.sum(flux_ha[valid])
    sum_br = np.sum(flux_br[valid])

    K = sum_ha / sum_br

    return K, n_valid


def apply_sub_operation(base_data, operand, offset):
    """
    Core arithmetic: result = base - operand + offset

    - base_data: ndarray from input (2D)
    - operand: scalar or ndarray (2D same shape)
    - offset: scalar
    - For integer types: compute in float64, clamp to dtype range, round, cast back.
    - For floats: compute in float64, cast back to original float dtype.
    """
    if base_data is None or base_data.ndim != 2:
        raise ValueError("Expected 2D primary image in input.")

    if isinstance(operand, np.ndarray):
        if operand.shape != base_data.shape:
            raise ValueError(
                f"Operand image shape {operand.shape} does not match input {base_data.shape}."
            )
        op = operand
    else:
        op = float(operand)

    if np.issubdtype(base_data.dtype, np.floating):
        work = base_data.astype(np.float64)
        work = work - op + offset

        if base_data.dtype == np.float32:
            return work.astype(np.float32)
        if base_data.dtype == np.float64:
            return work.astype(np.float64)
        return work.astype(base_data.dtype)

    info = np.iinfo(base_data.dtype)
    work = base_data.astype(np.float64)
    work = work - op + offset

    np.clip(work, info.min, info.max, out=work)
    work = np.rint(work)
    return work.astype(base_data.dtype)


def process_file(infile, outfile, operand_spec, file_index, offset,
                 continuum_K=None):
    """Load, process, and save single FITS file."""
    with fits.open(infile, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"File '{infile}' has no primary image data.")
        if hdul[0].data.ndim != 2:
            raise ValueError(f"File '{infile}' is not a 2D image.")

        data = hdul[0].data
        header = hdul[0].header

        # Get operand (scalar or file path)
        operand_raw = batch_utils.get_operand_for_file(operand_spec, file_index)

        # Resolve to numpy array (handles both scalar constants and FITS files)
        operand = batch_utils.resolve_operand_value(
            operand_raw, data.shape, data.dtype
        )

        # Scale operand for continuum subtraction
        if continuum_K is not None and isinstance(operand, np.ndarray):
            operand = operand.astype(np.float64) * continuum_K

        new_data = apply_sub_operation(data, operand, offset)

        if continuum_K is not None:
            header["HISTORY"] = f"Continuum subtracted (K={continuum_K:.6f})"

        hdul[0].data = new_data
        hdul[0].header = header
        hdul.writeto(outfile, overwrite=True)


def main():
    (input_pattern, output_pattern, operand_str, offset,
     continuum, continuum_snr) = parse_args(sys.argv)

    try:
        io_pairs = batch_utils.build_io_file_lists(input_pattern, output_pattern)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    if not io_pairs:
        sys.stderr.write("Error: no files to process.\n")
        sys.exit(1)

    try:
        operand_spec = batch_utils.build_operand_spec(operand_str, len(io_pairs))
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

    # Compute continuum scale factor
    continuum_K = None
    if continuum:
        # Load first input and first operand to compute K
        first_input = io_pairs[0][0]
        first_operand_raw = batch_utils.get_operand_for_file(operand_spec, 0)

        if isinstance(first_operand_raw, (int, float)):
            sys.stderr.write("Error: --continuum requires FITS file as operand, "
                             "not a constant.\n")
            sys.exit(1)

        sys.stderr.write(f"Continuum: detecting stars (snr={continuum_snr})...\n")

        with fits.open(first_input, memmap=False) as hdul:
            data_ha = hdul[0].data
            if data_ha.ndim != 2:
                sys.stderr.write("Error: continuum requires 2D input.\n")
                sys.exit(1)

        with fits.open(first_operand_raw, memmap=False) as hdul:
            data_broad = hdul[0].data
            if data_broad.ndim != 2:
                sys.stderr.write("Error: continuum requires 2D operand.\n")
                sys.exit(1)

        continuum_K, n_stars = compute_continuum_scale(
            data_ha, data_broad, snr=continuum_snr)

        print(f"Continuum scale: K={continuum_K:.6f} "
              f"(from {n_stars} matched stars)")
        print(f"  Result = input - {continuum_K:.6f} * operand"
              + (f" + {offset}" if offset != 0 else ""))

    total = len(io_pairs)
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            process_file(infile, outfile, operand_spec, i - 1, offset,
                         continuum_K)
            sys.stderr.write(f"\rProcessed {i} / {total} files")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")

    sys.stderr.write("\n")


if __name__ == "__main__":
    main()
