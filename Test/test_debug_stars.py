#!/usr/bin/env python3
"""Save a debug FITS with circles around detected stars."""
import sys, os
import numpy as np
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import star_utils

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} input.fit [output_debug.fit] [--snr N]")
    sys.exit(1)

# Parse args
args = sys.argv[1:]
snr = 10.0
if "--snr" in args:
    idx = args.index("--snr")
    snr = float(args[idx + 1])
    args = args[:idx] + args[idx + 2:]

infile = args[0]
outfile = args[1] if len(args) > 1 else os.path.splitext(infile)[0] + "_debug.fit"

with fits.open(infile, memmap=False) as hdul:
    data = hdul[0].data
    header = hdul[0].header

if data.ndim == 3:
    data = data[0]

print(f"Image: {infile}, shape={data.shape}, dtype={data.dtype}")

catalog = star_utils.detect_stars(data, snr=snr)
print(f"Detected: {catalog.n} stars (snr={snr}), median FWHM={np.median(catalog.fwhm):.2f}")

debug_img = star_utils.debug_mark_stars(data, catalog)

hdu = fits.PrimaryHDU(debug_img.astype(np.float32), header=header)
hdu.writeto(outfile, overwrite=True)
print(f"Saved: {outfile}")
