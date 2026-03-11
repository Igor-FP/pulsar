#!/usr/bin/env python3
"""Diagnostic: test hot pixel rejection on a single FITS file."""
import sys, os
import numpy as np
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import star_utils

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} file.fit")
    sys.exit(1)

with fits.open(sys.argv[1], memmap=False) as hdul:
    data = hdul[0].data
    if data.ndim == 3:
        data = data[0]

print(f"Image: {sys.argv[1]}, shape={data.shape}, dtype={data.dtype}")

# Run SEP
import sep
data64 = star_utils._ensure_sep_data(data)
bkg_img, bkg_rms = star_utils._estimate_background(data64)
data_sub = data64 - bkg_img

catalog = sep.extract(data_sub, thresh=10.0, err=bkg_rms)
print(f"SEP extracted: {len(catalog)} objects")

# Pass 1 only (hard filters)
n = len(catalog)
fwhm_a, fwhm_b = star_utils._fwhm_from_ab(catalog['a'], catalog['b'])
fwhm_worst = np.maximum(fwhm_a, fwhm_b)

mask1 = np.ones(n, dtype=bool)
mask1 &= (catalog['flag'] == 0)
margin = 10
h, w = data.shape
mask1 &= (catalog['x'] >= margin) & (catalog['x'] < w - margin)
mask1 &= (catalog['y'] >= margin) & (catalog['y'] < h - margin)
mask1 &= (fwhm_worst >= 1.0)
print(f"After pass 1 (hard filters): {np.sum(mask1)} objects")
print(f"  FWHM range: {fwhm_worst[mask1].min():.2f} — {fwhm_worst[mask1].max():.2f}")
print(f"  FWHM median: {np.median(fwhm_worst[mask1]):.2f}")

# Test hot pixel rejection on first N detections
N_SHOW = int(sys.argv[2]) if len(sys.argv) > 2 else 50
print(f"\nTesting _is_hot_pixel on first {N_SHOW} detections:")
indices = np.where(mask1)[0][:N_SHOW]
for i in indices:
    ix = int(catalog['xpeak'][i])
    iy = int(catalog['ypeak'][i])
    fwhm = fwhm_worst[i]
    peak = catalog['peak'][i]

    # Manual check
    if iy < 2 or iy >= h - 2 or ix < 2 or ix >= w - 2:
        print(f"  #{i} xpeak={ix} ypeak={iy} fwhm={fwhm:.2f} peak={peak:.0f} — EDGE")
        continue

    patch = data_sub[iy - 2:iy + 3, ix - 2:ix + 3]
    center = patch[2, 2]
    ring1 = np.array([
        patch[1,1], patch[1,2], patch[1,3],
        patch[2,1],             patch[2,3],
        patch[3,1], patch[3,2], patch[3,3]
    ])
    ring2 = np.array([
        patch[0,0], patch[0,1], patch[0,2], patch[0,3], patch[0,4],
        patch[1,0],                                       patch[1,4],
        patch[2,0],                                       patch[2,4],
        patch[3,0],                                       patch[3,4],
        patch[4,0], patch[4,1], patch[4,2], patch[4,3], patch[4,4]
    ])

    med1 = np.median(ring1)
    mad1 = np.median(np.abs(ring1 - med1))
    sigma1 = mad1 * 1.4826 if mad1 > 0 else np.std(ring1)
    mean2, std2 = np.mean(ring2), np.std(ring2)

    n_bright = int(np.sum(ring1 > mean2 + 3.0 * std2))
    is_hot = star_utils._is_hot_pixel(data_sub, iy, ix)

    print(f"  #{i} fwhm={fwhm:.2f} peak={peak:.0f} center={center:.0f} "
          f"med1={med1:.1f} sigma1={sigma1:.1f} "
          f"mean2={mean2:.1f} std2={std2:.1f} "
          f"bright_nb={n_bright} "
          f"cond1={(center-med1):.0f}>={10*sigma1:.0f}? "
          f"hot={is_hot}")

# Run full filter
mask_full = star_utils._filter_stars(data_sub, catalog, data.shape)
print(f"\nAfter full _filter_stars: {np.sum(mask_full)} objects")
if np.any(mask_full):
    print(f"  FWHM range: {fwhm_worst[mask_full].min():.2f} — {fwhm_worst[mask_full].max():.2f}")
    print(f"  FWHM median: {np.median(fwhm_worst[mask_full]):.2f}")
