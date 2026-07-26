#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Star detection and measurement utilities based on SEP (Source Extractor as Python).

Provides star detection, FWHM estimation, photometry, and shape measurement
for astronomical FITS images. All functions work with 2D numpy arrays.

Requires: numpy, sep (pip install sep)
"""

import sys
import numpy as np


# ---------------------------------------------------------
# Constants
# ---------------------------------------------------------

# Gaussian FWHM = 2 * sqrt(2 * ln(2)) * sigma
_FWHM_FACTOR = 2.0 * np.sqrt(2.0 * np.log(2.0))  # ~2.3548

# Star filtering defaults
_MIN_FWHM = 1.0        # pixels, below = hot pixel / cosmic ray
_BORDER_MARGIN = 10     # pixels from edge to exclude
_OUTLIER_SIGMA = 3.0    # sigma-clipping for FWHM and elongation outliers


# ---------------------------------------------------------
# SEP configuration
# ---------------------------------------------------------

# sep caps its internal "object pixel buffer" at 300000 active pixels over the
# detection threshold. A large extended bright source (a galaxy or rich
# nebulosity) can put more connected pixels over threshold than that, making
# sep.extract raise "internal pixel buffer full" — which would otherwise abort
# the whole caller. Raise the ceiling once (sep is imported lazily inside the
# functions below, so it is applied on first use, not at module import).
_PIXSTACK_LIMIT = 10_000_000
_pixstack_set = False


def _ensure_pixstack(sep):
    """Raise sep's object-pixel buffer once (idempotent)."""
    global _pixstack_set
    if not _pixstack_set:
        try:
            sep.set_extract_pixstack(_PIXSTACK_LIMIT)
        except Exception:
            pass  # older sep lacks the setter; the try/except around extract still guards
        _pixstack_set = True


# sep's default matched-filter convolution, combined with a per-pixel err map,
# dominates sep.extract on dense/nebulous fields (measured ~5x: 85s -> 17s on a
# 26 MP Sagittarius frame, same source count). The filter only helps pull the
# faintest sources out of the noise; every consumer here (staralign, bestof,
# rgbbalance, sub) ranks by flux and uses BRIGHT stars, so it is unnecessary.
# Disable it for speed while keeping the accurate per-pixel err thresholding.
_EXTRACT_FILTER_KERNEL = None


def _is_pixstack_overflow(exc):
    """True if exc is sep's object-pixel-buffer ('pixstack') overflow."""
    msg = str(exc).lower()
    return ("pixel buffer" in msg or "pixstack" in msg
            or "object pixels" in msg)


# ---------------------------------------------------------
# StarCatalog — result container
# ---------------------------------------------------------

class StarCatalog:
    """
    Result of star detection. All attributes are numpy arrays of length n.

    Attributes
    ----------
    x, y         : sub-pixel positions (0-based, SEP convention)
    flux         : brightness in ADU above background
    peak         : peak pixel value above background
    fwhm         : FWHM in pixels (worst axis = max of fwhm_a, fwhm_b)
    fwhm_a       : FWHM along major axis
    fwhm_b       : FWHM along minor axis
    elongation   : a/b ratio (1.0 = perfectly round)
    theta        : position angle in radians
    flag         : SEP extraction flags (0 = clean)
    n            : number of stars
    """
    __slots__ = ('x', 'y', 'flux', 'peak', 'fwhm', 'fwhm_a', 'fwhm_b',
                 'elongation', 'theta', 'flag', 'n')

    def __init__(self, x, y, flux, peak, fwhm_a, fwhm_b, elongation, theta, flag):
        self.x = np.asarray(x, dtype=np.float64)
        self.y = np.asarray(y, dtype=np.float64)
        self.flux = np.asarray(flux, dtype=np.float64)
        self.peak = np.asarray(peak, dtype=np.float64)
        self.fwhm_a = np.asarray(fwhm_a, dtype=np.float64)
        self.fwhm_b = np.asarray(fwhm_b, dtype=np.float64)
        self.fwhm = np.maximum(self.fwhm_a, self.fwhm_b)
        self.elongation = np.asarray(elongation, dtype=np.float64)
        self.theta = np.asarray(theta, dtype=np.float64)
        self.flag = np.asarray(flag, dtype=np.int32)
        self.n = len(self.x)

    def __len__(self):
        return self.n

    def __repr__(self):
        if self.n == 0:
            return "StarCatalog(0 stars)"
        med = np.median(self.fwhm)
        return f"StarCatalog({self.n} stars, median_fwhm={med:.2f})"

    def sorted_by_flux(self):
        """Return a new StarCatalog sorted by decreasing flux."""
        if self.n == 0:
            return self
        idx = np.argsort(-self.flux)
        return StarCatalog(
            self.x[idx], self.y[idx], self.flux[idx], self.peak[idx],
            self.fwhm_a[idx], self.fwhm_b[idx], self.elongation[idx],
            self.theta[idx], self.flag[idx]
        )

    @staticmethod
    def empty():
        """Return an empty catalog."""
        e = np.array([], dtype=np.float64)
        ei = np.array([], dtype=np.int32)
        return StarCatalog(e, e, e, e, e, e, e, e, ei)


# ---------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------

def _ensure_sep_data(data):
    """Convert data to C-contiguous float64 for SEP. Raises ValueError if not 2D."""
    if data.ndim != 2:
        raise ValueError(f"Expected 2D array, got {data.ndim}D")
    # SEP requires native byte order, C-contiguous
    result = np.ascontiguousarray(data, dtype=np.float64)
    if not result.dtype.isnative:
        result = result.byteswap().newbyteorder()
    return result


def _estimate_background(data, bw=64, bh=64):
    """
    Estimate spatially varying background and its RMS.

    Returns (bkg_image, bkg_rms) as 2D float64 arrays.
    """
    import sep
    bkg = sep.Background(data, bw=bw, bh=bh, fw=3, fh=3)
    return bkg.back(), bkg.rms()


def _fwhm_from_ab(a, b):
    """
    Convert SEP semi-axis lengths to FWHM values.

    Returns (fwhm_a, fwhm_b) — FWHM along major and minor axes.
    """
    fwhm_a = _FWHM_FACTOR * a
    fwhm_b = _FWHM_FACTOR * b
    return fwhm_a, fwhm_b


def _is_hot_pixel(data, iy, ix):
    """
    Check if a pixel is a hot pixel by examining its PSF footprint.

    A real star has a PSF: its neighbors (ring 1) are elevated above the
    surrounding background (ring 2). A hot pixel is isolated: its neighbors
    are at the same level as the background.

    Test:
    1. Center pixel must be > 10*sigma above ring 1 (bright isolated point)
    2. Ring 1 (8 neighbors) must NOT be > 2*sigma above ring 2 (no PSF halo)
    If both true -> hot pixel.

    Ring 1: 8 pixels in 3x3 around center (excluding center)
    Ring 2: 16 pixels in 5x5 border (excluding 3x3 core)
    """
    h, w = data.shape
    if iy < 2 or iy >= h - 2 or ix < 2 or ix >= w - 2:
        return False  # too close to edge for 5x5 check

    # Extract 5x5 patch
    patch = data[iy - 2:iy + 3, ix - 2:ix + 3]

    center = patch[2, 2]

    # Ring 1: 3x3 excluding center
    ring1 = np.array([
        patch[1, 1], patch[1, 2], patch[1, 3],
        patch[2, 1],              patch[2, 3],
        patch[3, 1], patch[3, 2], patch[3, 3]
    ])

    # Ring 2: 5x5 border (all pixels not in 3x3 core)
    ring2 = np.array([
        patch[0, 0], patch[0, 1], patch[0, 2], patch[0, 3], patch[0, 4],
        patch[1, 0],                                         patch[1, 4],
        patch[2, 0],                                         patch[2, 4],
        patch[3, 0],                                         patch[3, 4],
        patch[4, 0], patch[4, 1], patch[4, 2], patch[4, 3], patch[4, 4]
    ])

    # Condition 1: center is isolated bright point
    # Use median instead of mean — robust to hot pixel neighbors in ring1
    med1 = np.median(ring1)
    mad1 = np.median(np.abs(ring1 - med1))
    sigma1 = mad1 * 1.4826  # MAD to sigma
    if sigma1 <= 0:
        sigma1 = np.std(ring1)
    if sigma1 <= 0:
        return False
    if (center - med1) < 10.0 * sigma1:
        return False  # not an isolated bright point

    # Condition 2: ring 1 has no PSF halo (not elevated above ring 2)
    # A real star's PSF elevates multiple neighbors above the background.
    # A hot pixel is isolated: its neighbors are at background level.
    # Count ring1 pixels individually above 3*sigma of ring2.
    mean2 = np.mean(ring2)
    std2 = np.std(ring2)
    if std2 <= 0:
        return False

    n_bright_neighbors = np.sum(ring1 > mean2 + 3.0 * std2)
    if n_bright_neighbors >= 4:
        return False  # majority of neighbors bright -> real star PSF

    return True  # isolated bright pixel, no halo -> hot pixel


def _reject_hot_pixels(data, catalog, mask):
    """
    Post-filter: reject hot pixels from SEP catalog using PSF footprint test.

    For each detected source, checks if the peak pixel is an isolated bright
    point with no PSF halo. This correctly handles optics with small FWHM
    (~1.6 px) where real stars still have elevated neighbors.
    """
    xpeak = catalog['xpeak']
    ypeak = catalog['ypeak']

    indices = np.where(mask)[0]
    for i in indices:
        ix = int(xpeak[i])
        iy = int(ypeak[i])
        if _is_hot_pixel(data, iy, ix):
            mask[i] = False

    return mask


def _filter_stars(data, catalog, shape):
    """
    Build a boolean mask to keep only stellar objects.

    Three-pass filtering:
    1. Hard filters: flagged, border, too small
    2. Hot pixel rejection: PSF footprint test (isolated bright points)
    3. Adaptive filters: reject FWHM and elongation outliers relative
       to the median of the population. This handles smeared frames
       (all stars elongated) and varying seeing conditions.
    """
    h, w = shape
    n = len(catalog)
    if n == 0:
        return np.array([], dtype=bool)

    x = catalog['x']
    y = catalog['y']
    a = catalog['a']
    b = catalog['b']
    flag = catalog['flag']

    fwhm_a, fwhm_b = _fwhm_from_ab(a, b)
    fwhm_worst = np.maximum(fwhm_a, fwhm_b)
    elongation = np.where(b > 0, a / b, 99.0)

    # Pass 1: hard filters (non-negotiable)
    mask = np.ones(n, dtype=bool)
    mask &= (flag == 0)
    mask &= (x >= _BORDER_MARGIN) & (x < w - _BORDER_MARGIN)
    mask &= (y >= _BORDER_MARGIN) & (y < h - _BORDER_MARGIN)
    mask &= (fwhm_worst >= _MIN_FWHM)

    if np.sum(mask) < 3:
        return mask

    # Pass 2: reject hot pixels by PSF footprint analysis
    mask = _reject_hot_pixels(data, catalog, mask)

    if np.sum(mask) < 3:
        return mask

    # Pass 3: adaptive sigma-clipping on FWHM and elongation
    # Galaxies/nebula knots have larger FWHM than stars;
    # satellite trails have higher elongation.
    # On a smeared frame all stars share similar elongation,
    # so median+sigma naturally adapts.
    med_fwhm = np.median(fwhm_worst[mask])
    mad_fwhm = np.median(np.abs(fwhm_worst[mask] - med_fwhm))
    sigma_fwhm = mad_fwhm * 1.4826  # MAD to sigma
    sigma_fwhm = max(sigma_fwhm, 0.5)  # floor to avoid zero spread

    med_elong = np.median(elongation[mask])
    mad_elong = np.median(np.abs(elongation[mask] - med_elong))
    sigma_elong = mad_elong * 1.4826
    sigma_elong = max(sigma_elong, 0.1)

    mask &= (fwhm_worst <= med_fwhm + _OUTLIER_SIGMA * sigma_fwhm)
    mask &= (elongation <= med_elong + _OUTLIER_SIGMA * sigma_elong)

    return mask


def _build_catalog(catalog, mask, flux_values):
    """Build StarCatalog from SEP catalog, mask, and flux array."""
    if not np.any(mask):
        return StarCatalog.empty()

    c = catalog[mask]
    fwhm_a, fwhm_b = _fwhm_from_ab(c['a'], c['b'])
    elongation = np.where(c['b'] > 0, c['a'] / c['b'], 99.0)

    return StarCatalog(
        x=c['x'], y=c['y'],
        flux=flux_values[mask],
        peak=c['peak'],
        fwhm_a=fwhm_a, fwhm_b=fwhm_b,
        elongation=elongation,
        theta=c['theta'],
        flag=c['flag']
    )


# ---------------------------------------------------------
# Public API
# ---------------------------------------------------------

def estimate_fwhm(data, snr=38.0, bw=64, bh=64):
    """
    Estimate median FWHM of stars on the frame (worst axis).

    Uses a high SNR threshold to select only bright, well-measured stars.
    Returns FWHM in pixels, or NaN if no stars detected.

    Parameters
    ----------
    data : 2D numpy array
    snr  : minimum SNR for detection (default 38.0, higher = fewer but better stars)
    bw, bh : background mesh size in pixels
    """
    import sep
    _ensure_pixstack(sep)

    data = _ensure_sep_data(data)
    bkg_image, bkg_rms = _estimate_background(data, bw, bh)
    data_sub = data - bkg_image

    try:
        catalog = sep.extract(data_sub, thresh=snr, err=bkg_rms,
                              filter_kernel=_EXTRACT_FILTER_KERNEL)
    except Exception as exc:
        if _is_pixstack_overflow(exc):
            print("WARNING: sep pixel buffer overflow (large extended source); "
                  "no FWHM estimate", file=sys.stderr)
            return float('nan')
        raise
    mask = _filter_stars(data_sub, catalog, data.shape)

    if not np.any(mask):
        print("WARNING: no stars detected for FWHM estimation", file=sys.stderr)
        return float('nan')

    a = catalog['a'][mask]
    b = catalog['b'][mask]
    fwhm_a, fwhm_b = _fwhm_from_ab(a, b)
    fwhm_worst = np.maximum(fwhm_a, fwhm_b)

    return float(np.median(fwhm_worst))


def detect_stars(data, snr=5.0, photometry=False, bw=64, bh=64):
    """
    Detect stars and measure their properties.

    Parameters
    ----------
    data : 2D numpy array
    snr  : minimum SNR threshold for detection (default 5.0)
    photometry : if True, measure flux via aperture photometry with local
                 background annulus (accurate in nebulae/gradients, slower).
                 If False, use isophotal flux from SEP (fast, uses global
                 background model).
    bw, bh : background mesh size in pixels

    Returns
    -------
    StarCatalog with sub-pixel positions, flux, FWHM, elongation, etc.
    """
    import sep
    _ensure_pixstack(sep)

    data = _ensure_sep_data(data)
    bkg_image, bkg_rms = _estimate_background(data, bw, bh)
    data_sub = data - bkg_image

    try:
        catalog = sep.extract(data_sub, thresh=snr, err=bkg_rms,
                              filter_kernel=_EXTRACT_FILTER_KERNEL)
    except Exception as exc:
        if _is_pixstack_overflow(exc):
            # A galaxy-/nebula-dominated frame is not a star field: degrade to
            # "no stars" instead of crashing the caller.
            print("WARNING: sep pixel buffer overflow (large extended source); "
                  "treating frame as no stars", file=sys.stderr)
            return StarCatalog.empty()
        raise
    mask = _filter_stars(data_sub, catalog, data.shape)

    if not np.any(mask):
        return StarCatalog.empty()

    if not photometry:
        # Fast path: isophotal flux from SEP (above global background model)
        return _build_catalog(catalog, mask, catalog['flux'])

    # Accurate path: aperture photometry with local background annulus
    # Determine aperture radius from median FWHM of detected stars
    a_filt = catalog['a'][mask]
    b_filt = catalog['b'][mask]
    fwhm_a, fwhm_b = _fwhm_from_ab(a_filt, b_filt)
    median_fwhm = float(np.median(np.maximum(fwhm_a, fwhm_b)))

    r_aperture = 2.5 * median_fwhm / 2.0  # radius = 1.25 * FWHM
    r_inner = 2.0 * median_fwhm            # inner annulus radius
    r_outer = 3.0 * median_fwhm            # outer annulus radius

    # Ensure minimum sizes
    r_aperture = max(r_aperture, 2.0)
    r_inner = max(r_inner, r_aperture + 2.0)
    r_outer = max(r_outer, r_inner + 3.0)

    x_all = catalog['x']
    y_all = catalog['y']

    # Aperture sum on background-subtracted data
    ap_flux, ap_fluxerr, ap_flag = sep.sum_circle(
        data_sub, x_all, y_all, r_aperture
    )

    # Local background from annulus on ORIGINAL data (not bg-subtracted)
    ann_flux, ann_fluxerr, ann_flag = sep.sum_circann(
        data, x_all, y_all, r_inner, r_outer
    )
    ann_area = np.pi * (r_outer**2 - r_inner**2)
    ap_area = np.pi * r_aperture**2
    local_bg_per_pixel = ann_flux / ann_area

    # Net flux = aperture on original - local background * aperture area
    # Equivalent to: (original - local_bg) summed over aperture
    net_flux = sep.sum_circle(data, x_all, y_all, r_aperture)[0] \
               - local_bg_per_pixel * ap_area

    return _build_catalog(catalog, mask, net_flux)


def measure_flux_at(data, x, y, r_aperture, r_inner, r_outer, bw=64, bh=64):
    """
    Measure aperture photometry flux at given positions on a 2D channel.

    Uses local background annulus for accurate measurement in nebulae/gradients.

    Parameters
    ----------
    data : 2D numpy array (single channel)
    x, y : arrays of star positions (from detection on another channel)
    r_aperture : aperture radius in pixels
    r_inner, r_outer : background annulus radii
    bw, bh : background mesh size

    Returns
    -------
    net_flux : array of net fluxes (aperture - local background)
    """
    import sep

    data = _ensure_sep_data(data)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    # Aperture sum on raw data
    ap_flux = sep.sum_circle(data, x, y, r_aperture)[0]

    # Local background from annulus
    ann_flux = sep.sum_circann(data, x, y, r_inner, r_outer)[0]
    ann_area = np.pi * (r_outer**2 - r_inner**2)
    ap_area = np.pi * r_aperture**2
    local_bg_per_pixel = ann_flux / ann_area

    net_flux = ap_flux - local_bg_per_pixel * ap_area
    return net_flux


def debug_mark_stars(data, catalog, radius_scale=4.0):
    """
    Draw circles around detected stars on a copy of the image.

    Each circle has 1px width and diameter = radius_scale * 2 * FWHM.
    Circle pixels are brightened by 10% of the image full range above
    the local pixel value, so they are visible but don't obscure the image.

    Parameters
    ----------
    data : 2D numpy array (original image)
    catalog : StarCatalog from detect_stars()
    radius_scale : circle radius multiplier (default 4.0, so diameter = 8*FWHM)

    Returns
    -------
    2D numpy array (float64) with circles drawn
    """
    out = data.astype(np.float64).copy()
    if catalog.n == 0:
        return out

    full_range = float(np.max(out) - np.min(out))
    bump = 0.1 * full_range

    h, w = out.shape

    for i in range(catalog.n):
        cx = catalog.x[i]
        cy = catalog.y[i]
        r = catalog.fwhm[i] * radius_scale

        if r < 2.0:
            r = 2.0

        # Draw circle: iterate angles, oversample for smooth 1px line
        n_points = max(int(2 * np.pi * r * 4), 64)
        angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
        px = np.round(cx + r * np.cos(angles)).astype(int)
        py = np.round(cy + r * np.sin(angles)).astype(int)

        # Clip to image bounds and brighten
        valid = (px >= 0) & (px < w) & (py >= 0) & (py < h)
        out[py[valid], px[valid]] += bump

    return out
