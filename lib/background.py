# -*- coding: utf-8 -*-
"""
Background estimation library for astronomical images.

Cell-based background estimation with sigma-clipped medians,
median filtering, and polynomial surface fitting.

Usage:
    from background import estimate_background

    # Get background model (2D array same size as input)
    model = estimate_background(data2d)

    # Subtract background (preserves mean level)
    flattened = data2d - model + np.median(model)

    # With mask-center to exclude bright central object
    model = estimate_background(data2d, mask_center_d=0.6)
"""

import warnings
import numpy as np
from scipy.ndimage import median_filter

# Suppress warnings from nanmedian/nanstd on all-NaN slices (border cells)
warnings.filterwarnings('ignore', message='.*All-NaN slice.*', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*Degrees of freedom <= 0.*', category=RuntimeWarning)


def sigma_clipped_median(values, k=1.7, max_iter=5):
    """Iterative sigma-clipped median. Ignores zeros."""
    v = values[values != 0.0]
    if v.size == 0:
        return 0.0
    for _ in range(max_iter):
        med = np.median(v)
        std = np.std(v)
        if std == 0:
            break
        mask = np.abs(v - med) <= k * std
        if np.all(mask):
            break
        v = v[mask]
        if v.size == 0:
            return 0.0
    return float(np.median(v))


def _build_cell_grid(data, cell_size, clip_k, max_iter=5):
    """
    Divide image into cells, compute sigma-clipped median per cell.
    Vectorized: all cells processed simultaneously (releases GIL for threading).
    """
    h, w = data.shape
    ny = (h + cell_size - 1) // cell_size
    nx = (w + cell_size - 1) // cell_size
    n_cells = ny * nx

    # Pad to exact multiple of cell_size (zeros → will become NaN)
    pad_h = ny * cell_size - h
    pad_w = nx * cell_size - w
    if pad_h > 0 or pad_w > 0:
        padded = np.pad(data, ((0, pad_h), (0, pad_w)),
                        mode='constant', constant_values=0)
    else:
        padded = data

    # Reshape into cells: (ny*nx, cell_size²)
    cells = (padded.reshape(ny, cell_size, nx, cell_size)
                   .transpose(0, 2, 1, 3)
                   .reshape(n_cells, cell_size * cell_size))

    # Work in float64, zeros → NaN (ignored by nanmedian/nanstd)
    work = cells.astype(np.float64)
    work[work == 0.0] = np.nan

    # Iterative sigma-clipped median (vectorized across all cells)
    for _ in range(max_iter):
        med = np.nanmedian(work, axis=1)   # (n_cells,)
        std = np.nanstd(work, axis=1)      # (n_cells,)

        # Cells with std==0 are fully converged
        active = std > 0
        if not np.any(active):
            break

        # Deviation from median
        dev = np.abs(work - med[:, np.newaxis])
        threshold = (clip_k * std)[:, np.newaxis]

        # New outliers (only in active cells)
        new_outliers = dev > threshold
        new_outliers[~active] = False

        # Check if any cell still has new outliers
        if not np.any(new_outliers & np.isfinite(work)):
            break

        work[new_outliers] = np.nan

    # Final median (all-NaN slices → NaN → replaced with 0 below)
    result = np.nanmedian(work, axis=1)
    result = np.nan_to_num(result, nan=0.0)

    return result.reshape(ny, nx)


def _median_filter_ignore_zeros(grid, size):
    """Median filter on a grid, ignoring zero cells via sentinel trick."""
    zero_mask = (grid == 0.0)
    if not np.any(zero_mask):
        return median_filter(grid, size=size, mode='nearest')

    nz = grid[~zero_mask]
    if nz.size == 0:
        return grid.copy()

    sentinel = float(np.max(nz)) + abs(float(np.max(nz))) + 1.0
    work = np.where(zero_mask, sentinel, grid)
    med = median_filter(work, size=size, mode='nearest')
    med[med >= sentinel] = 0.0
    return med


def _apply_mask_center(grid, mask_center_d):
    """Zero out cells inside a central ellipse of size D * grid dimensions."""
    ny, nx = grid.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    ry, rx = mask_center_d * ny / 2.0, mask_center_d * nx / 2.0

    if ry < 1 or rx < 1:
        return grid

    gy, gx = np.mgrid[0:ny, 0:nx]
    ellipse = ((gy - cy) / ry) ** 2 + ((gx - cx) / rx) ** 2
    out = grid.copy()
    out[ellipse <= 1.0] = 0.0
    return out


def _fit_poly2d(data, order, coords=None):
    """
    Fit 2D polynomial surface to non-zero data points.
    coords: optional (norm_x, norm_y) for grid points.
    Returns (coeffs, terms) or (None, None) on failure.
    """
    h, w = data.shape

    if coords is not None:
        norm_x, norm_y = coords
        mask = (data != 0.0)
        N = int(mask.sum())
        if N == 0:
            return None, None
        z = data[mask]
        x_flat = norm_x[mask]
        y_flat = norm_y[mask]
    else:
        yy, xx = np.mgrid[0:h, 0:w]
        mask = (data != 0.0)
        N = int(mask.sum())
        if N == 0:
            return None, None
        z = data[mask]
        x = xx[mask].astype(np.float64)
        y = yy[mask].astype(np.float64)
        x_flat = (2.0 * (x / (w - 1.0)) - 1.0) if w > 1 else np.zeros_like(x)
        y_flat = (2.0 * (y / (h - 1.0)) - 1.0) if h > 1 else np.zeros_like(y)

    def num_terms(k):
        return (k + 1) * (k + 2) // 2

    used_order = min(order, 20)
    while used_order > 0 and num_terms(used_order) > N:
        used_order -= 1

    terms = [(i, j) for i in range(used_order + 1) for j in range(used_order + 1 - i)]
    cols = [(x_flat ** i) * (y_flat ** j) for (i, j) in terms]
    A = np.vstack(cols).T

    coeffs, _, _, _ = np.linalg.lstsq(A, z, rcond=None)
    return coeffs, terms


def _render_poly2d(height, width, coeffs, terms):
    """Evaluate 2D polynomial surface at every pixel."""
    yy, xx = np.mgrid[0:height, 0:width]

    x = (2.0 * (xx / (width - 1.0)) - 1.0) if width > 1 else np.zeros_like(xx, dtype=np.float64)
    y = (2.0 * (yy / (height - 1.0)) - 1.0) if height > 1 else np.zeros_like(yy, dtype=np.float64)

    model = np.zeros((height, width), dtype=np.float64)
    for c, (i, j) in zip(coeffs, terms):
        if c != 0.0:
            model += c * (x ** i) * (y ** j)

    return model


def estimate_background(data, cell_size=64, clip_k=1.7, poly_order=3,
                        median1=3, median2=5, border=2, mask_center_d=None):
    """
    Estimate background surface for a 2D image.

    Algorithm:
      1. Divide into cells, compute sigma-clipped median per cell
      2. Apply median filter on coarse grid (zeros excluded)
      3. Extend grid by border cells
      4. Apply second median filter on extended grid
      5. Fit polynomial surface (zeros excluded from fit)
      6. Evaluate polynomial at full resolution

    Parameters
    ----------
    data : 2D numpy array (float64 recommended)
    cell_size : int - cell size in pixels (default 64)
    clip_k : float - sigma clipping threshold (default 1.7)
    poly_order : int - polynomial order (default 3)
    median1 : int - first median filter size on grid (default 3)
    median2 : int - second median filter size on extended grid (default 5)
    border : int - border extension in cells (default 2)
    mask_center_d : float or None - if set, exclude central ellipse
                    of size D*W x D*H from estimation (default None)

    Returns
    -------
    model : 2D numpy array (float64), same shape as input.
            Background surface. To flatten: result = data - model + np.median(model)
    """
    h, w = data.shape
    data64 = np.asarray(data, dtype=np.float64)

    # Step 1: coarse grid
    grid = _build_cell_grid(data64, cell_size, clip_k)

    # Step 1b: mask center
    if mask_center_d is not None:
        grid = _apply_mask_center(grid, mask_center_d)

    # Step 2: first median filter
    grid = _median_filter_ignore_zeros(grid, median1)

    # Step 3: extend grid
    grid_ext = np.pad(grid, border, mode='edge')
    eny, enx = grid_ext.shape

    # Step 4: second median filter
    grid_ext = _median_filter_ignore_zeros(grid_ext, median2)

    # Step 5: polynomial fit
    gy_ext, gx_ext = np.mgrid[0:eny, 0:enx]
    pixel_y = ((gy_ext - border) + 0.5) * cell_size
    pixel_x = ((gx_ext - border) + 0.5) * cell_size

    norm_x = (2.0 * pixel_x / (w - 1.0) - 1.0) if w > 1 else np.zeros_like(pixel_x)
    norm_y = (2.0 * pixel_y / (h - 1.0) - 1.0) if h > 1 else np.zeros_like(pixel_y)

    coeffs, terms = _fit_poly2d(grid_ext, poly_order, coords=(norm_x, norm_y))
    if coeffs is None:
        return np.zeros((h, w), dtype=np.float64)

    # Step 6: render at full resolution
    return _render_poly2d(h, w, coeffs, terms)


def estimate_background_poly(data, cell_size=64, clip_k=1.7, poly_order=3,
                             median1=3, median2=5, border=2, mask_center_d=None):
    """
    Estimate background polynomial coefficients (for deferred rendering).

    Same algorithm as estimate_background, but returns (coeffs, terms)
    instead of the rendered model. Use render_background() to render later.
    Returns (None, None) if fit fails.
    """
    h, w = data.shape
    data64 = np.asarray(data, dtype=np.float64)

    grid = _build_cell_grid(data64, cell_size, clip_k)
    if mask_center_d is not None:
        grid = _apply_mask_center(grid, mask_center_d)
    grid = _median_filter_ignore_zeros(grid, median1)
    grid_ext = np.pad(grid, border, mode='edge')
    eny, enx = grid_ext.shape
    grid_ext = _median_filter_ignore_zeros(grid_ext, median2)

    gy_ext, gx_ext = np.mgrid[0:eny, 0:enx]
    pixel_y = ((gy_ext - border) + 0.5) * cell_size
    pixel_x = ((gx_ext - border) + 0.5) * cell_size
    norm_x = (2.0 * pixel_x / (w - 1.0) - 1.0) if w > 1 else np.zeros_like(pixel_x)
    norm_y = (2.0 * pixel_y / (h - 1.0) - 1.0) if h > 1 else np.zeros_like(pixel_y)

    return _fit_poly2d(grid_ext, poly_order, coords=(norm_x, norm_y))


def render_background(height, width, coeffs, terms):
    """Render background model from cached polynomial coefficients."""
    if coeffs is None:
        return np.zeros((height, width), dtype=np.float64)
    return _render_poly2d(height, width, coeffs, terms)
