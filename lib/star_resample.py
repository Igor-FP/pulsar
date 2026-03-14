#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pixel resampling for star alignment.

Applies a TPS (or projective/similarity) transformation to remap
an image from target geometry to reference geometry.

V1: uses scipy.ndimage.map_coordinates(order=3) with post-hoc clamping.
Future: Mitchell-Netravali, Lanczos-3, per-pixel adaptive filters (see todo.md).

TPS coordinate map optimization: evaluate TPS on a coarse grid (every Nth
pixel), then upscale to full resolution via bicubic interpolation. This
reduces TPS evaluations from H*W (~26M) to (H/N)*(W/N) (~25K for N=32).

Pure functions, no file I/O.
"""

import numpy as np
from scipy.ndimage import map_coordinates

# Coarse grid step for TPS coordinate map (pixels).
# TPS distortion is smooth, so even step=64 gives sub-pixel accuracy.
# step=32 is conservative for safety.
_TPS_GRID_STEP = 32


def _build_coordinate_map_tps_exact(tps_x, tps_y, output_shape):
    """Build coordinate map by evaluating TPS on EVERY pixel.

    DO NOT DELETE — reference implementation for accuracy validation.
    Correct but very slow: ~440s for 6248x4176 with ~340 control points.
    Use _build_coordinate_map_tps() (coarse grid + upscale) instead.
    """
    h, w = output_shape
    yy, xx = np.mgrid[0:h, 0:w]
    grid_pts = np.column_stack([xx.ravel(), yy.ravel()])
    tgt_x = tps_x(grid_pts).reshape(h, w)
    tgt_y = tps_y(grid_pts).reshape(h, w)
    return tgt_y, tgt_x


def _build_coordinate_map_tps(tps_x, tps_y, output_shape):
    """Build inverse coordinate map using TPS model via coarse grid + upscale.

    Evaluates TPS on a coarse grid (every _TPS_GRID_STEP pixels),
    then upscales coordinate maps to full resolution using bicubic
    interpolation. TPS distortion is smooth by definition (minimizes
    bending energy), so this introduces negligible error.

    Parameters
    ----------
    tps_x, tps_y : RBFInterpolator — maps reference coords → target coords
                   (inverse direction TPS)
    output_shape : tuple (H, W)

    Returns
    -------
    coords_y : ndarray (H, W) — target y coordinate for each output pixel
    coords_x : ndarray (H, W) — target x coordinate for each output pixel
    """
    h, w = output_shape
    step = _TPS_GRID_STEP

    # Coarse grid coordinates (always include edges)
    gy = np.arange(0, h, step, dtype=np.float64)
    gx = np.arange(0, w, step, dtype=np.float64)
    # Ensure last row/col is included
    if gy[-1] != h - 1:
        gy = np.append(gy, h - 1)
    if gx[-1] != w - 1:
        gx = np.append(gx, w - 1)

    # Build coarse grid points
    gxx, gyy = np.meshgrid(gx, gy)  # (ngy, ngx)
    coarse_pts = np.column_stack([gxx.ravel(), gyy.ravel()])

    # Evaluate TPS on coarse grid
    coarse_tx = tps_x(coarse_pts).reshape(len(gy), len(gx))
    coarse_ty = tps_y(coarse_pts).reshape(len(gy), len(gx))

    # Upscale to full resolution using bicubic interpolation
    # map_coordinates needs coordinates in the coarse grid's index space
    # Map full pixel coords → fractional coarse grid indices
    # gy/gx are irregularly spaced (last step may be shorter),
    # so use np.interp for mapping
    full_y = np.arange(h, dtype=np.float64)
    full_x = np.arange(w, dtype=np.float64)

    # Map pixel coordinates to coarse grid indices
    coarse_iy = np.interp(full_y, gy, np.arange(len(gy)))
    coarse_ix = np.interp(full_x, gx, np.arange(len(gx)))

    # Build full-res index arrays
    iy_2d, ix_2d = np.meshgrid(coarse_iy, coarse_ix, indexing='ij')

    # Interpolate coarse maps to full resolution
    coords_x = map_coordinates(coarse_tx, [iy_2d, ix_2d],
                               order=3, mode='nearest')
    coords_y = map_coordinates(coarse_ty, [iy_2d, ix_2d],
                               order=3, mode='nearest')

    return coords_y, coords_x


def _build_coordinate_map_projective(H_inv, output_shape):
    """Build inverse coordinate map using projective (homography) transform.

    Parameters
    ----------
    H_inv : ndarray (3, 3) — inverse homography (ref → tgt)
    output_shape : tuple (H, W)

    Returns
    -------
    coords_y : ndarray (H, W)
    coords_x : ndarray (H, W)
    """
    h, w = output_shape
    yy, xx = np.mgrid[0:h, 0:w]

    # Apply H_inv: [x_tgt, y_tgt, w] = H_inv @ [x_ref, y_ref, 1]
    x_tgt = H_inv[0, 0] * xx + H_inv[0, 1] * yy + H_inv[0, 2]
    y_tgt = H_inv[1, 0] * xx + H_inv[1, 1] * yy + H_inv[1, 2]
    w_tgt = H_inv[2, 0] * xx + H_inv[2, 1] * yy + H_inv[2, 2]

    # Normalize
    w_tgt[w_tgt == 0] = 1e-10
    x_tgt /= w_tgt
    y_tgt /= w_tgt

    return y_tgt, x_tgt


def _build_coordinate_map_similarity(M_inv, output_shape):
    """Build inverse coordinate map using similarity transform.

    Parameters
    ----------
    M_inv : ndarray (2, 3) — inverse similarity (ref → tgt)
    output_shape : tuple (H, W)

    Returns
    -------
    coords_y : ndarray (H, W)
    coords_x : ndarray (H, W)
    """
    h, w = output_shape
    yy, xx = np.mgrid[0:h, 0:w]

    x_tgt = M_inv[0, 0] * xx + M_inv[0, 1] * yy + M_inv[0, 2]
    y_tgt = M_inv[1, 0] * xx + M_inv[1, 1] * yy + M_inv[1, 2]

    return y_tgt, x_tgt


def build_inverse_tps(pts_ref, pts_tgt, smoothing=0.25):
    """Build inverse TPS: reference → target (for resampling).

    The forward TPS (tgt→ref) is used for matching. For resampling we need
    the inverse direction (ref→tgt): "for this output pixel in reference
    frame, where is the corresponding pixel in the target frame?"

    Parameters
    ----------
    pts_ref : ndarray (N, 2) — reference positions
    pts_tgt : ndarray (N, 2) — target positions
    smoothing : float

    Returns
    -------
    inv_tps_x, inv_tps_y : RBFInterpolator — maps ref coords → target coords
    """
    from scipy.interpolate import RBFInterpolator

    pts_ref = np.asarray(pts_ref, dtype=np.float64)
    pts_tgt = np.asarray(pts_tgt, dtype=np.float64)

    inv_tps_x = RBFInterpolator(pts_ref, pts_tgt[:, 0],
                                kernel='thin_plate_spline',
                                smoothing=smoothing)
    inv_tps_y = RBFInterpolator(pts_ref, pts_tgt[:, 1],
                                kernel='thin_plate_spline',
                                smoothing=smoothing)
    return inv_tps_x, inv_tps_y


def resample_image(image, coord_map=None,
                   inv_tps=None, inv_transform=None, model=None,
                   output_shape=None, clamping=True, fill_value=0.0):
    """Resample image using precomputed coordinate map or transform.

    Provide ONE of:
    - coord_map: precomputed (coords_y, coords_x) arrays
    - inv_tps: tuple (tps_x, tps_y) for inverse TPS
    - inv_transform: ndarray for inverse projective/similarity + model

    Parameters
    ----------
    image : ndarray (H, W) — single-channel image
    coord_map : tuple (coords_y, coords_x) — precomputed coordinate maps
    inv_tps : tuple (tps_x, tps_y) — inverse TPS model
    inv_transform : ndarray — inverse transform matrix
    model : str — 'projective' or 'similarity' (for inv_transform)
    output_shape : tuple (H, W) — output dimensions (required if no coord_map)
    clamping : bool — clamp negative values to zero (prevents ringing)
    fill_value : float — value for out-of-bounds pixels

    Returns
    -------
    ndarray (H_out, W_out) — resampled image in float64
    """
    image = np.asarray(image, dtype=np.float64)

    if coord_map is not None:
        coords_y, coords_x = coord_map
    elif inv_tps is not None:
        if output_shape is None:
            output_shape = image.shape
        coords_y, coords_x = _build_coordinate_map_tps(
            inv_tps[0], inv_tps[1], output_shape)
    elif inv_transform is not None:
        if output_shape is None:
            output_shape = image.shape
        if model == 'projective':
            coords_y, coords_x = _build_coordinate_map_projective(
                inv_transform, output_shape)
        elif model == 'similarity':
            coords_y, coords_x = _build_coordinate_map_similarity(
                inv_transform, output_shape)
        else:
            raise ValueError(f"Unknown model: {model}")
    else:
        raise ValueError("Provide coord_map, inv_tps, or inv_transform")

    # Build out-of-bounds mask before interpolation
    h_in, w_in = image.shape
    oob = ((coords_x < 0) | (coords_x >= w_in - 1) |
           (coords_y < 0) | (coords_y >= h_in - 1))

    # scipy.ndimage.map_coordinates expects (ndim, ...) coordinate array
    # order: row (y), col (x)
    coordinates = np.array([coords_y, coords_x])

    result = map_coordinates(image, coordinates, order=3,
                             mode='constant', cval=fill_value)

    # Clamping: prevent ringing artifacts (negative values from bicubic)
    if clamping:
        result = np.maximum(result, 0.0)

    # Fill out-of-bounds with fill_value
    result[oob] = fill_value

    return result


def resample_frame(data, pairs_ref, pairs_tgt, output_shape,
                   smoothing=0.25, clamping=True, fill_value=0.0):
    """High-level: resample a full frame (2D or 3D) using matched star pairs.

    Builds inverse TPS, computes coordinate map once (via coarse grid),
    applies to all channels.

    Parameters
    ----------
    data : ndarray — (H, W) for mono or (C, H, W) for RGB
    pairs_ref : ndarray (N, 2) — reference star positions (x, y)
    pairs_tgt : ndarray (N, 2) — target star positions (x, y)
    output_shape : tuple (H, W) — reference frame dimensions
    smoothing : TPS smoothness
    clamping : bool
    fill_value : float

    Returns
    -------
    ndarray — resampled data, same ndim as input, shape matches output_shape
    """
    # Build inverse TPS (ref → tgt) for resampling
    inv_tps_x, inv_tps_y = build_inverse_tps(pairs_ref, pairs_tgt, smoothing)

    # Build coordinate map once (coarse grid + upscale)
    coords_y, coords_x = _build_coordinate_map_tps(
        inv_tps_x, inv_tps_y, output_shape)
    coord_map = (coords_y, coords_x)

    if data.ndim == 2:
        return resample_image(data, coord_map=coord_map,
                              clamping=clamping, fill_value=fill_value)
    elif data.ndim == 3:
        n_channels = data.shape[0]
        result = np.empty((n_channels,) + output_shape, dtype=np.float64)
        for ch in range(n_channels):
            result[ch] = resample_image(data[ch], coord_map=coord_map,
                                        clamping=clamping, fill_value=fill_value)
        return result
    else:
        raise ValueError(f"Expected 2D or 3D data, got {data.ndim}D")
