#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Thin Plate Spline registration for star alignment.

Iterative TPS refinement: starts from a projective/similarity transform,
discovers new star pairs by proximity, refits TPS with growing inlier set.

Uses scipy.interpolate.RBFInterpolator with thin_plate_spline kernel
(approximating mode with smoothing parameter).

Pure functions, no file I/O.

References:
- Bookstein 1989 — TPS for image registration
- PixInsight StarAlignment — successive approximation TPS
"""

import numpy as np
from scipy.interpolate import RBFInterpolator
from scipy.spatial import cKDTree


def thin_control_points(pts_ref, pts_tgt, residuals, image_shape,
                        max_points=500, grid=(3, 3)):
    """Thin control points to limit TPS complexity while keeping spatial coverage.

    Divides the frame into a grid of zones, selects best points (lowest residual)
    per zone, ensuring uniform spatial coverage.

    Parameters
    ----------
    pts_ref : ndarray (N, 2) — reference positions
    pts_tgt : ndarray (N, 2) — target positions
    residuals : ndarray (N,) — fit residuals per pair
    image_shape : tuple (H, W)
    max_points : maximum total points to keep
    grid : tuple (rows, cols) — zone grid

    Returns
    -------
    indices : ndarray — indices of selected points
    """
    n = len(pts_ref)
    if n <= max_points:
        return np.arange(n)

    rows, cols = grid
    h, w = image_shape
    zone_h = h / rows
    zone_w = w / cols

    per_zone = max_points // (rows * cols)
    # Reserve some for the best global points
    per_zone = max(per_zone, 3)

    selected = set()

    for r in range(rows):
        for c in range(cols):
            # Zone boundaries
            y_lo = r * zone_h
            y_hi = (r + 1) * zone_h
            x_lo = c * zone_w
            x_hi = (c + 1) * zone_w

            # Find points in this zone (use reference positions)
            mask = ((pts_ref[:, 0] >= x_lo) & (pts_ref[:, 0] < x_hi) &
                    (pts_ref[:, 1] >= y_lo) & (pts_ref[:, 1] < y_hi))
            zone_idx = np.where(mask)[0]

            if len(zone_idx) == 0:
                continue

            # Sort by residual (best first) and take per_zone
            order = np.argsort(residuals[zone_idx])
            keep = zone_idx[order[:per_zone]]
            selected.update(keep.tolist())

    # If we have room, fill with globally best remaining
    if len(selected) < max_points:
        remaining = np.array([i for i in range(n) if i not in selected])
        if len(remaining) > 0:
            order = np.argsort(residuals[remaining])
            extra = max_points - len(selected)
            selected.update(remaining[order[:extra]].tolist())

    return np.array(sorted(selected), dtype=np.int32)


def build_tps(pts_ref, pts_tgt, smoothing=0.25):
    """Build approximating TPS model from control point pairs.

    Parameters
    ----------
    pts_ref : ndarray (N, 2) — reference (destination) positions
    pts_tgt : ndarray (N, 2) — target (source) positions
    smoothing : float — smoothness parameter (0 = interpolating, higher = smoother)

    Returns
    -------
    tps_x : RBFInterpolator — maps target coords → reference x
    tps_y : RBFInterpolator — maps target coords → reference y
    """
    pts_tgt = np.asarray(pts_tgt, dtype=np.float64)
    pts_ref = np.asarray(pts_ref, dtype=np.float64)

    tps_x = RBFInterpolator(pts_tgt, pts_ref[:, 0],
                            kernel='thin_plate_spline',
                            smoothing=smoothing)
    tps_y = RBFInterpolator(pts_tgt, pts_ref[:, 1],
                            kernel='thin_plate_spline',
                            smoothing=smoothing)
    return tps_x, tps_y


def apply_tps(tps_x, tps_y, pts):
    """Apply TPS model to transform points.

    Parameters
    ----------
    tps_x, tps_y : RBFInterpolator
    pts : ndarray (N, 2) — points to transform

    Returns
    -------
    ndarray (N, 2) — transformed points
    """
    pts = np.asarray(pts, dtype=np.float64)
    x_out = tps_x(pts)
    y_out = tps_y(pts)
    return np.column_stack([x_out, y_out])


def _apply_initial_transform(pts, transform, model):
    """Apply projective or similarity transform to points.

    Parameters
    ----------
    pts : ndarray (N, 2)
    transform : ndarray — (2, 3) for similarity, (3, 3) for projective
    model : str — 'similarity' or 'projective'

    Returns
    -------
    ndarray (N, 2)
    """
    pts = np.asarray(pts, dtype=np.float64)
    ones = np.ones((len(pts), 1), dtype=np.float64)
    pts_h = np.hstack([pts, ones])

    if model == 'similarity':
        return pts_h @ transform.T
    else:
        result_h = pts_h @ transform.T
        result_h[:, 0] /= result_h[:, 2]
        result_h[:, 1] /= result_h[:, 2]
        return result_h[:, :2]


def tps_iterative(xy_ref, xy_tgt, initial_transform, kdtree_ref=None,
                  model='projective', image_shape=None,
                  smoothing=0.25, max_iterations=10,
                  convergence_threshold=0.01, quality_threshold=1.0,
                  search_radius=10.0, max_control_points=500):
    """Iterative TPS refinement via successive approximations.

    Starting from an initial projective/similarity transform (from RANSAC),
    iteratively discovers new star pairs and refines the TPS model.

    Algorithm:
    1. Apply current model to predict target → reference positions
    2. Match predicted positions to reference stars (nearest neighbor)
    3. Filter by distance threshold
    4. Refit TPS on the new (expanded) set of pairs
    5. Repeat until convergence

    Parameters
    ----------
    xy_ref : ndarray (N1, 2) — all reference star positions
    xy_tgt : ndarray (N2, 2) — all target star positions
    initial_transform : ndarray — from RANSAC
    kdtree_ref : cKDTree (optional, built if None)
    model : 'similarity' or 'projective' — type of initial_transform
    image_shape : tuple (H, W) — for control point thinning
    smoothing : TPS smoothness
    max_iterations : iteration limit
    convergence_threshold : min improvement in median residual (pixels)
    quality_threshold : max median residual for pass (pixels)
    search_radius : initial search radius for pair discovery (pixels)
    max_control_points : max control points for TPS

    Returns
    -------
    dict with keys:
        'tps' : tuple (tps_x, tps_y) or None
        'pairs' : ndarray (K, 2) — final matched pairs (ref_idx, tgt_idx)
        'residuals' : ndarray (K,) — final residuals
        'median_residual' : float
        'n_iterations' : int
        'converged' : bool
        'quality_ok' : bool — median_residual <= quality_threshold
    """
    xy_ref = np.asarray(xy_ref, dtype=np.float64)
    xy_tgt = np.asarray(xy_tgt, dtype=np.float64)

    if kdtree_ref is None:
        kdtree_ref = cKDTree(xy_ref)

    result = {
        'tps': None,
        'pairs': None,
        'residuals': None,
        'median_residual': float('inf'),
        'n_iterations': 0,
        'converged': False,
        'quality_ok': False,
    }

    # Step 1: apply initial transform to all target stars
    predicted = _apply_initial_transform(xy_tgt, initial_transform, model)

    prev_median = float('inf')
    current_tps = None

    for iteration in range(max_iterations):
        # Find nearest reference star for each predicted target position
        dists, indices = kdtree_ref.query(predicted)

        # Filter by search radius
        good = dists < search_radius
        if np.sum(good) < 6:
            break

        tgt_idx = np.where(good)[0]
        ref_idx = indices[good]

        # Remove duplicates: if multiple targets match the same reference,
        # keep the one with smallest distance
        seen = {}
        for i in range(len(ref_idx)):
            ri = int(ref_idx[i])
            ti = int(tgt_idx[i])
            d = dists[good][i]
            if ri not in seen or d < seen[ri][1]:
                seen[ri] = (ti, d)

        ref_matched = np.array(list(seen.keys()), dtype=np.int32)
        tgt_matched = np.array([seen[r][0] for r in ref_matched], dtype=np.int32)

        if len(ref_matched) < 6:
            break

        pts_ref = xy_ref[ref_matched]
        pts_tgt = xy_tgt[tgt_matched]

        # Thin control points if too many
        # Compute residuals for thinning priority
        residuals = np.sqrt(np.sum((predicted[tgt_matched] - pts_ref) ** 2, axis=1))

        if image_shape is not None and len(ref_matched) > max_control_points:
            keep = thin_control_points(pts_ref, pts_tgt, residuals,
                                       image_shape, max_control_points)
            pts_ref = pts_ref[keep]
            pts_tgt = pts_tgt[keep]
            ref_matched = ref_matched[keep]
            tgt_matched = tgt_matched[keep]
            residuals = residuals[keep]

        # Build TPS
        try:
            tps_x, tps_y = build_tps(pts_ref, pts_tgt, smoothing)
        except Exception:
            break

        current_tps = (tps_x, tps_y)

        # Recompute predictions using TPS
        predicted = apply_tps(tps_x, tps_y, xy_tgt)

        # Compute residuals on matched pairs
        new_residuals = np.sqrt(
            np.sum((predicted[tgt_matched] - pts_ref) ** 2, axis=1))
        median_res = float(np.median(new_residuals))

        result['tps'] = current_tps
        result['pairs'] = np.column_stack([ref_matched, tgt_matched])
        result['residuals'] = new_residuals
        result['median_residual'] = median_res
        result['n_iterations'] = iteration + 1

        # Check convergence
        improvement = prev_median - median_res
        if improvement < convergence_threshold and iteration > 0:
            result['converged'] = True
            break

        prev_median = median_res

        # Tighten search radius as model improves
        search_radius = max(3.0 * median_res, 3.0)

    result['quality_ok'] = result['median_residual'] <= quality_threshold
    return result
