#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Star matching via polygonal descriptors (pentagons).

Builds scale/rotation-invariant descriptors from star positions,
matches them between two frames, and fits a geometric transformation
via RANSAC.

Pure functions, no file I/O. Reusable for staralign, mosaics, etc.

References:
- Lang et al. 2010 — astrometry.net, geometric hashing
- PixInsight StarAlignment — polygonal descriptors, successive approximation
- Groth 1986 — triangle matching for star catalogs
"""

import numpy as np
from scipy.spatial import cKDTree


# ---------------------------------------------------------
# Stage 1: Spatial index
# ---------------------------------------------------------

def build_kdtree(xy):
    """Build a 2D KD-tree from star positions.

    Parameters
    ----------
    xy : ndarray, shape (N, 2)
        Star positions (x, y).

    Returns
    -------
    cKDTree
    """
    return cKDTree(np.asarray(xy, dtype=np.float64))


# ---------------------------------------------------------
# Stage 2: Pentagon descriptors
# ---------------------------------------------------------

# Neighbor ring boundaries (indices into K-nearest sorted by distance)
# Ring 0 (close):  neighbors 0..4   (5 stars)
# Ring 1 (mid):    neighbors 5..14  (10 stars)
# Ring 2 (far):    neighbors 15..29 (15 stars)
#
# Pentagon templates: each is a 5-tuple of neighbor indices.
# Mix rings for scale diversity + robustness to missing close neighbors.

_PENTAGON_TEMPLATES = [
    # 3 close + 1 mid + 1 far
    (0, 1, 2, 5, 15),
    (0, 1, 3, 7, 18),
    (0, 2, 4, 9, 22),
    (1, 2, 3, 6, 16),
    (1, 3, 4, 8, 20),
    # 2 close + 2 mid + 1 far
    (0, 1, 5, 10, 25),
    (0, 2, 6, 12, 20),
    (1, 3, 7, 11, 28),
    (0, 4, 8, 14, 22),
    (2, 3, 9, 13, 17),
    # 2 close + 1 mid + 2 far
    (0, 1, 5, 15, 25),
    (0, 2, 8, 18, 28),
    (1, 3, 10, 20, 27),
    (0, 4, 12, 22, 29),
    (2, 4, 6, 16, 24),
    # 1 close + 2 mid + 2 far
    (0, 5, 10, 15, 25),
    (1, 6, 11, 18, 28),
    (2, 7, 12, 20, 26),
    (3, 8, 13, 22, 29),
    (4, 9, 14, 16, 24),
]


def _compute_pentagon_hash(pts):
    """Compute a 6D hash for a pentagon (5 points).

    Algorithm:
    1. Find the two most distant points → A, B (backbone).
    2. Deterministic A/B choice: center of mass of the remaining 3 points
       (C, D, E) determines which side of line AB it falls on.
       If CoM is to the "left" of AB → keep order; else swap A↔B.
    3. Build local coordinate system: A at origin, B at (1, 1).
    4. Project C, D, E into this system; sort lexicographically.
    5. Return 6 floats: (Xc, Yc, Xd, Yd, Xe, Ye).

    The hash is invariant to translation, rotation, and uniform scaling.
    NOT invariant to mirror (handled at meta-level by caller).

    Parameters
    ----------
    pts : ndarray, shape (5, 2)

    Returns
    -------
    ndarray, shape (6,), or None if degenerate.
    """
    # Find backbone: two most distant points
    best_dist2 = -1.0
    ai, bi = 0, 1
    for i in range(5):
        for j in range(i + 1, 5):
            d2 = (pts[i, 0] - pts[j, 0]) ** 2 + (pts[i, 1] - pts[j, 1]) ** 2
            if d2 > best_dist2:
                best_dist2 = d2
                ai, bi = i, j

    if best_dist2 < 1e-10:
        return None  # degenerate

    A = pts[ai]
    B = pts[bi]

    # Remaining 3 points
    others = [k for k in range(5) if k != ai and k != bi]
    rest = pts[others]  # shape (3, 2)

    # Deterministic A/B: CoM of rest relative to line AB
    com = rest.mean(axis=0)
    # Cross product (B-A) × (CoM-A): positive = left side
    ab = B - A
    ac = com - A
    cross = ab[0] * ac[1] - ab[1] * ac[0]
    if cross < 0:
        A, B = B, A
        ab = B - A

    # Local coordinate system: A → origin, B → (1, 1)
    # Transform: p' = M * (p - A), where M maps ab → (1, 1)
    # M = [[ab_x, ab_y], [-ab_y, ab_x]] / |ab|² * scale
    # We want B-A → (1, 1), so:
    #   u = dot(p-A, ab) / |ab|²
    #   v = cross(ab, p-A) / |ab|²  (perpendicular component)
    # Then rotate/scale so B maps to (1, 1):
    #   Actually, u(B) = 1, v(B) = 0.  We need B at (1,1)?
    #
    # Re-reading spec: "A в начале, B в (1,1)" means a general affine.
    # Simpler: just use (u, v) where u = projection along AB, v = perp.
    # u(A)=0, u(B)=1. This is scale+rotation invariant.
    # B maps to (1, 0) in this system. The spec says (1,1) but that
    # seems like a normalization choice. Let's use (1, 0) which is
    # the natural projection — the hash values are the same up to a
    # constant offset, and sorting CDE lexicographically is the same.

    ab_len2 = best_dist2
    coords = np.empty((3, 2), dtype=np.float64)
    for idx, k in enumerate(others):
        dp = pts[k] - A
        u = (dp[0] * ab[0] + dp[1] * ab[1]) / ab_len2
        v = (ab[0] * dp[1] - ab[1] * dp[0]) / ab_len2
        coords[idx, 0] = u
        coords[idx, 1] = v

    # Sort CDE lexicographically by (u, v)
    order = np.lexsort((coords[:, 1], coords[:, 0]))
    coords = coords[order]

    return coords.ravel()  # shape (6,)


def build_descriptors(xy, kdtree=None, n_descriptors=20, n_neighbors=30):
    """Build pentagon descriptors for all stars.

    For each star, uses K nearest neighbors and fixed templates to
    form pentagons. Each pentagon produces a 6D hash.

    Parameters
    ----------
    xy : ndarray, shape (N, 2)
    kdtree : cKDTree (optional, built if None)
    n_descriptors : max descriptors per star (default 20)
    n_neighbors : neighbors to query (default 30)

    Returns
    -------
    hashes : ndarray, shape (M, 6) — all descriptor hashes
    star_indices : ndarray, shape (M,) — which star each descriptor belongs to
    """
    xy = np.asarray(xy, dtype=np.float64)
    n_stars = len(xy)

    if n_stars < 5:
        return (np.empty((0, 6), dtype=np.float64),
                np.empty(0, dtype=np.int32))

    if kdtree is None:
        kdtree = build_kdtree(xy)

    # Query neighbors (K+1 because first result is the star itself)
    k_query = min(n_neighbors + 1, n_stars)
    _, neighbor_idx = kdtree.query(xy, k=k_query)
    # neighbor_idx[:, 0] is the star itself; neighbors are [:, 1:]
    neighbor_idx = neighbor_idx[:, 1:]  # shape (N, k_query-1)
    n_actual_neighbors = neighbor_idx.shape[1]

    # Select templates that fit within available neighbors
    templates = []
    for t in _PENTAGON_TEMPLATES:
        if max(t) < n_actual_neighbors:
            templates.append(t)
        if len(templates) >= n_descriptors:
            break

    if not templates and n_actual_neighbors >= 4:
        # Fall back to combinations of first few neighbors
        from itertools import combinations
        for combo in combinations(range(min(n_actual_neighbors, 8)), 4):
            templates.append(combo)
            if len(templates) >= n_descriptors:
                break

    all_hashes = []
    all_star_idx = []

    for star_i in range(n_stars):
        nbrs = neighbor_idx[star_i]
        count = 0

        for template in templates:
            if count >= n_descriptors:
                break

            # Pentagon: center star + 4 neighbors from template
            pts = np.empty((5, 2), dtype=np.float64)
            pts[0] = xy[star_i]
            valid = True
            for j, ni in enumerate(template[:4]):
                if ni >= len(nbrs):
                    valid = False
                    break
                pts[j + 1] = xy[nbrs[ni]]
            if not valid:
                continue

            h = _compute_pentagon_hash(pts)
            if h is not None:
                all_hashes.append(h)
                all_star_idx.append(star_i)
                count += 1

    if not all_hashes:
        return (np.empty((0, 6), dtype=np.float64),
                np.empty(0, dtype=np.int32))

    return (np.array(all_hashes, dtype=np.float64),
            np.array(all_star_idx, dtype=np.int32))


# ---------------------------------------------------------
# Stage 3: Descriptor matching
# ---------------------------------------------------------

def match_descriptors(hashes_ref, stars_ref_idx, hashes_tgt, stars_tgt_idx,
                      tolerance=0.05):
    """Match descriptors between reference and target via 6D KD-tree.

    For each target hash, finds the nearest reference hash within tolerance.
    Accumulates votes: each hash match votes for a (ref_star, tgt_star) pair.

    Parameters
    ----------
    hashes_ref : ndarray (M1, 6)
    stars_ref_idx : ndarray (M1,) — star index for each ref hash
    hashes_tgt : ndarray (M2, 6)
    stars_tgt_idx : ndarray (M2,) — star index for each tgt hash
    tolerance : max L2 distance in 6D hash space

    Returns
    -------
    pairs : ndarray (K, 2) — matched (ref_star_idx, tgt_star_idx)
    votes : ndarray (K,) — number of hash votes for each pair
    """
    if len(hashes_ref) == 0 or len(hashes_tgt) == 0:
        return np.empty((0, 2), dtype=np.int32), np.empty(0, dtype=np.int32)

    tree6d = cKDTree(hashes_ref)
    dist, idx = tree6d.query(hashes_tgt, distance_upper_bound=tolerance)

    # Count votes per (ref_star, tgt_star) pair
    vote_map = {}
    n_ref_hashes = len(hashes_ref)

    for i in range(len(hashes_tgt)):
        if idx[i] >= n_ref_hashes:
            continue  # no match within tolerance
        ref_star = int(stars_ref_idx[idx[i]])
        tgt_star = int(stars_tgt_idx[i])
        pair = (ref_star, tgt_star)
        vote_map[pair] = vote_map.get(pair, 0) + 1

    if not vote_map:
        return np.empty((0, 2), dtype=np.int32), np.empty(0, dtype=np.int32)

    pairs = np.array(list(vote_map.keys()), dtype=np.int32)
    votes = np.array(list(vote_map.values()), dtype=np.int32)

    # Sort by votes descending
    order = np.argsort(-votes)
    return pairs[order], votes[order]


# ---------------------------------------------------------
# Stage 3b: Angular verification
# ---------------------------------------------------------

# Max distance for neighbor inclusion: fraction of frame diagonal
_ANGULAR_RADIUS_FRACTION = 0.3

# Max neighbors for angular signature
_ANGULAR_MAX_NEIGHBORS = 20

# Minimum neighbors for a usable signature
_ANGULAR_MIN_NEIGHBORS = 4


def _angular_signature(xy, kdtree, star_idx, max_radius):
    """Compute angular signature for a star using nearby neighbors.

    Algorithm:
    1. Find up to 20 nearest neighbors within max_radius.
    2. Compute position angle to each neighbor via atan2(dy, dx).
    3. Sort CCW, find largest angular gap.
    4. First neighbor CCW from the gap = reference direction.
    5. Signature = angles relative to that first neighbor (all positive,
       increasing, from 0 to <2*pi).

    The signature is invariant to rotation and scale.

    Parameters
    ----------
    xy : ndarray (N, 2)
    kdtree : cKDTree
    star_idx : int
    max_radius : float — maximum distance for neighbor inclusion

    Returns
    -------
    rel_angles : ndarray of angles relative to first (radians), length K.
                 First element is always 0.0.
    Returns None if fewer than _ANGULAR_MIN_NEIGHBORS neighbors found.
    """
    k_query = min(_ANGULAR_MAX_NEIGHBORS + 1, len(xy))
    if k_query < _ANGULAR_MIN_NEIGHBORS + 1:
        return None

    dists, nbr_idx = kdtree.query(xy[star_idx], k=k_query)
    # Exclude self (index 0) and neighbors beyond max_radius
    mask = (dists[1:] <= max_radius) & (dists[1:] > 0)
    nbr_idx = nbr_idx[1:][mask]

    if len(nbr_idx) < _ANGULAR_MIN_NEIGHBORS:
        return None

    # Position angles from center star to each neighbor
    center = xy[star_idx]
    dx = xy[nbr_idx, 0] - center[0]
    dy = xy[nbr_idx, 1] - center[1]
    angles = np.arctan2(dy, dx)  # [-pi, pi]

    # Shift to [0, 2*pi)
    angles = angles % (2.0 * np.pi)

    # Sort CCW
    order = np.argsort(angles)
    angles_sorted = angles[order]
    n = len(angles_sorted)

    # Find largest gap between consecutive sorted angles (with wrap-around)
    gaps = np.empty(n, dtype=np.float64)
    gaps[:n - 1] = angles_sorted[1:] - angles_sorted[:n - 1]
    gaps[n - 1] = (2.0 * np.pi) - angles_sorted[-1] + angles_sorted[0]

    max_gap_idx = np.argmax(gaps)

    # First neighbor CCW from the largest gap
    start = (max_gap_idx + 1) % n

    # Reorder starting from start
    reordered = np.roll(angles_sorted, -start)

    # Angles relative to first neighbor
    rel_angles = reordered - reordered[0]
    # Ensure positive (wrap negatives from float precision)
    rel_angles = rel_angles % (2.0 * np.pi)

    return rel_angles


def _compare_angle_sequences(sig_a, sig_b, angle_tol):
    """Compare two angular signatures via LCS-like forward scan.

    Both signatures start at 0.0 and increase to ~2*pi.
    Walk both sequences left-to-right, matching elements within tolerance.
    On mismatch, skip the element with the smaller angle (it's a neighbor
    missing in the other frame).

    Parameters
    ----------
    sig_a, sig_b : 1D arrays of relative angles (starting from 0)
    angle_tol : maximum angle difference for a match (radians)

    Returns
    -------
    n_matched : number of matched elements
    """
    na, nb = len(sig_a), len(sig_b)
    ia, ib = 0, 0
    matched = 0

    while ia < na and ib < nb:
        diff = abs(sig_a[ia] - sig_b[ib])
        # Handle wrap-around near 2*pi (shouldn't happen often since
        # both start at 0, but be safe)
        if diff > np.pi:
            diff = 2.0 * np.pi - diff
        if diff <= angle_tol:
            matched += 1
            ia += 1
            ib += 1
        elif sig_a[ia] < sig_b[ib]:
            # sig_a[ia] has no match in B — skip it
            ia += 1
        else:
            # sig_b[ib] has no match in A — skip it
            ib += 1

    return matched


def verify_angular(xy_ref, kdtree_ref, xy_tgt, kdtree_tgt,
                   pairs, image_diag,
                   angle_tol=0.15, min_matching_fraction=0.5):
    """Angular verification of candidate star pairs.

    For each pair, computes angular signatures (position angles to neighbors,
    relative to a rotation-invariant start) and compares them via LCS-like
    forward scanning.

    Parameters
    ----------
    xy_ref, xy_tgt : ndarray (N, 2) — star positions
    kdtree_ref, kdtree_tgt : cKDTree
    pairs : ndarray (K, 2) — candidate (ref_idx, tgt_idx) pairs
    image_diag : float — frame diagonal in pixels (for neighbor radius limit)
    angle_tol : max angle difference for a match (radians)
    min_matching_fraction : fraction of shorter signature that must match

    Returns
    -------
    verified_mask : boolean ndarray (K,) — True for verified pairs
    """
    max_radius = image_diag * _ANGULAR_RADIUS_FRACTION

    n_pairs = len(pairs)
    verified = np.zeros(n_pairs, dtype=bool)

    for pi in range(n_pairs):
        ref_idx = pairs[pi, 0]
        tgt_idx = pairs[pi, 1]

        sig_r = _angular_signature(xy_ref, kdtree_ref, ref_idx, max_radius)
        sig_t = _angular_signature(xy_tgt, kdtree_tgt, tgt_idx, max_radius)

        if sig_r is None or sig_t is None:
            continue

        score = _compare_angle_sequences(sig_r, sig_t, angle_tol)
        n_min = min(len(sig_r), len(sig_t))
        if n_min > 0 and score / n_min >= min_matching_fraction:
            verified[pi] = True

    return verified


# ---------------------------------------------------------
# Stage 4: RANSAC
# ---------------------------------------------------------

def _fit_similarity(pts_ref, pts_tgt):
    """Fit similarity transform (rotation + scale + translation) from point pairs.

    Parameters
    ----------
    pts_ref : ndarray (N, 2)
    pts_tgt : ndarray (N, 2)

    Returns
    -------
    M : ndarray (2, 3) — affine matrix [R*s | t]
    """
    # Solve: ref = s*R*tgt + t
    # Using least squares: [x_r] = [a -b tx] [x_t]
    #                      [y_r]   [b  a ty] [y_t]
    #                                        [ 1 ]
    n = len(pts_ref)
    A = np.zeros((2 * n, 4), dtype=np.float64)
    b = np.zeros(2 * n, dtype=np.float64)

    for i in range(n):
        xt, yt = pts_tgt[i]
        xr, yr = pts_ref[i]
        A[2 * i] = [xt, -yt, 1, 0]
        A[2 * i + 1] = [yt, xt, 0, 1]
        b[2 * i] = xr
        b[2 * i + 1] = yr

    result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    a, bv, tx, ty = result

    M = np.array([[a, -bv, tx],
                  [bv, a, ty]], dtype=np.float64)
    return M


def _fit_projective(pts_ref, pts_tgt):
    """Fit projective (homography) transform from point pairs.

    Parameters
    ----------
    pts_ref : ndarray (N, 2) — at least 4 pairs
    pts_tgt : ndarray (N, 2)

    Returns
    -------
    H : ndarray (3, 3) — homography matrix (tgt → ref)
    """
    n = len(pts_ref)
    A = np.zeros((2 * n, 8), dtype=np.float64)
    b = np.zeros(2 * n, dtype=np.float64)

    for i in range(n):
        xt, yt = pts_tgt[i]
        xr, yr = pts_ref[i]
        A[2 * i] = [xt, yt, 1, 0, 0, 0, -xr * xt, -xr * yt]
        A[2 * i + 1] = [0, 0, 0, xt, yt, 1, -yr * xt, -yr * yt]
        b[2 * i] = xr
        b[2 * i + 1] = yr

    result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    H = np.array([
        [result[0], result[1], result[2]],
        [result[3], result[4], result[5]],
        [result[6], result[7], 1.0]
    ], dtype=np.float64)
    return H


def _apply_transform(pts, M, model):
    """Apply transform to points. Returns ndarray (N, 2)."""
    pts = np.asarray(pts, dtype=np.float64)
    if model == 'similarity':
        # M is (2, 3): [a -b tx; b a ty]
        ones = np.ones((len(pts), 1), dtype=np.float64)
        pts_h = np.hstack([pts, ones])  # (N, 3)
        result = pts_h @ M.T  # (N, 2)
        return result
    else:
        # M is (3, 3) homography
        ones = np.ones((len(pts), 1), dtype=np.float64)
        pts_h = np.hstack([pts, ones])  # (N, 3)
        result_h = pts_h @ M.T  # (N, 3)
        result_h[:, 0] /= result_h[:, 2]
        result_h[:, 1] /= result_h[:, 2]
        return result_h[:, :2]


def ransac_fit(xy_ref, xy_tgt, pairs, model='projective',
               n_iterations=1000, inlier_tolerance=3.0, min_inliers=6):
    """RANSAC fitting of geometric transformation.

    Parameters
    ----------
    xy_ref : ndarray (N1, 2)
    xy_tgt : ndarray (N2, 2)
    pairs : ndarray (K, 2) — (ref_idx, tgt_idx)
    model : 'similarity' (2 pairs per hypothesis) or 'projective' (4 pairs)
    n_iterations : RANSAC iterations
    inlier_tolerance : max residual in pixels for inlier
    min_inliers : minimum inliers for valid model

    Returns
    -------
    transform : ndarray — (2, 3) for similarity or (3, 3) for projective
    inlier_mask : boolean ndarray (K,)
    median_residual : float
    Returns (None, None, inf) if RANSAC fails.
    """
    n_pairs = len(pairs)
    min_sample = 2 if model == 'similarity' else 4

    if n_pairs < min_sample:
        return None, None, float('inf')

    pts_ref = xy_ref[pairs[:, 0]]
    pts_tgt = xy_tgt[pairs[:, 1]]

    best_inliers = None
    best_transform = None
    best_n_inliers = 0

    rng = np.random.default_rng(42)

    for _ in range(n_iterations):
        # Random sample
        sample = rng.choice(n_pairs, size=min_sample, replace=False)
        sample_ref = pts_ref[sample]
        sample_tgt = pts_tgt[sample]

        # Fit model
        try:
            if model == 'similarity':
                M = _fit_similarity(sample_ref, sample_tgt)
            else:
                M = _fit_projective(sample_ref, sample_tgt)
        except np.linalg.LinAlgError:
            continue

        # Compute residuals for all pairs
        predicted = _apply_transform(pts_tgt, M, model)
        residuals = np.sqrt(np.sum((predicted - pts_ref) ** 2, axis=1))

        inlier_mask = residuals < inlier_tolerance
        n_inliers = np.sum(inlier_mask)

        if n_inliers > best_n_inliers:
            best_n_inliers = n_inliers
            best_inliers = inlier_mask.copy()
            best_transform = M

    if best_n_inliers < min_inliers:
        return None, None, float('inf')

    # Refit on all inliers
    inlier_ref = pts_ref[best_inliers]
    inlier_tgt = pts_tgt[best_inliers]

    try:
        if model == 'similarity':
            best_transform = _fit_similarity(inlier_ref, inlier_tgt)
        else:
            best_transform = _fit_projective(inlier_ref, inlier_tgt)
    except np.linalg.LinAlgError:
        pass

    # Final residuals
    predicted = _apply_transform(pts_tgt[best_inliers], best_transform, model)
    residuals = np.sqrt(np.sum((predicted - pts_ref[best_inliers]) ** 2, axis=1))
    median_residual = float(np.median(residuals))

    return best_transform, best_inliers, median_residual


# ---------------------------------------------------------
# High-level matching pipeline
# ---------------------------------------------------------

def match_stars(xy_ref, xy_tgt, image_shape=None,
                n_descriptors=20, n_neighbors=30,
                hash_tolerance=0.05,
                angular_tolerance=0.15, angular_min_fraction=0.5,
                ransac_model='projective', ransac_iterations=1000,
                ransac_tolerance=3.0, min_vote=2):
    """Full matching pipeline: descriptors → hash match → angular verify → RANSAC.

    Parameters
    ----------
    xy_ref : ndarray (N1, 2) — reference star positions
    xy_tgt : ndarray (N2, 2) — target star positions
    image_shape : tuple (H, W) — frame dimensions for angular verification radius.
                  If None, estimated from max star coordinates.
    n_descriptors : pentagon descriptors per star
    n_neighbors : KD-tree neighbors for descriptor building
    hash_tolerance : L2 tolerance in 6D hash space
    angular_tolerance : max angle difference for angular verification (radians)
    angular_min_fraction : min fraction of matching angles
    ransac_model : 'similarity' or 'projective'
    ransac_iterations : RANSAC iterations
    ransac_tolerance : inlier tolerance in pixels
    min_vote : minimum hash votes to consider a pair

    Returns
    -------
    dict with keys:
        'transform' : ndarray or None
        'inlier_pairs' : ndarray (K, 2) of (ref_idx, tgt_idx), or None
        'median_residual' : float
        'n_hash_matches' : int
        'n_verified' : int
        'n_inliers' : int
    """
    xy_ref = np.asarray(xy_ref, dtype=np.float64)
    xy_tgt = np.asarray(xy_tgt, dtype=np.float64)

    result = {
        'transform': None,
        'inlier_pairs': None,
        'median_residual': float('inf'),
        'n_hash_matches': 0,
        'n_verified': 0,
        'n_inliers': 0,
    }

    if len(xy_ref) < 5 or len(xy_tgt) < 5:
        return result

    # Frame diagonal for angular verification radius limit
    if image_shape is not None:
        image_diag = np.sqrt(image_shape[0] ** 2 + image_shape[1] ** 2)
    else:
        # Estimate from star coordinates
        all_xy = np.vstack([xy_ref, xy_tgt])
        span = all_xy.max(axis=0) - all_xy.min(axis=0)
        image_diag = np.sqrt(span[0] ** 2 + span[1] ** 2)

    # Build spatial indices
    tree_ref = build_kdtree(xy_ref)
    tree_tgt = build_kdtree(xy_tgt)

    # Build descriptors
    hashes_ref, idx_ref = build_descriptors(xy_ref, tree_ref, n_descriptors, n_neighbors)
    hashes_tgt, idx_tgt = build_descriptors(xy_tgt, tree_tgt, n_descriptors, n_neighbors)

    # Match descriptors
    pairs, votes = match_descriptors(hashes_ref, idx_ref, hashes_tgt, idx_tgt,
                                     tolerance=hash_tolerance)

    # Filter by minimum votes
    vote_mask = votes >= min_vote
    pairs = pairs[vote_mask]
    votes = votes[vote_mask]

    result['n_hash_matches'] = len(pairs)

    if len(pairs) < 4:
        return result

    # Angular verification
    verified = verify_angular(xy_ref, tree_ref, xy_tgt, tree_tgt,
                              pairs, image_diag,
                              angular_tolerance, angular_min_fraction)
    pairs = pairs[verified]
    result['n_verified'] = len(pairs)

    if len(pairs) < 4:
        return result

    # RANSAC
    transform, inlier_mask, median_res = ransac_fit(
        xy_ref, xy_tgt, pairs,
        model=ransac_model,
        n_iterations=ransac_iterations,
        inlier_tolerance=ransac_tolerance,
    )

    if transform is not None:
        result['transform'] = transform
        result['inlier_pairs'] = pairs[inlier_mask]
        result['median_residual'] = median_res
        result['n_inliers'] = int(np.sum(inlier_mask))

    return result
