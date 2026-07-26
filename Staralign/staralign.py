#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
staralign - Star-based image registration for astronomical frames.

Aligns target FITS frames to a reference frame using star matching
with polygonal (pentagon) descriptors, RANSAC, and TPS refinement.

Supports 2D (mono) and 3D (C×H×W, RGB) FITS data.
"""

import sys
import os
import time
import threading
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from astropy.io import fits

# Import shared utilities
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils
import star_utils
import star_match
import star_tps
import star_resample


def usage():
    sys.stderr.write(
        "Usage:\n"
        "  staralign.py input_spec output_spec [--ref reference.fit] [options]\n"
        "\n"
        "  input_spec   - target frames: wildcard (*.fit), numbered (img0001.fit),\n"
        "                 @list.txt, or single file\n"
        "  output_spec  - output files: single, numbered pattern, or wildcard\n"
        "  --ref FILE   - reference frame (default: the first input frame)\n"
        "\n"
        "Options:\n"
        "  --model {tps|projective}  Registration model (default: tps)\n"
        "  --smoothness F            TPS smoothness, 0=exact (default: 0.25)\n"
        "  --descriptors N           Pentagons per star (default: 20, range: 5-50)\n"
        "\n"
        "Matching parameters:\n"
        "  --hash-tol F    Pentagon hash tolerance in 6D space (default: 0.05)\n"
        "                  Increase to 0.1-0.15 for fewer stars or different filters\n"
        "  --angle-tol F   Angular verification tolerance in radians (default: 0.15)\n"
        "                  Increase to 0.25-0.3 for sparser fields\n"
        "  --min-vote N    Min hash votes to accept a star pair (default: 2)\n"
        "                  Set to 1 for difficult cases (few stars, cross-filter)\n"
        "  --tolerance F   RANSAC inlier tolerance in pixels (default: 3.0)\n"
        "                  Increase to 5-10 for large distortion or scale differences\n"
        "\n"
        "Star detection:\n"
        "  --snr F         Detection SNR threshold (default: 5.0, lower=more stars)\n"
        "  --max-stars N   Initial max stars for matching (default: 150)\n"
        "                  On failure, auto-retries with 2x and 3x stars\n"
        "  --no-retry      Disable auto-retry with more stars on failure\n"
        "\n"
        "Other:\n"
        "  --threads N     Worker threads (default: CPU count - 1)\n"
        "                  Forced to 1 in --debug mode\n"
        "  --debug         Print detailed per-frame diagnostics and timing\n"
        "  --no-mirror     Disable mirror fallback attempt\n"
        "  --rgb           Chromatic mode: align R,B -> G channels within EACH input\n"
        "                  file (per-channel), instead of aligning frames to a reference\n"
        "\n"
        "Modes (when --ref is omitted):\n"
        "  Default        - align all input frames to the FIRST frame (star-based)\n"
        "  Single file    - chromatic correction: align R,B -> G channels within it\n"
        "  --rgb          - chromatic correction on every input file (fixes atmospheric\n"
        "                   refraction + chromatic aberration on colour frames)\n"
    )
    sys.exit(1)


def parse_args(argv):
    args = argv[1:]

    # Extract options
    ref_file = None
    model = 'tps'
    smoothness = 0.25
    n_descriptors = 20
    tolerance = 3.0
    hash_tol = 0.05
    angle_tol = 0.15
    min_vote = 2
    snr = 5.0
    max_stars = 150
    n_threads = max(1, (os.cpu_count() or 2) - 1)
    debug = False
    no_mirror = False
    no_retry = False
    rgb = False

    # Parse named options
    filtered = []
    i = 0
    while i < len(args):
        if args[i] == '--ref' and i + 1 < len(args):
            ref_file = args[i + 1]
            i += 2
        elif args[i] == '--model' and i + 1 < len(args):
            model = args[i + 1]
            if model not in ('tps', 'projective'):
                sys.stderr.write(f"Error: --model must be 'tps' or 'projective'\n")
                sys.exit(1)
            i += 2
        elif args[i] == '--smoothness' and i + 1 < len(args):
            smoothness = float(args[i + 1])
            i += 2
        elif args[i] == '--descriptors' and i + 1 < len(args):
            n_descriptors = int(args[i + 1])
            i += 2
        elif args[i] == '--tolerance' and i + 1 < len(args):
            tolerance = float(args[i + 1])
            i += 2
        elif args[i] == '--hash-tol' and i + 1 < len(args):
            hash_tol = float(args[i + 1])
            i += 2
        elif args[i] == '--angle-tol' and i + 1 < len(args):
            angle_tol = float(args[i + 1])
            i += 2
        elif args[i] == '--min-vote' and i + 1 < len(args):
            min_vote = int(args[i + 1])
            i += 2
        elif args[i] == '--snr' and i + 1 < len(args):
            snr = float(args[i + 1])
            i += 2
        elif args[i] == '--max-stars' and i + 1 < len(args):
            max_stars = int(args[i + 1])
            i += 2
        elif args[i] == '--threads' and i + 1 < len(args):
            n_threads = int(args[i + 1])
            i += 2
        elif args[i] == '--debug':
            debug = True
            i += 1
        elif args[i] == '--no-mirror':
            no_mirror = True
            i += 1
        elif args[i] == '--no-retry':
            no_retry = True
            i += 1
        elif args[i] == '--rgb':
            rgb = True
            i += 1
        else:
            filtered.append(args[i])
            i += 1

    if len(filtered) < 2:
        usage()

    input_spec = filtered[0]
    output_spec = filtered[1]

    return {
        'input_spec': input_spec,
        'output_spec': output_spec,
        'ref_file': ref_file,
        'model': model,
        'smoothness': smoothness,
        'n_descriptors': n_descriptors,
        'tolerance': tolerance,
        'hash_tol': hash_tol,
        'angle_tol': angle_tol,
        'min_vote': min_vote,
        'snr': snr,
        'max_stars': max_stars,
        'debug': debug,
        'no_mirror': no_mirror,
        'no_retry': no_retry,
        'rgb': rgb,
        'n_threads': max(1, n_threads),
    }


def detect_and_extract(data, snr, max_stars, debug=False, label=""):
    """Detect stars and return positions as (N, 2) array.

    For 3D data, detects on first channel.
    Returns all stars sorted by flux (up to max_stars).
    """
    if data.ndim == 3:
        detect_data = data[0]
    else:
        detect_data = data

    catalog = star_utils.detect_stars(detect_data, snr=snr)

    if len(catalog) == 0:
        return np.empty((0, 2), dtype=np.float64), 0

    # Sort by flux descending, take top max_stars
    n_detected = len(catalog)
    order = np.argsort(-catalog.flux)
    if len(order) > max_stars:
        order = order[:max_stars]

    xy = np.column_stack([catalog.x[order], catalog.y[order]])

    if debug:
        sys.stderr.write(f"  {label}: {n_detected} stars detected, "
                         f"using top {len(xy)}\n")

    return xy, n_detected


def align_single(target_data, target_header, xy_ref, ref_shape,
                 tree_ref, hashes_ref, idx_ref,
                 opts, debug=False, xy_tgt=None):
    """Align a single target frame to reference.

    Parameters
    ----------
    xy_tgt : ndarray (N, 2), optional
        Pre-detected target star positions. If None, detects internally.

    Returns
    -------
    tuple (aligned_data, info_dict) or (None, info_dict) on failure
    """
    info = {
        'n_stars_tgt': 0,
        'n_hash_matches': 0,
        'n_verified': 0,
        'n_inliers': 0,
        'n_tps_pairs': 0,
        'median_residual': float('inf'),
        'tps_median_residual': float('inf'),
        'mirrored': False,
        'success': False,
        'fail_reason': '',
    }

    output_shape = ref_shape
    timings = {}

    # Detect stars on target (or use pre-detected)
    t0 = time.time()
    if xy_tgt is not None:
        pass  # already provided
    else:
        xy_tgt, _ = detect_and_extract(target_data, opts['snr'], opts['max_stars'],
                                       debug, "target")
    timings['detect'] = time.time() - t0
    info['n_stars_tgt'] = len(xy_tgt)

    if len(xy_tgt) < 5:
        info['fail_reason'] = f'too few stars ({len(xy_tgt)})'
        if debug:
            sys.stderr.write(f"  FAIL: {info['fail_reason']}\n")
        return None, info

    # Track last failure reason across try_match calls
    last_fail = [None]

    def try_match(xy_t, mirrored=False):
        """Try matching with given target positions."""
        image_diag = np.sqrt(ref_shape[0] ** 2 + ref_shape[1] ** 2)
        prefix = "[mirror] " if mirrored else ""

        # Build target descriptors
        t1 = time.time()
        tree_tgt = star_match.build_kdtree(xy_t)
        hashes_tgt, idx_tgt = star_match.build_descriptors(
            xy_t, tree_tgt, opts['n_descriptors'])
        t_desc = time.time() - t1

        # Match descriptors
        t1 = time.time()
        pairs, votes = star_match.match_descriptors(
            hashes_ref, idx_ref, hashes_tgt, idx_tgt,
            tolerance=opts['hash_tol'])
        t_hash = time.time() - t1

        # All matches (before vote filter) for diagnostics
        n_all_matches = len(pairs)

        # Filter by votes
        vote_mask = votes >= opts['min_vote']
        pairs = pairs[vote_mask]

        mv = opts['min_vote']
        if debug:
            sys.stderr.write(f"  {prefix}descriptors: {len(hashes_tgt)} "
                             f"({t_desc:.3f}s), hash matches: {n_all_matches} "
                             f"(vote>={mv}: {len(pairs)}) ({t_hash:.3f}s)\n")

        if len(pairs) < 4:
            last_fail[0] = (f'hash match: {n_all_matches} total, '
                            f'{len(pairs)} with vote>={mv} (need >=4)')
            if debug:
                # Show vote distribution for diagnostics
                if n_all_matches > 0:
                    top_votes = votes[:min(10, len(votes))]
                    sys.stderr.write(f"  {prefix}  top votes: {top_votes.tolist()}\n")
                sys.stderr.write(f"  {prefix}  >> STOP: not enough hash matches\n")
            return None

        # Angular verification
        t1 = time.time()
        verified = star_match.verify_angular(
            xy_ref, tree_ref, xy_t, tree_tgt,
            pairs, image_diag,
            angle_tol=opts['angle_tol'])
        t_angular = time.time() - t1
        n_pre_angular = len(pairs)
        pairs = pairs[verified]

        if debug:
            sys.stderr.write(f"  {prefix}angular verify: {n_pre_angular} -> "
                             f"{len(pairs)} ({t_angular:.3f}s)\n")

        if len(pairs) < 4:
            last_fail[0] = (f'angular verify: {n_pre_angular} -> {len(pairs)} '
                            f'(need >=4)')
            if debug:
                sys.stderr.write(f"  {prefix}  >> STOP: angular verification "
                                 f"rejected too many pairs\n")
            return None

        # RANSAC
        t1 = time.time()
        ransac_model = 'projective' if opts['model'] == 'tps' else opts['model']
        transform, inlier_mask, median_res = star_match.ransac_fit(
            xy_ref, xy_t, pairs,
            model=ransac_model,
            inlier_tolerance=opts['tolerance'],
        )
        t_ransac = time.time() - t1

        if transform is None:
            last_fail[0] = (f'RANSAC failed ({len(pairs)} pairs, '
                            f'tolerance={opts["tolerance"]:.1f}px)')
            if debug:
                sys.stderr.write(f"  {prefix}  >> STOP: RANSAC could not fit "
                                 f"model from {len(pairs)} pairs\n")
            return None

        n_inliers = int(np.sum(inlier_mask))
        if debug:
            sys.stderr.write(f"  {prefix}RANSAC: {n_inliers} inliers, "
                             f"residual {median_res:.3f} px ({t_ransac:.3f}s)\n")

        timings['descriptors'] = t_desc
        timings['hash_match'] = t_hash
        timings['angular'] = t_angular
        timings['ransac'] = t_ransac

        return {
            'transform': transform,
            'pairs': pairs[inlier_mask],
            'median_residual': median_res,
            'n_hash_matches': len(pairs),
            'n_verified': len(pairs),
            'n_inliers': n_inliers,
            'ransac_model': ransac_model,
            'xy_tgt': xy_t,
        }

    # Try normal orientation
    match_result = try_match(xy_tgt)

    # Try mirror fallback
    if match_result is None and not opts['no_mirror']:
        if debug:
            sys.stderr.write("  Trying mirror fallback...\n")
        # Mirror target stars horizontally
        if target_data.ndim == 3:
            w = target_data.shape[2]
        else:
            w = target_data.shape[1]
        xy_tgt_mirror = xy_tgt.copy()
        xy_tgt_mirror[:, 0] = w - 1 - xy_tgt_mirror[:, 0]
        match_result = try_match(xy_tgt_mirror, mirrored=True)
        if match_result is not None:
            info['mirrored'] = True
            # Mirror the actual image data
            if target_data.ndim == 3:
                target_data = target_data[:, :, ::-1].copy()
            else:
                target_data = target_data[:, ::-1].copy()

    if match_result is None:
        if last_fail[0]:
            info['fail_reason'] = last_fail[0]
        else:
            info['fail_reason'] = 'matching failed (unknown)'
        if debug:
            sys.stderr.write(f"  FAIL: {info['fail_reason']}\n")
        return None, info

    info['n_hash_matches'] = match_result['n_hash_matches']
    info['n_verified'] = match_result['n_verified']
    info['n_inliers'] = match_result['n_inliers']
    info['median_residual'] = match_result['median_residual']

    xy_t = match_result['xy_tgt']
    inlier_pairs = match_result['pairs']
    pts_ref = xy_ref[inlier_pairs[:, 0]]
    pts_tgt = xy_t[inlier_pairs[:, 1]]

    # TPS refinement or direct projective resampling
    if opts['model'] == 'tps':
        t0 = time.time()
        tps_result = star_tps.tps_iterative(
            xy_ref, xy_t,
            match_result['transform'],
            kdtree_ref=tree_ref,
            model=match_result['ransac_model'],
            image_shape=ref_shape,
            smoothing=opts['smoothness'],
        )
        timings['tps_refine'] = time.time() - t0

        info['n_tps_pairs'] = len(tps_result['pairs']) if tps_result['pairs'] is not None else 0
        info['tps_median_residual'] = tps_result['median_residual']

        if debug:
            sys.stderr.write(f"  TPS: {info['n_tps_pairs']} pairs, "
                             f"{tps_result['n_iterations']} iterations, "
                             f"residual {tps_result['median_residual']:.3f} px, "
                             f"converged: {tps_result['converged']}, "
                             f"quality: {'OK' if tps_result['quality_ok'] else 'FAIL'} "
                             f"({timings['tps_refine']:.3f}s)\n")

        if tps_result['tps'] is not None and tps_result['pairs'] is not None and tps_result['quality_ok']:
            # Use TPS pairs for resampling
            final_ref = xy_ref[tps_result['pairs'][:, 0]]
            final_tgt = xy_t[tps_result['pairs'][:, 1]]

            t0 = time.time()
            aligned = star_resample.resample_frame(
                target_data.astype(np.float64),
                final_ref, final_tgt,
                output_shape,
                smoothing=opts['smoothness'],
            )
            timings['resample'] = time.time() - t0
        else:
            # TPS failed — do NOT fall back to raw projective (likely garbage)
            info['fail_reason'] = (
                f'TPS refinement failed '
                f'({info["n_tps_pairs"]} pairs, '
                f'residual {tps_result["median_residual"]:.3f} px)')
            if debug:
                sys.stderr.write(f"  FAIL: {info['fail_reason']}\n")
            return None, info
    else:
        # Direct projective/similarity resampling
        t0 = time.time()
        aligned = _resample_projective(
            target_data, pts_ref, pts_tgt,
            match_result['transform'],
            match_result['ransac_model'],
            output_shape)
        timings['resample'] = time.time() - t0
        info['tps_median_residual'] = info['median_residual']

    if debug:
        sys.stderr.write(f"  resample: ({timings['resample']:.3f}s)\n")
        parts = ' | '.join(f"{k}: {v:.3f}s" for k, v in timings.items())
        total_t = sum(timings.values())
        sys.stderr.write(f"  TIMING: {parts} | TOTAL: {total_t:.3f}s\n")

    info['success'] = True
    info['timings'] = timings
    return aligned, info


def _resample_projective(target_data, pts_ref, pts_tgt, transform, model,
                         output_shape):
    """Resample using projective/similarity transform (no TPS)."""
    try:
        if model == 'projective':
            inv_transform = np.linalg.inv(transform)
        else:
            # Invert similarity (2,3) → build 3x3, invert, extract 2x3
            M3 = np.eye(3, dtype=np.float64)
            M3[:2, :] = transform
            inv_M3 = np.linalg.inv(M3)
            inv_transform = inv_M3[:2, :]
    except np.linalg.LinAlgError:
        return None

    data_f = target_data.astype(np.float64)

    if data_f.ndim == 2:
        return star_resample.resample_image(
            data_f, inv_transform=inv_transform, model=model,
            output_shape=output_shape)
    else:
        result = np.empty((data_f.shape[0],) + output_shape, dtype=np.float64)
        for ch in range(data_f.shape[0]):
            result[ch] = star_resample.resample_image(
                data_f[ch], inv_transform=inv_transform, model=model,
                output_shape=output_shape)
        return result


def _get_retry_levels(base_max_stars, n_detected):
    """Build list of star counts to try (increasing).

    Uses multipliers [1, 2, 2.5, 3] to cover the search space densely
    enough to find the sweet spot between too-few and too-many stars.
    Each level is capped by n_detected. Deduplicates.
    """
    levels = []
    for mult in [1, 2, 2.5, 3]:
        n = min(int(base_max_stars * mult), n_detected)
        if not levels or n > levels[-1]:
            levels.append(n)
    return levels


def _process_frame(infile, outfile, ref_configs, retry_levels, ref_shape,
                   detect_limit, opts):
    """Process a single frame: load, align, save. Thread-safe.

    Returns dict with result info (no stderr output in parallel mode).
    """
    t0 = time.time()
    result = {
        'infile': infile,
        'outfile': outfile,
        'basename': os.path.basename(infile),
        'success': False,
        'retry_level': 0,
        'fail_reason': '',
        'info': None,
        'dt': 0.0,
    }

    try:
        hdu_tgt = fits.open(infile, memmap=False)
        tgt_data = hdu_tgt[0].data.copy()
        tgt_header = hdu_tgt[0].header.copy()
        orig_dtype = tgt_data.dtype
        hdu_tgt.close()

        # Detect target stars once (full catalog for retry)
        full_tgt_xy, _ = detect_and_extract(
            tgt_data, opts['snr'], detect_limit)

        # Try each retry level until success
        aligned = None
        info = None
        level_idx = 0
        for level_idx, rc in enumerate(ref_configs):
            n_stars = retry_levels[level_idx]
            tgt_xy = full_tgt_xy[:min(n_stars, len(full_tgt_xy))]

            aligned, info = align_single(
                tgt_data, tgt_header, rc['xy'], ref_shape,
                rc['tree'], rc['hashes'], rc['idx'],
                opts, debug=False, xy_tgt=tgt_xy)

            if aligned is not None:
                break

        if aligned is not None:
            # Convert back to original dtype
            if np.issubdtype(orig_dtype, np.integer):
                iinfo = np.iinfo(orig_dtype)
                aligned = np.clip(aligned, iinfo.min, iinfo.max)
                aligned = np.rint(aligned).astype(orig_dtype)
            elif np.issubdtype(orig_dtype, np.floating):
                aligned = np.nan_to_num(aligned, nan=0.0,
                                        posinf=0.0, neginf=0.0)
                aligned = aligned.astype(orig_dtype)

            # Update header
            tgt_header['HISTORY'] = (
                f'staralign: {info["n_inliers"]} inliers, '
                f'residual {info["tps_median_residual"]:.3f} px'
            )
            if info['mirrored']:
                tgt_header['HISTORY'] = 'staralign: mirrored horizontally'

            fits.PrimaryHDU(aligned, header=tgt_header).writeto(
                outfile, overwrite=True)

            result['success'] = True
            result['retry_level'] = level_idx
            result['info'] = info
        else:
            result['fail_reason'] = info.get('fail_reason', '?') if info else '?'
            result['info'] = info

    except Exception as e:
        result['fail_reason'] = str(e)

    result['dt'] = time.time() - t0
    return result


def _chromatic_align(opts):
    """Chromatic aberration correction: align R and B channels to G.

    Splits an RGB FITS into 3 mono channels, aligns R (ch0) and B (ch2)
    to G (ch1) using the standard star matching + TPS pipeline, then
    recombines into a corrected RGB output.
    """
    detect_limit = max(opts['max_stars'] * 3, 600)

    # Build IO pairs (may be single or batch)
    io_pairs = batch_utils.build_io_file_lists(opts['input_spec'],
                                               opts['output_spec'])
    total = len(io_pairs)

    n_threads = 1 if opts['debug'] else opts['n_threads']
    sys.stderr.write(f"Chromatic correction: {total} file(s)\n")

    n_success = 0
    n_fail = 0
    failed_files = []
    t_start = time.time()
    channel_names = ['R', 'G', 'B']

    for file_i, (infile, outfile) in enumerate(io_pairs, start=1):
        basename = os.path.basename(infile)
        t0 = time.time()
        try:
            hdu = fits.open(infile, memmap=False)
            data = hdu[0].data.copy()
            header = hdu[0].header.copy()
            orig_dtype = data.dtype
            hdu.close()

            if data.ndim != 3 or data.shape[0] < 3:
                sys.stderr.write(
                    f"[{file_i}/{total}] {basename}: FAIL "
                    f"not RGB ({data.ndim}D, {data.shape[0] if data.ndim==3 else 1} ch)"
                    f" -- chromatic mode needs a colour frame; to align frames pass"
                    f" --ref FRAME.fit (or give several frames)\n")
                n_fail += 1
                failed_files.append((basename, 'not RGB'))
                continue

            # G channel is the reference
            g_data = data[1].astype(np.float64)
            ref_shape = g_data.shape  # (H, W)

            if opts['debug']:
                sys.stderr.write(f"\n--- [{file_i}/{total}] {basename} ---\n")
                sys.stderr.write(f"  Frame {ref_shape[1]}x{ref_shape[0]}\n")

            # Detect stars on G channel
            g_xy, n_g_detected = detect_and_extract(
                g_data, opts['snr'], detect_limit,
                opts['debug'], "G-channel")

            if len(g_xy) < 5:
                sys.stderr.write(
                    f"[{file_i}/{total}] {basename}: FAIL "
                    f"too few stars on G ({len(g_xy)})\n")
                n_fail += 1
                failed_files.append((basename, f'too few stars on G ({len(g_xy)})'))
                continue

            # Build retry levels
            if opts['no_retry']:
                retry_levels = [min(opts['max_stars'], len(g_xy))]
            else:
                retry_levels = _get_retry_levels(opts['max_stars'], len(g_xy))

            # Pre-build G descriptors for each retry level
            ref_configs = []
            for n in retry_levels:
                xy = g_xy[:n]
                tree = star_match.build_kdtree(xy)
                hashes, idx = star_match.build_descriptors(
                    xy, tree, opts['n_descriptors'])
                ref_configs.append({
                    'xy': xy, 'tree': tree, 'hashes': hashes, 'idx': idx,
                })

            if opts['debug'] and len(retry_levels) > 1:
                sys.stderr.write(f"  G: {len(g_xy)} stars, "
                                 f"retry levels: {retry_levels}\n")

            # Result: start with original data
            result = data.astype(np.float64)
            ch_ok = True
            ch_results = []

            # Align R (ch0) and B (ch2) to G
            for ch_idx in [0, 2]:
                ch_name = channel_names[ch_idx]
                ch_data = data[ch_idx].astype(np.float64)

                if opts['debug']:
                    sys.stderr.write(f"\n  === Aligning {ch_name} to G ===\n")

                # Detect stars on this channel
                ch_xy, _ = detect_and_extract(
                    ch_data, opts['snr'], detect_limit,
                    opts['debug'], f"{ch_name}-channel")

                # Try each retry level
                aligned_ch = None
                info = None
                level_idx = 0
                for level_idx, rc in enumerate(ref_configs):
                    n_stars = retry_levels[level_idx]
                    tgt_xy = ch_xy[:min(n_stars, len(ch_xy))]

                    if level_idx > 0 and opts['debug']:
                        sys.stderr.write(
                            f"  {ch_name} auto-retry with {n_stars} stars...\n")

                    aligned_ch, info = align_single(
                        ch_data, header, rc['xy'], ref_shape,
                        rc['tree'], rc['hashes'], rc['idx'],
                        opts, opts['debug'], xy_tgt=tgt_xy)

                    if aligned_ch is not None:
                        break

                if aligned_ch is not None:
                    result[ch_idx] = aligned_ch
                    retry_tag = f" retry {level_idx}" if level_idx > 0 else ""
                    ch_results.append(
                        f"{ch_name}: {info['n_inliers']} inliers, "
                        f"{info['median_residual']:.2f} px"
                        f"{retry_tag}")
                else:
                    ch_ok = False
                    reason = info.get('fail_reason', '?') if info else '?'
                    ch_results.append(f"{ch_name}: FAIL {reason}")

            if ch_ok:
                # Convert back to original dtype
                if np.issubdtype(orig_dtype, np.integer):
                    iinfo = np.iinfo(orig_dtype)
                    result = np.clip(result, iinfo.min, iinfo.max)
                    result = np.rint(result).astype(orig_dtype)
                elif np.issubdtype(orig_dtype, np.floating):
                    result = np.nan_to_num(result, nan=0.0,
                                           posinf=0.0, neginf=0.0)
                    result = result.astype(orig_dtype)

                header['HISTORY'] = 'staralign: chromatic correction (R,B aligned to G)'
                fits.PrimaryHDU(result, header=header).writeto(
                    outfile, overwrite=True)
                n_success += 1

                dt = time.time() - t0
                detail = '; '.join(ch_results)
                sys.stderr.write(
                    f"[{file_i}/{total}] {basename}: OK ({detail})\n")
            else:
                n_fail += 1
                detail = '; '.join(ch_results)
                failed_files.append((basename, detail))
                sys.stderr.write(
                    f"[{file_i}/{total}] {basename}: FAIL ({detail})\n")

        except Exception as e:
            n_fail += 1
            failed_files.append((basename, str(e)))
            sys.stderr.write(f"\nError processing '{infile}': {e}\n")
            if opts['debug']:
                import traceback
                traceback.print_exc()

    dt_total = time.time() - t_start
    n_processed = n_success + n_fail

    # Summary
    parts = [f"{n_success} OK", f"{n_fail} failed"]
    sys.stderr.write(f"\nDone: {', '.join(parts)}\n")

    if n_processed > 0:
        fps = n_processed / dt_total
        mm = int(dt_total) // 60
        ss = dt_total - mm * 60
        sys.stderr.write(f"Time: {mm}m {ss:.0f}s for {n_processed} frames "
                         f"({fps:.2f} frames/sec)\n")

    if failed_files:
        sys.stderr.write(f"\nFailed files ({len(failed_files)}):\n")
        for fn, reason in failed_files:
            sys.stderr.write(f"  {fn}  —  {reason}\n")

    if n_fail > 0:
        sys.exit(1)


def main():
    opts = parse_args(sys.argv)

    # No --ref: default is star-based frame alignment to the FIRST frame.
    # Per-channel (chromatic) alignment happens only with an explicit --rgb, or
    # when there is a single input file.
    if opts['ref_file'] is None:
        try:
            files = batch_utils.expand_input_spec(opts['input_spec'])
        except Exception as e:
            sys.stderr.write(f"Error: {e}\n")
            sys.exit(1)
        if opts['rgb'] or len(files) == 1:
            _chromatic_align(opts)
            return
        opts['ref_file'] = files[0]
        sys.stderr.write(
            f"No --ref given: aligning frames to the first as reference "
            f"({os.path.basename(files[0])}).\n")
        # fall through to frame alignment below

    # Force single-threaded in debug mode
    n_threads = 1 if opts['debug'] else opts['n_threads']

    # Detection limit: detect enough stars for retry levels
    detect_limit = max(opts['max_stars'] * 3, 600)

    # Load reference
    sys.stderr.write(f"Loading reference: {opts['ref_file']}\n")
    hdu_ref = fits.open(opts['ref_file'], memmap=False)
    ref_data = hdu_ref[0].data.copy()
    hdu_ref.close()

    if ref_data.ndim == 3:
        ref_shape = ref_data.shape[1:]  # (H, W)
    else:
        ref_shape = ref_data.shape  # (H, W)

    # Detect reference stars (full catalog for retry levels)
    sys.stderr.write("Detecting reference stars...\n")
    full_ref_xy, n_ref_detected = detect_and_extract(
        ref_data, opts['snr'], detect_limit, opts['debug'], "reference")

    if len(full_ref_xy) < 5:
        sys.stderr.write(f"Error: only {len(full_ref_xy)} stars on reference, need >= 5\n")
        sys.exit(1)

    # Build retry levels
    if opts['no_retry']:
        retry_levels = [min(opts['max_stars'], len(full_ref_xy))]
    else:
        retry_levels = _get_retry_levels(opts['max_stars'], len(full_ref_xy))

    # Pre-build reference data for each retry level
    ref_configs = []
    for n in retry_levels:
        xy = full_ref_xy[:n]
        tree = star_match.build_kdtree(xy)
        hashes, idx = star_match.build_descriptors(xy, tree, opts['n_descriptors'])
        ref_configs.append({
            'xy': xy, 'tree': tree, 'hashes': hashes, 'idx': idx,
        })

    if len(retry_levels) > 1:
        sys.stderr.write(f"Reference: {len(full_ref_xy)} stars, "
                         f"frame {ref_shape[1]}x{ref_shape[0]}, "
                         f"retry levels: {retry_levels}\n")
    else:
        sys.stderr.write(f"Reference: {ref_configs[0]['xy'].shape[0]} stars, "
                         f"frame {ref_shape[1]}x{ref_shape[0]}\n")

    sys.stderr.write(f"Reference descriptors: {len(ref_configs[0]['hashes'])}\n")

    # Build IO pairs
    io_pairs = batch_utils.build_io_file_lists(opts['input_spec'],
                                               opts['output_spec'])
    total = len(io_pairs)
    sys.stderr.write(f"Processing {total} target frame(s)"
                     f" with {n_threads} thread(s)...\n")

    n_success = 0
    n_fail = 0
    n_retried = 0
    failed_files = []
    t_start = time.time()

    if n_threads > 1:
        # ---- Parallel mode ----
        lock = threading.Lock()
        counter = [0]  # completed count

        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            futures = {}
            for idx, (infile, outfile) in enumerate(io_pairs):
                f = executor.submit(
                    _process_frame, infile, outfile,
                    ref_configs, retry_levels, ref_shape,
                    detect_limit, opts)
                futures[f] = idx

            for future in as_completed(futures):
                result = future.result()
                with lock:
                    counter[0] += 1
                    i = counter[0]

                if result['success']:
                    n_success += 1
                    info = result['info']
                    lvl = result['retry_level']
                    if lvl > 0:
                        n_retried += 1
                    retry_tag = f" [retry {lvl}]" if lvl > 0 else ""
                    sys.stderr.write(
                        f"[{i}/{total}] {result['basename']}: OK "
                        f"{info['n_inliers']} inliers, "
                        f"{info['n_tps_pairs']} tps, "
                        f"{info['median_residual']:.2f} px"
                        f"{retry_tag}\n")
                    sys.stderr.flush()
                else:
                    n_fail += 1
                    failed_files.append(
                        (result['basename'], result['fail_reason']))
                    sys.stderr.write(
                        f"[{i}/{total}] {result['basename']}: FAIL "
                        f"{result['fail_reason']}\n")
                    sys.stderr.flush()

    else:
        # ---- Sequential mode (or debug) ----
        for i, (infile, outfile) in enumerate(io_pairs, start=1):
            t0 = time.time()

            if opts['debug']:
                sys.stderr.write(
                    f"\n--- [{i}/{total}] {os.path.basename(infile)} ---\n")

            try:
                t_load = time.time()
                hdu_tgt = fits.open(infile, memmap=False)
                tgt_data = hdu_tgt[0].data.copy()
                tgt_header = hdu_tgt[0].header.copy()
                orig_dtype = tgt_data.dtype
                hdu_tgt.close()
                t_load = time.time() - t_load

                if opts['debug']:
                    sys.stderr.write(f"  load FITS: ({t_load:.3f}s)\n")

                # Detect target stars once (full catalog for retry)
                full_tgt_xy, _ = detect_and_extract(
                    tgt_data, opts['snr'], detect_limit,
                    opts['debug'], "target")

                # Try each retry level until success
                aligned = None
                info = None
                level_idx = 0
                for level_idx, rc in enumerate(ref_configs):
                    n_stars = retry_levels[level_idx]
                    tgt_xy = full_tgt_xy[:min(n_stars, len(full_tgt_xy))]

                    if level_idx > 0 and opts['debug']:
                        sys.stderr.write(
                            f"  Auto-retry with {n_stars} stars "
                            f"(ref={rc['xy'].shape[0]}, "
                            f"tgt={len(tgt_xy)})...\n")

                    aligned, info = align_single(
                        tgt_data, tgt_header, rc['xy'], ref_shape,
                        rc['tree'], rc['hashes'], rc['idx'],
                        opts, opts['debug'], xy_tgt=tgt_xy)

                    if aligned is not None:
                        if level_idx > 0:
                            n_retried += 1
                            if opts['debug']:
                                sys.stderr.write(
                                    f"  Auto-retry succeeded at level "
                                    f"{level_idx} ({n_stars} stars)\n")
                        break

                if aligned is not None:
                    # Convert back to original dtype
                    if np.issubdtype(orig_dtype, np.integer):
                        iinfo = np.iinfo(orig_dtype)
                        aligned = np.clip(aligned, iinfo.min, iinfo.max)
                        aligned = np.rint(aligned).astype(orig_dtype)
                    elif np.issubdtype(orig_dtype, np.floating):
                        aligned = np.nan_to_num(aligned, nan=0.0,
                                                posinf=0.0, neginf=0.0)
                        aligned = aligned.astype(orig_dtype)

                    # Update header
                    tgt_header['HISTORY'] = (
                        f'staralign: {info["n_inliers"]} inliers, '
                        f'residual {info["tps_median_residual"]:.3f} px'
                    )
                    if info['mirrored']:
                        tgt_header['HISTORY'] = (
                            'staralign: mirrored horizontally')

                    t_save = time.time()
                    fits.PrimaryHDU(aligned, header=tgt_header).writeto(
                        outfile, overwrite=True)
                    t_save = time.time() - t_save
                    n_success += 1

                    dt = time.time() - t0
                    retry_tag = (f" [retry {level_idx}]"
                                 if level_idx > 0 else "")
                    if not opts['debug']:
                        sys.stderr.write(
                            f"\r[{i}/{total}] "
                            f"{os.path.basename(infile)}: OK "
                            f"{info['n_inliers']} inliers, "
                            f"{info['n_tps_pairs']} tps, "
                            f"{info['median_residual']:.2f} px"
                            f"{retry_tag}")
                        sys.stderr.flush()
                    else:
                        sys.stderr.write(
                            f"  save FITS: ({t_save:.3f}s)\n")
                        sys.stderr.write(
                            f"  TOTAL: {dt:.1f}s{retry_tag}\n")
                else:
                    n_fail += 1
                    reason = (info.get('fail_reason', '?')
                              if info else '?')
                    failed_files.append(
                        (os.path.basename(infile), reason))
                    if not opts['debug']:
                        sys.stderr.write(
                            f"\r[{i}/{total}] "
                            f"{os.path.basename(infile)}: FAIL "
                            f"{reason}")
                        sys.stderr.flush()

            except Exception as e:
                n_fail += 1
                failed_files.append(
                    (os.path.basename(infile), str(e)))
                sys.stderr.write(
                    f"\nError processing '{infile}': {e}\n")
                if opts['debug']:
                    import traceback
                    traceback.print_exc()

    dt_total = time.time() - t_start
    n_processed = n_success + n_fail

    # Summary line
    parts = [f"{n_success} OK", f"{n_fail} failed"]
    if n_retried > 0:
        parts.append(f"{n_retried} retried")
    sys.stderr.write(f"\n\nDone: {', '.join(parts)}\n")

    # Timing
    if n_processed > 0:
        fps = n_processed / dt_total
        mm = int(dt_total) // 60
        ss = dt_total - mm * 60
        sys.stderr.write(f"Time: {mm}m {ss:.0f}s for {n_processed} frames "
                         f"({fps:.2f} frames/sec")
        if n_threads > 1:
            sys.stderr.write(f", {n_threads} threads")
        sys.stderr.write(")\n")

    if failed_files:
        # Group by failure reason
        from collections import Counter
        reason_counts = Counter(reason for _, reason in failed_files)
        sys.stderr.write(f"\nFailure summary:\n")
        for reason, cnt in reason_counts.most_common():
            sys.stderr.write(f"  {cnt:4d}x  {reason}\n")

        sys.stderr.write(f"\nFailed files ({len(failed_files)}):\n")
        for fn, reason in failed_files:
            sys.stderr.write(f"  {fn}  —  {reason}\n")

    if n_fail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
