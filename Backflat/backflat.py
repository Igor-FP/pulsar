#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""backflat - Interactive masked diffusion background flattening for RGB FITS."""

import hashlib
import json
import math
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import textwrap
import time
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager

import numpy as np
from astropy.io import fits
from scipy import ndimage
from scipy.signal import fftconvolve

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


REFERENCE_DIAGONAL = 7515.0
DEFAULTS = {
    "median1": 40.0, "edge": 0.025, "blur1": 40.0, "median2": 60.0,
    "lake_mask": 100.0, "lake_blur": 250.0, "blur": 60.0,
    "center_plateau": 0.5, "center_mask": 150.0, "center_blur": 225.0,
    "margin": None, "diffusion_scale": 4.0, "median2_scale": 4.0,
}
PREPARATION_KEYS = ("median1", "edge", "blur1", "median2", "median_method", "median2_scale")


def ascii_text(value):
    """Keep diagnostics safe even when a filename or exception is not ASCII."""
    return str(value).encode("ascii", errors="backslashreplace").decode("ascii")


def report(message, error=False):
    stream = sys.stderr if error else sys.stdout
    stream.write(ascii_text(message) + "\n")
    stream.flush()


class OperationLog:
    """ASCII console timings and a lightweight status shared with the GUI."""

    def __init__(self):
        self.current = None
        self.last_message = "Ready"
        self.last_report = 0.0
        self.timings = {}

    def write(self, message):
        self.last_message = ascii_text(message)
        report("[{}] {}".format(time.strftime("%H:%M:%S"), self.last_message))

    def update(self, detail):
        if self.current is not None:
            label, start, _ = self.current
            self.current = (label, start, ascii_text(detail))
            now = time.perf_counter()
            if now-self.last_report >= 5:
                self.write(self.status())
                self.last_report = now

    def status(self):
        current = self.current
        if current is None:
            return self.last_message
        label, start, detail = current
        return "{}: {:.0f}s{}".format(label, time.perf_counter()-start,
                                     "; "+detail if detail else "")


@contextmanager
def operation(log, label):
    """Optional tracing keeps standalone numerical functions silent."""
    start = time.perf_counter()
    if log is not None:
        log.current = (label, start, "")
        log.last_report = start
        log.write(label + "...")
    try:
        yield
    except BaseException:
        if log is not None:
            log.current = None
            log.write("{} failed after {:.2f} s".format(label, time.perf_counter()-start))
        raise
    else:
        if log is not None:
            elapsed = time.perf_counter()-start
            log.timings[label] = elapsed
            log.current = None
            log.write("{} completed in {:.2f} s".format(label, elapsed))


def read_rgb(path):
    """Read physical FITS values without assuming a normalized intensity range."""
    with fits.open(path, memmap=False) as hdus:
        data = hdus[0].data
        header = hdus[0].header.copy()
        if data is None or data.ndim != 3 or data.shape[0] != 3:
            raise ValueError("Expected RGB FITS with shape (3, H, W): " + path)
        if min(data.shape[1:]) < 3:
            raise ValueError("Image dimensions must be at least 3 pixels.")
        if data.dtype.kind not in "iuf":
            raise ValueError("Unsupported FITS sample type: " + str(data.dtype))
        dtype = data.dtype
        image = np.array(data, dtype=np.float64, order="C")
    invalid = ~np.isfinite(image)
    if np.any(invalid):
        report("Warning: replaced {} non-finite input samples with zero in {}".format(
            np.count_nonzero(invalid), path), error=True)
        image[invalid] = 0.0
    return image, header, dtype


def clean_header(header):
    """Leave observation/WCS cards intact; regenerate image storage cards."""
    result = header.copy()
    for key in ("BSCALE", "BZERO", "BLANK", "DATAMIN", "DATAMAX", "CHECKSUM", "DATASUM"):
        result.remove(key, ignore_missing=True, remove_all=True)
    return result


def commit_fits(temporary, path, overwrite):
    """Publish a complete FITS file without replacing an unapproved output."""
    if overwrite:
        os.replace(temporary, path)
    else:
        try:
            if os.name == "nt":
                os.rename(temporary, path)  # Windows rename never replaces a file.
            else:
                os.link(temporary, path)  # Atomic no-clobber publication on POSIX.
        except FileExistsError:
            raise FileExistsError("Output already exists: {}\n"
                                  "Use --overwrite (-y) to replace existing outputs.".format(path)) from None


def atomic_fits(path, data, header, overwrite=False):
    """Replace a completed artifact only after its FITS write succeeds."""
    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".backflat-", suffix=".fit",
                                      dir=os.path.dirname(path))
    os.close(fd)
    try:
        fits.PrimaryHDU(data, header=header).writeto(temporary, overwrite=True)
        commit_fits(temporary, path, overwrite)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def finite_float(data, dtype=np.float32):
    """Sanitize after conversion as well, without clipping valid negatives."""
    with np.errstate(over="ignore", invalid="ignore"):
        result = np.asarray(data, dtype=dtype).copy()
    invalid = ~np.isfinite(result)
    if np.any(invalid):
        report("Warning: replaced {} unrepresentable output samples with zero.".format(
            np.count_nonzero(invalid)), error=True)
        result[invalid] = 0.0
    return result


def scaled_size(diameter, shape):
    return float(diameter) * math.hypot(*shape[-2:]) / REFERENCE_DIAGONAL


def validate_params(params):
    for key in DEFAULTS:
        value = params[key]
        if key == "margin" and value is None:
            continue
        if not math.isfinite(value) or value < 0:
            raise ValueError("--{} must be finite and non-negative.".format(key.replace("_", "-")))
    if not 0 <= params["edge"] < 0.5:
        raise ValueError("--edge must be in [0, 0.5); the actual strips must fit the image.")
    if not 0 < params["center_plateau"] < 1:
        raise ValueError("--center-plateau must be strictly between 0 and 1.")
    if params["center_blur"] <= 0 or params["center_mask"] <= 0:
        raise ValueError("Central smoothing is mandatory: --center-blur and --center-mask must be > 0.")
    for key in ("diffusion_scale", "median2_scale"):
        if params[key] not in (1, 2, 4):
            raise ValueError("--{} must be 1, 2 or 4.".format(key.replace("_", "-")))


def gaussian_plane(plane, sigma):
    """A unit-sum Gaussian, truncated at 4 sigma, with reflected edges.

    FFT is only a convolution accelerator, not a background model or a grid fit.
    Both paths implement the same sampled kernel and half-sample reflection.
    """
    plane = np.asarray(plane, dtype=np.float64)
    if sigma <= 0:
        return plane.copy()
    if sigma < 8 or plane.size < 65536:
        return ndimage.gaussian_filter(plane, sigma, mode="reflect", truncate=4.0)
    radius = int(4.0*sigma + 0.5)
    coordinates = np.arange(-radius, radius+1, dtype=np.float64)
    kernel = np.exp(-0.5*(coordinates/sigma)**2)
    kernel /= kernel.sum()
    work = plane
    for axis in (0, 1):
        padding = [(0, 0), (0, 0)]
        padding[axis] = (radius, radius)
        shape = [1, 1]
        shape[axis] = kernel.size
        padded = np.pad(work, padding, mode="symmetric")
        work = fftconvolve(padded, kernel.reshape(shape), mode="valid", axes=(axis,))
    return work


def gaussian(image, diameter, reference_shape=None):
    """User diameters are FWHM, expressed at the reference diagonal."""
    shape = image.shape if reference_shape is None else reference_shape
    sigma = scaled_size(diameter, shape) / math.sqrt(8.0*math.log(2.0))
    if image.ndim == 2:
        return gaussian_plane(image, sigma)
    result = np.empty(image.shape, dtype=np.float64)
    for c in range(3):
        result[c] = gaussian_plane(image[c], sigma)
    return result


def median_aperture(image, diameter, method="fast", progress=None,
                    reference_shape=None, pixel_size=1):
    """True circular footprint in float64; neither backend quantizes values."""
    shape = image.shape if reference_shape is None else reference_shape
    radius = scaled_size(diameter, shape)/(2.0*pixel_size)
    extent = int(math.floor(radius))
    if extent < 1:
        return np.array(image, dtype=np.float64, copy=True)
    y, x = np.ogrid[-extent:extent+1, -extent:extent+1]
    footprint = x*x+y*y <= radius*radius
    if method == "fast":
        try:
            import diplib as dip
        except ImportError:
            raise RuntimeError("Fast circular median requires diplib. Install it with pip install diplib==3.6.1 or use --median-mode exact.")
        kernel = dip.Kernel(dip.Image(footprint))
        result = np.empty(image.shape, dtype=np.float64)
        for c in range(3):
            # DIPlib's built-in mirror convention differs from SciPy's reflect.
            # Explicit padding makes both paths use half-sample reflection.
            padded = np.pad(np.asarray(image[c], dtype=np.float64), extent, mode="symmetric")
            # PyDIP 3.6.1 holds the GIL inside each call. Short strips and an
            # explicit yield keep the GUI worker cooperative. The halo covers
            # the entire footprint, so strip boundaries do not change results.
            for y0 in range(0, image.shape[1], 128):
                y1 = min(image.shape[1], y0+128)
                ranked = dip.MedianFilter(dip.Image(padded[y0:y1+2*extent]), kernel)
                result[c, y0:y1] = np.asarray(ranked)[extent:-extent, extent:-extent]
                if progress is not None:
                    progress.update("channel {}/3, rows {}/{}".format(c+1, y1, image.shape[1]))
                time.sleep(0)
        return result
    if method != "exact":
        raise ValueError("Unsupported median method: " + method)
    result = np.empty(image.shape, dtype=np.float64)
    # Bound the footprint workspace independently of image size.
    # Tile halos reproduce the same reflected image as a full-frame filter.
    for c in range(3):
        h, w = image.shape[1:]
        for y0 in range(0, h, 256):
            y1 = min(h, y0+256)
            lo, hi = max(0, y0-extent), min(h, y1+extent)
            filtered = ndimage.median_filter(image[c, lo:hi], footprint=footprint, mode="reflect")
            result[c, y0:y1] = filtered[y0-lo:y1-lo]
            if progress is not None:
                progress.update("channel {}/3, rows {}/{}".format(c+1, y1, h))
    return result


def area_reduce(plane, factor):
    """Float64 area average; reflect partial border blocks without clipping."""
    h, w = plane.shape
    padded = np.pad(np.asarray(plane, dtype=np.float64),
                    ((0, (-h) % factor), (0, (-w) % factor)), mode="symmetric")
    return padded.reshape(padded.shape[0]//factor, factor,
                          padded.shape[1]//factor, factor).mean(axis=(1, 3))


def restore_plane(plane, factor, shape):
    """Bilinear interpolation of block centers onto the original pixel centers."""
    h, w = shape
    result = np.empty(shape, dtype=np.float64)
    xs = (np.arange(w, dtype=np.float64)+0.5)/factor-0.5
    for start in range(0, h, 128):
        stop = min(h, start+128)
        ys = (np.arange(start, stop, dtype=np.float64)+0.5)/factor-0.5
        yy, xx = np.meshgrid(ys, xs, indexing="ij")
        result[start:stop] = ndimage.map_coordinates(plane, [yy, xx], order=1,
                                                    mode="nearest", prefilter=False)
    return result


def prepare(image, starless, params, return_stats=False, progress=None):
    """Pure stages 2-5; star removal is resolved by the I/O layer beforehand."""
    validate_params(params)
    if image.shape != starless.shape or image.ndim != 3 or image.shape[0] != 3:
        raise ValueError("Image and starless input must have the same (3, H, W) shape.")
    method = params.get("median_method", "fast")
    stats = {"median_backend": "DIPlib" if method == "fast" else "SciPy"}
    with operation(progress, "Median 1 (D={:g})".format(params["median1"])):
        work = median_aperture(starless, params["median1"], method, progress)
    with operation(progress, "Mirror edges"):
        work = mirror_edges(work, params["edge"])
    with operation(progress, "Preparation Gaussian"):
        work = gaussian(work, params["blur1"])
    # Reduce only an already blurred image, with at least 12 pixels per disk
    # and 2 pixels per preparation sigma after reduction. No value quantization.
    support = min(scaled_size(params["median2"], image.shape)/12,
                  scaled_size(params["blur1"], image.shape)/(2*math.sqrt(8*math.log(2))))
    factor = min(int(params["median2_scale"]), 2**int(math.floor(math.log2(max(1.0, support)))))
    if method == "exact":
        factor = 1
    stats["median2_scale"] = factor
    with operation(progress, "Median 2 (D={:g}, scale={})".format(params["median2"], factor)):
        if factor == 1:
            result = median_aperture(work, params["median2"], method, progress)
        else:
            coarse = np.stack([area_reduce(plane, factor) for plane in work])
            coarse = median_aperture(coarse, params["median2"], method, progress,
                                      reference_shape=image.shape, pixel_size=factor)
            result = np.empty_like(work)
            for c in range(3):
                result[c] = restore_plane(coarse[c], factor, image.shape[1:])
    return (result, stats) if return_stats else result


def diffuse_holes(prepared, mask, progress=None, scale=1, sigma_mask=None):
    """Grow RGB and opacity from fixed shores with small repeated Gaussians.

    Masked original samples never participate. Newly reached pixels inherit a
    normalized blur of already supported values; opacity is densified each step.
    Original unmasked samples are reinstated on every iteration. This is the
    specified progressive Gaussian fill, not a claim of an exact Laplace solve.
    """
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != prepared.shape[1:]:
        raise ValueError("Mask and image dimensions differ.")
    if not np.any(mask):
        return prepared.copy(), {"iterations": 0, "max_distance": 0.0, "diffusion_scale": 1}
    if np.all(mask):
        raise ValueError("The mask covers the entire image; no background shores remain.")
    h, w = mask.shape
    # Work on the hole bounding box plus all Gaussian support. No resampling.
    if sigma_mask is None:
        sigma_mask = max(0.5, scaled_size(20.0, mask.shape)/math.sqrt(8*math.log(2)))
    sigma_image = 4*sigma_mask
    halo = int(math.ceil(4*sigma_image))+2
    rows, columns = np.nonzero(mask)
    y0, y1 = max(0, rows.min()-halo), min(h, rows.max()+halo+1)
    x0, x1 = max(0, columns.min()-halo), min(w, columns.max()+halo+1)
    holes = mask[y0:y1, x0:x1]
    original = prepared[:, y0:y1, x0:x1]
    # Keep at least two samples per sigma for the advancing coverage front.
    factor = min(int(scale), 2**max(0, int(math.floor(math.log2(sigma_mask/2)))))
    if factor > 1:
        if progress is not None:
            progress.write("{}x reduced diffusion; full-resolution shores retained".format(factor))
        rh, rw = holes.shape
        # Only known sky contributes, including cells straddling a hole's edge.
        weights = area_reduce((~holes).astype(np.float64), factor)
        coarse = np.empty((3, *weights.shape), dtype=np.float64)
        for c in range(3):
            values = area_reduce(np.where(holes, 0.0, original[c]), factor)
            np.divide(values, weights, out=coarse[c], where=weights > 0)
            coarse[c, weights == 0] = 0.0
        filled, stats = diffuse_holes(coarse, weights == 0, progress=progress,
                                      scale=1, sigma_mask=sigma_mask/factor)
        # Interpolate at the original pixel centers; never replace real shores.
        result = prepared.copy()
        for c in range(3):
            restored = restore_plane(filled[c], factor, (rh, rw))
            target = result[c, y0:y1, x0:x1]
            target[holes] = restored[holes]
        stats["diffusion_scale"] = factor
        stats["max_distance"] *= factor
        return result, stats
    alpha = (~holes).astype(np.float64)
    work = np.where(holes[None], 0.0, original)
    distance = ndimage.distance_transform_edt(mask)
    max_distance = float(distance.max())
    del distance, rows, columns
    limit = max(32, int(math.ceil(max_distance/sigma_mask))*8+32)
    for iteration in range(1, limit+1):
        if progress is not None:
            progress.update("iteration {}".format(iteration))
        coverage = np.maximum(gaussian_plane(alpha, sigma_image), 0.0)
        grown = np.maximum(alpha, np.clip(2*gaussian_plane(alpha, sigma_mask), 0.0, 1.0))
        grown[~holes] = 1.0
        # Avoid numerical FFT leakage seeding far-away pixels with invented color.
        supported = coverage > 1e-8
        grown[~supported & holes] = alpha[~supported & holes]
        increment = np.maximum(grown-alpha, 0.0)
        for c in range(3):
            numerator = gaussian_plane(work[c]*alpha, sigma_image)
            blurred = np.divide(numerator, coverage, out=np.zeros_like(numerator), where=supported)
            combined = work[c]*alpha + blurred*increment
            work[c] = np.divide(combined, grown, out=np.zeros_like(combined), where=grown > 0)
            work[c, ~holes] = original[c, ~holes]
        alpha = grown
        if np.all(alpha[holes] >= 1.0-1e-6):
            break
    else:
        raise RuntimeError("Diffusion did not close the mask; inspect the mask and its shores.")
    result = prepared.copy()
    result[:, y0:y1, x0:x1] = work
    return result, {"iterations": iteration, "max_distance": max_distance,
                    "diffusion_scale": 1}


def center_mask(shape, plateau, diameter):
    """Rectangular distance-to-edge plateau, blurred and black-point corrected."""
    h, w = shape
    y = np.minimum(np.arange(h), np.arange(h)[::-1])[:, None]
    x = np.minimum(np.arange(w), np.arange(w)[::-1])[None, :]
    distance = np.minimum(y, x)
    mask = (distance > plateau*float(distance.max())).astype(np.float64)
    mask = gaussian(mask, diameter)
    border = max(mask[0].max(), mask[-1].max(), mask[:, 0].max(), mask[:, -1].max())
    peak = float(mask.max())
    if peak <= border + 1e-12:
        raise ValueError("Center mask has no plateau at this image size; reduce --center-mask.")
    mask = np.clip((mask-border)/(1.0-border), 0.0, 1.0)
    mask[0] = mask[-1] = 0.0
    mask[:, 0] = mask[:, -1] = 0.0
    return mask


def lake_blending_mask(mask, diameter):
    """Blur the object mask, then linearly stretch its full range to [0, 1]."""
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask) or np.all(mask):
        # Preserve empty/full masks rather than stretching convolution roundoff.
        return mask.astype(np.float64)
    blurred = gaussian(mask.astype(np.float64), diameter)
    low, high = float(blurred.min()), float(blurred.max())
    if high > low:
        blurred = (blurred-low)/(high-low)
    return np.clip(blurred, 0.0, 1.0)


def build_fmf(prepared, mask, params, return_stats=False, progress=None):
    """Pure stages 7-9: fill, lake blur, global blur, mandatory center blur."""
    validate_params(params)
    with operation(progress, "Fill masked regions"):
        filled, stats = diffuse_holes(prepared, mask, progress=progress,
                                      scale=params["diffusion_scale"])
    with operation(progress, "Blur and auto-level lake mask"):
        lake_mask = lake_blending_mask(mask, params["lake_mask"])
    with operation(progress, "Lake Gaussian and blend"):
        for c in range(3):
            lake = gaussian(filled[c], params["lake_blur"])
            filled[c] += lake_mask*(lake-filled[c])
    del lake_mask
    with operation(progress, "Global Gaussian"):
        work = gaussian(filled, params["blur"])
    del filled
    with operation(progress, "Center mask"):
        central = center_mask(mask.shape, params["center_plateau"], params["center_mask"])
    with operation(progress, "Center Gaussian and blend"):
        for c in range(3):
            stronger = gaussian(work[c], params["center_blur"])
            work[c] += central*(stronger-work[c])
    if not np.isfinite(work).all():
        raise ValueError("Background computation produced non-finite values.")
    stats["channel_means"] = work.mean(axis=(1, 2)).tolist()
    stats["neutral_level"] = float(np.mean(stats["channel_means"], dtype=np.float64))
    stats["masked_fraction"] = float(np.mean(mask))
    return (work, stats) if return_stats else work


def mirror_edges(image, fraction):
    """Reflect adjacent interior strips: left/right, then top/bottom."""
    out = np.array(image, dtype=np.float64, copy=True)
    h, w = out.shape[-2:]
    band = int(round(float(fraction) * math.hypot(h, w)))
    if band == 0:
        return out
    if band < 0 or 2 * band > min(h, w):
        raise ValueError("Edge strips are too wide for this image; reduce --edge.")
    out[..., :band] = out[..., band:2 * band][..., ::-1]
    out[..., w-band:] = out[..., w-2*band:w-band][..., ::-1]
    out[..., :band, :] = out[..., band:2*band, :][..., ::-1, :]
    out[..., h-band:, :] = out[..., h-2*band:h-band, :][..., ::-1, :]
    return out


def subtract(image, fmf):
    """Subtract RGB background and add its common mean level, retaining negatives."""
    if image.shape != fmf.shape or image.ndim != 3 or image.shape[0] != 3:
        raise ValueError("Image and background must have matching (3, H, W) shapes.")
    result = np.empty(image.shape, dtype=np.float64)
    channel_means = np.mean(fmf, axis=(1, 2), dtype=np.float64)
    neutral_level = np.mean(channel_means, dtype=np.float64)
    for c in range(3):
        np.subtract(image[c], fmf[c], out=result[c], dtype=np.float64)
        result[c] += neutral_level
    return result


def expand_mask(mask, margin):
    """Add a Euclidean margin once to the raw, editable mask."""
    mask = np.asarray(mask, dtype=bool)
    if margin <= 0 or not np.any(mask):
        return mask.copy()
    return ndimage.distance_transform_edt(~mask) <= margin


def read_mask(path, shape):
    """Read full-scale opacity; RAWMASK avoids compounding the saved margin."""
    with fits.open(path, memmap=False) as hdus:
        raw = hdus["RAWMASK"].data if "RAWMASK" in hdus else hdus[0].data
        if raw is None or raw.shape != tuple(shape):
            raise ValueError("Mask dimensions do not match the image: " + path)
        # Any nonzero opacity excludes the pixel: err on the side of more sky.
        mask = np.isfinite(raw) & (raw > 0)
        applied_margin = float(hdus[0].header.get("BFMARGIN", 0.0))
        has_raw = "RAWMASK" in hdus
    return mask, applied_margin, has_raw


def save_mask(path, raw_mask, effective_mask, header, margin, overwrite=False):
    """Archive a makemask-compatible uint16 primary image plus editable source."""
    hdr = clean_header(header)
    hdr["BFMARGIN"] = (float(margin), "Applied mask margin in image pixels")
    hdr.add_history("backflat: primary mask includes margin; RAWMASK is editable.")
    # A monochrome mask has no third image axis.
    hdr.remove("NAXIS3", ignore_missing=True)
    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".backflat-mask-", suffix=".fit",
                                      dir=os.path.dirname(path))
    os.close(fd)
    try:
        hdus = fits.HDUList([
            fits.PrimaryHDU(effective_mask.astype(np.uint16) * 65535, header=hdr),
            fits.ImageHDU(raw_mask.astype(np.uint16) * 65535, name="RAWMASK"),
        ])
        hdus.writeto(temporary, overwrite=True)
        hdus.close()
        commit_fits(temporary, path, overwrite)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def external_run(command, cwd=None, timeout=None):
    """Capture external messages; never inherit non-ASCII console output."""
    result = subprocess.run(command, cwd=cwd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=timeout)
    output = result.stdout.decode("utf-8", errors="replace")
    if result.returncode:
        raise RuntimeError("External tool failed (exit {}):\n{}".format(
            result.returncode, output[-8000:]))
    return output


def discover_engine(kind, executable=None):
    """Validate the selected optional engine without installing/activating it."""
    names = [executable] if executable else (["rc-astro"] if kind == "sxt" else ["starnet2"])
    path = next((shutil.which(name) for name in names if shutil.which(name)), None)
    website = "https://www.rc-astro.com/" if kind == "sxt" else "https://starnetastro.com/cli-tools/"
    if path is None:
        raise RuntimeError("{} executable not found. Install it and its models yourself ({}) "
                           "or use --starless FILE.".format(kind, website))
    if kind == "sxt":
        info = json.loads(external_run([path, "sxt", "--json"], timeout=60))
        if int(info.get("schemaVersion", 0)) < 3:
            raise RuntimeError("SxT CLI schema >= 3 required; update RC-Astro or use --starless.")
        if not info.get("license", {}).get("valid", False):
            raise RuntimeError("SxT license unavailable: {}. See {} or use --starless FILE.".format(
                info.get("license", {}).get("message", "not activated"), website))
        return path, "sxt:{}:ML{}".format(info.get("cliVersion"), info.get("mlVersion"))
    version = external_run([path, "--version"], timeout=30).strip()
    help_text = external_run([path, "--help"], timeout=30)
    if not all(flag in help_text for flag in ("--linear", "--input", "--output")):
        raise RuntimeError("StarNet2 with FITS and --linear support is required (2.6+). "
                           "Legacy StarNet++ is not supported. See " + website)
    return path, version


def cached_starless(input_path, image, header, engine, executable, cache_dir):
    """Disk-cache only the external star-removal step, keyed by source/engine."""
    path = shutil.which(executable or ("rc-astro" if engine == "sxt" else "starnet2"))
    if path is None:
        # Reuse the actionable dependency diagnostic from the normal probe.
        discover_engine(engine, executable)
    binary = os.stat(path)
    digest = hashlib.sha256()
    with open(input_path, "rb") as source:
        for chunk in iter(lambda: source.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    identity = "{}:{}:{}:{}:backflat-starless-v2".format(
        engine, os.path.realpath(path), binary.st_size, binary.st_mtime_ns)
    if engine == "sxt":
        # Older SxT caches retain the CLI's inverted FITS row order.
        identity += ":sxt-flip-y-v1"
    digest.update(identity.encode("utf-8"))
    key = digest.hexdigest()
    os.makedirs(cache_dir, exist_ok=True)
    cached = os.path.join(cache_dir, key + ".fit")
    if os.path.isfile(cached):
        data, cache_header, _ = read_rgb(cached)
        if data.shape != image.shape:
            raise ValueError("Cached starless shape mismatch: " + cached)
        report("Reusing starless cache: " + cached)
        return data, cache_header.get("BFENGINE", engine + " cached")
    path, version = discover_engine(engine, executable)
    report("Removing stars with {} ({})...".format(engine, version))
    with tempfile.TemporaryDirectory(prefix="backflat-", dir=cache_dir) as temporary:
        source = os.path.join(temporary, "input.fit")
        target = os.path.join(temporary, "starless.fit")
        # Use one scale for RGB, preserve its inverse, and never stretch here.
        # Engines interpret float FITS in normalized physical units.
        scale = max(float(np.max(image)), float(np.max(np.abs(image))), 1.0)
        fits.PrimaryHDU(image / scale, header=clean_header(header)).writeto(source)
        if engine == "sxt":
            command = [path, "sxt", source, "-o", target, "--json"]
        else:
            if min(image.shape[1:]) < 512:
                raise ValueError("StarNet2 requires at least 512x512; use --starless FILE.")
            command = [path, "--input", source, "--output", target, "--linear"]
        output = external_run(command, cwd=os.path.dirname(path))
        for line in output.splitlines():
            if "warning" in line.lower():
                report("External warning: " + line, error=True)
        if not os.path.isfile(target):
            raise RuntimeError("Star-removal tool did not create output:\n" + output[-4000:])
        data, _, _ = read_rgb(target)
        if data.shape != image.shape:
            raise ValueError("Star-removal tool changed the image dimensions.")
        if engine == "sxt":
            # RC-Astro CLI reverses FITS rows. Restore input coordinates once,
            # before caching or applying masks; channel and X order stay intact.
            data = data[:, ::-1, :].copy()
            report("Restored SxT FITS orientation: Ynew = H - 1 - Yold.")
        data *= scale
        cache_header = clean_header(header)
        cache_header["BFENGINE"] = ascii_text(version)
        cache_header.add_history(ascii_text("backflat starless cache: " + version))
        if engine == "sxt":
            cache_header.add_history("backflat: corrected RC-Astro CLI Y flip; input orientation restored")
            version += "; Y flip corrected"
            cache_header["BFENGINE"] = ascii_text(version)
        atomic_fits(cached, data, cache_header, overwrite=True)
    return data, version


def cleanup_caches(directories, protected_paths):
    """Remove recognized Backflat cache files, then each empty cache directory.

    Never recursively delete a user-selected --cache-dir or follow a reparse
    point. Paths and FITS provenance both identify files owned by Backflat.
    """
    success = True
    for directory, resolved in directories.items():
        if not os.path.lexists(directory):
            continue
        try:
            info = os.lstat(directory)
            reparse = getattr(info, "st_file_attributes", 0) & getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0)
            if (os.path.islink(directory) or reparse or
                    os.path.normcase(os.path.realpath(directory)) != resolved):
                raise OSError("Cache path is a link or its resolved location changed: " + directory)
            if not stat.S_ISDIR(info.st_mode):
                raise OSError("Cache path is not a directory: " + directory)
            with operation(OperationLog(), "Clean starless cache"):
                with os.scandir(directory) as entries:
                    candidates = [entry.path for entry in entries
                                  if re.fullmatch(r"[0-9a-f]{64}\.fit", entry.name)
                                  and entry.is_file(follow_symlinks=False)]
                removed = 0
                for path in candidates:
                    if os.path.normcase(os.path.realpath(path)) in protected_paths:
                        continue
                    try:
                        header = fits.getheader(path)
                    except (OSError, ValueError):
                        continue
                    history = header.get("HISTORY", [])
                    if not any(str(line).startswith("backflat starless cache: ") for line in history):
                        continue
                    os.remove(path)
                    removed += 1
                with os.scandir(directory) as entries:
                    remaining = next(entries, None) is not None
                if remaining:
                    report("Cache files removed: {}; directory retained with other files: {}".format(
                        removed, directory))
                else:
                    os.rmdir(directory)
                    report("Removed starless cache directory: " + directory)
        except OSError as exc:
            report("Error cleaning starless cache: " + str(exc), error=True)
            success = False
    return success


def block_mean(image, factor):
    """Area-average a display level, retaining partial blocks at the borders."""
    if factor <= 1:
        return image
    h, w = image.shape[-2:]
    ys = np.arange(0, h, factor)
    xs = np.arange(0, w, factor)
    counts_y = np.minimum(factor, h - ys)
    counts_x = np.minimum(factor, w - xs)
    result = np.empty((3, len(ys), len(xs)), dtype=np.float32)
    for c in range(3):
        sums = np.add.reduceat(np.add.reduceat(image[c], ys, axis=0), xs, axis=1)
        result[c] = sums / counts_y[:, None] / counts_x[None, :]
    return result


def block_any(mask, factor):
    """Retain thin mask contours when drawing a reduced display level."""
    return np.logical_or.reduceat(np.logical_or.reduceat(
        mask, np.arange(0, mask.shape[0], factor), axis=0),
        np.arange(0, mask.shape[1], factor), axis=1)


def display_rgb(image, midtone, bounds):
    """A shared RGB display stretch only; never used by the scientific stages."""
    low, high = bounds
    x = np.clip((image.astype(np.float32) - low) / (high - low), 0.0, 1.0)
    x = (midtone - 1.0) * x / ((2.0 * midtone - 1.0) * x - midtone)
    return np.ascontiguousarray(np.moveaxis(x * 255.0, 0, -1), dtype=np.uint8)


class Session:
    """Own the arrays and stage invalidation; no dependency on pygame."""

    def __init__(self, image, starless, header, params, mask_path, raw_mask=None, progress=None,
                 overwrite=False):
        self.image = image
        self.starless = starless
        self.header = header
        self.params = dict(params)
        self.mask_path = mask_path
        self.overwrite = overwrite
        self.mask_saved = False
        self.progress = progress if progress is not None else OperationLog()
        self.raw_mask = (np.zeros(image.shape[1:], dtype=bool) if raw_mask is None
                         else np.array(raw_mask, dtype=bool, copy=True))
        self.prepared = None
        self.fmf = None
        self.result = None
        self.effective_mask = None
        self.applied = None
        self.prepare_stats = {}
        self.stats = {}
        self.revision = 0

    def margin(self, params=None):
        value = (self.params if params is None else params)["margin"]
        return (0.01 * math.hypot(*self.raw_mask.shape) if value is None
                else scaled_size(value, self.raw_mask.shape))

    def archive(self):
        with operation(self.progress, "Archive mask"):
            margin = self.margin()
            effective = expand_mask(self.raw_mask, margin)
            save_mask(self.mask_path, self.raw_mask, effective, self.header, margin,
                      overwrite=self.overwrite or self.mask_saved)
            self.mask_saved = True
        return effective

    def apply(self, model=True):
        """Commit arrays together; a failed calculation keeps the last preview."""
        params = dict(self.params)
        start = time.perf_counter()
        self.progress.timings = {}
        redo_prepare = self.applied is None or any(
            self.applied.get(k) != params.get(k) for k in PREPARATION_KEYS)
        if redo_prepare:
            prepared, preparation_stats = prepare(self.image, self.starless, params,
                                                   return_stats=True, progress=self.progress)
        else:
            prepared, preparation_stats = self.prepared, self.prepare_stats
            self.progress.write("Reusing prepared reference in memory")
        if model:
            with operation(self.progress, "Expand object mask"):
                effective = expand_mask(self.raw_mask, self.margin(params))
            fmf, stats = build_fmf(prepared, effective, params, return_stats=True,
                                   progress=self.progress)
            with operation(self.progress, "Subtract RGB background and add neutral level"):
                result = subtract(self.image, fmf)
        else:
            effective, fmf, result, stats = None, None, None, {}
        self.prepared, self.fmf, self.result = prepared, fmf, result
        self.effective_mask = effective
        self.stats = stats
        self.prepare_stats = preparation_stats
        self.stats.update(preparation_stats)
        self.stats["elapsed"] = time.perf_counter() - start
        self.stats["timings"] = dict(self.progress.timings)
        self.progress.write("Computation completed in {:.2f} s".format(self.stats["elapsed"]))
        self.applied = params
        self.revision += 1
        return self.stats


class MaskEditor:
    """Pygame-only UI, with explicit Apply and one undo record per stroke."""

    def __init__(self, session, loaded_mask=False):
        os.environ.setdefault("PYGAME_HIDE_SUPPORT_PROMPT", "1")
        try:
            import pygame
        except ImportError:
            raise RuntimeError("GUI requires pygame. Install pygame or use --no-gui --mask FILE.")
        self.pg = pygame
        self.session = session
        pygame.init()
        info = pygame.display.Info()
        self.width = max(640, min(1600, info.current_w - 100))
        self.height = max(640, min(1000, info.current_h - 100))
        self.screen = pygame.display.set_mode((self.width, self.height), pygame.RESIZABLE)
        pygame.display.set_caption("backflat - Masked Diffusion Background")
        pygame.key.set_repeat(150, 30)
        self.font = pygame.font.SysFont("consolas,courier,monospace", 14)
        self.status_h = 28
        self.panel_w = 310
        self.h, self.w = session.raw_mask.shape
        self.cx, self.cy = self.w / 2.0, self.h / 2.0
        self.zoom = max(0.02, min(self.view_w / self.w, self.view_h / self.h, 1.0))
        self.brush = max(5, int(round(math.hypot(self.h, self.w) * 0.025)))
        self.mtf = 0.05
        self.mode = "parameters" if loaded_mask else "mask"
        self.view = "prepared"
        self.show_mask = True
        self.selected = 0
        self.keys = list(DEFAULTS)
        self.edit_text = None
        self.undo = []
        self.stroke = None
        self.last_point = None
        self.paint_button = None
        self.pan_start = None
        self.held = set()
        self.need_render = True
        self.message = "Mark diffuse objects with generous margins. Tab: parameters."
        self.levels = {}
        self.surfaces = {}
        self.mask_dirty = True
        self.margin_outline = None
        self.display_revision = -1
        self.mask_revision = 0
        self.future = None
        self.pool = ThreadPoolExecutor(max_workers=1)
        self.closing = False
        # Keep the same black/white levels across all four views and channels.
        sample = block_mean(session.image, max(1, int(math.ceil(max(self.h, self.w) / 1024))))
        self.bounds = tuple(float(x) for x in np.percentile(sample, [0.01, 99.99]))
        if self.bounds[1] <= self.bounds[0]:
            self.bounds = (self.bounds[0], self.bounds[0] + 1.0)

    @property
    def view_w(self):
        return max(100, self.width - self.panel_w)

    @property
    def view_h(self):
        return self.height - self.status_h

    def image_to_screen(self, x, y):
        return (self.view_w / 2.0 + (x - self.cx) * self.zoom,
                self.view_h / 2.0 + (y - self.cy) * self.zoom)

    def screen_to_image(self, x, y):
        return (self.cx + (x - self.view_w / 2.0) / self.zoom,
                self.cy + (y - self.view_h / 2.0) / self.zoom)

    def set_zoom(self, zoom, position):
        x, y = self.screen_to_image(*position)
        self.zoom = min(40.0, max(0.02, zoom))
        self.cx = x - (position[0] - self.view_w / 2.0) / self.zoom
        self.cy = y - (position[1] - self.view_h / 2.0) / self.zoom
        self.need_render = True

    def paint(self, position, erase=False):
        x, y = self.screen_to_image(*position)
        if self.last_point is None:
            self.last_point = (x, y)
        previous_x, previous_y = self.last_point
        distance = math.hypot(x - previous_x, y - previous_y)
        count = max(1, int(math.ceil(distance / max(1.0, self.brush / 4.0))))
        radius = self.brush / 2.0
        for t in np.linspace(0.0, 1.0, count + 1):
            px = previous_x + t * (x - previous_x)
            py = previous_y + t * (y - previous_y)
            x0, x1 = max(0, int(math.floor(px-radius))), min(self.w, int(math.ceil(px+radius+1)))
            y0, y1 = max(0, int(math.floor(py-radius))), min(self.h, int(math.ceil(py+radius+1)))
            if x0 >= x1 or y0 >= y1:
                continue
            yy, xx = np.ogrid[y0:y1, x0:x1]
            disk = (xx-px)**2 + (yy-py)**2 <= radius**2
            self.session.raw_mask[y0:y1, x0:x1][disk] = not erase
        self.last_point = (x, y)
        self.mask_dirty = True
        self.mask_revision += 1
        self.need_render = True

    def finish_stroke(self):
        if self.stroke is not None:
            self.undo.append(self.stroke)
            self.undo = self.undo[-20:]
        self.stroke = self.last_point = self.paint_button = None

    def undo_stroke(self):
        if self.undo:
            packed = self.undo.pop()
            self.session.raw_mask[:] = np.unpackbits(packed, count=self.h*self.w).reshape(self.h, self.w)
            self.mask_dirty = True
            self.mask_revision += 1
            self.need_render = True

    def start_apply(self):
        if self.future is not None:
            return
        self.finish_stroke()
        try:
            self.commit_edit()
            validate_params(self.session.params)
            self.session.archive()
        except (ValueError, OSError) as exc:
            self.message = ascii_text(exc)
            self.closing = False
            return
        self.message = "Computing... Display navigation remains available."
        self.future = self.pool.submit(self.session.apply)
        self.need_render = True

    def commit_edit(self):
        if self.edit_text is None:
            return
        key = self.keys[self.selected]
        value = None if key == "margin" and self.edit_text.strip().lower() == "auto" else float(self.edit_text)
        candidate = dict(self.session.params)
        candidate[key] = value
        validate_params(candidate)
        self.session.params = candidate
        self.edit_text = None
        if key == "margin":
            self.mask_dirty = True

    def adjust_value(self, direction, fine=False):
        key = self.keys[self.selected]
        value = self.session.params[key]
        if value is None:
            value = 0.01 * REFERENCE_DIAGONAL
        step = (0.001 if key == "edge" else 0.05 if key == "center_plateau" else 5.0)
        if fine:
            step /= 5.0
        candidate = dict(self.session.params)
        candidate[key] = max(0.0, round(value + direction*step, 6))
        if key in ("diffusion_scale", "median2_scale"):
            values = (1.0, 2.0, 4.0)
            candidate[key] = values[min(2, max(0, values.index(value)+direction))]
        try:
            validate_params(candidate)
        except ValueError:
            return
        self.session.params = candidate
        if key == "margin":
            self.mask_dirty = True
        self.need_render = True

    def surface(self):
        if self.display_revision != self.session.revision:
            self.levels.clear()
            self.surfaces.clear()
            self.display_revision = self.session.revision
        data = getattr(self.session, {"original": "image", "prepared": "prepared",
                                     "background": "fmf", "result": "result"}[self.view])
        if data is None:
            data = self.session.prepared
        factor = max(1, 2 ** int(max(0, math.floor(math.log2(1.0 / self.zoom)))))
        key = (self.view, factor)
        if key not in self.levels:
            # Limit retained display pyramids; full-resolution science arrays stay untouched.
            if len(self.levels) >= 5:
                self.levels.clear()
                self.surfaces.clear()
            self.levels[key] = block_mean(data, factor)
        surface_key = (key, self.mtf)
        if surface_key not in self.surfaces:
            self.surfaces.clear()
            rgb = display_rgb(self.levels[key], self.mtf, self.bounds)
            self.surfaces[surface_key] = self.pg.surfarray.make_surface(rgb.transpose(1, 0, 2))
        return self.surfaces[surface_key], factor

    def render(self):
        pg = self.pg
        self.screen.fill((24, 24, 24))
        self.screen.set_clip(pg.Rect(0, 0, self.view_w, self.view_h))
        surf, factor = self.surface()
        x0 = max(0, int(math.floor((self.cx-self.view_w/2/self.zoom)/factor)))
        y0 = max(0, int(math.floor((self.cy-self.view_h/2/self.zoom)/factor)))
        x1 = min(surf.get_width(), int(math.ceil((self.cx+self.view_w/2/self.zoom)/factor)))
        y1 = min(surf.get_height(), int(math.ceil((self.cy+self.view_h/2/self.zoom)/factor)))
        if x1 > x0 and y1 > y0:
            target_size = (max(1, round((x1-x0)*factor*self.zoom)),
                           max(1, round((y1-y0)*factor*self.zoom)))
            position = tuple(round(v) for v in self.image_to_screen(x0*factor, y0*factor))
            sub = surf.subsurface(pg.Rect(x0, y0, x1-x0, y1-y0))
            self.screen.blit(pg.transform.scale(sub, target_size), position)
            if self.show_mask:
                # Mask contours need nearest-neighbor geometry, unlike intensity previews.
                raw = self.session.raw_mask[y0*factor:min(self.h,y1*factor):factor,
                                            x0*factor:min(self.w,x1*factor):factor]
                overlay = pg.Surface((raw.shape[1], raw.shape[0]), pg.SRCALPHA)
                overlay.fill((255, 255, 255, 0))
                alpha = pg.surfarray.pixels_alpha(overlay)
                alpha[:] = raw.T.astype(np.uint8) * 51
                del alpha
                self.screen.blit(pg.transform.scale(overlay, target_size), position)
                if self.mask_dirty and self.stroke is None:
                    effective = expand_mask(self.session.raw_mask, self.session.margin())
                    self.margin_outline = effective & ~ndimage.binary_erosion(effective)
                    self.mask_dirty = False
                if self.margin_outline is not None:
                    border = block_any(self.margin_outline[y0*factor:min(self.h,y1*factor),
                                                           x0*factor:min(self.w,x1*factor)], factor)
                    outline = pg.Surface((border.shape[1], border.shape[0]), pg.SRCALPHA)
                    outline.fill((255, 210, 50, 0))
                    alpha = pg.surfarray.pixels_alpha(outline)
                    alpha[:] = border.T.astype(np.uint8) * 230
                    del alpha
                    self.screen.blit(pg.transform.scale(outline, target_size), position)
        mx, my = pg.mouse.get_pos()
        if self.mode == "mask" and mx < self.view_w and my < self.view_h:
            pg.draw.circle(self.screen, (240, 240, 240), (mx, my),
                           max(1, round(self.brush*self.zoom/2)), 1)
        self.screen.set_clip(None)
        pg.draw.rect(self.screen, (38, 38, 42), (self.view_w, 0, self.panel_w, self.height))
        labels = ["BACKFLAT / " + self.mode.upper(),
                  "Tab: mask / parameters", "1 original  2 prepared",
                  "3 background  4 result", "Home/End: MTF  +/-: zoom",
                  "Wheel: brush  Shift: x10", "LMB paint / RMB erase",
                  "Ctrl+Z: undo stroke  M: mask", "B/Enter: Apply  S: save mask",
                  "Q/Esc/close: apply, save, exit"]
        for i, text in enumerate(labels):
            self.screen.blit(self.font.render(text, True, (210, 210, 215)), (self.view_w+8, 8+i*20))
        for i, key in enumerate(self.keys):
            y = 218+i*20
            selected = i == self.selected
            if selected:
                pg.draw.rect(self.screen, (62, 65, 75), (self.view_w+4, y-2, self.panel_w-8, 19))
            value = self.session.params[key]
            shown = self.edit_text if selected and self.edit_text is not None else ("auto" if value is None else "{:g}".format(value))
            self.screen.blit(self.font.render("{}: {}".format(key.replace("_", "-"), shown), True,
                                             (255, 230, 150) if selected else (205, 205, 210)), (self.view_w+10, y))
        button_y = 218+len(self.keys)*20+10
        self.apply_rect = pg.Rect(self.view_w+10, button_y, self.panel_w-20, 32)
        pg.draw.rect(self.screen, (55, 90, 120), self.apply_rect)
        self.screen.blit(self.font.render("APPLY" if self.future is None else "COMPUTING...", True,
                                         (255, 255, 255)), (self.view_w+22, button_y+7))
        self.screen.blit(self.font.render("Click value; type or Left/Right", True, (175, 175, 175)),
                         (self.view_w+6, button_y+42))
        # Wrap messages instead of hiding the part that explains a failure.
        columns = max(1, (self.panel_w-14)//self.font.size("M")[0])
        for i, line in enumerate(textwrap.wrap(self.message[:240], width=columns,
                                                break_on_hyphens=False)):
            self.screen.blit(self.font.render(line, True, (230, 185, 120)),
                             (self.view_w+7, button_y+70+i*17))
        pg.draw.rect(self.screen, (0, 0, 0), (0, self.view_h, self.width, self.status_h))
        status = "{}  zoom={:.3f}x  MTF={:.4f}  brush diameter={} px  mask={:.2f}%".format(
            self.view, self.zoom, self.mtf, self.brush, 100*np.mean(self.session.raw_mask))
        if self.future is not None:
            status = self.session.progress.status()
        self.screen.blit(self.font.render(status, True, (220, 220, 220)), (6, self.view_h+5))
        pg.display.flip()
        self.need_render = False

    def event(self, event):
        pg = self.pg
        busy = self.future is not None
        if event.type == pg.QUIT:
            self.closing = True
            if not busy:
                self.start_apply()
        elif event.type == pg.VIDEORESIZE:
            self.width, self.height = max(640, event.w), max(640, event.h)
            self.screen = pg.display.set_mode((self.width, self.height), pg.RESIZABLE)
        elif event.type == pg.KEYUP:
            self.held.discard(event.key)
        elif event.type == pg.MOUSEWHEEL:
            step = 10 if pg.key.get_mods() & pg.KMOD_SHIFT else 1
            self.brush = max(1, min(max(self.w, self.h)*2, self.brush+event.y*step))
        elif event.type == pg.MOUSEBUTTONDOWN:
            x, y = event.pos
            if not busy and event.button == 1 and self.apply_rect.collidepoint(event.pos):
                self.start_apply()
            elif not busy and event.button == 1 and x >= self.view_w:
                index = (y-218)//20
                if 0 <= index < len(self.keys):
                    self.commit_edit()
                    self.selected = index
                    self.edit_text = ""
            elif x < self.view_w and y < self.view_h:
                if event.button == 2:
                    self.pan_start = event.pos
                elif not busy and self.mode == "mask" and event.button in (1, 3):
                    self.stroke = np.packbits(self.session.raw_mask)
                    self.paint_button = event.button
                    self.paint(event.pos, erase=event.button == 3)
        elif event.type == pg.MOUSEBUTTONUP:
            if event.button == self.paint_button:
                self.finish_stroke()
            if event.button == 2:
                self.pan_start = None
        elif event.type == pg.MOUSEMOTION:
            if self.pan_start is not None:
                self.cx -= event.rel[0]/self.zoom
                self.cy -= event.rel[1]/self.zoom
            elif self.paint_button is not None and not busy:
                self.paint(event.pos, erase=self.paint_button == 3)
        elif event.type == pg.KEYDOWN:
            repeatable = (pg.K_LEFT, pg.K_RIGHT, pg.K_UP, pg.K_DOWN, pg.K_BACKSPACE)
            if event.key in self.held and event.key not in repeatable:
                return
            self.held.add(event.key)
            if self.edit_text is not None and not busy:
                if event.key in (pg.K_RETURN, pg.K_KP_ENTER):
                    self.commit_edit()
                elif event.key == pg.K_ESCAPE:
                    self.edit_text = None
                elif event.key == pg.K_BACKSPACE:
                    self.edit_text = self.edit_text[:-1]
                elif event.unicode and event.unicode in "0123456789.eE+-auto":
                    self.edit_text += event.unicode
                self.need_render = True
                return
            if event.key in (pg.K_q, pg.K_ESCAPE):
                self.closing = True
                if not busy:
                    self.start_apply()
            elif event.key in (pg.K_HOME, pg.K_END):
                self.mtf = min(0.95, max(0.001, self.mtf*(1/1.4 if event.key == pg.K_HOME else 1.4)))
            elif event.key in (pg.K_PLUS, pg.K_EQUALS, pg.K_KP_PLUS, pg.K_MINUS, pg.K_KP_MINUS):
                factor = 1/1.25 if event.key in (pg.K_MINUS, pg.K_KP_MINUS) else 1.25
                self.set_zoom(self.zoom*factor, pg.mouse.get_pos())
            elif event.key in (pg.K_1, pg.K_2, pg.K_3, pg.K_4):
                self.view = {pg.K_1: "original", pg.K_2: "prepared", pg.K_3: "background", pg.K_4: "result"}[event.key]
            elif event.key == pg.K_m:
                self.show_mask = not self.show_mask
            elif not busy:
                if event.key == pg.K_TAB:
                    self.finish_stroke()
                    self.mode = "parameters" if self.mode == "mask" else "mask"
                elif event.key in (pg.K_b, pg.K_RETURN, pg.K_KP_ENTER):
                    self.mode = "parameters"
                    self.start_apply()
                elif event.key == pg.K_s:
                    self.session.archive()
                    self.message = "Mask saved: " + ascii_text(self.session.mask_path)
                elif event.key == pg.K_z and event.mod & pg.KMOD_CTRL:
                    self.undo_stroke()
                elif event.key in (pg.K_UP, pg.K_DOWN):
                    self.selected = (self.selected + (1 if event.key == pg.K_DOWN else -1)) % len(self.keys)
                elif event.key in (pg.K_LEFT, pg.K_RIGHT):
                    self.adjust_value(1 if event.key == pg.K_RIGHT else -1, bool(event.mod & pg.KMOD_SHIFT))
        self.need_render = True

    def run(self):
        clock = self.pg.time.Clock()
        self.render()
        try:
            while True:
                for event in self.pg.event.get():
                    try:
                        self.event(event)
                    except (ValueError, OSError) as exc:
                        self.message = ascii_text(exc)
                        self.need_render = True
                if self.future is not None and self.future.done():
                    future, self.future = self.future, None
                    try:
                        stats = future.result()
                        self.message = "Applied in {:.2f} s. Mask archived.".format(stats["elapsed"])
                        self.view = "result"
                        if self.closing:
                            return
                    except Exception as exc:
                        report("Background calculation failed: " + str(exc), error=True)
                        self.message = ascii_text(exc)
                        self.closing = False
                    self.need_render = True
                if self.future is not None:
                    message = self.session.progress.status()
                    if message != self.message:
                        self.message = message
                        self.need_render = True
                if self.need_render:
                    self.render()
                clock.tick(60)
        finally:
            self.pool.shutdown(wait=True)
            self.pg.quit()


def usage(stream=None):
    stream = sys.stdout if stream is None else stream
    stream.write(
        "backflat - Masked diffusion background flattening for RGB FITS\n\n"
        "Usage:\n"
        "  backflat.py input_spec output_spec --starless SPEC [options]\n"
        "  backflat.py input_spec output_spec --sxt [options]\n"
        "  backflat.py input_spec output_spec --starnet [options]\n\n"
        "Input: RGB primary FITS (3,H,W); single, =single, sequence, wildcard, @list.\n"
        "Output: output_spec is the corrected float32 FITS; negative values retained.\n"
        "The corrected image, background model and uint16 mask are always saved.\n"
        "Correction: out[c] = input[c] - background[c] + K.\n"
        "K = (mean(background[R]) + mean(background[G]) + mean(background[B])) / 3;\n"
        "the same neutral level is added to all three channels.\n"
        "Gaussian diameters are FWHM.\n"
        "Pixel diameters are specified at diagonal 7515 px and scale with the image.\n\n"
        "  --starless SPEC      Matching starless FITS, or matching batch input spec\n"
        "  --sxt                Use an installed, licensed RC-Astro CLI\n"
        "  --starnet            Use installed StarNet2 with FITS/linear support (2.6+)\n"
        "  --sxt-exe FILE       RC-Astro executable (otherwise search PATH)\n"
        "  --starnet-exe FILE   StarNet2 executable (otherwise search PATH)\n"
        "  --cache-dir DIR      Temporary starless cache (default: .backflat-cache\n"
        "                       beside the mask output); cleaned when the run ends\n"
        "  --mask SPEC          Input mono mask; positive values exclude pixels\n"
        "  --no-gui             Require an existing mask, calculate and save\n"
        "  --out-back SPEC      Background output name (default: background.fit)\n"
        "  --out-mask SPEC      Mask output name (default: back_mask.fit)\n"
        "  -y, --overwrite      Allow replacing existing output files\n"
        "  --median1 D          Circular median diameter (40)\n"
        "  --edge F             Mirror-strip width as diagonal fraction (0.025)\n"
        "  --blur1 D            Preparation Gaussian FWHM (40)\n"
        "  --median2 D          Second circular median diameter (60)\n"
        "  --median-mode MODE   fast (DIPlib, default) or exact (SciPy reference)\n"
        "                       Both use float64 circular medians without quantization\n"
        "  --lake-mask D        Lake-mask Gaussian FWHM (100), then auto-level to [0,1]\n"
        "  --lake-blur D        Filled-lake Gaussian FWHM (250)\n"
        "  --blur D             Global Gaussian FWHM (60)\n"
        "  --center-plateau F   Distance threshold, fraction of max edge distance (0.5)\n"
        "  --center-mask D      Center-mask Gaussian FWHM (150, must be >0)\n"
        "  --center-blur D      Central Gaussian FWHM (225, must be >0)\n"
        "  --margin D           Mask margin, reference pixels (default: 1% diagonal)\n"
        "  --diffusion-scale N  Maximum spatial reduction: 1 (full), 2 or 4 (default)\n"
        "  --median2-scale N    Maximum spatial reduction after blur: 1, 2 or 4 (default)\n"
        "                       Small kernels/images reduce less; exact medians stay full-size\n"
        "  -h, --help           Show this help\n\n"
        "Default background and mask paths: beside a single corrected output;\n"
        "for a batch, in <output_stem>_backflat/ for each corrected output.\n"
        "--out-back / --out-mask override these paths (batch: numbered patterns\n"
        "or directories). Explicit relative paths are relative to the working directory.\n"
        "The mask retains raw strokes. An existing mask at its output path is\n"
        "reused automatically; --mask overrides the source, not the destination.\n"
        "Existing outputs stop the run before processing unless --overwrite (-y) is set.\n"
        "This includes an existing mask archive being reopened for editing.\n"
        "No automatic object segmentation.\n\n"
        "Cache: only external star removal is cached on disk, until this run ends.\n"
        "After saving, error or Ctrl+C: remove Backflat cache files and the empty\n"
        "directory. Other files/subdirectories and input/output images are retained.\n"
        "With --starless no cache is created; an old selected cache is still cleaned.\n"
        "The next --sxt/--starnet run removes stars again. Force-kill cannot clean up.\n\n"
        "Accuracy: fast uses spatial reduction after blur; intensity is never quantized.\n"
        "Use --median2-scale 1 --diffusion-scale 1 for no spatial reduction.\n"
        "--median-mode exact selects SciPy; diffusion scale is still separate.\n"
        "--starless must match the input orientation; --sxt corrects the CLI Y flip.\n\n"
        "GUI: LMB paint, RMB erase; wheel brush diameter (+/-1, Shift +/-10).\n"
        "+/- zoom at cursor; middle drag pan; Home/End MTF; M toggle mask.\n"
        "Ctrl+Z undo stroke; Tab mask/parameters; 1/2/3/4 original/prepared/bg/result.\n"
        "Click parameter to type; Up/Down select; Left/Right adjust (Shift fine).\n"
        "B/Enter/Apply recompute; S save mask; Q/Esc/close apply, save and exit.\n"
        "Examples:\n"
        "  backflat input.fit corrected.fit --starless starless.fit\n"
        "  backflat input.fit corrected.fit --sxt --out-back sky.fit --out-mask objects.fit\n"
        "  backflat input.fit corrected.fit --starless starless.fit --no-gui -y\n"
    )


def parse_args(argv):
    if any(arg in ("-h", "--help") for arg in argv[1:]):
        usage()
        return None
    params = dict(DEFAULTS)
    options = {"starless": None, "mask": None, "out_back": None, "out_mask": None,
               "sxt_exe": None,
               "starnet_exe": None, "cache_dir": None, "sxt": False,
               "starnet": False, "no_gui": False, "overwrite": False, "median_mode": "fast"}
    positional = []
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--":
            positional.extend(argv[i+1:])
            break
        if arg in ("--sxt", "--starnet", "--no-gui", "--overwrite", "-y"):
            options["overwrite" if arg == "-y" else arg[2:].replace("-", "_")] = True
            i += 1
            continue
        if arg.startswith("--"):
            key = arg[2:].replace("-", "_")
            if key not in params and key not in options:
                raise ValueError("Unknown option: " + arg)
            if i+1 >= len(argv):
                raise ValueError("Missing value for " + arg)
            if key in params:
                try:
                    params[key] = float(argv[i+1])
                except ValueError:
                    raise ValueError("Expected a number after " + arg)
            else:
                options[key] = argv[i+1]
            i += 2
            continue
        positional.append(arg)
        i += 1
    if len(positional) != 2:
        usage(sys.stderr)
        raise SystemExit(1)
    if sum((bool(options["starless"]), options["sxt"], options["starnet"])) != 1:
        raise ValueError("Choose exactly one of --starless, --sxt or --starnet.")
    validate_params(params)
    if options["median_mode"] not in ("fast", "exact"):
        raise ValueError("--median-mode must be fast or exact.")
    params["median_method"] = options["median_mode"]
    return positional[0], positional[1], params, options


def matching_files(spec, count, label, broadcast=False):
    if spec is None:
        return [None]*count
    files = batch_utils.expand_input_spec(spec)
    if broadcast and len(files) == 1:
        return files*count
    if len(files) != count:
        raise ValueError("{} requires {} file(s); got {}.".format(label, count, len(files)))
    return files


def work_items(input_spec, output_spec, options):
    """Use shared expansion rules for all three mandatory output files."""
    inputs = batch_utils.expand_input_spec(input_spec)
    pairs = batch_utils.build_io_file_lists_from_list(inputs, output_spec)
    count = len(pairs)
    stars = matching_files(options["starless"], count, "--starless")
    masks = matching_files(options["mask"], count, "--mask", broadcast=True)
    backgrounds = ([p[1] for p in batch_utils.build_io_file_lists_from_list(inputs, options["out_back"])]
                   if options["out_back"] else [None]*count)
    archives = ([p[1] for p in batch_utils.build_io_file_lists_from_list(inputs, options["out_mask"])]
                if options["out_mask"] else [None]*count)
    items = []
    destinations = []
    for i, (source, output) in enumerate(pairs):
        output = os.path.abspath(output)
        folder = os.path.dirname(output)
        if count > 1:
            folder = os.path.join(folder, os.path.splitext(os.path.basename(output))[0]+"_backflat")
        archive = os.path.abspath(archives[i] or os.path.join(folder, "back_mask.fit"))
        background = os.path.abspath(backgrounds[i] or os.path.join(folder, "background.fit"))
        cache = (os.path.abspath(options["cache_dir"]) if options["cache_dir"] else
                 os.path.join(os.path.dirname(archive), ".backflat-cache"))
        items.append({"input": source, "output": output, "starless": stars[i],
                      "mask": masks[i] or (archive if os.path.isfile(archive) else None),
                      "archive": archive, "cache": cache, "bg": background})
        destinations.extend([output, archive, background])
    canonical = lambda p: os.path.normcase(os.path.realpath(p))
    targets = [canonical(p) for p in destinations]
    if len(set(targets)) != len(targets):
        raise ValueError("Output, background and mask archive paths must not collide.")
    source_paths = {canonical(p) for p in inputs + [s for s in stars if s]}
    for item in items:
        for path in [item["output"], item["archive"], item["bg"]]:
            if canonical(path) in source_paths:
                raise ValueError("An output would overwrite an input image: " + path)
        for path in [item["output"], item["bg"]]:
            if any(canonical(path) == canonical(mask) for mask in masks if mask):
                raise ValueError("An image output would overwrite the source mask: " + path)
        # Updating this image's own source mask is allowed; never overwrite a
        # mask that a later batch item still needs to read.
        for j, mask in enumerate(masks):
            if mask and canonical(item["archive"]) == canonical(mask) and items[j] is not item:
                raise ValueError("A mask output would overwrite another image's source mask: " + mask)
    existing = []
    for item in items:
        for key, label in (("output", "corrected"), ("bg", "background"), ("archive", "mask")):
            path = item[key]
            if os.path.isdir(path):
                raise ValueError("Output path is a directory: " + path)
            if not options["overwrite"] and os.path.lexists(path):
                existing.append("  {}: {}".format(label, path))
    if existing:
        raise FileExistsError("Output files already exist:\n" + "\n".join(existing) +
                              "\nUse --overwrite (-y) to replace existing outputs, or choose new output paths.")
    for item in items:
        if options["no_gui"] and item["mask"] is None:
            raise ValueError("--no-gui requires --mask FILE or an existing archive: " + item["archive"])
    return items


def output_header(header, params, item, engine_info, stats, background=False):
    result = clean_header(header)
    result.add_history("backflat: Masked Diffusion Background; Gaussian diameters=FWHM")
    result.add_history("backflat: reference diagonal=7515; circular median; float64 work")
    result.add_history("backflat: diffusion mask/image FWHM=20/80; alpha densification=2")
    result.add_history("backflat: lake mask Gaussian then linear min/max auto-level to [0,1]")
    for key, value in params.items():
        result.add_history(ascii_text("backflat: {}={}".format(key, "auto(1% diagonal)" if value is None else value)))
    result.add_history(ascii_text("backflat: mask=" + os.path.abspath(item["archive"])))
    result.add_history(ascii_text("backflat: starless=" + engine_info))
    result.add_history("backflat: iterations={}; means={}".format(
        stats["iterations"], ",".join("{:.12g}".format(v) for v in stats["channel_means"])))
    result.add_history("backflat: K=mean(Rmean,Gmean,Bmean)={:.12g}; common RGB level".format(
        stats["neutral_level"]))
    result.add_history("backflat: actual diffusion reduction={}".format(stats["diffusion_scale"]))
    result.add_history("backflat: actual median2 reduction={}".format(stats["median2_scale"]))
    result.add_history("backflat: median backend={}; no intensity quantization".format(
        stats["median_backend"]))
    result.add_history("backflat: " + ("background model, float32" if background else
                                      "out[c]=image[c]-background[c]+K; float32"))
    return result


def process_item(item, params, options):
    progress = OperationLog()
    report("Corrected output: " + item["output"])
    report("Background output: " + item["bg"])
    report("Mask output: " + item["archive"])
    with operation(progress, "Read input FITS"):
        image, header, _ = read_rgb(item["input"])
    if item["starless"]:
        with operation(progress, "Read starless FITS"):
            starless, _, _ = read_rgb(item["starless"])
        engine_info = os.path.abspath(item["starless"])
    else:
        engine = "sxt" if options["sxt"] else "starnet"
        with operation(progress, "Star removal or cache"):
            starless, engine_info = cached_starless(item["input"], image, header, engine,
                                                   options[engine+"_exe"], item["cache"])
    if starless.shape != image.shape:
        raise ValueError("The starless image must have the same shape and registration as the input.")
    raw = None
    effective_params = dict(params)
    if item["mask"]:
        raw, saved_margin, has_raw = read_mask(item["mask"], image.shape[1:])
        if has_raw and params["margin"] is None:
            effective_params["margin"] = saved_margin*REFERENCE_DIAGONAL/math.hypot(*raw.shape)
        report("Loaded mask: " + item["mask"])
    session = Session(image, starless, header, effective_params, item["archive"], raw, progress,
                      overwrite=options["overwrite"])
    report("Preparing background reference...")
    if options["no_gui"]:
        session.archive()
        session.apply()
    else:
        session.apply(model=False)
        MaskEditor(session, loaded_mask=raw is not None).run()
    # Archive even an empty mask: it records the explicit choice to exclude nothing.
    session.archive()
    report("Saved mask: " + item["archive"])
    hdr = output_header(header, session.params, item, engine_info, session.stats, background=True)
    with operation(progress, "Write background FITS"):
        atomic_fits(item["bg"], finite_float(session.fmf), hdr, overwrite=options["overwrite"])
    report("Saved background: " + item["bg"])
    hdr = output_header(header, session.params, item, engine_info, session.stats)
    with operation(progress, "Write corrected FITS"):
        atomic_fits(item["output"], finite_float(session.result), hdr, overwrite=options["overwrite"])
    report("Saved {} (mask {:.2f}%, diffusion {} iterations, compute {:.2f} s)".format(
        item["output"], 100*session.stats["masked_fraction"], session.stats["iterations"], session.stats["elapsed"]))


def main(argv=None):
    argv = sys.argv if argv is None else argv
    cache_dirs = {}
    protected_paths = set()
    exit_code = 0
    try:
        parsed = parse_args(argv)
        if parsed is None:
            return 0
        input_spec, output_spec, params, options = parsed
        items = work_items(input_spec, output_spec, options)
        protected_paths = {os.path.normcase(os.path.realpath(item[key]))
                           for item in items
                           for key in ("input", "starless", "mask", "output", "bg", "archive")
                           if item[key]}
        failed = 0
        for i, item in enumerate(items, 1):
            cache_dir = os.path.abspath(item["cache"])
            cache_dirs.setdefault(cache_dir, os.path.normcase(os.path.realpath(cache_dir)))
            report("[{}/{}] {}".format(i, len(items), item["input"]))
            try:
                process_item(item, params, options)
            except Exception as exc:
                failed += 1
                report("Error: " + str(exc), error=True)
        exit_code = 1 if failed else 0
    except (ValueError, OSError) as exc:
        report("Error: " + str(exc), error=True)
        exit_code = 1
    except KeyboardInterrupt:
        report("Interrupted. Existing mask archives are retained.", error=True)
        exit_code = 130
    finally:
        if not cleanup_caches(cache_dirs, protected_paths) and exit_code == 0:
            exit_code = 1
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
