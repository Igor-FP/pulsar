"""Numerical and FITS-contract tests for backflat; no external models needed."""

import contextlib
import io
import os
import tempfile
import unittest
from unittest import mock
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy import ndimage

import backflat


class BackgroundContracts(unittest.TestCase):
    def test_gaussian_fft_matches_spatial_filter(self):
        image = np.random.default_rng(21).normal(1000, 200, (257, 263))
        actual = backflat.gaussian_plane(image, 8.25)
        expected = ndimage.gaussian_filter(image, 8.25, mode="reflect", truncate=4.0)
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=1e-12)

    def test_fwhm_not_radius_or_kernel_size(self):
        image = np.zeros((301, 301))
        image[150, 150] = 1
        diameter = 20*backflat.REFERENCE_DIAGONAL/np.hypot(301, 301)
        blurred = backflat.gaussian(image, diameter)
        self.assertAlmostEqual(blurred[150, 160]/blurred[150, 150], 0.5, places=10)

    def test_fast_circular_median_matches_reference_exactly(self):
        rng = np.random.default_rng(99)
        image = rng.normal(1000, 300, (3, 283, 53))
        image[0] -= 1800
        diameter = 8*backflat.REFERENCE_DIAGONAL/np.hypot(283, 53)
        exact = backflat.median_aperture(image, diameter, "exact")
        fast = backflat.median_aperture(image, diameter, "fast")
        np.testing.assert_array_equal(fast, exact)

    def test_faint_background_is_preserved_without_quantization(self):
        y, x = np.mgrid[:37, :43]
        background = 1000+0.1*x/43+0.03*np.sin(y)
        starless = np.stack([background, background+0.01, background-0.02])
        original = starless.copy()
        original[:, 18, 21] = 1e8
        params = dict(backflat.DEFAULTS)
        params.update(median1=8*backflat.REFERENCE_DIAGONAL/np.hypot(37, 43),
                      median2=0, blur1=0, edge=0, median_method="fast")
        with_stars, stats = backflat.prepare(original, starless, params, return_stats=True)
        without_stars = backflat.prepare(starless, starless, params)
        np.testing.assert_array_equal(with_stars, without_stars)
        exact = backflat.median_aperture(starless, params["median1"], "exact")
        np.testing.assert_array_equal(with_stars, exact)
        self.assertEqual(stats["median_backend"], "DIPlib")

    def test_large_and_fractional_disks_preserve_float64_samples(self):
        y, x = np.mgrid[:75, :83]
        # Sub-float32 differences would disappear if the backend downcast data.
        plane = 1e6 + 1e-3*(x/83 + 0.3*np.sin(y))
        image = np.stack([plane, -plane, plane + 0.0002])
        for radius in (1.0, 10.0, 30.0, 30.35):
            with self.subTest(radius=radius):
                diameter = radius*2*backflat.REFERENCE_DIAGONAL/np.hypot(75, 83)
                actual = backflat.median_aperture(image, diameter, "fast")
                expected = backflat.median_aperture(image, diameter, "exact")
                np.testing.assert_array_equal(actual, expected)

    def test_diffusion_ignores_masked_signal_and_preserves_shores(self):
        y, x = np.mgrid[:75, :101]
        plane = 20+0.2*x+0.3*y
        image = np.stack([plane, 2*plane+10, 0.5*plane-80])
        holes = ((x-48)**2+(y-35)**2 < 20**2) | ((x-3)**2+(y-60)**2 < 10**2)
        changed = image.copy()
        changed[:, holes] = 1e6
        first, stats = backflat.diffuse_holes(image, holes)
        second, _ = backflat.diffuse_holes(changed, holes)
        np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(first[:, ~holes], image[:, ~holes])
        self.assertTrue(np.isfinite(first).all())
        self.assertGreater(stats["iterations"], 1)
        # The same spatial operator acts on each channel, including additive levels.
        np.testing.assert_allclose(first[1], 2*first[0]+10, atol=1e-10)
        np.testing.assert_allclose(first[2], 0.5*first[0]-80, atol=1e-10)

    def test_diffusion_constant_and_empty_mask(self):
        image = np.full((3, 61, 81), -1200.75)
        mask = np.zeros((61, 81), dtype=bool)
        out, stats = backflat.diffuse_holes(image, mask)
        np.testing.assert_array_equal(out, image)
        self.assertEqual(stats["iterations"], 0)
        mask[10:40, 10:70] = True
        out, _ = backflat.diffuse_holes(image, mask)
        np.testing.assert_allclose(out, image, atol=1e-9)
        with self.assertRaises(ValueError):
            backflat.diffuse_holes(image, np.ones_like(mask))

    def test_reduced_diffusion_keeps_shores_and_excludes_masked_values(self):
        y, x = np.mgrid[:157, :213]
        plane = 1000 + .001*x + .002*y
        image = np.stack([plane, 2*plane+5, -plane+3])
        # Odd dimensions, a border-touching hole, and a one-pixel hole exercise
        # partial blocks and reconstruction of small masked structures.
        mask = ((x-108)**2+(y-78)**2 < 40**2) | ((x-3)**2+(y-7)**2 < 12**2)
        mask[133, 187] = True
        altered = image.copy()
        altered[:, mask] = 1e10
        first, stats = backflat.diffuse_holes(image, mask, scale=4, sigma_mask=8.5)
        second, _ = backflat.diffuse_holes(altered, mask, scale=4, sigma_mask=8.5)
        self.assertEqual(stats["diffusion_scale"], 4)
        np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(first[:, ~mask], image[:, ~mask])
        np.testing.assert_allclose(first[1], 2*first[0]+5, atol=1e-9, rtol=0)
        np.testing.assert_allclose(first[2], -first[0]+3, atol=1e-9, rtol=0)
        constant = np.full_like(image, -1200.125)
        actual, _ = backflat.diffuse_holes(constant, mask, scale=4, sigma_mask=8.5)
        np.testing.assert_allclose(actual, constant, atol=1e-9, rtol=0)

    def test_spatial_reduction_preserves_pixel_centers_and_float64_precision(self):
        y, x = np.mgrid[:99, :127]
        plane = 1e6 + 1e-4*x + 3e-4*y
        coarse = backflat.area_reduce(plane, 4)
        restored = backflat.restore_plane(coarse, 4, plane.shape)
        self.assertEqual(restored.dtype, np.dtype('float64'))
        # Reflected partial edge blocks deliberately differ from extrapolating
        # a plane beyond the image. Interior samples must retain exact centers.
        np.testing.assert_allclose(restored[8:-8, 8:-8], plane[8:-8, 8:-8],
                                   atol=1e-9, rtol=0)

    def test_median_reduction_requires_preblur_and_exact_mode_stays_full_size(self):
        image = np.empty((3, 73, 91), dtype=np.float64)
        image[:] = np.array([-1200.125, 1e6+.001, 7000.75])[:, None, None]
        unit = backflat.REFERENCE_DIAGONAL/np.hypot(73, 91)
        params = dict(backflat.DEFAULTS, median1=0, edge=0, blur1=40*unit,
                      median2=60*unit, median2_scale=4, median_method="fast")
        actual, stats = backflat.prepare(image, image, params, return_stats=True)
        self.assertEqual(stats["median2_scale"], 4)
        np.testing.assert_allclose(actual, image, atol=1e-8, rtol=0)
        for changes in ({"blur1": 0}, {"median_method": "exact"}):
            _, stats = backflat.prepare(image, image, dict(params, **changes), return_stats=True)
            self.assertEqual(stats["median2_scale"], 1)

    def test_center_mask_protects_all_image_edges(self):
        mask = backflat.center_mask((111, 167), 0.5, 150)
        self.assertEqual(float(mask[0].max()), 0)
        self.assertEqual(float(mask[-1].max()), 0)
        self.assertEqual(float(mask[:, 0].max()), 0)
        self.assertEqual(float(mask[:, -1].max()), 0)
        self.assertGreater(mask[55, 83], 0.99)
        self.assertGreater(mask[55, 50], mask[10, 50])

    def test_blurred_lake_mask_uses_full_opacity_range(self):
        mask = np.zeros((61, 81), dtype=bool)
        mask[27:34, 37:44] = True
        diameter = 30*backflat.REFERENCE_DIAGONAL/np.hypot(61, 81)
        blurred = backflat.gaussian(mask.astype(np.float64), diameter)
        self.assertLess(blurred.max(), 0.1)
        actual = backflat.lake_blending_mask(mask, diameter)
        self.assertEqual(float(actual.min()), 0.0)
        self.assertEqual(float(actual.max()), 1.0)
        np.testing.assert_allclose(actual, (blurred-blurred.min())/np.ptp(blurred))
        for value in (False, True):
            constant = np.full(mask.shape, value)
            np.testing.assert_array_equal(backflat.lake_blending_mask(constant, diameter),
                                          constant.astype(np.float64))

    def test_full_pipeline_neutralizes_constant_rgb_background(self):
        image = np.empty((3, 81, 111))
        image[:] = np.array([-1500, 2800, 150000])[:, None, None]
        params = dict(backflat.DEFAULTS)
        raw_mask = np.zeros(image.shape[1:], dtype=bool)
        raw_mask[25:60, 30:85] = True
        prepared = backflat.prepare(image, image, params)
        fmf = backflat.build_fmf(prepared, raw_mask, params)
        output = backflat.subtract(image, fmf)
        np.testing.assert_allclose(fmf, image, atol=1e-8)
        np.testing.assert_allclose(output, np.mean([-1500, 2800, 150000]), atol=1e-8)

    def test_overall_mean_preserved_and_negative_values_retained(self):
        image = np.arange(3*8*11, dtype=np.float64).reshape(3, 8, 11)
        fmf = np.empty_like(image)
        fmf[0] = 3*image[0]
        fmf[1] = 7*image[1]+200
        fmf[2] = -2*image[2]-100
        output = backflat.subtract(image, fmf)
        self.assertAlmostEqual(output.mean(), image.mean())
        np.testing.assert_allclose(output + fmf - fmf.mean(), image)
        self.assertLess(output.min(), 0)

    def test_colored_background_removed_without_changing_object_color(self):
        y, x = np.mgrid[-3:4, -4:5]
        fmf = np.stack([100+2*x, 140+3*y, 100-x-y]).astype(np.float64)
        image = fmf.copy()
        signal = np.array([30.0, 60.0, 120.0])
        image[:, 3, 4] += signal
        output = backflat.subtract(image, fmf)
        sky = np.ones(x.shape, dtype=bool)
        sky[3, 4] = False
        level = 340.0/3
        np.testing.assert_allclose(output[:, sky], level, atol=1e-12)
        np.testing.assert_allclose(output[:, 3, 4]-level, signal, atol=1e-12)

    def test_mirror_corners_use_both_reflections(self):
        plane = np.arange(8*10).reshape(8, 10)
        image = np.stack([plane, plane+100, plane+200]).astype(np.float64)
        result = backflat.mirror_edges(image, 2/np.hypot(8, 10))
        ys = [3, 2, 2, 3, 4, 5, 5, 4]
        xs = [3, 2, 2, 3, 4, 5, 6, 7, 7, 6]
        np.testing.assert_array_equal(result, image[:, ys][:, :, xs])
        np.testing.assert_array_equal(image[0], plane)

    def test_mask_archive_does_not_accumulate_margin(self):
        raw = np.zeros((31, 41), dtype=bool)
        raw[15, 20] = True
        effective = backflat.expand_mask(raw, 3)
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "back_mask.fit")
            backflat.save_mask(path, raw, effective, fits.Header(), 3)
            loaded, margin, has_raw = backflat.read_mask(path, raw.shape)
            self.assertTrue(has_raw)
            self.assertEqual(margin, 3)
            np.testing.assert_array_equal(loaded, raw)
            np.testing.assert_array_equal(backflat.expand_mask(loaded, margin), effective)
            with fits.open(path, memmap=False) as hdus:
                self.assertEqual(hdus[0].data.dtype, np.dtype("uint16"))
                np.testing.assert_array_equal(hdus[0].data, effective.astype(np.uint16)*65535)
                self.assertEqual(hdus["RAWMASK"].data.dtype, np.dtype("uint16"))
                np.testing.assert_array_equal(hdus["RAWMASK"].data, raw.astype(np.uint16)*65535)
                for hdu in hdus:
                    self.assertEqual(hdu.header["BITPIX"], 16)
                    self.assertEqual(hdu.header["BZERO"], 32768)

    def test_fits_preserves_values_and_observation_metadata(self):
        header = fits.Header()
        header["EXPTIME"] = 300.0
        header["FILTER"] = "RGB"
        header["CRPIX1"] = 11.5
        image = np.linspace(-1200, 40000, 3*9*13).reshape(3, 9, 13)
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "output.fit")
            backflat.atomic_fits(path, backflat.finite_float(image), header)
            loaded, actual, _ = backflat.read_rgb(path)
            np.testing.assert_allclose(loaded, image, rtol=1e-6)
            self.assertEqual(actual["FILTER"], "RGB")
            self.assertEqual(actual["EXPTIME"], 300.0)
            self.assertEqual(actual["CRPIX1"], 11.5)

    def test_ascii_diagnostics_include_non_ascii_paths(self):
        message = "Error in " + chr(0x0444) + ".fit " + chr(0x2192) + " failed"
        sink = io.StringIO()
        with contextlib.redirect_stdout(sink):
            backflat.report(message)
        self.assertTrue(sink.getvalue().isascii())
        self.assertIn("\\u0444", sink.getvalue())

    def test_float32_overflow_is_sanitized_after_cast(self):
        source = np.array([-2.5, np.nan, np.inf, -np.inf, 1e200])
        with contextlib.redirect_stderr(io.StringIO()):
            output = backflat.finite_float(source)
        self.assertTrue(np.isfinite(output).all())
        np.testing.assert_array_equal(output, [-2.5, 0, 0, 0, 0])

    def test_display_block_mean_retains_partial_border_blocks(self):
        image = np.arange(3*5*7, dtype=np.float64).reshape(3, 5, 7)
        output = backflat.block_mean(image, 3)
        self.assertEqual(output.shape, (3, 2, 3))
        np.testing.assert_allclose(output[:, 1, 2], image[:, 3:, 6:].mean(axis=(1, 2)))
        np.testing.assert_allclose(output[:, 0, 0], image[:, :3, :3].mean(axis=(1, 2)))


class CommandLineContracts(unittest.TestCase):
    def test_headless_from_other_directory_with_ascii_redirect_and_unicode_path(self):
        script = str(Path(backflat.__file__).resolve())
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            name = chr(0x0444)+".fit"
            image = np.full((3, 51, 63), 3500, dtype=np.uint16)
            image[1] = 4200
            header = fits.Header()
            header["EXPTIME"] = 30
            fits.PrimaryHDU(image, header).writeto(folder/name)
            fits.PrimaryHDU(image).writeto(folder/"starless.fit")
            mask = np.zeros((51, 63), dtype=np.uint8)
            mask[15:25, 20:40] = 255
            fits.PrimaryHDU(mask).writeto(folder/"source_mask.fit")
            command = [sys.executable, "-B", script, name, "results/output.fit", "--starless", "starless.fit",
                       "--mask", "source_mask.fit", "--no-gui", "--starless-no-flip"]
            env = dict(os.environ, PYTHONIOENCODING="ascii:strict")
            with (folder/"log.txt").open("wb") as log:
                process = subprocess.run(command, cwd=directory, env=env, stdout=log, stderr=log)
            log = (folder/"log.txt").read_bytes()
            self.assertEqual(process.returncode, 0, log.decode("ascii"))
            self.assertTrue(log.isascii())
            self.assertIn(b"Read input FITS completed in", log)
            self.assertIn(b"Median 1", log)
            self.assertIn(b"Fill masked regions completed in", log)
            self.assertIn(b"Write corrected FITS completed in", log)
            for output in ["output.fit", "background.fit"]:
                with fits.open(folder/"results"/output, memmap=False) as hdus:
                    self.assertEqual(hdus[0].data.dtype.kind, "f")
                    self.assertEqual(hdus[0].data.dtype.itemsize, 4)
                    self.assertTrue(np.isfinite(hdus[0].data).all())
                    expected = image if output == "background.fit" else image.mean()
                    np.testing.assert_allclose(hdus[0].data, expected, atol=0.01)
                    self.assertEqual(hdus[0].header["EXPTIME"], 30)
                    self.assertIn("median_method=fast", str(hdus[0].header["HISTORY"]))
                    self.assertIn("K=mean(Rmean,Gmean,Bmean)", str(hdus[0].header["HISTORY"]))
                    if output == "output.fit":
                        self.assertIn("out[c]=image[c]-background[c]+K", str(hdus[0].header["HISTORY"]))
            with fits.open(folder/"results"/"back_mask.fit", memmap=False) as hdus:
                self.assertEqual(hdus[0].data.dtype, np.dtype("uint16"))
            self.assertIn(b"Saved background:", log)
            self.assertIn(b"Saved mask:", log)
            self.assertFalse((folder/"back_mask.fit").exists())

    def test_named_outputs_require_overwrite_and_reload_chosen_mask(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            image = np.full((3, 9, 13), 3500, dtype=np.uint16)
            fits.PrimaryHDU(image).writeto(folder/"input.fit")
            fits.PrimaryHDU(image).writeto(folder/"stars.fit")
            mask = np.zeros((9, 13), dtype=np.uint16)
            mask[4, 6] = 65535
            fits.PrimaryHDU(mask).writeto(folder/"source.fit")
            output = folder/"result"/"corrected.fit"
            background = folder/"models"/"sky.fit"
            archive = folder/"masks"/"objects.fit"
            args = ["backflat", str(folder/"input.fit"), str(output), "--starless", str(folder/"stars.fit"),
                    "--no-gui", "--out-back", str(background), "--out-mask", str(archive), "--starless-no-flip"]
            initial = args + ["--mask", str(folder/"source.fit")]
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(backflat.main(initial), 0)
            paths = (output, background, archive)
            saved = [path.read_bytes() for path in paths]
            self.assertFalse((output.parent/"back_mask.fit").exists())
            self.assertFalse((output.parent/"background.fit").exists())
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err), \
                 mock.patch.object(backflat, "process_item") as process:
                self.assertEqual(backflat.main(args), 1)
                process.assert_not_called()
            self.assertIn("Output files already exist", err.getvalue())
            self.assertIn("--overwrite", err.getvalue())
            for path in paths:
                self.assertIn(str(path), err.getvalue())
            self.assertEqual([path.read_bytes() for path in paths], saved)
            # Both spellings allow replacement and resume the selected archive.
            for i, flag in enumerate(("--overwrite", "-y")):
                with self.subTest(flag=flag):
                    fits.PrimaryHDU(image+100*(i+1)).writeto(folder/"input.fit", overwrite=True)
                    out = io.StringIO()
                    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
                        self.assertEqual(backflat.main(args+[flag]), 0)
                    self.assertIn("Loaded mask: " + str(archive), out.getvalue())
                    with fits.open(output, memmap=False) as hdus:
                        np.testing.assert_allclose(hdus[0].data, image+100*(i+1), atol=0.01)
                    with fits.open(archive, memmap=False) as hdus:
                        np.testing.assert_array_equal(hdus["RAWMASK"].data, mask)

    def test_any_existing_output_stops_before_processing(self):
        for name in ("corrected.fit", "background.fit", "back_mask.fit"):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                fits.PrimaryHDU(np.zeros((3, 5, 7))).writeto(folder/"input.fit")
                fits.PrimaryHDU(np.zeros((3, 5, 7))).writeto(folder/"stars.fit")
                existing = folder/name
                existing.write_bytes(b"keep this file")
                args = ["backflat", str(folder/"input.fit"), str(folder/"corrected.fit"),
                        "--starless", str(folder/"stars.fit")]
                error = io.StringIO()
                with contextlib.redirect_stderr(error), mock.patch.object(backflat, "process_item") as process:
                    self.assertEqual(backflat.main(args), 1)
                    process.assert_not_called()
                self.assertIn(str(existing), error.getvalue())
                self.assertEqual(existing.read_bytes(), b"keep this file")

    def test_batch_specs_keep_per_image_archives(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            for i in (1, 2):
                fits.PrimaryHDU(np.zeros((3, 5, 7))).writeto(folder/("input%04d.fit" % i))
                fits.PrimaryHDU(np.zeros((3, 5, 7))).writeto(folder/("stars%04d.fit" % i))
            parsed = backflat.parse_args(["backflat", str(folder/"input0001.fit"), str(folder/"out0001.fit"),
                                         "--starless", str(folder/"stars0001.fit")])
            items = backflat.work_items(parsed[0], parsed[1], parsed[3])
            self.assertEqual(len(items), 2)
            self.assertNotEqual(items[0]["archive"], items[1]["archive"])
            self.assertEqual(Path(items[0]["archive"]).name, "back_mask.fit")
            for i, item in enumerate(items, 1):
                expected = folder/("out%04d_backflat" % i)
                self.assertEqual(Path(item["archive"]), expected/"back_mask.fit")
                self.assertEqual(Path(item["bg"]), expected/"background.fit")
            # A collision in the second item must stop the entire batch up front.
            late = Path(items[1]["bg"])
            late.parent.mkdir()
            late.write_bytes(b"existing background")
            with self.assertRaisesRegex(FileExistsError, "background.fit"):
                backflat.work_items(parsed[0], parsed[1], parsed[3])
            parsed[3]["out_back"] = str(folder/"sky0001.fit")
            parsed[3]["out_mask"] = str(folder/"mask0001.fit")
            items = backflat.work_items(parsed[0], parsed[1], parsed[3])
            for i, item in enumerate(items, 1):
                self.assertEqual(Path(item["bg"]), folder/("sky%04d.fit" % i))
                self.assertEqual(Path(item["archive"]), folder/("mask%04d.fit" % i))
            # The overwrite flag cannot permit source/output or output/output collisions.
            parsed[3]["overwrite"] = True
            parsed[3]["out_back"] = str(folder/"input0001.fit")
            with self.assertRaisesRegex(ValueError, "overwrite an input image"):
                backflat.work_items(parsed[0], parsed[1], parsed[3])
            parsed[3]["out_back"] = parsed[3]["out_mask"]
            with self.assertRaisesRegex(ValueError, "must not collide"):
                backflat.work_items(parsed[0], parsed[1], parsed[3])

    def test_files_appearing_after_preflight_are_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"output.fit"
            path.write_bytes(b"created by another process")
            mask = np.zeros((5, 7), dtype=bool)
            with self.assertRaisesRegex(FileExistsError, "--overwrite"):
                backflat.atomic_fits(str(path), np.zeros((3, 5, 7)), fits.Header())
            with self.assertRaisesRegex(FileExistsError, "--overwrite"):
                backflat.save_mask(str(path), mask, mask, fits.Header(), 0)
            self.assertEqual(path.read_bytes(), b"created by another process")
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_cache_reuse_does_not_require_license_probe(self):
        with tempfile.TemporaryDirectory() as directory:
            image = np.arange(3*8*9, dtype=np.float64).reshape(3, 8, 9)/256
            source = os.path.join(directory, "input.fit")
            executable = os.path.join(directory, "engine.exe")
            Path(executable).write_bytes(b"test engine identity")
            fits.PrimaryHDU(image).writeto(source)
            def execute(command, **kwargs):
                output = command[command.index("-o")+1]
                fits.PrimaryHDU(image[:, ::-1, :]*0.5).writeto(output)
                return ""
            cache = os.path.join(directory, "cache")
            # Both older cache conventions must be ignored by the raw-row cache.
            os.makedirs(cache)
            binary = os.stat(executable)
            digest = backflat.hashlib.sha256(Path(source).read_bytes())
            identity = "{}:{}:{}:{}:backflat-starless-v2".format(
                "sxt", os.path.realpath(executable), binary.st_size, binary.st_mtime_ns)
            for suffix in ("", ":sxt-flip-y-v1"):
                old_digest = digest.copy()
                old_digest.update((identity+suffix).encode("utf-8"))
                fits.PrimaryHDU(np.full_like(image, -123)).writeto(
                    os.path.join(cache, old_digest.hexdigest()+".fit"))
            with mock.patch.object(backflat.shutil, "which", return_value=executable), \
                 mock.patch.object(backflat, "discover_engine", return_value=(executable, "test:1")) as probe, \
                 mock.patch.object(backflat, "external_run", side_effect=execute) as run, \
                 contextlib.redirect_stdout(io.StringIO()):
                first, _ = backflat.cached_starless(source, image, fits.Header(), "sxt", None, cache)
                probe.side_effect = RuntimeError("offline")
                second, _ = backflat.cached_starless(source, image, fits.Header(), "sxt", None, cache)
                np.testing.assert_array_equal(first, image[:, ::-1, :]*0.5)
                np.testing.assert_array_equal(first, second)
                self.assertEqual(run.call_count, 1)
                self.assertEqual(probe.call_count, 1)


class StarlessOrientationContracts(unittest.TestCase):
    @staticmethod
    def pair(shape=(257, 385), seed=42):
        rng = np.random.default_rng(seed)
        y, x = np.mgrid[:shape[0], :shape[1]]
        texture = 150*ndimage.gaussian_filter(rng.normal(size=shape), 3)
        sky = 2000+0.5*x+0.3*y+texture
        starless = np.stack([sky, sky*1.2+30, sky*0.8-20])
        stars = np.zeros(shape)
        stars[rng.integers(0, shape[0], 150), rng.integers(0, shape[1], 150)] = 15000
        stars = ndimage.gaussian_filter(stars, 0.9)
        return starless+stars, starless

    def test_known_pair_both_directions_and_intensity_scaling(self):
        image, starless = self.pair()
        for flip in (False, True):
            for gain, offset in ((1, 0), (0.03, 1200), (7, -500)):
                with self.subTest(flip=flip, gain=gain):
                    supplied = (starless[:, ::-1, :] if flip else starless)*gain+offset
                    before = supplied.copy()
                    actual, decision = backflat.orient_starless(image, supplied)
                    self.assertEqual(decision["flip_y"], flip)
                    np.testing.assert_array_equal(actual, starless*gain+offset)
                    np.testing.assert_array_equal(supplied, before)
                    if not flip:
                        self.assertIs(actual, supplied)

    def test_odd_size_reduction_commutes_with_flip_and_removes_plane(self):
        image, _ = self.pair((1033, 1247))
        a = backflat.orientation_preview(image)
        b = backflat.orientation_preview(image[:, ::-1, :])
        self.assertLessEqual(max(a.shape), 1024)
        np.testing.assert_allclose(a[::-1], b, atol=1e-14, rtol=0)
        y, x = np.mgrid[:1033, :1247]
        planar = np.broadcast_to(4000+3*x-2*y, (3, 1033, 1247)).copy()
        with self.assertRaisesRegex(backflat.StarlessOrientationError, "insufficient"):
            backflat.orient_starless(planar, planar[:, ::-1, :])

    def test_flat_symmetric_nearly_symmetric_and_unrelated_are_ambiguous(self):
        image, starless = self.pair()
        symmetric = (starless+starless[:, ::-1, :])/2
        unrelated = self.pair(seed=180)[1]
        cases = [(image, symmetric), (symmetric, symmetric+0.001*(starless-symmetric)),
                 (np.ones_like(image), np.ones_like(image)), (image, unrelated)]
        for source, supplied in cases:
            with self.subTest(variation=float(np.std(supplied))):
                with self.assertRaisesRegex(backflat.StarlessOrientationError, "--starless-flip-y"):
                    backflat.orient_starless(source, supplied)

    def test_manual_overrides_are_explicit_and_mutually_exclusive(self):
        image = np.arange(3*9*13).reshape(3, 9, 13).astype(np.float64)
        for mode in ("no-flip", "flip-y"):
            actual, decision = backflat.orient_starless(image, image, mode)
            np.testing.assert_array_equal(actual, image[:, ::-1, :] if mode == "flip-y" else image)
            self.assertEqual(decision["scores"], [])
            self.assertEqual(decision["mode"], mode)
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            backflat.parse_args(["backflat", "input.fit", "output.fit", "--starless", "starless.fit",
                                 "--starless-flip-y", "--starless-no-flip"])

    def test_cli_records_decision_and_preserves_input_mask_coordinates(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            image, starless = self.pair()
            header = fits.Header({"CRPIX1": 90.5, "CRPIX2": 120.5})
            fits.PrimaryHDU(image, header).writeto(folder/"input.fit")
            mask = np.zeros(image.shape[1:], dtype=np.uint16)
            mask[30:45, 60:100] = 65535
            fits.PrimaryHDU(mask).writeto(folder/"mask.fit")
            backgrounds = []
            for flipped in (False, True):
                source = folder/("flipped.fit" if flipped else "aligned.fit")
                fits.PrimaryHDU(starless[:, ::-1, :] if flipped else starless).writeto(source)
                output = folder/str(flipped)/"output.fit"
                args = ["backflat", str(folder/"input.fit"), str(output), "--starless", str(source),
                        "--mask", str(folder/"mask.fit"), "--no-gui"]
                log, error = io.StringIO(), io.StringIO()
                with contextlib.redirect_stdout(log), contextlib.redirect_stderr(error):
                    self.assertEqual(backflat.main(args), 0, error.getvalue())
                self.assertEqual("Corrected starless Y orientation" in log.getvalue(), flipped)
                with fits.open(output, memmap=False) as hdus:
                    history = str(hdus[0].header["HISTORY"])
                    self.assertIn("applied="+("flip-y" if flipped else "as-is"), history)
                    self.assertIn("hp-corr-v1", history)
                    self.assertIn("corr as-is=", history)
                    self.assertEqual(hdus[0].header["CRPIX2"], 120.5)
                with fits.open(output.parent/"back_mask.fit", memmap=False) as hdus:
                    np.testing.assert_array_equal(hdus["RAWMASK"].data, mask)
                backgrounds.append(fits.getdata(output.parent/"background.fit", memmap=False))
            np.testing.assert_array_equal(*backgrounds)

    def test_all_engine_sources_and_cache_use_the_same_measured_orientation(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            image, starless = self.pair((513, 769))
            source, executable = folder/"input.fit", folder/"engine.exe"
            fits.PrimaryHDU(image).writeto(source)
            executable.write_bytes(b"simulated engine")
            scale = max(float(np.max(np.abs(image))), 1)
            for engine in ("sxt", "starnet"):
                for flip in (False, True):
                    raw = (starless[:, ::-1, :] if flip else starless)/scale
                    def execute(command, **kwargs):
                        output_flag = "-o" if engine == "sxt" else "--output"
                        fits.PrimaryHDU(raw).writeto(command[command.index(output_flag)+1])
                        return ""
                    with self.subTest(engine=engine, flip=flip), \
                         mock.patch.object(backflat.shutil, "which", return_value=str(executable)), \
                         mock.patch.object(backflat, "discover_engine", return_value=(str(executable), "test")) as probe, \
                         mock.patch.object(backflat, "external_run", side_effect=execute), \
                         contextlib.redirect_stdout(io.StringIO()):
                        cache = str(folder/(engine+str(flip)))
                        loaded, _ = backflat.cached_starless(str(source), image, fits.Header(), engine, None, cache)
                        actual, decision = backflat.orient_starless(image, loaded)
                        self.assertEqual(decision["flip_y"], flip)
                        np.testing.assert_allclose(actual, starless, atol=1e-12)
                        probe.side_effect = RuntimeError("must use cache")
                        loaded, _ = backflat.cached_starless(str(source), image, fits.Header(), engine, None, cache)
                        again, repeated = backflat.orient_starless(image, loaded)
                        np.testing.assert_array_equal(actual, again)
                        self.assertEqual(decision, repeated)

    def test_batch_stops_at_changed_or_ambiguous_orientation(self):
        for modes in ((True, True), (False, True, False), (False, "symmetric", False)):
            with self.subTest(modes=modes), tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                image, starless = self.pair()
                fits.PrimaryHDU(np.zeros(image.shape[1:], dtype=np.uint16)).writeto(folder/"mask.fit")
                for i, mode in enumerate(modes, 1):
                    supplied = ((starless+starless[:, ::-1, :])/2 if mode == "symmetric" else
                                starless[:, ::-1, :] if mode else starless)
                    fits.PrimaryHDU(image).writeto(folder/("input%04d.fit" % i))
                    fits.PrimaryHDU(supplied).writeto(folder/("stars%04d.fit" % i))
                args = ["backflat", str(folder/"input0001.fit"), str(folder/"output0001.fit"),
                        "--starless", str(folder/"stars0001.fit"), "--mask", str(folder/"mask.fit"), "--no-gui"]
                log, error = io.StringIO(), io.StringIO()
                with contextlib.redirect_stdout(log), contextlib.redirect_stderr(error):
                    status = backflat.main(args)
                self.assertEqual(status, 0 if len(modes) == 2 else 1, error.getvalue())
                self.assertTrue((folder/"output0001.fit").is_file())
                if status:
                    self.assertIn("ambiguous" if modes[1] == "symmetric" else "Inconsistent", error.getvalue())
                    self.assertFalse((folder/"output0002.fit").exists())
                    self.assertFalse((folder/"output0003.fit").exists())
                    self.assertFalse((folder/"output0002_backflat"/"back_mask.fit").exists())
                else:
                    self.assertTrue((folder/"output0002.fit").is_file())


class CacheCleanupContracts(unittest.TestCase):
    @staticmethod
    def cache_file(folder, digit="a", data=None):
        folder.mkdir(parents=True, exist_ok=True)
        path = folder/(digit*64+".fit")
        header = fits.Header()
        header.add_history("backflat starless cache: test")
        fits.PrimaryHDU(np.ones((3, 9, 13)) if data is None else data, header).writeto(path)
        return path

    def test_engine_cache_removed_after_gui_headless_error_and_interrupt(self):
        for outcome in ("gui", "headless", "error", "interrupt"):
            with self.subTest(outcome=outcome), tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                image = np.full((3, 9, 13), 500.0)
                fits.PrimaryHDU(image).writeto(folder/"input.fit")
                fits.PrimaryHDU(np.zeros((9, 13), dtype=np.uint16)).writeto(folder/"mask.fit")
                executable = folder/"engine.exe"
                executable.write_bytes(b"test engine")
                cache = folder/".backflat-cache"
                self.cache_file(cache)  # Also remove recognized leftovers from an earlier run.
                args = ["backflat", str(folder/"input.fit"), str(folder/"output.fit"),
                        "--sxt", "--sxt-exe", str(executable), "--mask", str(folder/"mask.fit"),
                        "--starless-no-flip"]
                if outcome != "gui":
                    args.append("--no-gui")

                def execute(command, **kwargs):
                    if outcome == "error":
                        raise RuntimeError("simulated engine error")
                    if outcome == "interrupt":
                        raise KeyboardInterrupt()
                    data = fits.getdata(command[2], memmap=False)
                    fits.PrimaryHDU(data[:, ::-1, :]).writeto(command[command.index("-o")+1])
                    return ""

                def editor(session, **kwargs):
                    def close_window():
                        self.assertTrue(cache.is_dir())
                        session.archive()
                        session.apply()
                    return mock.Mock(run=close_window)

                out, err = io.StringIO(), io.StringIO()
                with mock.patch.object(backflat.shutil, "which", return_value=str(executable)), \
                     mock.patch.object(backflat, "discover_engine", return_value=(str(executable), "test")), \
                     mock.patch.object(backflat, "external_run", side_effect=execute), \
                     mock.patch.object(backflat, "MaskEditor", side_effect=editor), \
                     contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                    status = backflat.main(args)
                self.assertEqual(status, {"gui": 0, "headless": 0, "error": 1, "interrupt": 130}[outcome], err.getvalue())
                self.assertFalse(cache.exists())
                self.assertIn("Removed starless cache directory:", out.getvalue())
                if status == 0:
                    for name in ("output.fit", "background.fit", "back_mask.fit"):
                        self.assertTrue((folder/name).is_file())

    def test_custom_cache_preserves_inputs_other_files_and_subdirectories(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            cache = folder/"custom"
            old = self.cache_file(cache)
            source = self.cache_file(cache, "b", np.full((3, 9, 13), 500.0))
            unrelated = cache/("c"*64+".fit")
            fits.PrimaryHDU(np.ones((3, 9, 13))).writeto(unrelated)
            (cache/"keep.txt").write_bytes(b"unrelated")
            (cache/"nested").mkdir()
            (cache/"nested"/"keep.txt").write_bytes(b"unrelated nested file")
            fits.PrimaryHDU(np.zeros((9, 13))).writeto(folder/"mask.fit")
            args = ["backflat", str(source), str(folder/"output.fit"), "--starless", str(source),
                    "--cache-dir", str(cache), "--mask", str(folder/"mask.fit"), "--no-gui", "--starless-no-flip"]
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(backflat.main(args), 0)
            self.assertFalse(old.exists())
            for path in (source, unrelated, cache/"keep.txt", cache/"nested"/"keep.txt"):
                self.assertTrue(path.is_file())

    def test_batch_cache_survives_until_last_item_and_preflight_failure_keeps_it(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            cache = folder/"cache"
            cached = self.cache_file(cache)
            for i in (1, 2):
                fits.PrimaryHDU(np.ones((3, 9, 13))).writeto(folder/("input%04d.fit" % i))
                fits.PrimaryHDU(np.ones((3, 9, 13))).writeto(folder/("stars%04d.fit" % i))
            args = ["backflat", str(folder/"input0001.fit"), str(folder/"output0001.fit"),
                    "--starless", str(folder/"stars0001.fit"), "--cache-dir", str(cache)]
            existing = folder/"output0002.fit"
            existing.write_bytes(b"existing output")
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(backflat.main(args), 1)
            self.assertTrue(cached.is_file())

            def process(*args):
                self.assertTrue(cached.is_file())

            with mock.patch.object(backflat, "process_item", side_effect=process) as run, \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(backflat.main(args+["-y"]), 0)
                self.assertEqual(run.call_count, 2)
            self.assertFalse(cache.exists())

    def test_cleanup_refuses_changed_location_and_reports_permission_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            cache = folder/".backflat-cache"
            path = self.cache_file(cache)
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertFalse(backflat.cleanup_caches({str(cache): str(folder/"wrong")}, set()))
            self.assertTrue(path.is_file())
            fits.PrimaryHDU(np.ones((3, 9, 13))).writeto(folder/"input.fit")
            args = ["backflat", str(folder/"input.fit"), str(folder/"output.fit"),
                    "--starless", str(folder/"input.fit")]
            error = io.StringIO()
            with mock.patch.object(backflat, "process_item"), \
                 mock.patch.object(backflat.os, "remove", side_effect=PermissionError("locked: " + str(path))), \
                 contextlib.redirect_stderr(error), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(backflat.main(args), 1)
            self.assertIn("Error cleaning starless cache:", error.getvalue())
            self.assertIn(str(path), error.getvalue())
            self.assertTrue(path.is_file())


class GuiContracts(unittest.TestCase):
    def test_paint_undo_zoom_apply_and_save_in_dummy_display(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.dict(os.environ, {
            "SDL_VIDEODRIVER": "dummy", "SDL_AUDIODRIVER": "dummy", "PYGAME_HIDE_SUPPORT_PROMPT": "1"}):
            image = np.full((3, 91, 127), 500.0)
            session = backflat.Session(image, image, fits.Header(), dict(backflat.DEFAULTS),
                                       os.path.join(directory, "back_mask.fit"))
            session.apply(model=False)
            ui = backflat.MaskEditor(session)
            pg = ui.pg
            try:
                ui.render()
                point = tuple(round(v) for v in ui.image_to_screen(60, 45))
                ui.event(pg.event.Event(pg.MOUSEBUTTONDOWN, pos=point, button=1))
                ui.event(pg.event.Event(pg.MOUSEBUTTONUP, pos=point, button=1))
                self.assertTrue(session.raw_mask[45, 60])
                self.assertEqual(len(ui.undo), 1)
                ui.undo_stroke()
                self.assertFalse(session.raw_mask.any())
                before = ui.screen_to_image(*point)
                ui.set_zoom(ui.zoom*2, point)
                np.testing.assert_allclose(ui.screen_to_image(*point), before)
                session.raw_mask[30:60, 35:80] = True
                ui.start_apply()
                ui.future.result(timeout=15)
                self.assertTrue(Path(session.mask_path).exists())
                np.testing.assert_allclose(session.result, image, atol=1e-9)
                ui.render()
            finally:
                ui.pool.shutdown(wait=True)
                pg.quit()


if __name__ == "__main__":
    unittest.main()
