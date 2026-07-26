```
██████╗ ██╗   ██╗██╗     ███████╗ █████╗ ██████╗
██╔══██╗██║   ██║██║     ██╔════╝██╔══██╗██╔══██╗
██████╔╝██║   ██║██║     ███████╗███████║██████╔╝
██╔═══╝ ██║   ██║██║     ╚════██║██╔══██║██╔══██╗
██║     ╚██████╔╝███████╗███████║██║  ██║██║  ██║
╚═╝      ╚═════╝ ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝
```
### **P**ortable **U**tility **L**ibrary for **S**hell **A**strophotography **R**outines

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Platform: Windows | Linux](https://img.shields.io/badge/platform-Windows%20|%20Linux-lightgrey.svg)]()

---

## About

**PULSAR** is a general-purpose software system for reduction and processing of astronomical CCD/CMOS imaging data. Each tool performs a single well-defined operation on FITS images; tools compose into automated pipelines covering the full workflow from raw frames to calibrated, combined, and enhanced results.

The project is inspired by **[IRIS](http://www.astrosurf.com/buil/us/iris/iris.htm)** by Christian Buil. PULSAR brings that philosophy to the modern command line: each tool does one thing well, tools compose into pipelines, and everything is scriptable.

### Key Capabilities

- **Star-based alignment** — robust geometric matching using pentagon descriptors with RANSAC outlier rejection and Thin Plate Spline refinement. Handles rotation, scale, and optical distortion. Works across filters and varying star populations.

- **Optimal stacking** — weighted combination where each frame contributes according to its measured SNR and FWHM. Outlier rejection uses smooth sigma-fade instead of hard clipping, preserving faint signal near bright objects. Iterative refinement with background-compensated deviation maps.

- **Background estimation** — polynomial surface fitted to a sigma-clipped cell grid with median filtering and border extension. Accurately models sky gradients, vignetting, and light pollution while ignoring stars and extended objects. Central object masking for galaxies and nebulae.

- **Nonlinear stretch and color processing** — Midtone Transfer Function with automatic black/white level detection, multi-curve blending, and color preservation through HSL or luminance-ratio methods. LRGB composition with automatic RGB balancing by star photometry.

- **Calibration pipeline** — automatic dark/flat matching by exposure time, filter, and date. Optimized dark coefficient fitting. Hot pixel detection and correction. Master frame creation from raw calibration data. Observatory automation: best-flat selection by trial application quality, equipment maintenance interval tracking for master flats validation & selection.

- **Astrometric solving** — WCS coordinate solutions via astrometry.net, reprojection to tangent plane, subpixel alignment. Distortion-only rectification without derotation for mosaics and surveys. Supports RGB FITS input. Works on Windows through WSL.

### Why Command Line?

- **Automation** — scriptable pipelines from raw to final result, no GUI interaction needed
- **Reproducibility** — same command, same result, always
- **Scalability** — process hundreds of frames in batch with parallel execution
- **Integration** — works with schedulers, remote access, and any workflow

---

## Features

### Calibration
- Dark/bias subtraction with coefficient optimization
- Flat division with automatic filter-based selection
- Hot pixel cosmetic correction (from coordinate list or sigma-based detection)
- Automatic master dark and master flat creation

### Alignment
- Star-based registration with geometric descriptor matching, RANSAC, and TPS
- WCS astrometric solving and tangent-plane reprojection (or reuse of an embedded/prior WCS to skip re-solving)
- FFT-based alignment with rotation and scale search
- Chromatic aberration correction (R/B channel alignment to G)

### Stacking
- Optimal weighted stacking with SNR and FWHM-based weights
- Smooth outlier rejection (sigma-fade) with iterative deviation refinement
- Simple summation and median combining
- Brightness normalization between frames

### Image Processing
- Nonlinear stretch (MTF) with auto black/white levels and multi-curve blending
- Background field flattening with polynomial fitting
- LRGB composition with automatic RGB balancing
- Max-SNR super-luminance (inverse-variance blend of L with synthetic R+G+B)
- Color preservation through HSL and luminance-ratio methods, with optional background desaturation
- RGB color balance by star photometry
- Bayer demosaicing, pixel binning, cropping

### Conversion
- Camera RAW to FITS (Canon CR2/CR3) with EXIF mapping and Bayer CFA preservation
- FITS to TIFF, JPEG, PNG with auto stretch for screen preview; JPEG/PNG preserve full FITS header as metadata
- TIFF to FITS with header recovery

### Utilities
- Frame sorting by time with session splitting
- Best frame selection by FWHM
- AstroBin acquisition session CSV generator

- Mosaics and sky surveys (in development)
- Object detection — comet/asteroid identification, nova patrol (in development)

---

## Installation

### Step 1 — Clone the repository

```bash
git clone https://github.com/Igor-FP/pulsar.git
cd pulsar
```

### Step 2 — Install dependencies

#### Windows — one-click installer (recommended)

Double-click **`setup.bat`** in the project root. It will:
- Check Python version (3.6+)
- Install pip if missing
- Install all Python dependencies
- Add commands to PATH

That's it — everything is ready to use.

#### Manual install (Linux / macOS / advanced)

```bash
pip install -r requirements.txt
```

On Linux/macOS, run scripts directly: `python Add/add.py --help`

<details>
<summary>Selective install (minimal)</summary>

```bash
pip install numpy astropy              # required — core functionality
pip install scipy                      # autocalibrate, autoflat, autosolve, fft_align, crop, staralign
pip install sep                        # staralign, bestof, rgbbalance, stack (star detection)
pip install psutil                     # stack (optional, adaptive memory management)
pip install Pillow                     # fits2tiff, tiff2fits
pip install rawpy exifread             # raw2fits (CR2 fallback reader)
pip install reproject                  # autosolve (WCS reprojection)
pip install opencv-python              # debayer (--method vng)
```
</details>

### Step 3 — Add commands to PATH (manual, Windows only)

> **Note:** If you used `setup.bat`, this step is already done.

To add commands to PATH without the full installer:
```batch
Commands\setup.bat
```

### Step 4 — astrometry.net (optional, for autosolve.py)

Only needed if you use astrometric solving. On Windows requires WSL (Windows 10+):

```bash
# Install WSL (PowerShell as Administrator):
wsl --install

# In WSL:
sudo apt update
sudo apt install astrometry.net

# Download index files for your FOV (http://data.astrometry.net/4200/):
#   4216 = 44'-2°    4213 = 4°-5.6°    4210 = 11°-16°
#   4215 = 2°-2.8°   4212 = 5.6°-8°    4209 = 16°-22°
#   4214 = 2.8°-4°   4211 = 8°-11°     4207 = 30°-40°
cd /usr/share/astrometry
sudo wget http://data.astrometry.net/4200/index-4212.fits
```

See [SCRIPTS.md](SCRIPTS-english.md) for the full index file reference.

---

## Quick Start

### Usage Examples

```bash
# Create master dark from a folder of darks
python MakeDark/makedark.py /path/to/darks

# Create master flat for all filters
python MakeFlat/makeflat.py /path/to/flats

# Automatic calibration of image series
python Autocalibrate/autocalibrate.py raw*.fit calibrated/ darks/ flats/

# Median combine
python Med/med.py calibrated*.fit stacked.fit

# Astrometric solving and alignment
python Autosolve/autosolve.py --rectify --align *.fit aligned/
```

On Windows after running `setup.bat`, commands are available directly:
```batch
makedark C:\Darks
makeflat C:\Flats
autocalibrate raw*.fit calibrated\ darks\ flats\
med calibrated*.fit stacked.fit
autosolve --rectify --align *.fit aligned\
```

---

## Script Reference

| Script | Purpose |
|--------|---------|
| **add.py** | Add images/constants |
| **sub.py** | Subtract images/constants |
| **mul.py** | Multiply images/constants |
| **div.py** | Divide images/constants |
| **arith.py** | Universal arithmetic |
| **sum.py** | Stack summation |
| **med.py** | Median combining |
| **stack.py** | Optimal weighted stacking with smooth outlier rejection |
| **calibrate.py** | Calibration (dark/bias/flat/cosme) |
| **autocalibrate.py** | Auto-calibration with dark/flat matching |
| **normalize.py** | Brightness normalization (regression) |
| **ngain.py** | Gain normalization (multiply to target median) |
| **noffset.py** | Offset normalization (add to target median) |
| **autoflat.py** | Background field flattening |
| **cosme.py** | Hot pixel correction |
| **make_cosme.py** | Hot pixel list generation |
| **makedark.py** | Master dark creation |
| **makeflat.py** | Master flat creation |
| **darkopt.py** | Optimized dark subtraction |
| **sortfits.py** | Time-based sorting |
| **autosolve.py** | Astrometry and reprojection |
| **fits2tiff.py** | FITS to TIFF/JPEG/PNG conversion |
| **tiff2fits.py** | TIFF to FITS conversion |
| **raw2fits.py** | Camera RAW to FITS conversion (currently Canon CR2/CR3) |
| **staralign.py** | Star-based image registration |
| **fft_align.py** | FFT-based alignment |
| **absession.py** | AstroBin acquisition session CSV |
| **binxy.py** | Software 2×2 / 4×4 pixel binning |
| **crop.py** | Crop FITS images (by size/center or margins) |
| **debayer.py** | Demosaic Bayer-pattern FITS to RGB |
| **hotfix.py** | Remove single hot (and cold) pixels |
| **mtf.py** | Nonlinear brightness stretch with auto levels |
| **lrgb.py** | Combine luminance channel with RGB color |
| **rgbbalance.py** | RGB color balance and brightness normalization |

Full documentation: **[SCRIPTS.md](SCRIPTS-english.md)** (English) | **[SCRIPTS.md](SCRIPTS.md)** (Russian)

---

## Input Formats

All scripts support flexible file specification:

```bash
# Single file
python Add/add.py image.fit output.fit 100

# Numbered sequence (auto-discovers all files)
python Add/add.py light0001.fit output0001.fit dark.fit

# Wildcard mask
python Med/med.py *.fit combined.fit

# List file
python Calibrate/calibrate.py @list.txt output0001.fit -d dark.fit

# Output to directory (preserves filenames)
python Autosolve/autosolve.py input*.fit output_dir/
```

---

## Platforms

| Platform | Status |
|----------|--------|
| **Windows 10/11** | Full support (batch wrappers included) |
| **Windows 7** | Supported (except autosolve.py — requires WSL, available from Windows 10) |
| **Linux** | Full support (run scripts with python) |
| **macOS** | Should work (not tested) |

On Windows, astrometry (autosolve.py) uses WSL with astrometry.net installed.

---

## Project Structure

```
PULSAR/
├── Commands/          # Batch wrappers for Windows PATH
├── lib/               # Shared modules (batch_utils.py)
├── Add/               # add.py
├── Sub/               # sub.py
├── Mul/               # mul.py
├── Div/               # div.py
├── Arith/             # arith.py
├── Sum/               # sum.py
├── Med/               # med.py
├── Stack/             # stack.py
├── Calibrate/         # calibrate.py
├── Autocalibrate/     # autocalibrate.py
├── Normalize/         # normalize.py
├── NGain/             # ngain.py
├── NOffset/           # noffset.py
├── Autoflat/          # autoflat.py
├── Cosme/             # cosme.py
├── MakeCosme/         # make_cosme.py
├── MakeDark/          # makedark.py
├── MakeFlat/          # makeflat.py
├── DarkOpt/           # darkopt.py
├── SortFits/          # sortfits.py
├── Autosolve/         # autosolve.py
├── Fits2tiff/         # fits2tiff.py
├── Tiff2fits/         # tiff2fits.py
├── Raw2fits/          # raw2fits.py
├── Staralign/         # staralign.py
├── FFT_Align/         # fft_align.py
├── Absession/         # absession.py
├── Binxy/             # binxy.py
├── Crop/              # crop.py
├── Debayer/           # debayer.py
├── Hotfix/            # hotfix.py
├── Mtf/               # mtf.py
├── Lrgb/              # lrgb.py
├── Rgbbalance/        # rgbbalance.py
├── Samples*/          # Test data
├── setup.bat          # Windows one-click installer
├── SCRIPTS.md         # Detailed documentation (Russian)
├── CLAUDE.md          # AI assistant instructions
├── requirements.txt   # Python dependencies
└── README.md          # This file
```

---

## Development

### Adding a New Script

1. Create directory: `NewScript/`
2. Create script: `NewScript/newscript.py`
3. Create test: `NewScript/run.bat`
4. Create wrapper: `Commands/newscript.bat`
5. Update documentation: `SCRIPTS.md`

### Code Guidelines

- Use `batch_utils.py` for file handling
- Support all FITS data types (int8-64, float32-64)
- Never write NaN/Inf to output files
- Preserve FITS headers, update only necessary fields

---

## License

MIT License

---

## Author

**Igor Chekalin**

---

## Acknowledgments

- **Christian Buil** — for [IRIS](http://www.astrosurf.com/buil/us/iris/iris.htm), the inspiration for this project
- **Dustin Lang** and the **astrometry.net** team — for the open-source astrometric solver
- The amateur astrophotography community

---

## Related Projects

- [IRIS](http://www.astrosurf.com/buil/us/iris/iris.htm) — GUI astronomical image processing
- [Siril](https://siril.org/) — free astronomical image processing
- [PixInsight](https://pixinsight.com/) — professional processing platform
- [Astrometry.net](http://astrometry.net/) — astrometric solving service

---

PixInsight is a registered trademark of Pleiades Astrophoto S.L. This project is not affiliated with or endorsed by Pleiades Astrophoto.
