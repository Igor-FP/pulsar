<!-- ============== VARIANT 1 - Classic ============== -->
<!--
```
 ____  _   _ _     ____    _    ____
|  _ \| | | | |   / ___|  / \  |  _ \
| |_) | | | | |   \___ \ / _ \ | |_) |
|  __/| |_| | |___ ___) / ___ \|  _ <
|_|    \___/|_____|____/_/   \_\_| \_\

Portable Utility Library for Shell Astrophotography Routines
```
-->

<!-- ============== VARIANT 2 - Box with stars ============== -->
<!--
```
╔═══════════════════════════════════════════════════════════════════╗
║  ★  P U L S A R  ★                                                ║
║  Portable Utility Library for Shell Astrophotography Routines    ║
╚═══════════════════════════════════════════════════════════════════╝
```
-->

<!-- ============== VARIANT 3 - Minimalist ============== -->
<!--
```
    ____  __  ____   ___   ___  ____
   / __ \/ / / / /  / __\ / _ \/ __ \
  / /_/ / /_/ / /__/\__ \/ /_\ / /_/ /
 / .___/\____/____/\___//_/ \_\____/
/_/  Portable Utility Library for Shell Astrophotography Routines
```
-->

<!-- ============== VARIANT 4 - Cosmic ============== -->
<!--
```
       *  .  *       .       *    .        *   .
    .    *       P U L S A R        .  *      .
  *    .    *                   .        *
 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Portable Utility Library for Shell Astrophotography Routines
 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
       .    *    .        *   .      *     .
```
-->

<!-- ============== VARIANT 5 - Block (big) ============== -->
<!--
```
██████╗ ██╗   ██╗██╗     ███████╗ █████╗ ██████╗
██╔══██╗██║   ██║██║     ██╔════╝██╔══██╗██╔══██╗
██████╔╝██║   ██║██║     ███████╗███████║██████╔╝
██╔═══╝ ██║   ██║██║     ╚════██║██╔══██║██╔══██╗
██║     ╚██████╔╝███████╗███████║██║  ██║██║  ██║
╚═╝      ╚═════╝ ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝

Portable Utility Library for Shell Astrophotography Routines
```
-->

<!-- ============== VARIANT 6 - Classic with letter meanings ============== -->
<!--
```
 ____  _   _ _     ____    _    ____
|  _ \| | | | |   / ___|  / \  |  _ \
| |_) | | | | |   \___ \ / _ \ | |_) |
|  __/| |_| | |___ ___) / ___ \|  _ <
|_|    \___/|_____|____/_/   \_\_| \_\

[P]ortable [U]tility [L]ibrary for [S]hell [A]strophotography [R]outines
```
-->

<!-- ============== CHOSEN VARIANT (uncomment one above and paste here) ============== -->

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

## Vision and Goals

**PULSAR** is a suite of command-line tools for complete automation of amateur astrophotography processing — from raw frame calibration to final results.

### Inspiration

This project is ideologically inspired by **[IRIS](http://www.astrosurf.com/buil/us/iris/iris.htm)** by Christian Buil — a legendary astronomical image processing software. PULSAR brings the IRIS philosophy to the modern command line, enabling full automation through scripts and pipelines.

### Project Goals

- **Full automation** — from RAW to final result with zero manual intervention
- **Calibration** — automatic dark/flat selection based on metadata (EXPTIME, FILTER, date)
- **Stacking** — summation and median combining with alignment
- **Astrometry** — automatic WCS solving and reprojection to a common coordinate grid
- **Object detection** — comet and asteroid identification, nova patrol (in development)
- **Mosaics and surveys** — panorama stitching and sky surveys (in development)

### Why Command Line?

- **Automation**: easily integrates into scripts and pipelines
- **Reproducibility**: one command — one result, always
- **Scalability**: process thousands of frames with a single command
- **Integration**: works with any scheduler, CI/CD, remote access

---

## Features

### Calibration
- Dark/bias subtraction with coefficient optimization
- Flat division with automatic filter-based selection
- Hot pixel cosmetic correction
- Automatic master dark and master flat creation

### Arithmetic
- Add, subtract, multiply, divide images
- Numeric constants supported as operands
- Works with all data types (int8-64, float32-64)

### Stacking
- Summation with exposure time tracking
- Median combining (parallel tiled processing)
- Brightness normalization between frames

### Astrometry and Alignment
- WCS solving via astrometry.net
- Reprojection to tangent plane (TAN/gnomonic projection)
- FFT alignment with subpixel accuracy
- Rotation and scale correction

### Utilities
- Time-based sorting with session splitting
- Hot pixel list generation
- Flat field gradient correction

---

## Quick Start

### Installation

```bash
git clone https://github.com/ichekailin/PULSAR.git
cd PULSAR

# Windows: add Commands to PATH
Commands\setup.bat

# Linux: run scripts directly with python
python Add/add.py --help
```

### Dependencies

```bash
pip install numpy astropy scipy reproject
```

For astrometry (autosolve.py), [astrometry.net](http://astrometry.net/) is required (via WSL on Windows, native on Linux).

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

On Windows with setup.bat, commands are available directly:
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
| **calibrate.py** | Calibration (dark/bias/flat/cosme) |
| **autocalibrate.py** | Auto-calibration with dark/flat matching |
| **normalize.py** | Brightness normalization (regression) |
| **ngain.py** | Gain normalization (multiply to target median) |
| **noffset.py** | Offset normalization (add to target median) |
| **autoflat.py** | Flat field gradient correction |
| **cosme.py** | Hot pixel correction |
| **make_cosme.py** | Hot pixel list generation |
| **makedark.py** | Master dark creation |
| **makeflat.py** | Master flat creation |
| **darkopt.py** | Optimized dark subtraction |
| **sortfits.py** | Time-based sorting |
| **autosolve.py** | Astrometry and reprojection |
| **fft_align.py** | FFT-based alignment |

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
├── FFT_Align/         # fft_align.py
├── Samples*/          # Test data
├── SCRIPTS.md         # Detailed documentation (Russian)
├── CLAUDE.md          # AI assistant instructions
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
