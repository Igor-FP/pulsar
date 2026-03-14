# PULSAR - Script Reference

A toolkit for batch processing of astronomical FITS images.

---

## Script Quick Reference

| Script | Purpose |
|--------|---------|
| **add.py** | Addition: `result = input + operand + offset` |
| **sub.py** | Subtraction: `result = input - operand + offset` |
| **mul.py** | Multiplication: `result = input * operand * scale` |
| **div.py** | Division: `result = (input / operand) * scale` |
| **arith.py** | Universal arithmetic (add/sub/mul/div in one script) |
| **sum.py** | Stack summation with EXPTIME update |
| **med.py** | Median combining (tiled parallel mode) |
| **calibrate.py** | Calibration: dark, bias, flat, cosmetic correction |
| **autocalibrate.py** | Auto-calibration from dark/flat trees with EXPTIME/FILTER matching |
| **normalize.py** | Brightness normalization between frames (linear regression) |
| **ngain.py** | Gain normalization (scale median to target value) |
| **noffset.py** | Offset normalization (shift median to target value) |
| **autoflat.py** | Flat field gradient correction (polynomial fitting) |
| **cosme.py** | Hot pixel correction from coordinate list |
| **make_cosme.py** | Hot pixel list generation from dark frame |
| **makedark.py** | Master dark and cosme.lst creation from raw darks |
| **makeflat.py** | Master flat creation per filter from raw flats |
| **newflat.py** | Add entry to maintenance log (for autocalibrate --flatlog) |
| **darkopt.py** | Optimized dark subtraction (K coefficient fitting) |
| **sortfits.py** | FITS sorting by time, session splitting |
| **autosolve.py** | Astrometric solving (WCS), reprojection, alignment |
| **fits2tiff.py** | FITS to TIFF conversion (8/16/32-bit) |
| **tiff2fits.py** | TIFF to FITS conversion (with header recovery) |
| **raw2fits.py** | Camera RAW to FITS conversion (raw Bayer CFA, currently Canon CR2/CR3) |
| **fft_align.py** | FFT-based frame alignment (rotation, scale, shift) |
| **absession.py** | Generate AstroBin acquisition session CSV |
| **binxy.py** | Software 2×2 / 4×4 pixel binning |
| **crop.py** | Crop FITS images (by size/center or margins) |
| **debayer.py** | Demosaic Bayer-pattern FITS to RGB |
| **hotfix.py** | Remove single hot (and cold) pixels |
| **mtf.py** | Midtone Transfer Function |
| **rgbbalance.py** | RGB color balance and brightness normalization |
| **bestof.py** | Select best frames by FWHM (seeing quality) |
| **rgb.py** | Merge/split RGB channels (3 mono ↔ 1 RGB FITS) |
| **staralign.py** | Star-based image registration (Automatic, Thin Plate Spline) |
| **xisf2fits.py** | Convert XISF (PixInsight) files to FITS |

---

## General Rules for All Scripts

### Input Data Formats (input_spec)

- **Single file**: `image.fit`
- **Numbered sequence**: `light0001.fit` → automatically discovers light0002.fit, light0003.fit, ...
- **Wildcard mask**: `*.fit`, `light_*.fit`
- **List file**: `@list.txt` or `list.txt` (one path per line)

### Numeric Constants as Operands

Arithmetic scripts (add, sub, mul, div, arith) support numeric constants as operands. A constant is interpreted as a "virtual file" filled with that value.

**Allowed constant formats** (strict parsing):
- Integers: `123`, `-456`, `+789`
- Decimals: `123.456`, `-0.5`, `+1.0`

**NOT allowed**:
- Exponential notation: `1e5`, `1E-3`
- Special values: `inf`, `-inf`, `nan`
- Digit separators: `1_000`, `1,000`
- Incomplete numbers: `.5`, `5.`

**Important**: At least one argument must be a file — image dimensions and FITS header are taken from it.

**Use case**: Useful for calibration on thermally stabilized CCDs where BIAS frames can be replaced with zero or a known offset constant.

### Output Data Formats (output_spec)

- **Single file**: `output.fit`
- **Numbered pattern**: `out0001.fit` → auto-increment for each input file
- **Directory**: `output_dir/` → preserves input filenames

### Data Types

**Supported formats:**
- Integer: int8, int16, int32, int64 (signed/unsigned)
- Floating-point: float32, float64

**Output conversion:**
- Integer FITS: result is clamped to type range and rounded
- Floating-point FITS: saved as-is
- NaN, Inf, -Inf: replaced with 0 (never written to output files)

**When float64 is used:**
- Multiplication/division by coefficients ~1.0 (normalization, flat correction)
- Mean calculation (rounding error accumulation)
- Operations with fractional coefficients

**When float64 is NOT needed:**
- Addition/subtraction of integers or images
- Median (element selection, no arithmetic)
- Min/max operations

---

## Detailed Script Descriptions

---

### add.py

**Purpose**: Add a value or image to input frames.

**Formula**: `result = input + operand + offset`

**Syntax**:
```
add.py input_spec output_spec operand [offset]
```

**Parameters**:
- `input_spec` — input files
- `output_spec` — output files
- `operand` — numeric constant OR FITS file OR numbered FITS pattern
- `offset` — optional numeric offset (default 0)

**Examples**:
```bash
add light0001.fit cal0001.fit 100
add *.fit out0001.fit bias.fit
add light0001.fit result0001.fit dark0001.fit 500
add image.fit result.fit 1024           # add constant to all pixels
```

---

### sub.py

**Purpose**: Subtract a value or image from input frames.

**Formula**: `result = input - operand + offset`

**Syntax**:
```
sub.py input_spec output_spec operand [offset]
```

**Parameters**:
- `input_spec` — input files
- `output_spec` — output files
- `operand` — numeric constant OR FITS file (value to subtract)
- `offset` — optional offset (default 0)

**Options**:
- `--continuum [snr]` — Continuum subtraction mode. Detects stars in both input and operand, cross-matches them (tolerance 1.5 px), and scales the operand so star flux matches before subtraction. Stars subtract to zero, leaving only emission line signal (e.g. H-alpha). SNR threshold for star detection (default: 38).

**Formula with --continuum**: `result = input - K * operand + offset`, where K = sum(flux_input) / sum(flux_operand) computed from matched star photometry.

**Examples**:
```bash
sub light0001.fit dark_sub0001.fit master_dark.fit
sub raw.fit calibrated.fit dark.fit 100
sub light0001.fit cal0001.fit 0         # subtract zero (for thermally stable CCDs without bias)
sub light0001.fit cal0001.fit 1024      # subtract offset constant
sub ha.fit continuum_sub.fit red.fit --continuum
sub ha.fit output.fit broadband.fit --continuum 50
sub ha0001.fit out0001.fit red0001.fit 1000 --continuum   # with offset
```

---

### mul.py

**Purpose**: Multiply input frames by a value or image.

**Formula**: `result = input * operand * scale`

**Syntax**:
```
mul.py input_spec output_spec operand [scale]
```

**Parameters**:
- `input_spec` — input files
- `output_spec` — output files
- `operand` — multiplier (numeric constant or FITS)
- `scale` — additional multiplier (default 1.0)

**Examples**:
```bash
mul light0001.fit scaled0001.fit 2.5
mul image.fit result.fit gain_map.fit
mul light0001.fit doubled0001.fit 2     # multiply by constant
```

---

### div.py

**Purpose**: Divide input frames by a value or image.

**Formula**: `result = (input / operand) * scale`

**Notes**: Division by zero yields 0.

**Syntax**:
```
div.py input_spec output_spec operand [scale]
```

**Parameters**:
- `input_spec` — input files
- `output_spec` — output files
- `operand` — divisor (numeric constant or FITS)
- `scale` — result multiplier (default 1.0)

**Examples**:
```bash
div light0001.fit flat_div0001.fit master_flat.fit 5000
div image.fit normalized.fit 32768      # divide by constant (range normalization)
div light0001.fit norm0001.fit 1 5000   # divide by 1 with scaling (equivalent to mul)
```

---

### arith.py

**Purpose**: Universal arithmetic tool (combines add/sub/mul/div).

**Syntax**:
```
arith.py (add|sub|mul|div) input_spec output_spec operand [param]
```

**Operations**:
- `add` → `result = input + operand + offset`
- `sub` → `result = input - operand + offset`
- `mul` → `result = input * operand * scale`
- `div` → `result = (input / operand) * scale`

**Parameters**:
- `operand` — numeric constant OR FITS file OR numbered pattern
- `param` — for add/sub this is offset (default 0), for mul/div this is scale (default 1)

**Examples**:
```bash
arith add light0001.fit out0001.fit 100
arith sub light0001.fit cal0001.fit dark.fit 500
arith sub light0001.fit cal0001.fit 0           # subtract zero (bias replacement for thermally stable CCDs)
arith div image.fit result.fit flat.fit 5000
arith mul light0001.fit scaled0001.fit 2.5      # multiply by constant
```

---

### sum.py

**Purpose**: Sum a stack of FITS images.

**Features**:
- Updates EXPTIME (total exposure)
- Updates DATE-OBS (earliest time)
- For integer FITS: scales result to type maximum
- For floating-point FITS: computes average

**Syntax**:
```
sum.py input_spec output.fit
```

**Parameters**:
- `input_spec` — input files (sequence, mask, or single file)
- `output.fit` — output file

**Examples**:
```bash
sum light0001.fit stacked.fit
sum *.fit combined.fit
```

---

### med.py

**Purpose**: Median combining of a FITS image stack.

**Features**:
- Tiled parallel mode for large images
- Exact per-pixel median
- Uses all CPU cores

**Syntax**:
```
med.py input_spec output.fit [--tile N]
```

**Parameters**:
- `input_spec` — input files
- `output.fit` — output file
- `--tile N` — tile size in pixels (default 2048, 0 = no tiling)

**Examples**:
```bash
med dark0001.fit master_dark.fit
med flat*.fit master_flat.fit --tile 1024
med @list.txt median.fit --tile 0
```

---

### calibrate.py

**Purpose**: Image calibration (dark/bias subtraction, flat division, cosmetic correction).

**Formula**: `result = ((raw - bias) - dark * OPTIMIZ) * K / flat`

**Syntax**:
```
calibrate.py input_spec output_spec -d dark.fit [-b bias.fit] [-f flat.fit K] [-c cosme.lst] [-optimize|-o]
```

**Parameters**:
- `input_spec` — raw files
- `output_spec` — output files
- `-d dark.fit` — master dark (required)
- `-b bias.fit` — master bias (optional)
- `-f flat.fit K` — master flat and multiplier K (optional)
- `-c cosme.lst` — hot pixel list (optional)
- `-optimize` or `-o` — optimize dark coefficient for each frame

**Examples**:
```bash
calibrate raw0001.fit cal0001.fit -d dark10s.fit
calibrate raw0001.fit cal0001.fit -d dark.fit -b bias.fit -f flat.fit 5000
calibrate raw0001.fit cal0001.fit -d dark.fit -f flat.fit 5000 -c cosme.lst -o
```

---

### autocalibrate.py

**Purpose**: Automatic calibration with dark/flat matching by EXPTIME and FILTER.

**Features**:
- Searches darks by EXPTIME and JD (closest in time)
- Searches flats by FILTER and JD (closest, no more than N days in the future)
- Automatically finds cosme.lst in darks folder (by name: dark600s.fit → cosme600s.lst)
- MULT multiplier is computed automatically as flat frame median

**Syntax**:
```
autocalibrate.py [options] rawfiles.fit out_path dark_path flat_path
```

**Parameters**:
- `rawfiles.fit` — raw files (file, list, mask)
- `out_path` — output files base name
- `dark_path` — directory with master darks and cosme*.lst
- `flat_path` — directory with master flats

**Options**:
- `--bestflat` — auto-select best flat for each session (by minimum σ after division)
- `--debug` — save previews to ./debug/ for visual flat selection verification
- `--flat-future-days N` — max days in the future for flat search (default 2)
- `--flatlog FILE` — CSV flat update log for strict interval-based matching

**Calibration formula**: `result = ((raw - dark) * MULT) / flat`

Order: dark subtraction → cosmetic correction → flat division

**Output file format**: `<base>_exp<EXPTIME>_<FILTER>_<N>.fit`

**--bestflat mode**:
- Groups files by filter and session (noon-to-noon by local time)
- For each group creates preview: 4×4 downscale after dark subtraction
- For each flat: 16px median blur → downscale → division → gradient removal → σ
- Selects flat with minimum σ (fewer dust/vignetting artifacts)

**--flatlog mode**:
- CSV format: `DATETIME_UTC,CAMERA_ID[,COMMENT]`
- Camera ID is matched as substring against INSTRUME
- Flat selection logic:
  1. First searches in current interval (between log entries)
  2. If multiple flats in interval — takes closest (no newer than +N days)
  3. If nothing with +N days but interval has a flat — takes it (interval = priority)
  4. If interval is empty — searches globally (by +N days, as without flatlog)
  5. L fallback at the end if exact filter not found

**Examples**:
```bash
# Basic usage
autocalibrate raw*.fit calibrated darks/ flats/

# With best flat auto-selection
autocalibrate --bestflat raw*.fit calibrated darks/ flats/

# With debug previews
autocalibrate --bestflat --debug raw*.fit calibrated darks/ flats/

# With camera maintenance log
autocalibrate --flatlog maintenance.csv raw*.fit cal darks/ flats/

# Combined options
autocalibrate --bestflat --flat-future-days 3 --flatlog maint.csv raw*.fit out darks/ flats/
```

---

### newflat.py

**Purpose**: Add entry to camera maintenance log (for --flatlog in autocalibrate).

**Use case**: Marks the moment when flat frames become invalid (sensor cleaning, filter change, etc.)

**Syntax**:
```
newflat.py --camera CAMERA_ID --log FILE [--date DATETIME] [--comment TEXT]
```

**Parameters**:
- `--camera` — camera identifier (must be a substring of INSTRUME in FITS)
- `--log` — path to CSV log file
- `--date` — optional: date/time in ISO format (default: current UTC)
- `--comment` — optional: comment

**Log format**:
```csv
# Maintenance log
DATETIME_UTC,CAMERA_ID,COMMENT
2024-05-18T14:30:00,2600MM,Cleaned sensor
2024-12-21T10:00:00,2600MM,Changed dust cover
```

**Examples**:
```bash
# Add entry with current time
newflat --camera 2600MM --log maintenance.csv

# With specified date
newflat --camera 2600MM --log maintenance.csv --date 2024-05-18T14:30:00

# With comment
newflat --camera "ASI2600" --log maint.csv --comment "Sensor cleaning"
```

---

### normalize.py

**Purpose**: Brightness normalization of images relative to a reference frame.

**Model**: `I = B * R + C` → normalized result: `(I - C) / B`

**Syntax**:
```
normalize.py input_spec output_spec [basefile.fit] [method]
```

**Parameters**:
- `input_spec` — input files
- `output_spec` — output files
- `basefile.fit` — reference frame (default: first input)
- `method` — normalization method:
  - `1` — linear regression (default)
  - `2` — robust regression (sigma-clipping)
  - `3` — global iterative normalization of all frames

**Examples**:
```bash
normalize light0001.fit norm0001.fit
normalize light0001.fit norm0001.fit reference.fit 2
normalize *.fit normalized0001.fit 3
```

---

### ngain.py

**Purpose**: Gain normalization — scales each frame's median to target value.

**Formula**: `result = input * (target_median / current_median)`

**Equivalent**: NGAIN in Iris.

**Features**:
- Supports all data types: signed/unsigned int 8/16/32/64, float 32/64
- All computations in float64
- Proper conversion back with range boundary control
- Zero median images are copied unchanged

**Syntax**:
```
ngain.py input_spec output_spec target_median
```

**Parameters**:
- `input_spec` — input files
- `output_spec` — output files
- `target_median` — target median value for all output images

**Examples**:
```bash
ngain flat0001.fit norm_flat0001.fit 10000
ngain *.fit normalized0001.fit 5000
ngain @list.txt out0001.fit 32768
```

---

### noffset.py

**Purpose**: Offset normalization — shifts each frame's median to target value by adding a constant.

**Formula**: `result = input + (target_median - current_median)`

**Equivalent**: NOFFSET in Iris.

**Features**:
- Supports all data types: signed/unsigned int 8/16/32/64, float 32/64
- All computations in float64
- Proper conversion back with range boundary control
- Preserves relative pixel brightness (unlike ngain)

**Syntax**:
```
noffset.py input_spec output_spec target_median
```

**Parameters**:
- `input_spec` — input files
- `output_spec` — output files
- `target_median` — target median value for all output images

**Examples**:
```bash
noffset bias0001.fit norm_bias0001.fit 1000
noffset *.fit normalized0001.fit 5000
noffset light0001.fit out0001.fit 10000
```

---

### autoflat.py

**Purpose**: Flat field gradient correction (polynomial background fitting).

**Algorithm**:
1. Zero pixel mask expansion
2. Median filtering
3. Min-binning for background extraction
4. Polynomial surface fitting
5. Correction: `result = input - model + offset`

**Syntax**:
```
autoflat.py [-d] input_spec output_spec [poly_order]
```

**Parameters**:
- `-d` — debug mode (saves intermediate images)
- `input_spec` — input files
- `output_spec` — output files
- `poly_order` — polynomial order (default 1: plane)

**Examples**:
```bash
autoflat flat0001.fit corrected0001.fit
autoflat -d flat.fit debug_flat.fit 2
```

---

### cosme.py

**Purpose**: Hot pixel correction from coordinate list.

**Algorithm**: Replaces hot pixels with mean value of neighbors (3×3).

**Syntax**:
```
cosme.py input_spec output_spec cosme.lst [-t]
```

**Parameters**:
- `input_spec` — input files
- `output_spec` — output files
- `cosme.lst` — hot pixel list (format: `P x y`)
- `-t` — test mode: creates mask instead of correction

**cosme.lst format**:
```
P 123 456
P 789 101
# comment
```

**Examples**:
```bash
cosme light0001.fit clean0001.fit cosme.lst
cosme image.fit mask.fit cosme.lst -t
```

---

### make_cosme.py

**Purpose**: Generate hot pixel list from master dark.

**Alias**: `find_hot` (the `find_hot` command calls this same script)

**Algorithm**: Finds 10000 brightest pixels in the image.

**Syntax**:
```
make_cosme.py input.fit cosme.lst
```

**Parameters**:
- `input.fit` — master dark
- `cosme.lst` — output list file

**Examples**:
```bash
make_cosme master_dark.fit cosme.lst
```

---

### makedark.py

**Purpose**: Automatic master dark and cosmetic correction list creation from raw darks.

**Meta-script**: Uses sub.py, med.py, make_cosme.py.

**Algorithm**:
1. Scans input files, selects those with `IMAGETYP='Dark Frame'`
2. Groups by exposure time (`EXPTIME` or `EXPOSURE`)
3. For each group:
   - Subtracts bias from each dark
   - Median combines
   - Generates hot pixel list

**Syntax**:
```
makedark.py input_spec [bias_spec]
```

**Parameters**:
- `input_spec` — directory OR mask OR dark sequence
- `bias_spec` — optional: numeric constant OR FITS file OR bias file list/mask (default 0)

**Output files** (to current directory):
- `dark<exp>.fit` — master dark for each exposure
- `cosme<exp>.lst` — hot pixel list
- `bias.fit` — master bias (if bias file list was specified)

**Name format**: `dark300s.fit` for >= 1s, `dark500ms.fit` for < 1s

**Examples**:
```bash
makedark /path/to/darks                   # all darks from folder, bias=0
makedark dark0001.fit                     # dark sequence, bias=0
makedark *.fit 1024                       # mask, bias=constant
makedark /path/to/darks master_bias.fit   # with ready master bias
makedark /path/to/darks bias0001.fit      # will create bias.fit from sequence
```

---

### makeflat.py

**Purpose**: Automatic master flat creation for each filter from raw flats.

**Meta-script**: Uses sub.py, ngain.py, med.py, cosme.py, makedark.py.

**Algorithm**:
1. Scans input files, selects those with `IMAGETYP='Flat Frame'`
2. Groups by filter (`FILTER`)
3. Validates: all files in group must have same exposure
4. Searches for `dark<exp>.fit` and `cosme<exp>.lst` (first in current, then in input directory)
5. If darks not found — runs `makedark` automatically
6. For each filter group:
   - Subtracts corresponding dark
   - Normalizes (ngain) to target_median
   - Median combines
   - Applies cosmetic correction

**Syntax**:
```
makeflat.py input_spec [target_median]
```

**Parameters**:
- `input_spec` — directory OR mask OR flat sequence
- `target_median` — optional: target median for normalization (default 5000)

**Output files** (to current directory):
- `flat_<filter>.fit` — master flat for each filter

**Filter codes**:
| FILTER | File code |
|--------|-----------|
| OIII | flat_o.fit |
| SII | flat_s.fit |
| Ha | flat_h.fit |
| L | flat_l.fit |
| R | flat_r.fit |
| G | flat_g.fit |
| B | flat_b.fit |
| (other) | flat_<full_name>.fit |

**Examples**:
```bash
makeflat /path/to/flats                   # all flats from folder, median=5000
makeflat /path/to/flats 10000             # with different target_median
makeflat flat*.fit                        # file mask
makeflat @list.txt 8000                   # from list
```

---

### darkopt.py

**Purpose**: Optimized dark subtraction with K coefficient fitting.

**Formula**: `result = input - K * dark`

**K fitting algorithm**:
- Splits image into 64×64 tiles
- Excludes tiles with zeros or saturation
- Selects darkest tile
- Computes K using least squares method

**Syntax**:
```
darkopt.py input_spec output_spec master_dark.fit [cosme.lst]
```

**Parameters**:
- `input_spec` — input files
- `output_spec` — output files
- `master_dark.fit` — master dark
- `cosme.lst` — optional hot pixel list

**Examples**:
```bash
darkopt light0001.fit opt0001.fit master_dark.fit
darkopt light0001.fit opt0001.fit dark.fit cosme.lst
```

---

### sortfits.py

**Purpose**: Sort FITS files by observation time, split into sessions.

**Syntax**:
```
sortfits.py input_spec output_pattern [-s|--sessions] [--gap-hours H]
```

**Parameters**:
- `input_spec` — input files
- `output_pattern` — output name pattern
- `-s`, `--sessions` — session mode (output names: `<base>_Sssss_Fffff.fit`)
- `--gap-hours H` — gap in hours for session splitting (default 1.0)
- `--auto` — auto-naming mode: output names based on FITS headers: `{OBJECT}_exp{TIME}_{FILTER}[_S{SSSS}]_{NNNN}.fit`
- `--name NAME` — override OBJECT header value (use with --auto)
- `--group-num` — restart numbering per group (use with --auto)

**Examples**:
```bash
sortfits light0001.fit sorted0001.fit
sortfits *.fit session.fit --sessions --gap-hours 2
sortfits *.fit sorted/ --auto
sortfits *.fit sorted/ --auto --name NGC1097
sortfits *.fit sorted/ --auto --group-num
```

---

### autosolve.py

**Purpose**: Astrometric solving (WCS), reprojection to tangent plane (TAN), subpixel alignment.

**Features**:
- Uses astrometry.net (solve-field) via WSL
- Reprojection to gnomonic (TAN) projection
- Subpixel alignment via FFT
- WCS refit from .corr files

**Syntax**:
```
autosolve.py [options] input_spec output_spec
```

**Main options**:
- `--no-solve` — skip astrometric solving
- `--rectify` — reproject to tangent plane (TAN)
- `--rect-center-ra` — projection center RA (degrees)
- `--rect-center-dec` — projection center Dec (degrees)
- `--individual` — separate projection center for each file
- `--align` — subpixel alignment
- `--ref file.fit` — reference frame for alignment
- `--scale-low`, `--scale-high` — scale range (arcsec/pixel)
- `--radius` — search radius (degrees)

**Reprojection (--rectify)**:
Transforms image to gnomonic (TAN) projection — projection onto a plane tangent to the celestial sphere. By default center is taken from first file; can be set explicitly via `--rect-center-ra` and `--rect-center-dec`.

**Examples**:
```bash
autosolve light0001.fit solved0001.fit
autosolve --rectify --align --ref ref.fit light0001.fit aligned0001.fit
autosolve --rectify --rect-center-ra 180.5 --rect-center-dec 45.2 *.fit out0001.fit
```

**Installing astrometry.net (Windows, via WSL)**:

```bash
# 1. Install WSL (from PowerShell as Administrator):
wsl --install

# 2. In WSL, install astrometry.net:
sudo apt update
sudo apt install astrometry.net

# 3. Download index files (choose based on your field of view):
#    http://data.astrometry.net/4200/
#    Files index-42XX.fits — number determines FOV range:
#      4219 = 00-10'     (narrow field, small sensor / long focal length)
#      4218 = 10-22'
#      4217 = 22-44'
#      4216 = 44'-2°
#      4215 = 2°-2.8°
#      4214 = 2.8°-4°
#      4213 = 4°-5.6°
#      4212 = 5.6°-8°    (typical full-frame + medium focal length)
#      4211 = 8°-11°
#      4210 = 11°-16°
#      4209 = 16°-22°    (wide-angle lenses)
#      4208 = 22°-30'
#      4207 = 30°-40°    (ultra wide-angle)
#    Download those covering your telescope/lens FOV.

# Example downloads:
cd /usr/share/astrometry
sudo wget http://data.astrometry.net/4200/index-4212.fits
sudo wget http://data.astrometry.net/4200/index-4213.fits

# 4. Verify installation:
solve-field --help
```

**Python dependencies**: scipy, reproject (`pip install scipy reproject`)

---

### fits2tiff.py

**Purpose**: Convert FITS images to TIFF format.

**Features**:
- Selectable bit depth: 8-bit (stretch), 16-bit (clamp), 32-bit float (normalized to [0,1])
- Auto bit depth detection from FITS data type
- Wildcard output pattern: `*.tif` preserves input filenames
- 32-bit float normalized to [0.0, 1.0] for Photoshop compatibility

**Dependencies**: Pillow (`pip install Pillow`)

**Syntax**:
```
fits2tiff.py [--bits 8|16|32] [--flip] input_spec output_spec
```

**Parameters**:
- `input_spec` — input files (standard formats: file, mask, sequence, @list)
- `output_spec` — output `.tif`/`.tiff` file, numbered pattern (`out0001.tif`), or wildcard (`*.tif`)
- `--bits 8` — linear stretch [min,max] → [0,255], uint8
- `--bits 16` — clamp to [0,65535], uint16
- `--bits 32` — normalize [min,max] → [0.0,1.0], float32
- (no --bits) — auto: int8/uint8→8, int16/uint16→16, everything else→32
- `--flip` — vertical flip (FITS bottom-left → TIFF top-left)

**Examples**:
```bash
# Single file (auto bit depth)
fits2tiff image.fit image.tif

# Batch conversion preserving filenames
fits2tiff *.fit *.tif

# 8-bit for previews
fits2tiff --bits 8 light0001.fit preview0001.tif

# Numbered output
fits2tiff --bits 16 *.fit out0001.tif
```

---

### tiff2fits.py

**Purpose**: Convert TIFF images back to FITS format with header recovery.

**Features**:
- Automatic header recovery from existing output FITS file being overwritten
- Explicit header source via `--header`
- Supports TIFF 8-bit (L), 16-bit (I;16), 32-bit float (F)
- RGB TIFF automatically converted to grayscale

**Dependencies**: Pillow (`pip install Pillow`)

**Header logic** (priority order):
1. `--header source.fit` — headers from specified file
2. Output file already exists — headers read from it before overwriting
3. Neither — minimal header (dimensions and data type only)

**Syntax**:
```
tiff2fits.py [--flip] [--header source.fit] input_spec output_spec
```

**Parameters**:
- `input_spec` — input TIFF files (file, mask, numbered sequence, @list)
- `output_spec` — output FITS file, numbered pattern (`out0001.fit`), or wildcard (`*.fit`)
- `--header F` — take FITS headers from file F (applied to all outputs)
- `--flip` — vertical flip (reverses fits2tiff --flip)

**Examples**:
```bash
# Reverse conversion — headers taken from existing .fit files being overwritten
tiff2fits *.tif *.fit

# With explicit header source
tiff2fits edited.tif result.fit --header original.fit

# Batch conversion with numbering
tiff2fits img0001.tif out0001.fit --header reference.fit
```

**Typical workflow** (Photoshop):
```bash
# 1. Export to TIFF
fits2tiff *.fit *.tif

# 2. Edit in Photoshop, save .tif

# 3. Import back — headers picked up from existing .fit files
tiff2fits *.tif *.fit
```

---

### raw2fits.py

**Purpose**: Convert Camera RAW files to FITS with raw Bayer sensor data. Currently supports Canon CR2/CR3.

**Features**:
- Default mode: raw CFA Bayer mosaic (2D monochrome uint16, original sensor values)
- Optional `--debaer`: linear demosaicing to 3-channel RGB FITS (uint16)
- Full EXIF to FITS header mapping (camera, exposure, ISO, date, temperature)
- Built-in CR3 ISOBMFF parser (no external tools needed for Canon's new format)
- Canon MakerNote parsing for camera body temperature (CCD-TEMP)
- Auto dark/light frame detection based on bright pixel fraction
- Wildcard output pattern: `*.fit` preserves input filenames

**Dependencies**: rawpy (`pip install rawpy`), exifread (`pip install exifread`)

**FITS header mapping**:
| Source | FITS keyword | Description |
|--------|-------------|-------------|
| Make + Model | INSTRUME | Camera make and model |
| ExposureTime | EXPTIME, EXPOSURE | Exposure time in seconds |
| ISOSpeedRatings | GAIN | ISO speed |
| DateTimeOriginal | DATE-OBS | Observation date/time (ISO format) |
| (from DATE-OBS) | JD | Julian Date |
| Canon ShotInfo | CCD-TEMP | Camera body temperature (Celsius) |
| raw_pattern | BAYERPAT | Bayer pattern (e.g. RGGB) |
| (auto-detect) | IMAGETYP | "Dark Frame" or "Light Frame" |
| — | FILTER | "CFA" (raw) or "RGB" (debaer) |
| black_level_per_channel | BLKLVL_R/G/B/G2 | Per-channel black levels |

**Syntax**:
```
raw2fits.py [--debaer] [--debug] input_spec output_spec
```

**Parameters**:
- `input_spec` — input RAW files (file, wildcard `*.cr2`, numbered sequence `IMG0001.CR3`, @list.txt)
- `output_spec` — output FITS file, numbered pattern (`out0001.fit`), or wildcard (`*.fit`)
- `--debaer` — demosaic (debayer) to linear RGB instead of raw CFA
- `--debug` — write detailed diagnostic log (`<input_name>.log`) for each file

**Examples**:
```bash
# Single file
raw2fits IMG_0001.CR2 test.fit

# Batch — preserve filenames
raw2fits *.cr3 *.fit

# Numbered output
raw2fits IMG0001.CR3 out0001.fit

# Demosaiced RGB
raw2fits --debaer *.cr3 *.fit

# Debug mode (creates .log files)
raw2fits --debug IMG_0001.CR3 test.fit
```

---

### fft_align.py

**Purpose**: FFT-based frame alignment (rotation, scale, subpixel shift).

**Features**:
- Phase correlation for shift detection
- Pyramidal angle and scale search
- Post-processing: local affine correction
- Parallel processing

**Syntax**:
```
fft_align.py reference input_spec output_spec [options]
```

**Main options**:
- `--superfine` — high-precision pyramidal mode
- `--post-correction` — local affine correction
- `--max-angle N` — maximum search angle (default 5°)
- `--scale-delta N` — scale range ±N (default 0.01)
- `--flux` — flux preservation mode (linear interpolation)

**Examples**:
```bash
fft_align ref.fit light0001.fit aligned0001.fit
fft_align ref.fit *.fit out0001.fit --superfine --post-correction
fft_align ref.fit light0001.fit aligned0001.fit --flux --max-angle 2
```

---

### absession.py

**Purpose**: Generate AstroBin-compatible CSV acquisition session list from FITS files.

**Features**:
- Recursive directory scanning for FITS files
- Two modes: FITS header reading (default) or filename parsing (`--parsename`, no dependencies)
- Session date from DATE-OBS header (or file mtime as fallback); pre-noon files count as previous evening
- Grouping by date, filter, and exposure time
- CSV output compatible with AstroBin's CSV import format
- Summary statistics per filter

**Dependencies**: astropy (default mode), none (`--parsename` mode)

**AstroBin Filter IDs**:

> **IMPORTANT:** The filter IDs in the script are EXAMPLES. You MUST edit the
> `ASTROBIN_FILTER_IDS` dictionary in `absession.py` to match YOUR filters
> on AstroBin. Search for your filter at:
> `https://app.astrobin.com/equipment/explorer/filter?page=1`
> The ID is the number in the filter page URL, e.g.:
> `.../filter/4388/filter-name-...` → ID = 4388

| Filter | ID (example) |
|--------|--------------|
| L | 5652 |
| R | 5656 |
| G | 5646 |
| B | 5642 |
| Ha | 4388 |
| SII | 4396 |
| OIII | 4392 |

**Syntax**:
```
absession.py [options] [directory]
```

**Parameters**:
- `directory` — path to scan (default: current directory)
- `--parsename` — extract filter and exposure from filenames (no dependencies)
- `--out FILE` — write CSV to file (default: stdout)
- `--flat` — do not recurse into subdirectories

**Filename convention (`--parsename` mode)**:
```
Name_Type_FILTER_SeqNum_Binning_EXPOSURE_Temp.fit

Example: NGC253_L_R_54205_Bin1x1_120s_-10C.fit
         |      | | |     |      |    |
         |      | | |     |      |    +-- temperature (ignored)
         |      | | |     |      +------ exposure: 120 seconds
         |      | | |     +------------- binning (ignored)
         |      | | +------------------- sequence number (ignored)
         |      | +--------------------- FILTER: R, G, B, L, Ha, SII, OIII
         |      +----------------------- exposure type (ignored)
         +------------------------------ object name

Mosaic variant (panel number in object name, e.g. "M8_02"):
M8_02_L_G_11788_Bin1x1_120s_-10C.fit
|  |  | | |     |      |    |
|  |  | | |     |      |    +-- temperature (ignored)
|  |  | | |     |      +------ exposure: 120 seconds
|  |  | | |     +------------- binning (ignored)
|  |  | | +------------------- sequence number (ignored)
|  |  | +--------------------- FILTER: G
|  |  +----------------------- exposure type (ignored)
|  +-------------------------- mosaic panel number (ignored)
+------------------------------ object name
```
Detection: 7 fields (excl. extension) = normal frame, 8 fields = mosaic.
Only FILTER and EXPOSURE fields are used; all others are ignored.
Configure this naming pattern in APT, N.I.N.A., SGPro, or your imaging capture software.

**Session date logic**: files with DATE-OBS (or mtime) before noon are counted as part of the previous evening's session.

**Examples**:
```bash
# Current directory (reads FITS headers)
absession

# Specific directory
absession /path/to/NGC3576

# Filename parsing mode (no astropy needed)
absession --parsename /path/to/data

# Save to file
absession --out sessions.csv /path/to/data

# Non-recursive
absession --flat .
```

---

### staralign.py

**Purpose**: Star-based image registration for astronomical frames.

Automatically detects stars, matches them across frames, and aligns images with sub-pixel accuracy. Handles scale differences, rotation, mirror, and non-linear distortions (TPS). Robust to cross-filter matching (e.g. R-reference vs Ha-narrowband).

**Dependencies**: scipy, sep

**Two operating modes:**

1. **Reference alignment** (with `--ref`) — aligns all target frames to a reference frame
2. **Chromatic correction** (without `--ref`) — aligns R and B channels to G within each RGB file. Corrects atmospheric refraction and chromatic aberration

**Syntax**:
```
staralign input_spec output_spec --ref reference.fit [options]
staralign input_spec output_spec                     (chromatic mode)
```

**Parameters**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ref FILE` | — | Reference frame (omit for chromatic mode) |
| `--model {tps\|projective}` | tps | Registration model |
| `--smoothness F` | 0.25 | TPS smoothness (0 = exact interpolation) |
| `--descriptors N` | 20 | Pentagons per star (5-50) |

**Matching parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hash-tol F` | 0.05 | Pentagon hash tolerance in 6D. Increase to 0.1-0.15 for cross-filter |
| `--angle-tol F` | 0.15 | Angular verification tolerance (radians). Increase to 0.25-0.3 for sparse fields |
| `--min-vote N` | 2 | Min hash votes per star pair. Set to 1 for difficult cases |
| `--tolerance F` | 3.0 | RANSAC inlier tolerance (pixels). Increase to 5-10 for large distortion |

**Star detection:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--snr F` | 5.0 | Detection SNR threshold (lower = more stars) |
| `--max-stars N` | 150 | Initial max stars. Auto-retries with 2x, 2.5x, 3x on failure |
| `--no-retry` | off | Disable auto-retry with more stars |

**Other:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threads N` | CPU-1 | Worker threads (forced to 1 in --debug mode) |
| `--debug` | off | Detailed per-frame diagnostics and timing |
| `--no-mirror` | off | Disable mirror fallback |

**Auto-retry**: on matching failure, automatically retries with more stars: levels [N, 2N, 2.5N, 3N], capped by detected count. Disable with `--no-retry`.

**Recommended settings for difficult cases:**

| Scenario | Recommended Settings |
|----------|---------------------|
| Same filter, same scope | defaults work |
| Cross-filter (R vs Ha) | `--max-stars 200` |
| Different focal lengths | `--tolerance 5` or higher |
| Very sparse field (<30 stars) | `--min-vote 1 --angle-tol 0.25` |
| Very different star counts | `--max-stars 300 --hash-tol 0.1` |

**Output**: each line shows `[N/M] filename: OK inliers, tps_pairs, residual_px` or `FAIL reason`.
- **inliers** — star pairs that survived RANSAC outlier rejection
- **tps_pairs** — star pairs used for final TPS warp (typically 30-300)
- **residual_px** — median positional error in pixels (lower = better)

Final summary: counts, timing, list of failed files.

**Examples**:
```bash
# Align all frames to reference
staralign *.fit out0001.fit --ref ref.fit

# Cross-filter: R-reference, Ha-targets
staralign Ha_*.fit aligned0001.fit --ref R_ref.fit --max-stars 200

# Multi-threaded
staralign *.fit out0001.fit --ref ref.fit --threads 4

# Chromatic correction in RGB file
staralign color.fit corrected.fit

# Batch chromatic correction
staralign *.fit out0001.fit
```

---

## Dependencies

**Required**:
- Python 3.6+
- numpy
- astropy

**For advanced features**:
- scipy (fft_align, autoflat, autosolve, staralign)
- reproject (autosolve reprojection)
- astrometry.net (autosolve solving)
- Pillow (fits2tiff, tiff2fits)
- rawpy (raw2fits)
- exifread (raw2fits)
- sep (bestof, rgbbalance --autostar, sub --continuum, staralign)
- xisf (xisf2fits)

---

## Installation

**Windows — one-click installer:** run `setup.bat` in the project root. It checks Python, installs dependencies, and adds commands to PATH.

**PATH only (no dependency install):** run `Commands\setup.bat` to add the Commands folder to PATH.

After that, scripts are available as commands: `add`, `sub`, `med`, `calibrate`, etc.

On Linux, run scripts directly with python:
```bash
python Add/add.py --help
```

---

### binxy.py

**Purpose**: Software 2×2 / 4×4 pixel binning.

Reduces image dimensions by summing adjacent pixels. Supports 2D (H×W) and 3D (N×H×W) FITS images.

**Syntax**:
```
binxy input_spec output_spec --2|--4 [--keep_bitness]
```

**Parameters**:
- `--2` — 2×2 binning (reduce dimensions by 2)
- `--4` — 4×4 binning (two passes of 2×2)
- `--keep_bitness` — integer: divide sum by 4, keep original dtype (no promotion)

**Behavior by dtype**:

| Input | Operation | Output |
|-------|-----------|--------|
| Integer | Sum 2×2 block | Same or promoted if overflow (8→16→32→64) |
| Integer + `--keep_bitness` | Sum 2×2 → `// 4` | Same as input |
| Float | Average 2×2 block (×0.25) | Same as input |

Updates headers: `XPIXSZ`, `YPIXSZ`, `XBINNING`, `YBINNING`. Odd dimensions cropped to even.

**Examples**:
```bash
# 2x2 binning
binxy img0001.fit out0001.fit --2

# 4x4 binning (sequence)
binxy img0001.fit out0001.fit --4

# 2x2, keep original dtype (integer divide by 4)
binxy img0001.fit out0001.fit --2 --keep_bitness
```

---

### crop.py

**Purpose**: Crop FITS images by size/center or margins.

Two modes: fixed output size with center point, or trim margins from edges. Supports auto-centering on the brightest object. Out-of-bounds regions are zero-padded.

**Syntax**:
```
crop input_spec output_spec [options]
```

**Mode 1 — Crop to size**:
- `--width WW` — output width in pixels
- `--height HH` — output height in pixels
- `--center XX YY` — center point (default: image center)
- `--autocenter [P]` — auto-detect center by brightness (P = threshold, default 50)
- `--positions FILE` — per-file centers from CSV

**Positions CSV format**:
```
filename,cx,cy
009A7220.fit,1960,1091
009A7221.fit,1960,1091
```
Extra columns (e.g. `radius`) are ignored.

**Mode 2 — Trim margins**:
- `--top N`, `--bottom N`, `--left N`, `--right N` — pixels to trim

**Examples**:
```bash
# Crop 1000x1000 from image center
crop img.fit out.fit --width 1000 --height 1000

# Crop around specific point
crop img.fit out.fit --width 1000 --height 1000 --center 3000 2000

# Auto-center on brightest object
crop img.fit out.fit --width 1000 --height 1000 --autocenter

# Trim margins
crop img0001.fit out0001.fit --top 100 --bottom 100 --left 200 --right 200
```

---

### debayer.py

**Purpose**: Demosaic Bayer-pattern FITS images to 3-channel RGB.

Input: 2D monochrome FITS with CFA Bayer pattern (H×W).
Output: 3-channel RGB FITS (3×H×W).

**Syntax**:
```
debayer input_spec output_spec [--pattern PAT] [--method METHOD]
```

**Parameters**:
- `--pattern PAT` — Bayer pattern: RGGB, BGGR, GRBG, GBRG (default: from BAYERPAT header or RGGB)
- `--method METHOD` — demosaicing algorithm

**Methods**:

| Method | Quality | Speed | Dependency |
|--------|---------|-------|------------|
| `bilinear` | Basic | Fast | None (numpy) |
| `vng` | High (VNG) | Moderate | OpenCV (`pip install opencv-python`) |

**Examples**:
```bash
# Default (bilinear, pattern from header)
debayer img0001.fit out0001.fit

# VNG with explicit pattern
debayer *.fit rgb/ --method vng --pattern RGGB
```

---

### hotfix.py

**Purpose**: Remove single hot (and cold) pixels from FITS images.

For each pixel, computes mean and std of its 8 neighbors. If the pixel exceeds `mean + N*sigma`, replaces it with the neighbor mean. Stars are preserved: their PSF spans multiple pixels, so neighbors are also bright. Fully vectorized (no Python pixel loops). Supports 2D and 3D FITS.

**Syntax**:
```
hotfix input_spec output_spec [options]
```

**Parameters**:
- `--sigma N` — detection threshold in standard deviations (default 5)
- `--floor N` — minimum noise level in ADU (default: auto)
- `--cold` — also fix cold (dead) pixels below `mean - N*sigma`
- `--debug` — diagnostic output for the first file

**Examples**:
```bash
# Default (5 sigma, hot only)
hotfix img0001.fit out0001.fit

# Stricter threshold
hotfix *.fit fixed/ --sigma 7

# Fix both hot and cold pixels
hotfix img0001.fit out0001.fit --sigma 4 --cold
```

---

### mtf.py

**Purpose**: Midtone Transfer Function (as in PixInsight PixelMath).

Applies a nonlinear tone curve: `mtf(x, m) = (1-m)*x / (m + x*(1-2*m))`

The curve passes through (0,0), (m, 0.5), (1,1) where m is the midtones balance.
Pixel values are normalized to [0,1] using dtype's full range (uint16: 0..65535).
Supports 2D and 3D FITS images.

**Syntax**:
```
mtf input_spec output_spec K [options]
```

**Parameters**:
- `K` — midtones balance (0..1). K<0.5 brightens, K=0.5 identity, K>0.5 darkens
- `--preserve_color` — color images: MTF on per-pixel average, scale R/G/B proportionally
- `--k2 K2` — second MTF parameter for dual-MTF blending
- `--blend B` — blend factor (0.0 = only K, 1.0 = only K2, default 0.5)

**Examples**:
```bash
# Brighten midtones
mtf img0001.fit out0001.fit 0.2

# Color-preserving stretch
mtf img0001.fit out0001.fit 0.25 --preserve_color

# Dual MTF blending
mtf img.fit out.fit 0.15 --k2 0.45 --blend 0.3
```

---

### rgbbalance.py

**Purpose**: RGB color balance and brightness normalization for color FITS.

Neutralizes background color, applies white balance (auto or manual), and normalizes brightness across a sequence of frames. Input: 3-channel RGB FITS (3×H×W).

**Syntax**:
```
rgbbalance input_spec output_spec [options]
```

**Parameters**:
- `--auto [file]` — auto white balance: scale R & B to match Green channel range. Reference from file (default: first input)
- `--autoeach` — auto white balance computed independently per file
- `--rgb R G B` — manual per-channel scaling coefficients
- `--kmin N` — percent of darkest pixels for black level estimation (default 5)
- `--kmax N` — percent of brightest pixels for white level estimation (default 5)
- `--autostar [file]` — star-based white balance: detects stars on G channel, cross-matches with R and B (tolerance 1.5 px), measures flux via aperture photometry, computes K from total star flux ratios. More accurate than `--auto` for narrowband/broadband mixing. Reference from file (default: first input).
- `--snr N` — SNR threshold for --autostar star detection (default: 38)

`--rgb`, `--auto`, `--autoeach`, `--autostar` are mutually exclusive. Without any — only background neutralization and inter-frame brightness normalization.

**Algorithm**:
1. Quantile medians for each channel (bottom kmin%, top kmax%)
2. Black neutralization: shift channels to average min_median
3. Color balance around black point (auto: `K_R = range_G/range_R`)
4. Brightness normalization between frames

**Examples**:
```bash
# Auto white balance
rgbbalance img0001.fit out0001.fit --auto

# Auto white balance, specific reference
rgbbalance img0001.fit out0001.fit --auto ref.fit

# Manual coefficients
rgbbalance img0001.fit out0001.fit --rgb 1.1 2.0 0.5

# Only background + brightness normalization
rgbbalance img0001.fit out0001.fit

# Star-based white balance
rgbbalance img.fit out.fit --autostar
rgbbalance img0001.fit out0001.fit --autostar ref.fit --snr 50
```

---

### bestof.py

**Purpose**: Select best frames by FWHM (seeing quality).

Analyses star FWHM across a set of FITS images using SEP-based star detection. Can copy best N% frames, filter by FWHM threshold, and/or write a CSV report.

**Dependencies**: sep (`pip install sep`)

**Syntax**:
```
bestof input_spec [output_spec] [options]
```

**Selection modes** (at least one required unless --csv only):
- `--best P` — copy best P percent of frames (lowest FWHM), 1-99
- `--threshold T` — copy frames with FWHM <= T pixels

**Options**:
- `--csv FILE` — write CSV report: fwhm,filename (sorted by FWHM)
- `--snr N` — SNR threshold for star detection (default: 10)

**Examples**:
```bash
bestof light0001.fit best0001.fit --best 70
bestof *.fit selected/ --threshold 3.5
bestof *.fit --csv report.csv
bestof *.fit best/ --best 80 --csv report.csv
```

---

### rgb.py

**Purpose**: Merge/split RGB channels in FITS files.

Merge: combines 3 monochrome 2D FITS (R, G, B) into one 3-channel RGB FITS (3, H, W).
Split: extracts R, G, B channels from RGB FITS (3, H, W) into 3 monochrome FITS. Sets FILTER header to R/G/B for each channel.

**Syntax**:
```
rgb --merge inR inG inB --out outRGB
rgb --split inRGB --out outR outG outB
```

**Parameters**:
All file arguments support standard specs: single file, wildcard, numbered sequence, @list.txt. For merge, all three inputs must expand to the same number of files.

**Examples**:
```bash
# Merge single files
rgb --merge r.fit g.fit b.fit --out rgb.fit

# Merge sequences
rgb --merge r0001.fit g0001.fit b0001.fit --out rgb0001.fit

# Split single file
rgb --split rgb.fit --out r.fit g.fit b.fit

# Split sequence
rgb --split rgb0001.fit --out r0001.fit g0001.fit b0001.fit
```

---

### xisf2fits.py

**Purpose**: Convert XISF (PixInsight) files to FITS format.

Preserves pixel data as-is (float32/float64/uint16, no type conversion). Restores FITS header keywords from XISF FITSKeyword metadata. Handles RGB (HWC→CHW axis conversion) and mono images.

**Dependencies**: xisf (`pip install xisf`)

**Syntax**:
```
xisf2fits input_spec output_spec
```

**Parameters**:
- `input_spec` — single .xisf file, wildcard (*.xisf), or @list.txt
- `output_spec` — single .fit file, numbered pattern (out0001.fit), or directory (outdir/)

**Examples**:
```bash
xisf2fits image.xisf image.fit
xisf2fits *.xisf converted/
xisf2fits @list.txt out0001.fit
```

---

## Testing

Each script folder contains `run.bat` for testing:

```bash
cd Add
run.bat
```

Test data is located in `Samples1/`, `Samples2/`, `1/` folders.
