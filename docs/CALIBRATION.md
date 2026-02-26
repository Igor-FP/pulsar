# Calibration Scripts Guide

This guide covers `calibrate.py` (manual calibration) and `autocalibrate.py` (automatic calibration with dark/flat matching).

---

## Quick Start

### calibrate.py — Manual Calibration

Simple calibration with explicit master frames:

```bash
# Basic: subtract dark only
calibrate raw001.fit out.fit -d master_dark.fit

# Full: dark + bias + flat with multiplier
calibrate raw001.fit out.fit -d dark.fit -b bias.fit -f flat.fit 25000

# With optimization and cosmetic correction
calibrate raw*.fit calibrated.fit -d dark.fit -f flat.fit 25000 -c cosme.lst -o
```

### autocalibrate.py — Automatic Calibration

Auto-matches darks and flats from directory trees:

```bash
# Basic usage
autocalibrate lights/*.fit output darks/ flats/

# With best flat selection per session
autocalibrate --bestflat lights/*.fit output darks/ flats/

# With maintenance log for flat intervals
autocalibrate --flatlog maintenance.csv lights/*.fit output darks/ flats/

# Debug mode (saves preview images)
autocalibrate --bestflat --debug lights/*.fit output darks/ flats/
```

---

## Calibration Formula

Both scripts use the same formula:

```
RESULT = (RAW - BIAS - DARK × K) × MULT / FLAT
```

Where:
- **RAW** — input light frame
- **BIAS** — master bias (optional, `calibrate.py` only)
- **DARK** — master dark frame
- **K** — dark multiplier (1.0 or optimized per-frame)
- **FLAT** — master flat field
- **MULT** — flat multiplier = median(FLAT)

---

## calibrate.py — Detailed Reference

### Synopsis

```
calibrate.py rawfiles.fit outfiles.fit -d dark.fit [-b bias.fit] [-f flat.fit K] [-c cosme.lst] [-o]
```

### Options

| Option | Description |
|--------|-------------|
| `-d dark.fit` | Master dark frame (required) |
| `-b bias.fit` | Master bias frame (optional) |
| `-f flat.fit K` | Master flat and multiplier K (e.g., `-f flat.fit 25000`) |
| `-c cosme.lst` | Cosmetic correction pixel list |
| `-o` or `-optimize` | Compute optimal dark multiplier per frame |

### Dark Optimization

With `-o` flag, the script computes optimal dark scale factor K for each frame:

1. Divides image into 256×256 tiles
2. Excludes tiles with saturated or zero pixels
3. For each valid tile, computes K that minimizes residual: `sum((R - B - K×D)²)`
4. Returns median K across all valid tiles

This compensates for temperature-dependent dark current variations.

### Examples

```bash
# Process numbered sequence
calibrate light0001.fit calib.fit -d dark.fit -f flat_L.fit 24500

# Process wildcard pattern with optimization
calibrate *.fit output/ -d dark.fit -f flat.fit 25000 -o

# Full calibration pipeline
calibrate @filelist.txt result.fit -d dark.fit -b bias.fit -f flat.fit 25000 -c cosme.lst -o
```

---

## autocalibrate.py — Detailed Reference

### Synopsis

```
autocalibrate.py [options] rawfiles.fit output_path dark_path flat_path
```

### Options

| Option | Description |
|--------|-------------|
| `--bestflat` | Auto-select best flat per observing session |
| `--debug` | Save debug images to `./debug/` folder |
| `--flat-future-days N` | Max days in future to accept flat (default: 2) |
| `--flatlog FILE` | CSV file with flat renewal timestamps |

### Output Naming

Output files are named automatically:
```
<base>_exp<EXPTIME>_<FILTER>_<N>.fit
```

Example: `output_exp120_L_001.fit`, `output_exp120_L_002.fit`, ...

---

## Directory Structure

### What is REQUIRED vs RECOMMENDED

**REQUIRED** (will not work without):
- All metadata is read from **FITS headers only** (EXPTIME, JD, FILTER)
- File names and folder structure are **NOT parsed** — you can put everything in one folder

**RECOMMENDED** (for convenience):
- Organize by date/camera/telescope for easier management
- Use meaningful file names for human readability
- Keep cosmetic files next to corresponding darks

### Dark Directory

Directories are scanned **recursively** — any structure works:

```
darks/                    # All these structures work equally well:
├── dark120s.fit          # Everything in root
├── dark300s.fit
├── cosme.lst             # Single cosme file for all
│
├── 2024_05/              # Or organized by date
│   ├── dark120s.fit
│   └── cosme120s.lst
│
└── camera_A/             # Or by camera
    └── dark60s.fit
```

**Cosmetic file discovery** (in order of priority):
1. Same folder, name derived from dark file name:
   - If name starts with "dark": replace "dark" prefix with "cosme"
     - `dark120s.fit` → `cosme120s.lst`
     - `Dark_2024.fit` → `cosme_2024.lst`
   - Otherwise: replace all occurrences of "dark" with "cosme"
     - `master_dark_60s.fit` → `master_cosme_60s.lst`
2. If not found: most recent `.lst` file anywhere in dark tree (fallback)

### Flat Directory

Same as darks — scanned **recursively**, structure doesn't matter:

```
flats/                    # All these work:
├── flat_L.fit            # Everything in root
├── flat_Ha.fit
│
├── 2024_05/              # Or by date
│   ├── flat_L.fit
│   └── flat_Ha.fit
│
└── Telescope_A/          # Or by setup
    └── 2024_06/
        └── flat_L.fit
```

Flats are matched by **FILTER header** and **JD header** — not by file names or paths.

---

## Required FITS Headers

### Light Frames (RAW)

| Header | Required | Description |
|--------|----------|-------------|
| `EXPTIME` or `EXPOSURE` | Yes | Exposure time in seconds |
| `JD` | Yes | Julian Date of observation |
| `FILTER` or `FILTERS` | Yes | Filter name (L, R, G, B, Ha, etc.) |
| `INSTRUME` | For flatlog | Camera identifier |

### Dark Frames

| Header | Required | Description |
|--------|----------|-------------|
| `EXPTIME` or `EXPOSURE` | Yes | Exposure time in seconds |
| `JD` | Yes | Julian Date when dark was taken |

### Flat Frames

| Header | Required | Description |
|--------|----------|-------------|
| `FILTER` or `FILTERS` | Yes | Filter name |
| `JD` | Yes | Julian Date when flat was taken |

---

## Dark Matching Algorithm

Darks are matched by exposure time:

1. **Exact match**: Find darks with same EXPTIME (±1e-6 tolerance)
   - If multiple, select the most recent (highest JD)
2. **Closest match**: If no exact match, find closest EXPTIME
   - Tie-breaker: most recent JD

```
Example: RAW has EXPTIME=120s

Available darks:
  dark_120s.fit (JD=2460400)  ← Selected (exact match, latest)
  dark_120s.fit (JD=2460350)
  dark_60s.fit  (JD=2460450)
```

---

## Flat Matching Algorithms

### Standard Mode (default)

Flats are matched by filter and date:

1. Find flats with matching FILTER
2. Filter by date: `flat_JD ≤ raw_JD + flat_future_days`
3. Select most recent (highest JD)
4. If no match, try L filter as fallback

```
Example: RAW taken 2024-05-20, filter=Ha, --flat-future-days=2

Valid flats: JD ≤ 2024-05-22
  flat_Ha_2024-05-18.fit  ← Selected (most recent valid)
  flat_Ha_2024-05-10.fit
  flat_Ha_2024-05-25.fit  ✗ (too far in future)
```

### With Maintenance Log (--flatlog)

When optical train is modified (camera rotated, filters changed, etc.), flats become invalid. The maintenance log defines intervals:

```csv
# maintenance.csv
# DATETIME_UTC,CAMERA_ID[,COMMENT]
2024-03-15T14:00:00,ASI2600MM,Initial setup
2024-05-18T12:30:00,ASI2600MM,Rotated camera
2024-07-01T10:00:00,ASI2600MM,Changed filter wheel
```

**Matching logic**:
1. Find maintenance interval containing the observation
2. Search for flats within that interval first
3. Apply `--flat-future-days` limit within interval
4. If interval has flats but none within limit, use latest from interval
5. If interval empty, fall back to global search

```
Example: RAW taken 2024-06-15, camera=ASI2600MM

Maintenance intervals:
  2024-03-15 — 2024-05-18  (interval A)
  2024-05-18 — 2024-07-01  (interval B) ← RAW is here
  2024-07-01 — ...         (interval C)

Only flats from 2024-05-18 to 2024-07-01 are considered.
```

### Best Flat Selection (--bestflat)

Automatically selects the best matching flat by analyzing image quality:

**Algorithm**:

1. **Prepare light previews** (per observing session):
   - Downscale 8×8
   - Apply median filter (size=16) to remove stars
   - Combine via pixel-wise median across session

2. **Prepare flat previews**:
   - Downscale 8×8
   - Normalize to median 5000

3. **Evaluate each candidate flat**:
   - Divide: `light_preview / flat_preview`
   - Apply high-pass filter (removes gradients)
   - Apply gaussian blur (suppresses noise)
   - Compute sigma (standard deviation)

4. **Select best**: The flat that produces the flattest result with fewest mid-frequency local deviations on trial application is selected for the entire session

**Sessions**: Observations are grouped by local noon-to-noon periods:
- 12:01 today → 12:00 tomorrow = today's session
- This groups all observations from one night together

**Debug output** (`--debug`):
```
debug/
└── 2024_05_20_L/
    ├── light_median.fit           # Combined light preview
    ├── hp_flat_L_20240518_sigma12.3_best.fit  # Best flat
    └── hp_flat_L_20240510_sigma15.7.fit       # Other candidates
```

---

## Cosmetic Correction

Hot pixel correction using coordinate lists.

### Format

```
# cosme.lst
# X Y coordinates (0-indexed)
1234 567
2345 678
3456 789
```

### Auto-discovery

For darks, cosmetic files are found automatically:
- `dark120s.fit` → `cosme120s.lst`
- If not found, most recent `.lst` in dark tree is used

### Correction Method

Hot pixels are replaced with median of surrounding pixels.

---

## Processing Order

Both scripts apply corrections in this order:

1. **Bias subtraction** (if provided)
2. **Dark subtraction** (with optional scaling)
3. **Cosmetic correction**
4. **Flat division** (with multiplier)
5. **NaN/Inf cleanup**
6. **Clamp to original dtype range**

---

## Tips and Best Practices

### Recommended Library Organization

This structure is **optional but recommended** for easier management:

```
CalibrationLibrary/
├── ASI2600MM/
│   ├── darks/
│   │   ├── 2024_05/
│   │   │   ├── dark120s.fit
│   │   │   └── cosme120s.lst
│   │   └── 2024_06/
│   │       └── dark120s.fit
│   └── flats/
│       ├── Telescope_A/
│       │   ├── 2024_05/
│       │   │   └── flat_L.fit
│       │   └── 2024_06/
│       │       └── flat_L.fit
│       └── Telescope_B/
│           └── flat_L.fit
└── maintenance.csv
```

Remember: **folder structure is for your convenience only** — the scripts read all metadata from FITS headers.

### Flat Multiplier

The multiplier (MULT) compensates for flat division:
- Should equal median(FLAT) to preserve image brightness
- `autocalibrate.py` computes this automatically
- For `calibrate.py`, check flat median: `fitshdr flat.fit | grep -i median`

### Checking Results

1. Compare histograms before/after calibration
2. Check for dust donuts (flat mismatch)
3. Verify hot pixels are corrected
4. Look for gradient residuals (dark mismatch)

---

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| "No usable darks found" | Missing EXPTIME/JD headers | Add headers to darks |
| "No usable flats found" | Missing FILTER/JD headers | Add headers to flats |
| Dust donuts visible | Wrong flat selected | Use `--bestflat` or check dates |
| Gradient in result | Dark mismatch | Use `-o` optimization or better dark |
| Hot pixels remain | No cosme file found | Create cosme.lst for dark |
| "Shape mismatch" | Different image dimensions | Ensure all frames same size |
