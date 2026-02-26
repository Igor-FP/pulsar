# Creating Master Calibration Frames

This guide covers creating master bias, dark, flat frames and cosmetic correction lists using PULSAR tools.

---

## Quick Start

### Master Bias

```bash
med bias*.fit master_bias.fit
```

### Master Darks + Cosme Lists

```bash
# Auto-create darks grouped by exposure time + generate cosme lists
makedark darks_folder/

# With bias subtraction
makedark darks_folder/ master_bias.fit
```

### Master Flats

```bash
# Auto-create flats grouped by filter (requires darks in current dir)
makeflat flats_folder/

# With custom normalization target
makeflat flats_folder/ 10000
```

### Hot Pixel List (manual)

```bash
# Create cosme list from existing master dark
make_cosme master_dark.fit cosme.lst
```

---

## Calibration Frame Types

| Frame | What it is | How to capture | Processing tool |
|-------|-----------|----------------|-----------------|
| **Bias** | Readout noise pattern | Shortest exposure, covered | `med` |
| **Dark** | Thermal noise (dark current) | Same exposure as lights, covered | `makedark` |
| **Flat** | Optical response (vignetting, dust) | Uniform light, 25-40% histogram | `makeflat` |
| **Cosme** | Hot pixel coordinate list | — | `makedark` (auto) or `make_cosme` |

### Capture Guidelines

**Bias frames:**
- Set the shortest possible exposure (camera minimum)
- Cover the telescope/lens completely
- Take 20-50 frames

**Dark frames:**
- Same exposure time as your light frames
- Same temperature as during imaging session
- Cover the telescope completely
- Take 20-50 frames per exposure time

**Flat frames:**
- Uniform illumination (twilight sky, light panel, or T-shirt over telescope)
- Exposure to reach 25-40% histogram level
  - Closer to sky background level reduces nonlinearity effects
  - Too low = excessive flat noise contribution
- Take 20-50 frames per filter
- Same optical configuration as light frames

---

## Tool Reference

### med.py — Median Combine

Universal tool for combining any FITS images via pixel-wise median.

```
med.py input_spec output.fit [--memory N]
```

| Option | Description |
|--------|-------------|
| `--memory N` | Max memory in GB for stack (default: 16) |

**Input formats:**

```bash
med bias*.fit master_bias.fit      # Wildcard pattern
med bias0001.fit master_bias.fit   # Numbered sequence (auto-discovers all)
med @biaslist.txt master_bias.fit  # List file
```

**When to use:**
- Creating master bias (always)
- Creating master darks when bias subtraction is not needed
- Any generic median stacking task

---

### makedark.py — Master Dark Creator

Automated dark frame processing with bias subtraction and hot pixel detection.

```
makedark.py input_spec [bias_spec]
```

| Argument | Description |
|----------|-------------|
| `input_spec` | Directory, file mask, or list of dark frames |
| `bias_spec` | Optional: bias file, file list, or numeric constant (default: 0) |

**Output files (in current directory):**

| File | Description |
|------|-------------|
| `dark<exp>.fit` | Master dark for each exposure time |
| `cosme<exp>.lst` | Hot pixel list for each exposure time |
| `bias.fit` | Master bias (if bias_spec was file list) |

**How it works:**

1. Scans input for files with `IMAGETYP='Dark Frame'`
2. Groups by exposure time (from `EXPTIME` header)
3. Subtracts bias from each frame (if provided)
4. Median combines each group
5. Detects hot pixels and writes cosme lists (10,000 hottest pixels)

**Examples:**

```bash
makedark darks/                    # Basic (no bias subtraction)
makedark darks/ master_bias.fit   # With master bias file
makedark darks/ bias*.fit          # With bias frames (creates bias.fit)
makedark darks/ 100                # With numeric bias constant
```

**Exposure time format:** `dark120s.fit` (seconds), `dark500ms.fit` (milliseconds)

---

### makeflat.py — Master Flat Creator

Automated flat frame processing with dark subtraction and normalization.

```
makeflat.py input_spec [target_median]
```

| Argument | Description |
|----------|-------------|
| `input_spec` | Directory, file mask, or list of flat frames |
| `target_median` | Normalization target (default: 5000) |

**Output files:** `flat_<filter>.fit` for each filter

**How it works:**

1. Scans input for files with `IMAGETYP='Flat Frame'`
2. Groups by filter (from `FILTER` header)
3. Validates exposure consistency within each group
4. Finds corresponding `dark<exp>.fit` and `cosme<exp>.lst`
5. Subtracts darks, applies cosmetic correction
6. Normalizes each frame to target median
7. Median combines each group

**Prerequisites:** For each flat exposure time, needs `dark<exp>.fit` and `cosme<exp>.lst` in current or input directory.

**Filter naming:** `flat_l.fit` (Luminance), `flat_r.fit`/`flat_g.fit`/`flat_b.fit` (RGB), `flat_h.fit` (H-alpha), `flat_o.fit` (OIII), `flat_s.fit` (SII)

---

### make_cosme.py — Hot Pixel List Creator

Creates cosmetic correction lists by finding the hottest pixels in an image.

**Aliases:** `make_cosme`, `find_hot`

```
make_cosme.py input.fit output.lst
```

Finds 10,000 hottest (brightest) pixels and writes coordinates to `.lst` file.

**Output format (IRIS-compatible):**

The `.lst` format is compatible with **IRIS** (© Christian Buil):

```
P 1234 567
P 2345 678
```

- Each line: `P X Y` (P = pixel marker, X Y = zero-indexed coordinates)
- Lines starting with `#` are comments
- CRLF line endings for Windows compatibility

**When to use:**
- `makedark` creates cosme lists automatically
- Use `make_cosme` manually for existing master darks without cosme files

---

## Required FITS Headers

| Frame Type | Header | Required | Description |
|------------|--------|----------|-------------|
| Dark | `IMAGETYP` | Yes | Must be "Dark Frame" |
| Dark | `EXPTIME` | Yes | Exposure time in seconds |
| Flat | `IMAGETYP` | Yes | Must be "Flat Frame" |
| Flat | `EXPTIME` | Yes | Exposure time in seconds |
| Flat | `FILTER` | Yes | Filter name |
| Bias | — | — | `med.py` works with any FITS files |

---

## Complete Workflow Example

```bash
cd /calibration/2024_05_20/

# Step 1: Create master bias
med bias*.fit bias.fit

# Step 2: Create master darks (also creates cosme lists)
makedark darks/ bias.fit
# Output: dark120s.fit, dark300s.fit, cosme120s.lst, cosme300s.lst

# Step 3: Create master flats
makeflat flats/
# Output: flat_l.fit, flat_r.fit, flat_ha.fit, etc.

# Step 4: Verify results
ls -la *.fit *.lst
fitshdr dark120s.fit
fitshdr flat_l.fit
```

---

## Tips and Best Practices

### Temperature Matching

Dark current depends on sensor temperature:
- Capture darks at the same temperature as lights
- Or use camera cooling to maintain consistent temperature
- Check `CCD-TEMP` header if available

### Sufficient Frame Count

More frames = cleaner master:
- **Minimum**: 15-20 frames
- **Recommended**: 30-50 frames
- **Diminishing returns** after ~50 frames

### Flat Exposure Level

- Target 25-40% of full well capacity
- Closer to sky background level minimizes nonlinearity effects
- Too low increases flat noise contribution
- Consistent exposure within filter group

### Storage Organization

Recommended structure:
```
CalibrationLibrary/
├── 2024_05_20/
│   ├── bias.fit
│   ├── dark120s.fit
│   ├── dark300s.fit
│   ├── cosme120s.lst
│   ├── cosme300s.lst
│   ├── flat_l.fit
│   └── flat_ha.fit
└── 2024_06_15/
    └── ...
```

---

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| "No dark frames found" | Missing IMAGETYP header | Add `IMAGETYP='Dark Frame'` to headers |
| "No flat frames found" | Missing IMAGETYP header | Add `IMAGETYP='Flat Frame'` to headers |
| "Inconsistent exposures" | Mixed exposure times in flat group | Separate flats by exposure time |
| "Dark not found for exposure" | Missing master dark | Run `makedark` first |
| Hot pixels remain | Cosme list not found | Check cosme file exists and naming |
| Gradient in flat | Uneven illumination | Reshoot flats with better light source |
