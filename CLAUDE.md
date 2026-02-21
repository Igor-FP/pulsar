# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

AstroBatch is a Python 3 suite of command-line tools for batch processing astronomical FITS images. It provides calibration, arithmetic, combination, and analysis operations for CCD/CMOS imaging data.

## Architecture

### Core Design Pattern

All tools share a consistent architecture:
1. Import `batch_utils.py` from `lib/` for file handling
2. Define `usage()` and `parse_args()` for CLI
3. Process input files using astropy.io.fits and numpy
4. Output FITS files with preserved/updated headers

### File Pattern System (Critical)

The `lib/batch_utils.py` module provides unified input/output handling that all tools must use:

**Argument types (checked in this order by `parse_argument()`):**

1. **Numeric constant** (strict format, checked FIRST):
   - Integers: `123`, `-456`, `+789`
   - Decimals: `123.456`, `-0.5`, `+1.0`
   - NOT accepted: `1e5`, `inf`, `nan`, `1_000`, `1,000`

2. **Wildcard mask**: `*.fit`, `light_*.fit`, `img???.fits`

3. **List file**: `@list.txt` or any `.txt`/`.lst` file

4. **Numbered sequence**: `light0001.fit` → auto-discovers light0002.fit, ...

5. **Single file**: `image.fit`

**Output patterns:**
- Single file: `output.fit`
- Numbered pattern: `output0001.fit` (auto-increments)

### Key Functions in batch_utils.py

```python
# Universal argument parser (constant or file list)
parse_argument(arg)            # Returns float OR List[str]
parse_numeric_constant(s)      # Returns float OR None (strict parsing)

# File-only functions
expand_input_spec(spec)        # Converts file spec to List[str]
build_io_file_lists(in, out)   # Creates input→output file pairs

# Operand handling (for arithmetic operations)
build_operand_spec(op, n)      # Parse operand, validate list length against n
get_operand_for_file(spec, i)  # Get operand for file at index i
resolve_operand_value(op, shape, dtype)  # Expand scalar to array or load FITS

# Validation
validate_has_file_input(*specs)  # Ensure at least one arg is a file (not constant)
```

## Key Tools

| Tool | Purpose |
|------|---------|
| add.py, arith.py | Image arithmetic (add/sub/mul/div) |
| sum.py | Stack summation with exposure time handling |
| med.py | Tiled parallel median combine |
| calibrate.py | Dark/bias/flat/cosmetic calibration pipeline |
| normalize.py | Linear regression normalization between frames |
| ngain.py | Normalize by gain (multiply to target median) |
| noffset.py | Normalize by offset (add to target median) |
| autoflat.py | Generate flat fields via polynomial fitting |
| autosolve.py | WCS solving and astrometric rectification (WSL-aware) |
| cosme.py | Hot pixel correction |
| makedark.py | Meta-script: create master darks + cosme lists |
| makeflat.py | Meta-script: create master flats per filter |
| sortfits.py | Organize FITS by metadata |
| fft_align.py | FFT-based image alignment |

## Dependencies

- Python 3.6+
- numpy, astropy, scipy
- Optional: reproject (WCS work), astrometry.net (autosolve.py)

## Running Tools

```bash
# Direct Python
python Add/add.py input.fit output.fit 100

# Via batch wrapper (after running Commands/setup.bat)
add input.fit output.fit 100
```

## Testing

Each tool directory contains `run.bat` test scripts. Test data is in `Samples1/`, `Samples2/`, and `1/` directories.

## Development Guidelines

### Project Structure for New Scripts

1. **Create directory**: `ScriptName/` with main script `scriptname.py`
2. **Create test file**: `ScriptName/run.bat` for local testing
3. **Create command wrapper**: `Commands/scriptname.bat`
4. **Update documentation**: Add to `SCRIPTS.md`

### Script Portability (Critical)

Scripts must run from ANY directory. Use relative path to import batch_utils:

```python
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils
```

### Input/Output Handling

Always use `batch_utils.py` for file handling - never write custom glob/expansion logic:

```python
# Build input→output pairs (handles all spec types automatically)
io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)

for infile, outfile in io_pairs:
    process_file(infile, outfile)
```

Supported input specs: single file, numbered sequence, wildcards, @list.txt

### Numeric Constant Support

Arithmetic scripts support numeric constants as operands. A constant is treated as a "virtual file" filled with that value.

**For arithmetic operations (add, sub, mul, div, arith):**

```python
# Get raw operand (float or file path)
operand_raw = batch_utils.get_operand_for_file(operand_spec, file_index)

# Resolve to array (handles both constants and files)
operand = batch_utils.resolve_operand_value(operand_raw, data.shape, data.dtype)
```

**Use case:** Calibration on thermostable CCDs where BIAS can be replaced with 0 or a known constant offset.

**For new scripts:** Decide whether numeric constants make sense:
- **Support them:** If the operation can logically use a constant (e.g., subtraction, division)
- **Reject them:** If a file is always required, use `expand_input_spec()` instead of `parse_argument()` — it will raise `FileNotFoundError` for numeric strings

**Validation:** When constants are supported, at least one argument must be a file (to get dimensions):
```python
batch_utils.validate_has_file_input(input_spec, operand_spec)  # Raises if all are constants
```

### Data Type Support (Mandatory)

All scripts MUST support the full range of FITS data types:

**Integer types:**
- Signed: int8, int16, int32, int64
- Unsigned: uint8, uint16, uint32, uint64

**Floating-point types:**
- float32, float64

### When to Use float64 (Performance Consideration)

**NOT every operation needs float64 conversion.** Choose based on operation type:

**float64 REQUIRED (precision loss otherwise):**
- Multiplication/division by coefficients close to 1.0 (normalization, gain correction)
- Computing mean/average (accumulation errors)
- Any operation with fractional coefficients
- Dark scaling with optimal K factor
- Flat field division

```python
# Normalization - MUST use float64
work = data.astype(np.float64)
result = work * (target_median / current_median)
```

**float64 NOT NEEDED (keep native dtype for speed):**
- Addition/subtraction of integers or whole images (bias subtraction, dark subtraction)
- Median computation (just selecting an element, no arithmetic)
- Min/max operations
- Comparison operations

```python
# Median combine - NO float64 needed, work with native dtype
stack = np.array([img1, img2, img3])  # Keep original dtype
result = np.median(stack, axis=0)     # Median preserves dtype sense
```

```python
# Simple subtraction - can stay in native dtype if no overflow risk
# But for safety with potential negative results, use int64 or float64
result = img1.astype(np.int64) - img2.astype(np.int64)
```

**float64 RECOMMENDED (borderline cases):**
- Sum of many images (overflow risk for int16/uint16)
- Operations that may produce negative intermediate values from unsigned types

**Decision rule:** If the operation involves division or multiplication by non-integer values, use float64. If it's pure addition/subtraction/median/min/max, consider staying native.

### Output Conversion

**Integer conversion (clamping + rounding):**
```python
if np.issubdtype(orig_dtype, np.integer):
    info = np.iinfo(orig_dtype)
    arr = np.clip(data, info.min, info.max)
    arr = np.rint(arr)  # Round to nearest integer
    return arr.astype(orig_dtype)
```

**Float conversion:**
```python
if np.issubdtype(orig_dtype, np.floating):
    return data.astype(orig_dtype)
```

### NaN/Inf Policy (Critical)

**NEVER write NaN or Inf values to output FITS files.**

Always sanitize before writing:

```python
# Option 1: Replace with zero
work = np.nan_to_num(work, nan=0.0, posinf=0.0, neginf=0.0)

# Option 2: Check and handle explicitly
bad = ~np.isfinite(work)
if np.any(bad):
    work[bad] = 0.0  # or appropriate fill value
```

For division operations, prevent NaN/Inf at source:
```python
# Safe division
safe_divisor = divisor.copy()
safe_divisor[safe_divisor == 0] = np.nan  # Mark zeros
result = dividend / safe_divisor
result = np.nan_to_num(result, nan=0.0)  # Replace NaN with 0
```

### Batch File Wrapper Format

Create `Commands/scriptname.bat`:

```batch
@echo off
set "SCRIPTTMP=%~dp0..\ScriptName\scriptname.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
```

### Other Guidelines

- Preserve FITS headers; update only necessary fields (EXPTIME, HISTORY, etc.)
- Tools are independent - do not import from other tool modules (except cosme.py for calibration)
- Add HISTORY entries documenting applied operations
- Print progress for batch processing (progress bar or counter)
- Use `memmap=False` when reading FITS to avoid issues with BZERO/BSCALE
