#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Shared helpers for FITS batch-processing tools.

Utilities for file handling, input expansion, operand resolution,
fast median computation, and optimized median combine operations.
"""

import os
import re
import glob
from typing import List, Tuple, Optional, Union


# ---------------------------------------------------------
# Numbered sequence detection
# ---------------------------------------------------------

_SEQ_RE = re.compile(r"^(.*?)(\d+)(\.(?:fit|fits))$", re.IGNORECASE)


def extract_index_from_filename(path: str) -> Optional[int]:
    """
    Extract numeric index from a file name like 'image0123.fit'.
    Returns integer index or None if no numeric part.
    """
    base = os.path.basename(path)
    m = _SEQ_RE.match(base)
    if not m:
        return None
    return int(m.group(2))


def find_numbered_sequence(first_path: str) -> List[Tuple[str, int, str]]:
    """
    Discover a contiguous numbered sequence starting from first_path.
    Example:
        find_numbered_sequence("light0001.fit") ->
        [
            ("light0001.fit", 1, ".fit"),
            ("light0002.fit", 2, ".fit"),
            ...
        ]

    Returns empty list if first_path does not match numbered pattern.
    """
    base = os.path.basename(first_path)
    m = _SEQ_RE.match(base)
    if not m:
        return []

    prefix, digits, ext = m.group(1), m.group(2), m.group(3)
    width = len(digits)
    start_index = int(digits)

    dir_ = os.path.dirname(os.path.abspath(first_path)) or "."
    seq: List[Tuple[str, int, str]] = []

    idx = start_index
    while True:
        name = f"{prefix}{str(idx).zfill(width)}{ext}"
        full = os.path.join(dir_, name)
        if not os.path.isfile(full):
            break
        seq.append((full, idx, ext.lower()))
        idx += 1

    return seq


# ---------------------------------------------------------
# Input / output mapping helpers
# ---------------------------------------------------------

def _expand_list_file(path: str) -> List[str]:
    """Internal: read paths from a simple text file."""
    lst: List[str] = []
    base_dir = os.path.dirname(os.path.abspath(path)) or "."

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if not os.path.isabs(s):
                s = os.path.join(base_dir, s)
            lst.append(s)
    if not lst:
        raise ValueError(f"No entries found in list file: {path}")
    return lst


def has_wildcards(s: str) -> bool:
    """Check if string contains glob wildcard characters."""
    return any(c in s for c in '*?[')


def expand_input_spec(spec: str) -> List[str]:
    """
    Expand input specification:
      - wildcard masks: *.fit, img???.fits, etc.
      - @list.txt / list.txt / list.lst
      - numbered patterns
      - single file

    Wildcards:
      If the spec contains '*' or '?', glob.glob() is used and all matching
      files are returned (sorted). If no files match, FileNotFoundError is raised.
    """
    spec = spec.strip()

    # Wildcard mask: first priority (as requested)
    if ("*" in spec) or ("?" in spec):
        matches = glob.glob(spec)
        files = [os.path.abspath(m) for m in matches if os.path.isfile(m)]
        files.sort()
        if not files:
            raise FileNotFoundError(f"No files match wildcard pattern: {spec}")
        return files

    # @list.txt
    if spec.startswith("@"):
        return _expand_list_file(spec[1:])

    # list.txt / list.lst
    if spec.lower().endswith((".txt", ".lst")) and os.path.isfile(spec):
        return _expand_list_file(spec)

    # Explicit single file: =filename (bypass sequence detection)
    if spec.startswith("="):
        single = spec[1:]
        if os.path.isfile(single):
            return [os.path.abspath(single)]
        raise FileNotFoundError(f"File not found: {single}")

    # Single file or numbered sequence
    if os.path.isfile(spec):
        # Try numbered sequence (need at least 2 files to be a sequence)
        seq = find_numbered_sequence(spec)
        if len(seq) > 1:
            return [name for (name, _, _) in seq]

        # Single file
        return [os.path.abspath(spec)]

    raise FileNotFoundError(f"Input specification not recognized or file not found: {spec}")


def build_io_file_lists(
    input_spec: str,
    output_spec: str
) -> List[Tuple[str, str]]:
    """
    Build list of (input_file, output_file) pairs.

    Input logic:
        follows expand_input_spec()

    Output logic:
        - If output_spec is a single file (no numbered pattern), and
          exactly one input file exists → single output file.
        - If output_spec has a numbered pattern → produce matching number
          of output files; indices increment from the pattern start index.

    Valid numbered pattern:
        prefix + digits + .fit(s)
        e.g. "out0001.fit" or "result_0050.fits"
    """

    inputs = expand_input_spec(input_spec)
    if not inputs:
        raise ValueError("No input files found.")

    # Check if output_spec is a directory (existing dir, or ends with separator)
    is_dir = os.path.isdir(output_spec) or output_spec.endswith(os.sep) or output_spec.endswith('/')
    if is_dir:
        # Output to directory, preserving input filenames
        out_dir = os.path.abspath(output_spec)
        os.makedirs(out_dir, exist_ok=True)
        return [(inp, os.path.join(out_dir, os.path.basename(inp))) for inp in inputs]

    # Single input -> single output
    if len(inputs) == 1:
        return [(inputs[0], output_spec)]

    # Multiple inputs => output must be numbered pattern
    base = os.path.basename(output_spec)
    m = _SEQ_RE.match(base)
    if not m:
        raise ValueError(
            "Output pattern must contain a numeric field when multiple "
            "input files are provided (e.g. out0001.fit), or specify a directory."
        )

    prefix, digits, ext = m.group(1), m.group(2), m.group(3)
    width = len(digits)
    start_index = int(digits)
    out_dir = os.path.dirname(os.path.abspath(output_spec)) or "."

    out: List[Tuple[str, str]] = []
    idx = start_index
    for inp in inputs:
        fname = f"{prefix}{str(idx).zfill(width)}{ext}"
        out.append((inp, os.path.join(out_dir, fname)))
        idx += 1

    return out


def build_io_file_lists_from_list(
    input_files: List[str],
    output_spec: str
) -> List[Tuple[str, str]]:
    """
    Build list of (input_file, output_file) pairs from an already-expanded
    list of input files. Same output logic as build_io_file_lists().
    """
    if not input_files:
        raise ValueError("No input files provided.")

    # Directory output
    is_dir = os.path.isdir(output_spec) or output_spec.endswith(os.sep) or output_spec.endswith('/')
    if is_dir:
        out_dir = os.path.abspath(output_spec)
        os.makedirs(out_dir, exist_ok=True)
        return [(inp, os.path.join(out_dir, os.path.basename(inp))) for inp in input_files]

    # Single input -> single output
    if len(input_files) == 1:
        return [(input_files[0], output_spec)]

    # Multiple inputs => numbered pattern
    base = os.path.basename(output_spec)
    m = _SEQ_RE.match(base)
    if not m:
        raise ValueError(
            "Output pattern must contain a numeric field when multiple "
            "input files are provided (e.g. out0001.fit), or specify a directory."
        )

    prefix, digits, ext = m.group(1), m.group(2), m.group(3)
    width = len(digits)
    start_index = int(digits)
    out_dir = os.path.dirname(os.path.abspath(output_spec)) or "."

    out: List[Tuple[str, str]] = []
    idx = start_index
    for inp in input_files:
        fname = f"{prefix}{str(idx).zfill(width)}{ext}"
        out.append((inp, os.path.join(out_dir, fname)))
        idx += 1

    return out


# ---------------------------------------------------------
# Numeric constant parsing (strict)
# ---------------------------------------------------------

# Strict numeric constant pattern:
#   - Optional sign: + or -
#   - Integer part: one or more digits
#   - Optional fractional part: dot + one or more digits
# Does NOT match: exponential (1e5), inf, nan, separators (1_000, 1,000)
_NUMERIC_CONST_RE = re.compile(r"^[+-]?\d+(?:\.\d+)?$")


def parse_numeric_constant(s: str) -> Optional[float]:
    """
    Parse a strict numeric constant.

    Accepted formats:
        123, -456, +789          (integers)
        123.456, -0.5, +1.0      (decimals with dot)

    NOT accepted:
        1e5, 1E-3                (exponential notation)
        inf, -inf, nan           (special values)
        1_000, 1,000             (digit separators)
        .5, 5.                   (missing integer/fractional part)

    Returns:
        float value if valid constant, None otherwise.
    """
    s = s.strip()
    if not _NUMERIC_CONST_RE.match(s):
        return None
    return float(s)


# ---------------------------------------------------------
# Operand helpers (for arith.py)
# ---------------------------------------------------------

def detect_scalar_operand(s: str) -> Optional[float]:
    """
    Try to interpret operand as a scalar (strict numeric constant).
    Returns float or None.

    Uses strict parsing - only plain numbers with optional sign and decimal point.
    Does not accept exponential notation, inf, nan, etc.
    """
    return parse_numeric_constant(s)


def parse_argument(arg: str) -> Union[float, List[str]]:
    """
    Universal argument parser. Checks in order:

    1. Numeric constant (strict format: 123, -45.67, +0.5)
       → returns float

    2. Wildcard mask (*.fit, img???.fits)
       → returns sorted list of matching file paths

    3. List file (@list.txt or file with .txt/.lst extension)
       → returns list of paths from file

    4. Numbered sequence (light0001.fit → finds light0002.fit, ...)
       → returns list of all files in sequence

    5. Single file
       → returns list with one element

    Raises FileNotFoundError if nothing matches.

    Returns:
        float           → numeric constant
        list[str]       → list of file paths
    """
    arg = arg.strip()

    # 1. Numeric constant (FIRST priority)
    const = parse_numeric_constant(arg)
    if const is not None:
        return const

    # 2. Wildcard mask
    if ("*" in arg) or ("?" in arg):
        matches = glob.glob(arg)
        files = [os.path.abspath(m) for m in matches if os.path.isfile(m)]
        files.sort()
        if not files:
            raise FileNotFoundError(f"No files match wildcard pattern: {arg}")
        return files

    # 3. List file (@list.txt or .txt/.lst extension)
    if arg.startswith("@"):
        return _expand_list_file(arg[1:])

    if arg.lower().endswith((".txt", ".lst")) and os.path.isfile(arg):
        return _expand_list_file(arg)

    # 4. Single file or numbered sequence
    if os.path.isfile(arg):
        seq = find_numbered_sequence(arg)
        if seq:
            return [name for (name, _, _) in seq]
        return [os.path.abspath(arg)]

    raise FileNotFoundError(f"Argument not recognized or file not found: {arg}")


def build_operand_spec(
    op_string: str,
    num_inputs: Optional[int] = None
) -> Union[float, List[str]]:
    """
    Interpret operand: scalar, list file, FITS sequence, or wildcard mask.

    Args:
        op_string: operand specification string
        num_inputs: optional number of input files for validation

    Returns:
        float           → scalar operand (applies to all inputs)
        list[str]       → list of FITS paths (must match num_inputs if provided)

    Raises:
        ValueError: if operand is a file list with wrong length
    """
    result = parse_argument(op_string)

    # Validate list length if num_inputs provided
    if isinstance(result, list) and num_inputs is not None:
        if len(result) == 1:
            # Single file operand applies to all inputs (broadcast)
            pass
        elif len(result) != num_inputs:
            raise ValueError(
                f"Operand file count ({len(result)}) does not match "
                f"input file count ({num_inputs}). "
                f"Use a single file or matching sequence."
            )

    return result


def get_operand_for_file(
    operand_spec: Union[float, List[str]],
    index: int
) -> Union[float, str]:
    """
    Given the operand spec (scalar or list),
    return the operand for file at position 'index'.

    For single-file operand list, returns the same file for all indices (broadcast).
    """
    if isinstance(operand_spec, (int, float)):
        return float(operand_spec)

    if len(operand_spec) == 1:
        # Broadcast single file to all inputs
        return operand_spec[0]

    if index < 0 or index >= len(operand_spec):
        raise IndexError(
            f"Operand index {index} out of range (have {len(operand_spec)} operand files)."
        )
    return operand_spec[index]


# ---------------------------------------------------------
# Operand resolution (creates numpy arrays)
# ---------------------------------------------------------

def resolve_operand_value(
    operand: Union[float, str],
    reference_shape: Tuple[int, int],
    reference_dtype: "np.dtype"
) -> "np.ndarray":
    """
    Resolve operand to a numpy array.

    Args:
        operand: scalar float OR path to FITS file
        reference_shape: (height, width) for scalar expansion
        reference_dtype: dtype for scalar array (used for type hints only,
                        actual computation should be in float64)

    Returns:
        numpy array of shape reference_shape

    If operand is a scalar, creates array filled with that value.
    If operand is a file path, reads the FITS data.

    Note: This function requires numpy and astropy.io.fits to be imported
    by the calling module. It returns data as-is from file, or as float64
    for scalar constants.
    """
    import numpy as np
    from astropy.io import fits

    if isinstance(operand, (int, float)):
        # Scalar constant → create filled array
        return np.full(reference_shape, float(operand), dtype=np.float64)

    # File path → read FITS
    if not os.path.isfile(operand):
        raise FileNotFoundError(f"Operand file not found: {operand}")

    with fits.open(operand, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"Operand file '{operand}' has no primary image data.")
        data = hdul[0].data

    if data.shape != reference_shape:
        raise ValueError(
            f"Operand image shape {data.shape} does not match "
            f"reference shape {reference_shape}."
        )

    return data


def validate_has_file_input(*specs: Union[float, List[str]]) -> None:
    """
    Validate that at least one argument is a file (not a numeric constant).

    This is required because we need at least one file to determine:
    - Image dimensions (shape)
    - Data type
    - FITS header for output

    Args:
        *specs: parsed argument specifications (from parse_argument)

    Raises:
        ValueError: if ALL arguments are numeric constants
    """
    for spec in specs:
        if isinstance(spec, list):
            return  # Found at least one file spec

    raise ValueError(
        "At least one argument must be a file (not a numeric constant). "
        "Cannot determine image dimensions and metadata from constants alone."
    )


# ---------------------------------------------------------
# Fast median computation
# ---------------------------------------------------------

def _partition_median(arr: "np.ndarray") -> float:
    """
    Compute exact median using numpy.partition (C implementation of introselect).
    O(n) time complexity.

    Args:
        arr: 1D numpy array

    Returns:
        Median value as float
    """
    import numpy as np

    n = len(arr)
    k = n // 2

    if n % 2 == 1:
        # Odd length: single middle element
        return float(np.partition(arr, k)[k])
    else:
        # Even length: average of two middle elements
        # partition with two indices is still O(n)
        part = np.partition(arr, [k - 1, k])
        return float(part[k - 1] + part[k]) / 2.0


def fast_median(
    data: "np.ndarray",
    sample_size: int = 100000,
    exact_threshold: int = 1000000
) -> float:
    """
    Compute median of array, using sampling for very large arrays.

    Uses numpy.partition (C-implemented introselect, O(n)) for exact median.
    For very large arrays, takes a random sample first for speed.

    Args:
        data: numpy array (any shape, will be flattened)
        sample_size: number of elements to sample for large arrays (default 100000)
        exact_threshold: use exact median below this size (default 1000000)

    Returns:
        Median value as float

    Performance:
        - Uses np.partition which is O(n) and implemented in C
        - For 16 megapixel image: ~50-100ms exact, ~5-10ms sampled
        - Sampled accuracy: typically within 0.05% of true median
    """
    import numpy as np

    flat = data.ravel()
    n = len(flat)

    if n == 0:
        return 0.0

    if n == 1:
        return float(flat[0])

    # For arrays up to threshold, compute exact median with np.partition
    if n <= exact_threshold:
        return _partition_median(flat)

    # For very large arrays, use random sampling
    # This gives excellent accuracy for image normalization
    rng = np.random.default_rng()
    sample_indices = rng.choice(n, size=min(sample_size, n), replace=False)
    sample = flat[sample_indices]

    return _partition_median(sample)


# ---------------------------------------------------------
# Fast median combine (optimized for I/O)
# ---------------------------------------------------------

def _load_file_data(args):
    """Worker: load entire file data. Returns (index, data_array)."""
    from astropy.io import fits
    import numpy as np

    idx, filepath, ref_shape = args
    with fits.open(filepath, memmap=False) as hdul:
        data = hdul[0].data
        if data is None:
            raise ValueError(f"No data in primary HDU: {filepath}")
        if data.shape != ref_shape:
            raise ValueError(f"Shape mismatch: {filepath} has {data.shape}, expected {ref_shape}")
        return idx, data.astype(np.float64)


def _load_file_strip(args):
    """Worker: load horizontal strip from file. Returns (index, strip_array)."""
    from astropy.io import fits
    import numpy as np

    idx, filepath, y0, y1, width = args
    with fits.open(filepath, memmap=False) as hdul:
        data = hdul[0].data
        if data is None:
            raise ValueError(f"No data in primary HDU: {filepath}")
        strip = data[y0:y1, :].astype(np.float64)
        return idx, strip


def fast_median_combine(
    file_list: "List[str]",
    ref_shape: "Tuple[int, int]",
    ref_dtype: "np.dtype",
    max_memory_gb: float = 16.0,
    max_workers: int = None,
    progress_callback=None
) -> "np.ndarray":
    """
    Fast median combine using parallel I/O and single-pass processing.

    Strategy:
    1. If all data fits in memory: load all files in parallel, compute median
    2. If too large: process in horizontal strips (each file read once per strip,
       but strips are sequential for better caching)

    Args:
        file_list: List of FITS file paths
        ref_shape: Expected image shape (height, width)
        ref_dtype: Original data type (for reference only, computation in float64)
        max_memory_gb: Maximum memory to use for stack (default 4 GB)
        max_workers: Max threads for parallel loading (default: CPU count)
        progress_callback: Optional callback(current, total, message) for progress

    Returns:
        Median image as float64 numpy array

    Mathematical correctness:
        - Uses np.partition for O(n) exact median
        - For even n: median = (v[n//2-1] + v[n//2]) / 2
        - All computations in float64 for precision
    """
    import numpy as np
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import os

    n_files = len(file_list)
    if n_files == 0:
        raise ValueError("No files to combine")

    height, width = ref_shape
    pixels_per_image = height * width
    bytes_per_pixel = 8  # float64

    # Calculate memory requirements
    stack_bytes = n_files * pixels_per_image * bytes_per_pixel
    max_bytes = max_memory_gb * 1024 * 1024 * 1024

    if max_workers is None:
        max_workers = os.cpu_count() or 4

    def report(current, total, msg):
        if progress_callback:
            progress_callback(current, total, msg)

    # Full-frame mode: load all into memory
    if stack_bytes <= max_bytes:
        report(0, n_files, "Loading files (full-frame)")

        stack = np.empty((n_files, height, width), dtype=np.float64)
        tasks = [(i, f, ref_shape) for i, f in enumerate(file_list)]

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_load_file_data, t) for t in tasks]
            done = 0
            for future in as_completed(futures):
                idx, data = future.result()
                stack[idx] = data
                done += 1
                report(done, n_files, "Loading")

        # Compute median using partition
        report(0, 1, "Computing median")
        if n_files % 2 == 1:
            k = n_files // 2
            median = np.partition(stack, k, axis=0)[k]
        else:
            k1, k2 = n_files // 2 - 1, n_files // 2
            part = np.partition(stack, [k1, k2], axis=0)
            median = (part[k1] + part[k2]) / 2.0

        report(1, 1, "Median complete")
        return median

    # Strip-based mode: process in horizontal strips
    # Calculate optimal strip height to fit in memory
    strip_bytes_per_row = n_files * width * bytes_per_pixel
    max_strip_rows = int(max_bytes / strip_bytes_per_row)
    strip_height = max(1, min(max_strip_rows, height))

    n_strips = (height + strip_height - 1) // strip_height
    report(0, n_strips, f"Processing in {n_strips} strips")

    median = np.empty((height, width), dtype=np.float64)

    for strip_idx in range(n_strips):
        y0 = strip_idx * strip_height
        y1 = min(y0 + strip_height, height)
        strip_h = y1 - y0

        # Load this strip from all files in parallel
        strip_stack = np.empty((n_files, strip_h, width), dtype=np.float64)
        tasks = [(i, f, y0, y1, width) for i, f in enumerate(file_list)]

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_load_file_strip, t) for t in tasks]
            for future in as_completed(futures):
                idx, strip_data = future.result()
                strip_stack[idx] = strip_data

        # Compute median for this strip
        if n_files % 2 == 1:
            k = n_files // 2
            median[y0:y1, :] = np.partition(strip_stack, k, axis=0)[k]
        else:
            k1, k2 = n_files // 2 - 1, n_files // 2
            part = np.partition(strip_stack, [k1, k2], axis=0)
            median[y0:y1, :] = (part[k1] + part[k2]) / 2.0

        report(strip_idx + 1, n_strips, f"Strip {strip_idx + 1}/{n_strips}")

    return median
