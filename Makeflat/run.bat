@echo off
setlocal

REM Test makeflat.py
REM Usage: place flat frames (and optionally dark frames) in a test directory

set SCRIPT=%~dp0makeflat.py

echo ============================================
echo Testing makeflat.py
echo ============================================

REM Example 1: Process flats from directory (will auto-run makedark if needed)
REM python %SCRIPT% C:\path\to\flats

REM Example 2: Process flats with custom normalization target
REM python %SCRIPT% C:\path\to\flats 10000

REM Example 3: Process flats from file mask
REM python %SCRIPT% flat*.fit

echo.
echo To test, uncomment one of the examples above and provide valid paths.
echo Output files (flat_*.fit) will be created in current directory.
echo.
echo Note: Requires dark frames in the same directory or pre-built
echo       dark*.fit and cosme*.lst files.

endlocal
