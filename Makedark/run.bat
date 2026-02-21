@echo off
setlocal

REM Test makedark.py
REM Usage: place dark frames in a test directory and run this script

set SCRIPT=%~dp0makedark.py

echo ============================================
echo Testing makedark.py
echo ============================================

REM Example 1: Process darks from directory with bias=0 (default)
REM python %SCRIPT% C:\path\to\darks

REM Example 2: Process darks with numeric bias constant
REM python %SCRIPT% C:\path\to\darks 1024

REM Example 3: Process darks with master bias file
REM python %SCRIPT% C:\path\to\darks master_bias.fit

REM Example 4: Process darks with bias file list (creates bias.fit)
REM python %SCRIPT% dark0001.fit bias0001.fit

echo.
echo To test, uncomment one of the examples above and provide valid paths.
echo Output files (darkXXXs.fit, cosmeXXXs.lst) will be created in current directory.

endlocal
