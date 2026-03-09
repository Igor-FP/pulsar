@echo off
setlocal EnableDelayedExpansion
chcp 65001 >nul 2>&1
title PULSAR Setup

echo.
echo  ============================================
echo   PULSAR - Setup
echo  ============================================
echo.

:: Track what was done
set "DID_PYTHON=0"
set "DID_PIP=0"
set "DID_DEPS=0"
set "DID_PATH=0"

:: ============================================
:: Step 1: Check Python
:: ============================================
echo [1/4] Checking Python...

python --version >nul 2>&1
if errorlevel 1 (
    echo.
    echo  ERROR: Python not found.
    echo  Download from: https://www.python.org/downloads/
    echo  IMPORTANT: Check "Add Python to PATH" during installation.
    echo.
    pause
    exit /b 1
)

:: Parse version
for /f "tokens=2 delims= " %%V in ('python --version 2^>^&1') do set "PYVER=%%V"
for /f "tokens=1,2 delims=." %%A in ("%PYVER%") do (
    set "PYMAJOR=%%A"
    set "PYMINOR=%%B"
)

if !PYMAJOR! LSS 3 (
    echo  ERROR: Python !PYVER! is too old. Python 3.6+ is required.
    echo  Download from: https://www.python.org/downloads/
    echo.
    pause
    exit /b 1
)
if !PYMAJOR! EQU 3 if !PYMINOR! LSS 6 (
    echo  ERROR: Python !PYVER! is too old. Python 3.6+ is required.
    echo  Download from: https://www.python.org/downloads/
    echo.
    pause
    exit /b 1
)

echo  Found Python !PYVER!
set "DID_PYTHON=1"

:: ============================================
:: Step 2: Check pip
:: ============================================
echo.
echo [2/4] Checking pip...

pip --version >nul 2>&1
if errorlevel 1 (
    echo  pip not found, installing...
    python -m ensurepip --default-pip >nul 2>&1
    pip --version >nul 2>&1
    if errorlevel 1 (
        echo.
        echo  ERROR: Could not install pip.
        echo  Try manually: python -m ensurepip --default-pip
        echo.
        pause
        exit /b 1
    )
    echo  pip installed successfully.
    set "DID_PIP=1"
) else (
    for /f "tokens=1,2 delims= " %%A in ('pip --version 2^>^&1') do set "PIPVER=%%B"
    echo  Found pip !PIPVER!
)

:: ============================================
:: Step 3: Install dependencies
:: ============================================
echo.
echo [3/4] Checking dependencies...

set "NEED_INSTALL=0"
python -c "import numpy" >nul 2>&1
if errorlevel 1 set "NEED_INSTALL=1"
python -c "import astropy" >nul 2>&1
if errorlevel 1 set "NEED_INSTALL=1"

if !NEED_INSTALL! EQU 1 (
    echo  Installing required packages...
    echo.
    pip install -r "%~dp0requirements.txt"
    if errorlevel 1 (
        echo.
        echo  WARNING: Some packages may have failed to install.
        echo  You can retry manually: pip install -r requirements.txt
        echo.
    ) else (
        echo.
        echo  Dependencies installed successfully.
    )
    set "DID_DEPS=1"
) else (
    echo  Core dependencies already installed.
)

:: ============================================
:: Step 4: Add Commands/ to PATH
:: ============================================
echo.
echo [4/4] Checking PATH...

set "CMDDIR=%~dp0Commands"
:: Remove trailing backslash if present
if "!CMDDIR:~-1!"=="\" set "CMDDIR=!CMDDIR:~0,-1!"

:: Check if already in PATH
echo !PATH! | find /I "%CMDDIR%" >nul
if not errorlevel 1 (
    echo  Commands/ already in PATH.
    goto :summary
)

echo  Commands/ is not in PATH.
echo.
echo  Add to PATH?
echo  [Y] System PATH (requires admin)
echo  [U] User PATH only
echo  [N] Skip
choice /c YUN /m "Choose: "

if errorlevel 3 goto :summary
if errorlevel 2 (
    set "TARGET=USER"
    set "REGKEY=HKCU\Environment"
) else (
    set "TARGET=SYSTEM"
    set "REGKEY=HKLM\SYSTEM\CurrentControlSet\Control\Session Manager\Environment"
)

:: If system PATH selected, ensure admin
if /I "!TARGET!"=="SYSTEM" (
    net session >nul 2>&1
    if errorlevel 1 (
        echo.
        echo  Admin rights required for System PATH.
        echo  Falling back to User PATH...
        set "TARGET=USER"
        set "REGKEY=HKCU\Environment"
    )
)

:: Read current registry PATH
set "OLDPATH="
for /f "tokens=2,*" %%A in ('reg query "!REGKEY!" /v Path 2^>nul') do set "OLDPATH=%%B"

:: Double-check it's not already there (registry PATH vs runtime PATH)
echo !OLDPATH! | find /I "%CMDDIR%" >nul
if not errorlevel 1 (
    echo  Already in !TARGET! PATH.
    goto :summary
)

:: Build and set new PATH
if defined OLDPATH (
    set "NEWPATH=!OLDPATH!;%CMDDIR%"
) else (
    set "NEWPATH=%CMDDIR%"
)

reg add "!REGKEY!" /v Path /t REG_EXPAND_SZ /d "!NEWPATH!" /f >nul
if errorlevel 1 (
    echo  ERROR: Failed to update PATH.
    goto :summary
)

:: Broadcast WM_SETTINGCHANGE
powershell -Command "[Environment]::SetEnvironmentVariable('Path',[Environment]::GetEnvironmentVariable('Path','User'),'User')" >nul 2>&1

echo  Added to !TARGET! PATH: %CMDDIR%
echo  Restart your console to apply.
set "DID_PATH=1"

:: ============================================
:: Step 5: Summary
:: ============================================
:summary
echo.
echo  ============================================
echo   Setup complete
echo  ============================================
echo.
if !DID_PYTHON! EQU 1 echo  [OK] Python !PYVER!
if !DID_PIP! EQU 1 (
    echo  [OK] pip installed
) else (
    echo  [OK] pip
)
if !DID_DEPS! EQU 1 (
    echo  [OK] Dependencies installed
) else (
    echo  [OK] Dependencies (already present)
)
if !DID_PATH! EQU 1 (
    echo  [OK] PATH updated
) else (
    echo  [--] PATH (unchanged)
)
echo.
echo  Optional: install astrometry.net for plate solving
echo  (see SCRIPTS.md, autosolve section)
echo.
pause
