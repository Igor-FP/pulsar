@echo off
setlocal EnableDelayedExpansion

:: Ask user which PATH to modify
echo Add current folder to PATH:
echo [Y] System (requires admin)
echo [N] Current user only
choice /c YN /m "Choose target: "

:: If user pressed N (2) -> USER, otherwise (Y=1) -> SYSTEM
if errorlevel 2 (
    set "TARGET=USER"
) else (
    set "TARGET=SYSTEM"
)

set "CURR=%~dp0"
:: Remove trailing backslash
if "%CURR:~-1%"=="\" set "CURR=%CURR:~0,-1%"

:: If system PATH selected, ensure admin or relaunch elevated
if /I "%TARGET%"=="SYSTEM" (
    :: Check admin rights
    net session >nul 2>&1
    if errorlevel 1 (
        echo Requesting administrative privileges...
        powershell -Command "Start-Process '%~f0' -Verb runAs"
        exit /b
    )
    set "REGKEY=HKLM\SYSTEM\CurrentControlSet\Control\Session Manager\Environment"
) else (
    set "REGKEY=HKCU\Environment"
)

:: Read current PATH (if missing, initialize as empty)
set "OLDPATH="
for /f "tokens=2,*" %%A in ('reg query "%REGKEY%" /v Path 2^>nul') do set "OLDPATH=%%B"

:: Check if current folder is already in PATH
echo !OLDPATH! | find /I "%CURR%" >nul
if not errorlevel 1 (
    echo Already in PATH:
    echo %CURR%
    pause
    exit /b
)

:: Build new PATH value
if defined OLDPATH (
    set "NEWPATH=!OLDPATH!;%CURR%"
) else (
    set "NEWPATH=%CURR%"
)

:: Update registry
reg add "%REGKEY%" /v Path /t REG_EXPAND_SZ /d "!NEWPATH!" /f >nul

:: Broadcast WM_SETTINGCHANGE so running apps pick up the new PATH
powershell -Command "[Environment]::SetEnvironmentVariable('Path',[Environment]::GetEnvironmentVariable('Path','User'),'User')" >nul 2>&1

echo.
echo === DONE ===
echo Added to %TARGET% PATH:
echo %CURR%
echo.
echo Restart your console to apply changes.
pause
