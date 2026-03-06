@echo off
set "SCRIPTTMP=%~dp0..\Tiff2fits\tiff2fits.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
