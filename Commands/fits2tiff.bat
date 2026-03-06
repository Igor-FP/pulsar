@echo off
set "SCRIPTTMP=%~dp0..\Fits2tiff\fits2tiff.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
