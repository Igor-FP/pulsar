@echo off
set "SCRIPTTMP=%~dp0..\Calibrate\calibrate.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
