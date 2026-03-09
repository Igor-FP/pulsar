@echo off
set "SCRIPTTMP=%~dp0..\Rgbbalance\rgbbalance.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
