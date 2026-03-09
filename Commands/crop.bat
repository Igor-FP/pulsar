@echo off
set "SCRIPTTMP=%~dp0..\Crop\crop.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
