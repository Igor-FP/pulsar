@echo off
set "SCRIPTTMP=%~dp0..\Mtf\mtf.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
