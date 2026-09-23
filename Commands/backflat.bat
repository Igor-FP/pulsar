@echo off
set "SCRIPTTMP=%~dp0..\Backflat\backflat.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
