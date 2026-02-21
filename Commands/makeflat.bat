@echo off
set "SCRIPTTMP=%~dp0..\Makeflat\makeflat.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
