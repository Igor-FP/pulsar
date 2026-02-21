@echo off
set "SCRIPTTMP=%~dp0..\Autoflat\autoflat.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
