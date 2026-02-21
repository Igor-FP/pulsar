@echo off
set "SCRIPTTMP=%~dp0..\Autosolve\autosolve.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
