@echo off
set "SCRIPTTMP=%~dp0..\Mul\mul.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
