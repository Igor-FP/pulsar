@echo off
set "SCRIPTTMP=%~dp0..\Sum\sum.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
