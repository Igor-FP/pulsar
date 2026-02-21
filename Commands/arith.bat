@echo off
set "SCRIPTTMP=%~dp0..\Arith\arith.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
