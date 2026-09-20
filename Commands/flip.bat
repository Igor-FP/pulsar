@echo off
set "SCRIPTTMP=%~dp0..\Flip\flip.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
