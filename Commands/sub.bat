@echo off
set "SCRIPTTMP=%~dp0..\Sub\sub.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
