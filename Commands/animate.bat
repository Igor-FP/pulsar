@echo off
set "SCRIPTTMP=%~dp0..\Animate\animate.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
