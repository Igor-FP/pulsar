@echo off
set "SCRIPTTMP=%~dp0..\Lrgb\lrgb.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
