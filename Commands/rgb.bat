@echo off
set "SCRIPTTMP=%~dp0..\Rgb\rgb.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
