@echo off
set "SCRIPTTMP=%~dp0..\Autocalibrate\autocalibrate.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*

