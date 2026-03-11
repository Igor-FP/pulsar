@echo off
set "SCRIPTTMP=%~dp0..\Xisf2fits\xisf2fits.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
