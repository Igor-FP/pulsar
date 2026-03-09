@echo off
set "SCRIPTTMP=%~dp0..\Raw2fits\raw2fits.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
