@echo off
set "SCRIPTTMP=%~dp0..\Newflat\newflat.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*

