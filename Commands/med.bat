@echo off
set "SCRIPTTMP=%~dp0..\Median\med.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
