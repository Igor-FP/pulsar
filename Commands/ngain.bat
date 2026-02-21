@echo off
set "SCRIPTTMP=%~dp0..\Ngain\ngain.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
