@echo off
set "SCRIPTTMP=%~dp0..\Binxy\binxy.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
