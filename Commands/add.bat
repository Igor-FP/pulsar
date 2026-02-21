@echo off
set "SCRIPTTMP=%~dp0..\Add\add.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
