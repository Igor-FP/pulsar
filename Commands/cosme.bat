@echo off
set "SCRIPTTMP=%~dp0..\Cosme\cosme.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
