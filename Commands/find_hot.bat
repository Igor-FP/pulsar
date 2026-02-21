@echo off
set "SCRIPTTMP=%~dp0..\Cosme\make_cosme.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
