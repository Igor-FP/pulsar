@echo off
set "SCRIPTTMP=%~dp0..\Stack\stack.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
