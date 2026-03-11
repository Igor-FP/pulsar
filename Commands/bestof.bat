@echo off
set "SCRIPTTMP=%~dp0..\Bestof\bestof.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
