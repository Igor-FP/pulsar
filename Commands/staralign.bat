@echo off
set "SCRIPTTMP=%~dp0..\Staralign\staralign.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
