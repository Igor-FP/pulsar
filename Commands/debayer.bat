@echo off
set "SCRIPTTMP=%~dp0..\Debayer\debayer.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
