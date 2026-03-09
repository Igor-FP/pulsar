@echo off
set "SCRIPTTMP=%~dp0..\Hotfix\hotfix.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
