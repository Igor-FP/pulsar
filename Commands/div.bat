@echo off
set "SCRIPTTMP=%~dp0..\Div\div.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
