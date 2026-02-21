@echo off
set "SCRIPTTMP=%~dp0..\Makedark\makedark.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
