@echo off
set "SCRIPTTMP=%~dp0..\MakeMask\makemask.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
