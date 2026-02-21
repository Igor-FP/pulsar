@echo off
set "SCRIPTTMP=%~dp0..\Noffset\noffset.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
