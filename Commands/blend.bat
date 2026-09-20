@echo off
set "SCRIPTTMP=%~dp0..\Blend\blend.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
