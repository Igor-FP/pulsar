@echo off
set "SCRIPTTMP=%~dp0..\Sortfits\Sortfits.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
