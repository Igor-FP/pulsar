@echo off
set "SCRIPTTMP=%~dp0..\CometAlign\cometalign.py"

echo Running %SCRIPTTMP%
python "%SCRIPTTMP%" %*
