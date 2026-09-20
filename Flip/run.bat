@echo off
REM Quick local test for flip.py
REM --y = reverse the Y coordinate (Ynew = H - Yold - 1), i.e. top <-> bottom.
python "%~dp0flip.py" ..\Samples1\*.fit flipped0001.fit --y
pause
