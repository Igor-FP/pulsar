@echo off
REM Quick local test for makemask.py
python "%~dp0makemask.py" ..\Samples1\*.fit out_mask0001.fit --black 30%% --white 99.5%% --grow 2 --invert
pause
