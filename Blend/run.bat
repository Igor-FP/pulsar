@echo off
REM Quick local test for blend.py
REM 1) build a mask from the frames, 2) fade the masked regions toward 0
REM    (operand = constant 0), softening the mask with an MTF first.
python "%~dp0..\MakeMask\makemask.py" ..\Samples1\*.fit mask0001.fit --black 40%% --white 99%%
python "%~dp0blend.py" ..\Samples1\*.fit out_blend0001.fit mask0001.fit 0 --mtf 0.3
pause
