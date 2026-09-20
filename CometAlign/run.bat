@echo off
REM Quick local test for cometalign.py
REM Interactive GUI: mark the comet nucleus on the first frame, press Tab to
REM switch to the last frame, mark it there, then Enter to process.
python "%~dp0cometalign.py" ..\Samples1\*.fit comet0001.fit

REM Headless (no GUI) - supply nucleus pixel coords on the first/last frames:
REM python "%~dp0cometalign.py" ..\Samples1\*.fit comet0001.fit --start 512 480 --stop 530 505
pause
