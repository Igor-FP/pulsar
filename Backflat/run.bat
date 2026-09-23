@echo off
REM Supply an RGB FITS image and a matching starless RGB FITS image.
REM Default outputs: output.fit, background.fit and back_mask.fit.
REM Example: run.bat input.fit output.fit --starless starless.fit
REM Custom names: add --out-back sky.fit --out-mask objects.fit
REM Existing outputs: add --overwrite or -y
REM Temporary .backflat-cache files and its empty directory are removed at exit.
REM Headless with saved default mask: add --no-gui --mask back_mask.fit -y
REM Numerical tests: python -B -m unittest discover -s Backflat -p test_backflat.py
python "%~dp0backflat.py" %*
