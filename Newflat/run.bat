@echo off
REM Test script for newflat.py

echo === Test 1: Add entry with current UTC time ===
python newflat.py --camera "2600MM" --log test_maintenance.csv

echo.
echo === Test 2: Add entry with specific date ===
python newflat.py --camera "2600MM" --log test_maintenance.csv --date "2024-05-18T14:30:00"

echo.
echo === Test 3: Add entry with comment ===
python newflat.py --camera "2600MM" --log test_maintenance.csv --date "2024-12-21T10:00:00" --comment "Changed dust cover, cleaned sensor"

echo.
echo === Contents of test log ===
type test_maintenance.csv

echo.
echo === Cleanup ===
del test_maintenance.csv

