@echo off
echo Starting MELTS Calculator...
echo.

REM Check if Python is installed
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo Error: Python is not installed or not in PATH
    echo Please install Python 3.8 or higher from https://www.python.org/
    pause
    exit /b 1
)

REM Check if required packages are installed
echo Checking dependencies...
python -c "import numpy, scipy, pandas, matplotlib" >nul 2>&1
if %errorlevel% neq 0 (
    echo Installing required packages...
    pip install -r requirements.txt
    if %errorlevel% neq 0 (
        echo Error: Failed to install packages
        pause
        exit /b 1
    )
)

REM Run the application
echo.
echo Launching MELTS Calculator...
python melts_calculator.py

if %errorlevel% neq 0 (
    echo.
    echo Error: Failed to start MELTS Calculator
    pause
)

exit /b 0