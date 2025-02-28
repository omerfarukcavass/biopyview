@echo off
set PROJECT_NAME=BioPyView	
set FOLDER_NAME="biopyview-main"
set PYTHON_VERSION=3.10.11

echo Installing %PROJECT_NAME%...

:: Step 1: Check if Python is installed
where python >nul 2>nul
if %errorlevel% neq 0 (
    echo Error: Python is not installed. Please install Python %PYTHON_VERSION% or later.
    exit /b 1
)

:: Step 2: Create virtual environment
echo Creating virtual environment...
python -m venv venv

:: Step 3: Activate virtual environment and install dependencies
call venv\Scripts\activate
echo Installing dependencies...
pip install --upgrade pip
pip install -r requirements.txt

echo Installation complete.

:: Step 4: Run the app immediately
set /p RUN_APP="Do you want to run the app now? (y/n): "
if /I "%RUN_APP%"=="y" (
    python main.py
)

echo Done!
pause
