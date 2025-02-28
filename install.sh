#!/bin/bash

# Set project name and GitHub repo
PROJECT_NAME="BioPyView"
FOLDER_NAME="biopyview-main"
PYTHON_VERSION="3.10.11"

echo "Installing $PROJECT_NAME..."

# Step 1: Check if Python 3.10.11+ is installed
PYTHON_CMD=$(command -v python3 || command -v python)
if [ -z "$PYTHON_CMD" ]; then
    echo "Error: Python is not installed. Please install Python 3.10.11+."
    exit 1
fi

PYTHON_VERSION_INSTALLED=$($PYTHON_CMD -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
if [[ "$PYTHON_VERSION_INSTALLED" < "$PYTHON_VERSION" ]]; then
    echo "Error: Python $PYTHON_VERSION or later is required. You have $PYTHON_VERSION_INSTALLED."
    exit 1
fi

# Step 2: Create virtual environment
echo "Creating virtual environment..."
$PYTHON_CMD -m venv venv

# Step 3: Activate virtual environment
source venv/bin/activate

# Step 4: Install dependencies
echo "Installing dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

echo "Installation complete."

# Step 5: Run the app immediately
read -p "Do you want to run the app now? (y/n): " RUN_APP
if [[ "$RUN_APP" == "y" ]]; then
    python main.py
fi

echo "Done!"
