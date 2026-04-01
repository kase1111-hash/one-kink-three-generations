@echo off
title MKGF Calculation Validator
color 0B

echo.
echo  ============================================================
echo   One Kink, Three Generations — Calculation Validator
echo  ============================================================
echo.
echo   This checks every key formula against known analytical
echo   results. No physics background required to read the output.
echo.
echo   Each test shows:
echo     - What is being checked (in plain English)
echo     - What the code computed
echo     - What the answer should be
echo     - PASS or FAIL
echo.
echo  ============================================================
echo.

:: Change to the directory where this batch file lives
cd /d "%~dp0"

:: Check Python
where python >nul 2>&1
if %errorlevel% neq 0 (
    echo  [ERROR] Python not found. Install from https://python.org
    pause
    exit /b 1
)

:: Check numpy
pip show numpy >nul 2>&1
if %errorlevel% neq 0 (
    echo  Installing numpy...
    pip install numpy
)

:: Run validator
python validate.py
set RESULT=%errorlevel%

echo.
if %RESULT% equ 0 (
    echo  ============================================================
    echo   ALL TESTS PASSED — the math checks out.
    echo  ============================================================
) else (
    echo  ============================================================
    echo   SOME TESTS FAILED — see details above.
    echo  ============================================================
)
echo.
pause
