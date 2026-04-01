@echo off
title MKGF Numerical Scripts — Full Run
color 0E

echo.
echo  ============================================================
echo   One Kink, Three Generations — Numerical Computations
echo  ============================================================
echo.
echo   Running all 10 scripts. Each produces terminal output + PNG.
echo   Estimated time: ~2 minutes (global fit takes ~60 seconds).
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

:: Check dependencies
echo  [1/12] Checking dependencies...
pip show numpy >nul 2>&1
if %errorlevel% neq 0 (
    echo         Installing requirements...
    pip install -r requirements.txt
) else (
    echo         numpy, scipy, matplotlib — OK
)
echo.

:: Create output folder
if not exist results mkdir results

echo  [2/12] Kink Profiles (Ch. 2-3)...
echo  -------------------------------------------------------
python kink_profiles.py
if exist triple_kink_profiles.png move /y triple_kink_profiles.png results\ >nul
if exist schrodinger_potential.png move /y schrodinger_potential.png results\ >nul
echo.

echo  [3/12] Overlap Integrals (Ch. 7)...
echo  -------------------------------------------------------
python overlap_integrals.py
if exist overlap_integrals.png move /y overlap_integrals.png results\ >nul
echo.

echo  [4/12] Kibble-Zurek Selection (Ch. 5)...
echo  -------------------------------------------------------
python kibble_zurek.py
if exist kibble_zurek.png move /y kibble_zurek.png results\ >nul
echo.

echo  [5/12] Global Chi-Squared Fit (Ch. 11)...
echo         (This one takes ~60 seconds — running optimizer)
echo  -------------------------------------------------------
python global_chi2_fit.py
if exist global_chi2_fit.png move /y global_chi2_fit.png results\ >nul
echo.

echo  [6/12] Parameter Scan (Ch. 6, 11)...
echo  -------------------------------------------------------
python parameter_scan.py --resolution 60
if exist parameter_scan.png move /y parameter_scan.png results\ >nul
echo.

echo  [7/12] V_ub Diameter Bound (Ch. 12)...
echo  -------------------------------------------------------
python vub_diameter_bound.py
if exist vub_diameter_bound.png move /y vub_diameter_bound.png results\ >nul
echo.

echo  [8/12] Radiative Stability (Ch. 10)...
echo  -------------------------------------------------------
python radiative_stability.py
if exist radiative_stability.png move /y radiative_stability.png results\ >nul
echo.

echo  [9/12] Seesaw Fit (Ch. 9)...
echo  -------------------------------------------------------
python seesaw_fit.py
if exist seesaw_fit.png move /y seesaw_fit.png results\ >nul
echo.

echo  [10/12] Baryogenesis (Ch. 19)...
echo  -------------------------------------------------------
python baryogenesis.py
if exist baryogenesis.png move /y baryogenesis.png results\ >nul
echo.

echo  [11/12] Axion Relic Density (Ch. 20)...
echo  -------------------------------------------------------
python axion_relic.py
if exist axion_relic.png move /y axion_relic.png results\ >nul
echo.

echo  [12/12] Assembling results...
echo  -------------------------------------------------------

:: Count figures
set COUNT=0
for %%f in (results\*.png) do set /a COUNT+=1

echo.
echo  ============================================================
echo   COMPLETE
echo  ============================================================
echo.
echo   Figures saved to results\ folder:
echo.
dir /b results\*.png 2>nul
echo.
echo   Total: %COUNT% figures generated.
echo.
echo   To validate the calculations, run:
echo     validate.bat
echo  ============================================================
echo.
pause
