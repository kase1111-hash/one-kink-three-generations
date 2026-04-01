# Numerical Scripts

Standalone Python implementations of the key calculations from *One Kink, Three Generations*.

## Available

| Script | Chapter | Description |
|--------|---------|-------------|
| `kink_profiles.py` | 2, 3 | Triple-kink BPS solution, zero-mode wave functions, Schrödinger potential |
| `overlap_integrals.py` | 7 | Yukawa matrix from zero-mode overlaps, CKM matrix extraction, mass hierarchies |
| `kibble_zurek.py` | 5 | KZ correlation length, freeze-out, N=3 selection probability, metastability |
| `global_chi2_fit.py` | 11 | Full χ² fit of 19 observables to 18 moduli torus parameters |
| `parameter_scan.py` | 6, 11 | 2D scan over (θ₁, θ₂) moduli space with χ² landscape and mass/mixing contours |
| `vub_diameter_bound.py` | 12 | Diameter bound on \|V_ub\|, 3.4σ tension, SU(5) Clebsch–Gordan correction |

## Coming Soon

| Script | Chapter | Description |
|--------|---------|-------------|
| `radiative_stability.py` | 10 | One-loop Coleman–Weinberg correction to sub-kink positions |
| `seesaw_fit.py` | 9 | PMNS matrix and neutrino mass predictions from geometric seesaw |
| `baryogenesis.py` | 19 | η_B estimate from kink-wall CP violation and sphaleron rate |
| `axion_relic.py` | 20 | Kink-axion relic density via misalignment mechanism |

## Quick Start

```bash
pip install -r requirements.txt

# Core pipeline (kink → overlaps → fit)
python kink_profiles.py
python overlap_integrals.py
python global_chi2_fit.py

# Exploration
python parameter_scan.py --alpha 4.3 --resolution 100
python vub_diameter_bound.py
python kibble_zurek.py
```

Each script produces terminal output with key numerical results and saves a PNG figure.

## Requirements

```
numpy >= 1.20
scipy >= 1.7
matplotlib >= 3.5
```
