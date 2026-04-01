# Numerical Scripts

Standalone Python implementations of the key calculations from *One Kink, Three Generations*.

## Scripts

| Script | Chapter | Description |
|--------|---------|-------------|
| `kink_profiles.py` | 2, 3 | Triple-kink BPS solution, zero-mode wave functions, Schrödinger potential |
| `overlap_integrals.py` | 7 | Yukawa matrix from zero-mode overlaps, CKM matrix, mass hierarchies |
| `kibble_zurek.py` | 5 | KZ correlation length, freeze-out, N=3 selection, metastability |
| `global_chi2_fit.py` | 11 | Full χ² fit of 19 observables to 18 moduli torus parameters |
| `parameter_scan.py` | 6, 11 | 2D scan over (θ₁, θ₂) moduli space with χ² landscape |
| `vub_diameter_bound.py` | 12 | Diameter bound on \|V_ub\|, 3.4σ tension, SU(5) Clebsch–Gordan correction |
| `radiative_stability.py` | 10 | Coleman–Weinberg corrections, sub-kink position shifts, mass hierarchy robustness |
| `seesaw_fit.py` | 9 | PMNS matrix, neutrino masses, Σm_ν and \|m_ee\| predictions |
| `baryogenesis.py` | 19 | Kink-wall CP source, sphaleron rate, η_B/η_obs ≈ 0.4–5 |
| `axion_relic.py` | 20 | Kink-axion mass, misalignment relic density, detection prospects |

## Quick Start

```bash
pip install -r requirements.txt

# Run everything and collect figures
run_all.bat

# Validate the math (23 independent checks)
validate.bat
```

Or run individual scripts:

```bash
python kink_profiles.py
python overlap_integrals.py --alpha 4.3 --theta1 0.35
python parameter_scan.py --resolution 100
python axion_relic.py
```

Each script produces terminal output with key results and saves a PNG figure.

## Validation

`validate.py` runs 23 independent checks across 7 test groups, comparing computed values to known analytical results:

| Group | Tests | Checks |
|-------|-------|--------|
| Kink profile | 4 | Boundary values, zero-crossing, kink width |
| Zero modes | 3 | Normalization, width, peak value |
| Overlaps | 5 | Gaussian overlaps at 4 separations, mass hierarchy |
| Kibble–Zurek | 3 | KZ exponent, correlation length, Euclidean action |
| CKM matrix | 4 | Unitarity, Cabibbo angle, Jarlskog invariant |
| Diameter bound | 2 | \|V_ub\| ceiling, required separation |
| Seesaw | 2 | Neutrino mass scale, atmospheric splitting |

No physics background required to read the output.

## Requirements

```
numpy >= 1.20
scipy >= 1.7
matplotlib >= 3.5
```
