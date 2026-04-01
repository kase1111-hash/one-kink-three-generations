# Numerical Scripts

Standalone Python implementations of the key calculations from *One Kink, Three Generations*.

## Planned Scripts

| Script | Chapter | Description |
|--------|---------|-------------|
| `overlap_integrals.py` | 7 | Yukawa matrix from zero-mode profile overlaps with localized Higgs |
| `global_chi2_fit.py` | 11 | Full χ² fit of 19 observables to moduli torus parameters |
| `parameter_scan.py` | 6, 11 | Scan over (θ₁, θ₂, α) moduli space with χ² landscape |
| `vub_diameter_bound.py` | 12 | Diameter bound calculation on S¹/Z₂ and |V_ub| prediction |
| `kibble_zurek.py` | 5 | KZ correlation length, defect density, and N=3 selection probability |
| `kink_profiles.py` | 2, 3 | Triple-kink BPS solution and zero-mode wave functions |
| `radiative_stability.py` | 10 | One-loop Coleman–Weinberg correction to sub-kink positions |
| `seesaw_fit.py` | 9 | PMNS matrix and neutrino mass predictions from geometric seesaw |
| `baryogenesis.py` | 19 | η_B estimate from kink-wall CP violation and sphaleron rate |
| `axion_relic.py` | 20 | Kink-axion relic density via misalignment mechanism |

## Requirements

```
numpy
scipy
matplotlib
```

## Status

Scripts are being extracted and cleaned from the interactive research sessions that produced the book. Each script will be self-contained, documented, and reproduce the corresponding chapter's key numerical results.

**Full version coming soon.**
