# One Kink, Three Generations

**The Standard Model from Extra-Dimensional Topology**

*Kase Branham — Independent Researcher*

ISBN 979-8-254-58282-3 | [Academia.edu](https://independent.academia.edu/KaseBranham) | [Amazon](https://www.amazon.com/dp/B0F9Q4MR53)

---

## The Idea

A single real scalar field in one compact extra dimension forms a topological defect — a triple-kink — that binds exactly three chiral fermion zero modes via the Jackiw–Rebbi mechanism. The Callias index theorem guarantees the count. Kibble–Zurek dynamics during the cosmological phase transition selects winding number N=3 over all alternatives.

From the positions, widths, and phases of three sub-kinks on the moduli torus T², the framework derives:

- The fermion mass hierarchy (five orders of magnitude from ν to t)
- The CKM and PMNS mixing matrices
- CP violation from a single geometric phase
- Neutrino masses via a geometric seesaw
- Strong CP suppression (geometric Nelson–Barr)
- Proton stability (Z₂ baryon parity)
- Electroweak baryogenesis from kink-wall dynamics
- A dark matter candidate (the kink-axion, mₐ ~ 6–60 μeV)

**18 geometric parameters. 22 observables. 19 fitted. 3 predicted.**

The single structural tension is |V_ub|, explained by a finite-diameter bound on the compact space.

## The Book

The book assembles the Multi-Kink Generation Framework (MKGF) — a research program spanning 25 papers — into a single unified treatment. Every chapter opens with an accessible "Before the Math" section requiring no equations beyond algebra, followed by the full technical derivation.

**Chapter List**

| # | Title | Topic |
|---|-------|-------|
| 1 | The Generation Problem | Why three? The question and its history |
| 2 | The Mathematical Toolkit | Kinks, index theorems, overlaps, moduli |
| 3 | Three Modes from One Defect | Triple-kink construction and zero modes |
| 4 | Gauge Fields and Anomalies | Gauge universality and anomaly cancellation |
| 5 | Why Three: Cosmic Selection | Kibble–Zurek selection of N=3 |
| 6 | The Moduli Torus | Parameter space T² and the flavor map |
| 7 | Flavor from Geometry | Overlap integrals, FCNCs, mass matrices |
| 8 | CP Violation Without New Phases | CKM phase from kink geometry |
| 9 | Leptons and the Seesaw | PMNS matrix and neutrino masses |
| 10 | Stability: Radiative and Topological | Loop corrections and metastability |
| 11 | The Global Fit: 21 of 22 | χ² fit to 19/20 observables + 3 predictions |
| 12 | The Diameter Bound | |V_ub| tension and proton decay |
| 13 | Strong CP: Geometric Nelson–Barr | θ̄ suppression from kink symmetry |
| 14 | The Two-Path Obstruction | Circle topology and phase control |
| 15 | The BPS Triple-Kink and Kink-Axion | BPS construction and axion dark matter |
| 16 | Gauge Embedding: The SU(5) Test | Adjoint embedding and Clebsch–Gordan fix |
| 17 | The Relaxion as a Kink | Hierarchy problem from kink rolling |
| 18 | The G' Sector and Back-Reaction | Dark confinement and barrier generation |
| 19 | Baryogenesis from Kink Walls | Baryon asymmetry from wall CP violation |
| 20 | Dark Matter: What Works, What Doesn't | Kink-axion, KK parity, G' baryons |
| 21 | The Cosmological Timeline | Planck scale to today in one narrative |
| 22 | Open Problems and the Road Ahead | What's solved, what's not, where to look |

## Repository Contents

```
one-kink-three-generations/
├── README.md
├── figures/               # Python scripts generating all book figures
│   ├── triple_kink_profile.py
│   ├── zero_mode_localization.py
│   ├── moduli_torus_scan.py
│   ├── overlap_hierarchy.py
│   ├── ckm_phase_geometry.py
│   ├── kibble_zurek_selection.py
│   └── ...
├── numerical/             # Global fit and parameter scans
│   ├── global_chi2_fit.py
│   ├── parameter_scan.py
│   ├── vub_diameter_bound.py
│   └── requirements.txt
├── paper-map.md           # Which paper → which chapter
├── errata.md              # Corrections (maintained post-publication)
└── LICENSE
```

## Paper Map

The book is based on the following research papers, all available on [Academia.edu](https://independent.academia.edu/KaseBranham):

| Paper | Topic | Chapter |
|-------|-------|---------|
| MKGF I | Chiral fermion localization, Callias index | Ch. 3 |
| MKGF II | Gauge universality, anomaly cancellation | Ch. 4 |
| MKGF III | Kibble–Zurek selection of N=3 | Ch. 5 |
| MKGF IV–V | Experimental bridge, validation tools | Ch. 11 |
| MKGF VI | Moduli torus T², two-path interference | Ch. 6 |
| MKGF VII | FCNC suppression from geometric separation | Ch. 7 |
| MKGF VIII–IX | CP violation, CKM structure | Ch. 8 |
| MKGF X | PMNS matrix and seesaw mechanism | Ch. 9 |
| MKGF XI–XII | Topological stability, Higgs localization | Ch. 10 |
| MKGF XIII | Global χ² fit (19/20 observables) | Ch. 11 |
| MKGF XIV | Signatures, V_ub, proton decay | Ch. 11–12 |
| MKGF XV | Geometric Nelson–Barr, strong CP | Ch. 13 |
| MKGF XVI | Two-path obstruction | Ch. 14 |
| GP I–II | SU(5) adjoint embedding, Clebsch–Gordan | Ch. 16 |
| Rel I–II | Relaxion as kink, G' sector | Ch. 17–18 |
| BPS | BPS triple-kink and kink-axion | Ch. 15 |
| Baryo | Baryogenesis from kink walls | Ch. 19 |
| KK | KK parity and dark matter | Ch. 20 |
| Timeline | Cosmological timeline | Ch. 21 |

## Methodology

This book was written using **natural language programming** — directing large language models (Anthropic Claude) through iterative prose to perform symbolic computation, numerical optimization, figure generation, and manuscript preparation. All substantive physical decisions — what questions to ask, which results to trust, when to change direction — were made by the author. The methodology is described in the Preface.

## Key Results

| Observable | MKGF | Experiment | Status |
|-----------|------|------------|--------|
| m_u/m_t | Geometric overlap | ~10⁻⁵ | ✓ Fitted |
| m_c/m_t | Geometric overlap | ~10⁻² | ✓ Fitted |
| \|V_us\| | Overlap ratio | 0.2245 ± 0.0008 | ✓ Fitted |
| \|V_cb\| | Overlap ratio | 0.0405 ± 0.0015 | ✓ Fitted |
| \|V_ub\| | Diameter bound | 0.00369 ± 0.00011 | 3.4σ tension |
| δ_CKM | Single geometric phase | 1.196 ± 0.045 | ✓ Fitted |
| Σm_ν | Prediction | ~59 meV | Testable (DESI, CMB-S4) |
| Ordering | Prediction | Normal | Testable (JUNO) |
| \|m_ee\| | Prediction | ~3.8 meV | Testable (nEXO) |

## Citation

```
@book{branham2026onekink,
  title     = {One Kink, Three Generations: The Standard Model 
               from Extra-Dimensional Topology},
  author    = {Branham, Kase},
  year      = {2026},
  isbn      = {979-8-254-58282-3},
  publisher = {Independently published}
}
```

## License

Code in this repository is released under the MIT License. The book text and figures are © 2026 Kase Branham, all rights reserved.

## Contact

- **Academia.edu:** [independent.academia.edu/KaseBranham](https://independent.academia.edu/KaseBranham)
- **GitHub:** [github.com/kase1111-hash](https://github.com/kase1111-hash)
- **Email:** Kase1111@gmail.com
