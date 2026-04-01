#!/usr/bin/env python3
"""
axion_relic.py — Kink-axion dark matter via the misalignment mechanism.

Chapter 20 of "One Kink, Three Generations"

The kink modulus has a pseudo-Goldstone boson — the kink-axion —
with mass m_a ~ 6–60 μeV determined by the PQ breaking scale f_a.
Its relic density is set by the misalignment mechanism:

    Ω_a h² ≈ 0.12 × (f_a / 10¹² GeV)^{7/6} × θ_i²

For f_a ~ 10¹² GeV and θ_i ~ O(1), the kink-axion can account
for all of dark matter.

Other dark matter candidates in the MKGF:
  - KK parity: BROKEN by the kink background → no stable KK WIMP
  - G' dark baryon: possible supplement at ~10 GeV

This script computes:
  - Kink-axion mass from f_a
  - Relic density from misalignment
  - Detection prospects (ADMX, ABRACADABRA)

Usage:
    python axion_relic.py

Author: Kase Branham — Independent Researcher
"""

import numpy as np
# numpy compatibility
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid
import matplotlib.pyplot as plt


# Constants
OMEGA_DM = 0.12      # observed dark matter density (Planck)
H_0 = 67.4           # Hubble constant (km/s/Mpc)
T_QCD = 0.15         # QCD phase transition temperature (GeV)
M_PLANCK = 1.22e19   # Planck mass (GeV)


def axion_mass(f_a):
    """
    Axion mass from the PQ breaking scale:
    
        m_a ≈ 6 μeV × (10¹² GeV / f_a)
    
    For the kink-axion: f_a ~ 10¹¹ – 10¹³ GeV
    → m_a ~ 0.6 – 60 μeV
    """
    return 6e-6 * (1e12 / f_a)  # eV


def misalignment_relic_density(f_a, theta_i=1.0):
    """
    Relic density from the vacuum misalignment mechanism:
    
        Ω_a h² ≈ 0.12 × (f_a / 10¹² GeV)^{7/6} × θ_i²
    
    The axion field starts displaced from its minimum by angle θ_i.
    After m_a > H, it oscillates coherently — these oscillations
    are cold dark matter.
    
    Parameters
    ----------
    f_a : float
        PQ breaking scale (GeV)
    theta_i : float
        Initial misalignment angle (0 to π)
    """
    return 0.12 * (f_a / 1e12)**(7.0/6) * theta_i**2


def oscillation_temperature(f_a):
    """
    Temperature when axion oscillations begin (m_a = 3H):
    
        T_osc ~ 1 GeV × (f_a / 10¹² GeV)^{-1/2}
    """
    return 1.0 * (f_a / 1e12)**(-0.5)  # GeV


def kk_parity_status():
    """
    Check KK parity as a dark matter mechanism.
    
    In symmetric orbifolds, KK parity (y → -y) would stabilize
    the lightest KK-odd particle. But the kink background φ(y)
    explicitly breaks this symmetry:
    
        φ(-y) ≠ ±φ(y) for the triple-kink
    
    Result: KK parity is broken → no stable KK WIMP.
    """
    return {
        'mechanism': 'KK parity',
        'status': 'BROKEN',
        'reason': 'Triple-kink breaks y → -y symmetry',
        'consequence': 'Lightest KK particle decays — not dark matter',
    }


def gprime_dark_baryon(M_gprime=10.0, n_f=2):
    """
    G' dark baryon as a WIMP-like supplement.
    
    The dark SU(3)' with n_f=2 flavors confines at Λ' ~ ΛQCD.
    The lightest dark baryon has mass ~ 3Λ' ~ few GeV.
    
    This is viable but not the primary DM candidate.
    """
    Lambda_prime = M_gprime / 3  # confinement scale
    return {
        'mechanism': "G' dark baryon",
        'mass': f'~{M_gprime:.0f} GeV',
        'status': 'POSSIBLE (supplement)',
        'sigma_SI': '~10⁻⁴⁷ cm² (below current limits)',
    }


def admx_sensitivity(f_a_range=(1e11, 1e13)):
    """
    ADMX haloscope sensitivity range.
    
    ADMX is sensitive to axions with:
        m_a ~ 2–40 μeV (f_a ~ 1.5×10¹¹ – 3×10¹² GeV)
    
    and coupling:
        g_aγγ ~ α/(2π f_a)
    """
    m_low = axion_mass(f_a_range[1])
    m_high = axion_mass(f_a_range[0])
    return m_low, m_high


def run_axion_analysis(f_a=1e12, theta_i=1.0):
    """Full kink-axion dark matter analysis."""
    
    m_a = axion_mass(f_a)
    omega = misalignment_relic_density(f_a, theta_i)
    T_osc = oscillation_temperature(f_a)
    kk = kk_parity_status()
    gp = gprime_dark_baryon()
    
    print("=" * 60)
    print("  DARK MATTER: KINK-AXION AND ALTERNATIVES")
    print("=" * 60)
    
    print(f"\n  Kink-axion parameters:")
    print(f"    f_a = {f_a:.1e} GeV (PQ breaking scale)")
    print(f"    m_a = {m_a*1e6:.1f} μeV = {m_a:.2e} eV")
    print(f"    θ_i = {theta_i:.2f} (initial misalignment)")
    print(f"    T_osc = {T_osc:.2f} GeV (oscillation onset)")
    
    print(f"\n  Relic density:")
    print(f"    Ω_a h² = {omega:.4f}")
    print(f"    Ω_DM h² = {OMEGA_DM} (observed)")
    print(f"    Ω_a/Ω_DM = {omega/OMEGA_DM:.2f}")
    
    # Find f_a that gives exactly Ω_DM
    f_a_exact = 1e12 * (OMEGA_DM / (0.12 * theta_i**2))**(6.0/7)
    m_a_exact = axion_mass(f_a_exact)
    print(f"\n  For Ω_a = Ω_DM (θ_i = {theta_i}):")
    print(f"    f_a = {f_a_exact:.2e} GeV")
    print(f"    m_a = {m_a_exact*1e6:.1f} μeV")
    
    print(f"\n  Alternative candidates:")
    print(f"    {kk['mechanism']}: {kk['status']}")
    print(f"      {kk['reason']}")
    print(f"      {kk['consequence']}")
    print(f"    {gp['mechanism']}: {gp['status']}")
    print(f"      Mass: {gp['mass']}")
    
    print(f"\n  Detection prospects:")
    m_low, m_high = admx_sensitivity()
    print(f"    ADMX range: {m_low*1e6:.1f}–{m_high*1e6:.0f} μeV")
    in_range = m_low <= m_a <= m_high
    print(f"    Kink-axion ({m_a*1e6:.1f} μeV): "
          f"{'IN RANGE' if in_range else 'outside current range'}")
    print(f"    ABRACADABRA: sensitive to f_a ~ 10¹²–10¹⁶ GeV")
    print(f"    CASPEr: sensitive to m_a ~ 10⁻¹² – 10⁻⁶ eV")
    
    # --- Visualization ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Kink-Axion Dark Matter', fontsize=14, fontweight='bold')
    
    # Panel 1: Relic density vs f_a
    ax = axes[0, 0]
    f_range = np.logspace(10, 14, 500)
    for ti, color, ls in [(0.5, '#7F77DD', '--'), (1.0, '#1D9E75', '-'),
                           (2.0, '#D4A850', '-.')]:
        omega_arr = [misalignment_relic_density(f, ti) for f in f_range]
        ax.loglog(f_range, omega_arr, color=color, linestyle=ls,
                 linewidth=2, label=fr'$\theta_i = {ti}$')
    ax.axhline(y=OMEGA_DM, color='red', linestyle='--', linewidth=1.5,
              label=fr'$\Omega_{{DM}}h^2 = {OMEGA_DM}$')
    ax.axvline(x=f_a, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel(r'$f_a$ (GeV)', fontsize=12)
    ax.set_ylabel(r'$\Omega_a h^2$', fontsize=12)
    ax.set_title('Relic density vs PQ scale')
    ax.legend(fontsize=8)
    ax.set_ylim(1e-4, 1e4)
    
    # Panel 2: Axion mass vs f_a
    ax = axes[0, 1]
    m_range = [axion_mass(f) for f in f_range]
    ax.loglog(f_range, np.array(m_range) * 1e6, 'b-', linewidth=2)
    ax.axhspan(m_low * 1e6, m_high * 1e6, alpha=0.15, color='green',
              label='ADMX range')
    ax.plot(f_a, m_a * 1e6, 'r*', markersize=15, label=f'MKGF benchmark')
    ax.set_xlabel(r'$f_a$ (GeV)', fontsize=12)
    ax.set_ylabel(r'$m_a$ ($\mu$eV)', fontsize=12)
    ax.set_title('Axion mass')
    ax.legend()
    
    # Panel 3: DM candidate scorecard
    ax = axes[1, 0]
    ax.axis('off')
    scorecard = [
        ['Candidate', 'Mass', 'Viable?', 'Testable?'],
        ['Kink-axion', f'{m_a*1e6:.0f} μeV', '✓ YES', 'ADMX, CASPEr'],
        ['KK WIMP', '~TeV', '✗ NO', '(KK parity broken)'],
        ["G' baryon", '~10 GeV', '? Maybe', 'Direct detection'],
        ['Breather', '~GeV', '? Maybe', 'Collider'],
    ]
    table = ax.table(cellText=scorecard, loc='center', cellLoc='center',
                    colWidths=[0.22, 0.18, 0.18, 0.32])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.8)
    for j in range(4):
        table[0, j].set_facecolor('#1D9E75')
        table[0, j].set_text_props(color='white', fontweight='bold')
    table[1, 2].set_facecolor('#E1F5EE')  # kink-axion: green
    table[2, 2].set_facecolor('#FCEBEB')  # KK: red
    ax.set_title('Dark Matter Scorecard', fontsize=12, fontweight='bold', pad=20)
    
    # Panel 4: θ_i tuning
    ax = axes[1, 1]
    theta_range = np.linspace(0.01, np.pi, 200)
    omega_theta = [misalignment_relic_density(f_a, t) for t in theta_range]
    ax.plot(theta_range, np.array(omega_theta) / OMEGA_DM, 'b-', linewidth=2)
    ax.axhline(y=1, color='red', linestyle='--', label=r'$\Omega_a = \Omega_{DM}$')
    ax.fill_between(theta_range, 0.5, 2.0, alpha=0.1, color='green',
                   label='Factor-of-2 range')
    theta_exact = np.sqrt(OMEGA_DM / (0.12 * (f_a/1e12)**(7/6)))
    ax.axvline(x=theta_exact, color='orange', linestyle=':', linewidth=2,
              label=fr'$\theta_i = {theta_exact:.2f}$')
    ax.set_xlabel(r'Initial misalignment $\theta_i$', fontsize=12)
    ax.set_ylabel(r'$\Omega_a / \Omega_{DM}$', fontsize=12)
    ax.set_title(fr'Misalignment tuning ($f_a = {f_a:.0e}$ GeV)')
    ax.legend(fontsize=8)
    ax.set_ylim(0, 5)
    
    plt.tight_layout()
    plt.savefig('axion_relic.png', dpi=150, bbox_inches='tight')
    print(f"\nSaved: axion_relic.png")
    plt.show()


if __name__ == '__main__':
    run_axion_analysis()
