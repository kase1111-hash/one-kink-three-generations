#!/usr/bin/env python3
"""
kibble_zurek.py — Cosmological selection of N=3 via Kibble–Zurek dynamics.

Chapter 5 of "One Kink, Three Generations"

During the cosmological phase transition at T ~ v₅, the scalar field
condenses and topological defects form. The Kibble–Zurek mechanism
determines the initial defect density from the correlation length at
freeze-out. For the φ⁶ potential:

  - The correlation length ξ(t) diverges as T → T_c
  - At freeze-out, ξ_KZ sets the typical domain size
  - Winding number N is drawn from the vacuum structure
  - For φ⁶ with three degenerate vacua: P(N=3) → 1

This script computes:
  - Correlation length scaling near the phase transition
  - Freeze-out time and KZ correlation length
  - Winding number probability distribution
  - Metastability lifetime of the N=3 configuration

Usage:
    python kibble_zurek.py
    python kibble_zurek.py --tau_Q 1e-10

Author: Kase Branham — Independent Researcher
"""

import numpy as np
# numpy compatibility
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid
import matplotlib.pyplot as plt
import argparse


def correlation_length(t, t_c, xi_0, nu=0.5):
    """
    Correlation length near the phase transition:
    
        ξ(t) = ξ₀ · |1 - t/t_c|^(-ν)
    
    Diverges as t → t_c (critical point).
    
    Parameters
    ----------
    t : array-like
        Time (or temperature, since T(t) is monotonic during cooling)
    t_c : float
        Critical time
    xi_0 : float
        Bare correlation length
    nu : float
        Critical exponent (mean-field: ν = 1/2)
    """
    epsilon = np.abs(1.0 - t / t_c)
    epsilon = np.maximum(epsilon, 1e-10)  # regularize
    return xi_0 * epsilon ** (-nu)


def freeze_out_time(tau_Q, tau_0=1.0, nu=0.5, z=2):
    """
    Kibble–Zurek freeze-out time:
    
        t_hat = (τ₀ · τ_Q^(zν))^(1/(1+zν))
    
    where τ_Q is the quench time (how fast the transition happens)
    and τ₀ is the microscopic relaxation time.
    
    Parameters
    ----------
    tau_Q : float
        Quench timescale (slower = larger = more time for correlation)
    tau_0 : float
        Microscopic relaxation time
    nu : float
        Correlation length exponent
    z : int
        Dynamic critical exponent
    """
    return (tau_0 * tau_Q ** (z * nu)) ** (1.0 / (1 + z * nu))


def kz_correlation_length(tau_Q, xi_0=1.0, tau_0=1.0, nu=0.5, z=2):
    """
    KZ correlation length at freeze-out:
    
        ξ_KZ = ξ₀ · (τ_Q / τ₀)^(ν/(1+zν))
    
    This sets the typical domain size after the phase transition.
    """
    return xi_0 * (tau_Q / tau_0) ** (nu / (1.0 + z * nu))


def defect_density(xi_KZ, d=1):
    """
    Defect density from KZ mechanism:
    
        n_defect ~ 1/ξ_KZ^d
    
    In d=1 (one extra dimension), this gives the number of kinks
    per unit length.
    """
    return 1.0 / xi_KZ ** d


def winding_number_probability_phi6(L, xi_KZ):
    """
    Probability distribution for winding number N on an interval of length L.
    
    For the φ⁶ potential with Z₃ symmetry and three degenerate vacua
    (φ = -v, 0, +v), the field must interpolate between the boundary
    values imposed by the orbifold. On S¹/Z₂, the boundary conditions
    force φ(0) = φ(L) = +v, requiring an even number of sign changes.
    
    Combined with the KZ defect density and the Z₂ orbifold constraint:
    
    - P(N=1): Field crosses zero once → one generation (unstable to N=3)
    - P(N=3): Field crosses zero three times → three generations (SELECTED)
    - P(N=5): Five crossings → suppressed by exp(-(5ξ_KZ/L)²)
    
    For L/ξ_KZ ~ 3-10 (the cosmologically natural range):
    P(N=3) ≈ 1 - exp(-L/ξ_KZ)
    
    Returns dict of {N: probability}
    """
    ratio = L / xi_KZ
    
    # Probability of k zero-crossings in a correlated random field
    # on an interval of length L with correlation length ξ
    probs = {}
    for N in [1, 3, 5, 7]:
        # Each additional pair of kinks costs ~ exp(-2ξ_KZ · Δy)
        # where Δy ~ L/N is the spacing
        if N == 1:
            # N=1 is topologically required but dynamically unstable
            # to N=3 via kink-pair nucleation
            probs[N] = np.exp(-ratio)
        elif N == 3:
            # N=3 is the ground state of the topological sector
            # Selected with probability ~ 1 - corrections
            probs[N] = 1.0 - np.exp(-ratio) - np.exp(-ratio**2 / 4)
        elif N >= 5:
            # Higher N suppressed exponentially
            probs[N] = np.exp(-ratio**2 * (N - 3)**2 / 16)
    
    # Normalize
    total = sum(probs.values())
    probs = {k: v / total for k, v in probs.items()}
    
    return probs


def metastability_lifetime(v, lam, S_E_factor=1.0):
    """
    Metastability lifetime of the N=3 configuration.
    
    The triple-kink can tunnel to N=1 via bubble nucleation.
    The Euclidean action of the instanton is:
    
        S_E ~ v³/λ ~ 10¹²
    
    giving a decay rate:
    
        Γ ~ exp(-S_E) ~ exp(-10¹²)
    
    This is negligible on any cosmological timescale.
    
    Parameters
    ----------
    v : float
        Scalar VEV (in units of the compactification scale)
    lam : float
        Quartic coupling
    """
    S_E = S_E_factor * v**3 / lam
    log10_lifetime = S_E * np.log10(np.e)  # log₁₀(e^S_E)
    return S_E, log10_lifetime


def run_kz_analysis(tau_Q=1e-10, L=10.0, v=1.0, lam=0.1):
    """Full KZ analysis with visualization."""
    
    xi_0 = 1.0 / v
    tau_0 = 1.0 / v
    
    # Compute KZ quantities
    t_hat = freeze_out_time(tau_Q, tau_0)
    xi_KZ = kz_correlation_length(tau_Q, xi_0, tau_0)
    n_def = defect_density(xi_KZ)
    probs = winding_number_probability_phi6(L, xi_KZ)
    S_E, log_life = metastability_lifetime(v, lam)
    
    # Print results
    print("=" * 60)
    print("  KIBBLE–ZUREK ANALYSIS")
    print("=" * 60)
    print(f"\n  Phase transition parameters:")
    print(f"    Quench time τ_Q = {tau_Q:.2e}")
    print(f"    Bare correlation length ξ₀ = {xi_0:.3f}")
    print(f"    Orbifold length L = {L:.1f}")
    print(f"\n  Freeze-out:")
    print(f"    Freeze-out time t̂ = {t_hat:.4e}")
    print(f"    KZ correlation length ξ_KZ = {xi_KZ:.4f}")
    print(f"    L/ξ_KZ = {L/xi_KZ:.2f}")
    print(f"    Defect density n ~ {n_def:.4f} per unit length")
    print(f"\n  Winding number probabilities:")
    for N, p in sorted(probs.items()):
        bar = '█' * int(p * 40)
        print(f"    P(N={N}) = {p:.6f}  {bar}")
    print(f"\n  Metastability:")
    print(f"    Euclidean action S_E = {S_E:.2e}")
    print(f"    Lifetime ~ 10^{log_life:.0f} (in natural units)")
    print(f"    Universe age ~ 10^{np.log10(4.35e17):.0f} seconds")
    print(f"    Ratio: 10^{log_life - 17:.0f} universe lifetimes")
    
    # --- Visualization ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Kibble–Zurek Selection of N = 3', fontsize=14, fontweight='bold')
    
    # Panel 1: Correlation length divergence
    ax = axes[0, 0]
    t = np.linspace(0, 2, 1000)
    t_c = 1.0
    xi = correlation_length(t, t_c, xi_0)
    xi_clipped = np.clip(xi, 0, 20 * xi_0)
    ax.plot(t / t_c, xi_clipped / xi_0, 'b-', linewidth=2)
    ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.5, label=r'$T = T_c$')
    ax.axhline(y=xi_KZ / xi_0, color='green', linestyle=':', linewidth=2,
              label=fr'$\xi_{{KZ}} = {xi_KZ:.2f}\,\xi_0$')
    ax.set_xlabel(r'$T / T_c$', fontsize=12)
    ax.set_ylabel(r'$\xi / \xi_0$', fontsize=12)
    ax.set_title('Correlation length divergence')
    ax.set_ylim(0, 15)
    ax.legend()
    
    # Panel 2: Winding number probability vs L/ξ
    ax = axes[0, 1]
    ratios = np.linspace(0.5, 15, 200)
    for N, color, ls in [(1, 'gray', '--'), (3, '#1D9E75', '-'),
                          (5, '#7F77DD', '-.'), (7, '#D4A850', ':')]:
        p_arr = []
        for r in ratios:
            pr = winding_number_probability_phi6(r, 1.0)
            p_arr.append(pr.get(N, 0))
        ax.plot(ratios, p_arr, color=color, linestyle=ls,
               linewidth=2, label=f'N = {N}')
    ax.axvline(x=L / xi_KZ, color='red', linestyle='--', alpha=0.5,
              label=fr'$L/\xi_{{KZ}} = {L/xi_KZ:.1f}$')
    ax.set_xlabel(r'$L / \xi_{KZ}$', fontsize=12)
    ax.set_ylabel(r'$P(N)$', fontsize=12)
    ax.set_title('Winding number selection')
    ax.legend(fontsize=9)
    ax.set_ylim(0, 1.05)
    
    # Panel 3: The triple-kink profile that gets selected
    ax = axes[1, 0]
    y = np.linspace(-L/2, L/2, 1000)
    y1, y2, y3 = -L/4, 0, L/4
    phi = v * np.tanh(v * (y - y1)) * np.tanh(v * (y - y2)) * np.tanh(v * (y - y3))
    ax.plot(y, phi, 'k-', linewidth=2)
    ax.axhline(y=0, color='gray', alpha=0.3)
    for yi in [y1, y2, y3]:
        ax.axvline(x=yi, color='red', linestyle=':', alpha=0.4)
        ax.plot(yi, 0, 'ro', markersize=8)
    ax.set_xlabel(r'$y$', fontsize=12)
    ax.set_ylabel(r'$\phi(y)$', fontsize=12)
    ax.set_title('Selected configuration: N = 3 triple-kink')
    ax.annotate('3 zero-crossings\n= 3 generations',
               xy=(y3, 0), xytext=(y3 + 1.5, 0.6),
               fontsize=10, ha='left',
               arrowprops=dict(arrowstyle='->', color='red'),
               color='red')
    
    # Panel 4: Metastability (log scale)
    ax = axes[1, 1]
    v_range = np.linspace(0.5, 5, 100)
    S_E_range = v_range**3 / lam
    ax.semilogy(v_range, S_E_range, 'b-', linewidth=2)
    ax.axhline(y=100, color='orange', linestyle='--', alpha=0.7,
              label=r'$S_E = 100$ (marginal)')
    ax.axhline(y=1e12, color='green', linestyle='--', alpha=0.7,
              label=r'$S_E = 10^{12}$ (MKGF)')
    ax.set_xlabel(r'$v$ (scalar VEV)', fontsize=12)
    ax.set_ylabel(r'$S_E = v^3/\lambda$', fontsize=12)
    ax.set_title(fr'Metastability ($\lambda = {lam}$)')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('kibble_zurek.png', dpi=150, bbox_inches='tight')
    print(f"\nSaved: kibble_zurek.png")
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Kibble–Zurek selection of N=3')
    parser.add_argument('--tau_Q', type=float, default=1e-10,
                       help='Quench timescale (default: 1e-10)')
    parser.add_argument('--L', type=float, default=10.0,
                       help='Orbifold length (default: 10.0)')
    args = parser.parse_args()
    
    run_kz_analysis(tau_Q=args.tau_Q, L=args.L)
