#!/usr/bin/env python3
"""
kink_profiles.py — Triple-kink BPS solution and zero-mode wave functions.

Chapters 2–3 of "One Kink, Three Generations"

Computes:
  - The BPS triple-kink scalar profile φ(y) for the φ⁶ potential
  - The three chiral zero-mode wave functions f_a(y) localized at each sub-kink
  - Visualization of both

The key physics: a scalar field with superpotential W(φ) = v²φ - φ³/3
produces a potential V(φ) = ½(W')² = ½(v² - φ²)² with vacua at ±v.
The triple-kink (winding number N=3) has three zero-crossings, each
trapping one chiral fermion via the Jackiw–Rebbi mechanism.

Usage:
    python kink_profiles.py
    python kink_profiles.py --alpha 4.3 --L 10.0

Author: Kase Branham — Independent Researcher
"""

import numpy as np
# numpy compatibility
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import argparse


def phi6_superpotential_derivative(phi, v=1.0):
    """W'(φ) = v² - φ² for the φ⁶ potential."""
    return v**2 - phi**2


def single_kink(y, y0=0.0, v=1.0):
    """
    Single BPS kink solution: φ(y) = v·tanh(v·(y - y0)).
    
    This is the exact solution of φ' = W'(φ) = v² - φ²
    interpolating from -v to +v.
    """
    return v * np.tanh(v * (y - y0))


def product_ansatz_triple_kink(y, y1, y2, y3, v=1.0):
    """
    Product ansatz for the triple-kink profile (Eq. 2.15 in the book):
    
        φ(y) ≈ v · tanh(v(y-y1)) · tanh(v(y-y2)) · tanh(v(y-y3)) / v²
    
    This is a heuristic approximation that captures the three zero-crossings
    at y1, y2, y3. The exact BPS triple-kink requires numerical integration
    of the full φ⁶ Bogomolny equation.
    
    Parameters
    ----------
    y : array-like
        Extra-dimensional coordinate
    y1, y2, y3 : float
        Positions of the three sub-kinks (zero-crossings)
    v : float
        Scalar VEV
    """
    k1 = np.tanh(v * (y - y1))
    k2 = np.tanh(v * (y - y2))
    k3 = np.tanh(v * (y - y3))
    return v * k1 * k2 * k3


def zero_mode_profile(y, ya, G, phi_func, v=1.0):
    """
    Jackiw–Rebbi zero-mode wave function localized at sub-kink position ya.
    
    The zero mode of the Dirac operator D = γ⁵(∂_y - Gφ(y)) is:
    
        f_a(y) ∝ exp(-G ∫₀ʸ φ(y') dy')
    
    Near sub-kink ya, this is approximately Gaussian:
    
        f_a(y) ≈ (Gv/π)^{1/4} · exp(-½ Gv (y - ya)²)
    
    Parameters
    ----------
    y : array-like
        Extra-dimensional coordinate
    ya : float
        Sub-kink position (zero-crossing of φ)
    G : float
        5D Yukawa coupling
    phi_func : callable
        The scalar field profile φ(y)
    v : float
        Scalar VEV
    """
    # Gaussian approximation (valid for well-separated sub-kinks)
    alpha = G * v  # compactness parameter
    width = 1.0 / np.sqrt(alpha)
    f = (alpha / np.pi) ** 0.25 * np.exp(-0.5 * alpha * (y - ya) ** 2)
    return f


def compute_kink_width(v=1.0):
    """
    Kink width δ = 1/v.
    
    The characteristic scale over which the scalar field transitions
    between vacua. Sub-kinks must be separated by several δ for the
    product ansatz and Gaussian zero-mode approximation to hold.
    """
    return 1.0 / v


def schrodinger_potential(y, G, phi_func, v=1.0):
    """
    The Schrödinger-like potential from squaring the Dirac operator:
    
        V_S(y) = G²φ(y)² - Gφ'(y)
    
    The zero modes correspond to the ground state of -∂²_y + V_S(y).
    The three quasi-degenerate bound states near each sub-kink
    become the three fermion generations.
    """
    dy = 1e-6
    phi = phi_func(y)
    dphi = (phi_func(y + dy) - phi_func(y - dy)) / (2 * dy)
    return G**2 * phi**2 - G * dphi


def plot_profiles(y1=-3.0, y2=0.0, y3=3.0, v=1.0, G=4.3, L=12.0):
    """Generate the triple-kink profile and zero-mode visualization."""
    
    alpha = G * v
    y = np.linspace(-L/2, L/2, 2000)
    
    # Scalar field profile
    phi = product_ansatz_triple_kink(y, y1, y2, y3, v)
    
    # Zero-mode profiles
    f1 = zero_mode_profile(y, y1, G, None, v)
    f2 = zero_mode_profile(y, y2, G, None, v)
    f3 = zero_mode_profile(y, y3, G, None, v)
    
    # Normalize
    for f in [f1, f2, f3]:
        norm = np.trapz(f**2, y)
        f /= np.sqrt(norm) if norm > 0 else 1.0
    
    # Schrödinger potential
    phi_func = lambda yy: product_ansatz_triple_kink(yy, y1, y2, y3, v)
    V_S = schrodinger_potential(y, G, phi_func, v)
    
    # --- Figure 1: Scalar profile + zero modes ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle('Triple-Kink Profile and Zero Modes', fontsize=14, fontweight='bold')
    
    # Top: scalar field
    ax1.plot(y, phi, 'k-', linewidth=2, label=r'$\phi(y)$ (triple-kink)')
    ax1.axhline(y=v, color='gray', linestyle='--', alpha=0.5, label=r'$\pm v$')
    ax1.axhline(y=-v, color='gray', linestyle='--', alpha=0.5)
    ax1.axhline(y=0, color='gray', linestyle='-', alpha=0.2)
    for yi, label in [(y1, r'$y_1$'), (y2, r'$y_2$'), (y3, r'$y_3$')]:
        ax1.axvline(x=yi, color='red', linestyle=':', alpha=0.4)
    ax1.set_ylabel(r'$\phi(y) / v$', fontsize=12)
    ax1.set_title(f'Scalar field profile (product ansatz, $\\alpha = Gv = {alpha:.1f}$)')
    ax1.legend(loc='upper left')
    ax1.set_ylim(-1.5*v, 1.5*v)
    
    # Bottom: zero modes
    colors = ['#D4A850', '#1D9E75', '#7F77DD']  # gold, teal, violet
    labels = ['Generation 1', 'Generation 2', 'Generation 3']
    for f, c, lbl in zip([f1, f2, f3], colors, labels):
        ax2.fill_between(y, 0, f**2, alpha=0.3, color=c)
        ax2.plot(y, f**2, color=c, linewidth=2, label=lbl)
    
    ax2.set_xlabel(r'Extra dimension $y$', fontsize=12)
    ax2.set_ylabel(r'$|f_a(y)|^2$', fontsize=12)
    ax2.set_title('Zero-mode probability densities (Jackiw–Rebbi)')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('triple_kink_profiles.png', dpi=150, bbox_inches='tight')
    print("Saved: triple_kink_profiles.png")
    
    # --- Figure 2: Schrödinger potential ---
    fig2, ax3 = plt.subplots(figsize=(10, 5))
    V_clipped = np.clip(V_S, -10*G, 10*G)
    ax3.plot(y, V_clipped, 'k-', linewidth=1.5)
    ax3.set_xlabel(r'Extra dimension $y$', fontsize=12)
    ax3.set_ylabel(r"$V_S(y) = G^2\phi^2 - G\,d\phi/dy$", fontsize=12)
    ax3.set_title(f'Schrödinger potential ($G = {G:.1f}$, $v = {v:.1f}$)')
    ax3.set_ylim(-2*G, 5*G)
    
    # Mark the three wells
    for yi, lbl in [(y1, 'Well 1'), (y2, 'Well 2'), (y3, 'Well 3')]:
        ax3.annotate(lbl, xy=(yi, -G), fontsize=9, ha='center',
                    color='red', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('schrodinger_potential.png', dpi=150, bbox_inches='tight')
    print("Saved: schrodinger_potential.png")
    
    # --- Print key parameters ---
    print(f"\nKey parameters:")
    print(f"  Compactness parameter α = Gv = {alpha:.2f}")
    print(f"  Kink width δ = 1/v = {1/v:.3f}")
    print(f"  Zero-mode width w = 1/√(Gv) = {1/np.sqrt(alpha):.3f}")
    print(f"  Sub-kink separations: Δy₁₂ = {y2-y1:.2f}, Δy₂₃ = {y3-y2:.2f}")
    print(f"  Separation in kink widths: {(y2-y1)*v:.1f}δ, {(y3-y2)*v:.1f}δ")
    
    # Overlap integrals (preview)
    O12 = np.trapz(f1 * f2, y)
    O23 = np.trapz(f2 * f3, y)
    O13 = np.trapz(f1 * f3, y)
    print(f"\nOverlap integrals (zero-mode × zero-mode):")
    print(f"  ⟨f₁|f₂⟩ = {O12:.6e}")
    print(f"  ⟨f₂|f₃⟩ = {O23:.6e}")
    print(f"  ⟨f₁|f₃⟩ = {O13:.6e}")
    print(f"  Hierarchy: O₁₂/O₂₃ = {O12/O23:.3f}" if O23 != 0 else "")
    
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Triple-kink BPS solution and zero-mode profiles')
    parser.add_argument('--alpha', type=float, default=4.3,
                       help='Compactness parameter α = Gv (default: 4.3)')
    parser.add_argument('--L', type=float, default=12.0,
                       help='Orbifold length (default: 12.0)')
    parser.add_argument('--sep', type=float, default=3.0,
                       help='Sub-kink separation (default: 3.0)')
    args = parser.parse_args()
    
    v = 1.0
    G = args.alpha / v
    sep = args.sep
    
    plot_profiles(y1=-sep, y2=0.0, y3=sep, v=v, G=G, L=args.L)
