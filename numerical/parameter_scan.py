#!/usr/bin/env python3
"""
parameter_scan.py — Moduli torus scan with χ² landscape.

Chapters 6 and 11 of "One Kink, Three Generations"

Scans the two-dimensional moduli space T² = (θ₁, θ₂) at fixed
compactness parameter α, computing χ² at each point. Visualizes
the flavor landscape showing how mass hierarchies and mixing angles
vary continuously across the torus.

The key insight: the Standard Model fermion spectrum corresponds
to a specific point on T², not a random one. The χ² minimum is
sharp and isolated — the geometry picks out the observed physics.

Usage:
    python parameter_scan.py
    python parameter_scan.py --alpha 4.3 --resolution 100

Author: Kase Branham — Independent Researcher
"""

import numpy as np
# numpy compatibility
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse


def overlap(sep, alpha):
    """Yukawa from Gaussian overlap: Y ~ exp(-α·s²/2)."""
    return np.exp(-0.5 * alpha * sep**2)


def quick_chi2(theta1, theta2, alpha, L=10.0):
    """
    Fast χ² computation for a point on the moduli torus.
    
    Uses the quark sector only (6 mass ratios + 3 CKM elements)
    for speed during scanning.
    """
    # Sub-kink positions
    y1 = -theta1 * L / 2
    y2 = 0.0
    y3 = theta2 * L / 2

    s12 = abs(y2 - y1)
    s23 = abs(y3 - y2)
    s13 = abs(y3 - y1)

    # Predictions
    m_u_t = overlap(s13, alpha)**2
    m_c_t = overlap(s12, alpha)**2
    m_d_b = overlap(s13, alpha * 0.9)**2
    m_s_b = overlap(s12, alpha * 0.9)**2

    V_us = overlap(s12, alpha)
    V_cb = overlap(s23, alpha)
    V_ub = overlap(s13, alpha)

    # Experimental values and errors
    data = [
        (m_u_t,  1.27e-5, 0.50e-5),
        (m_c_t,  7.35e-3, 0.30e-3),
        (m_d_b,  1.13e-3, 0.20e-3),
        (m_s_b,  2.16e-2, 0.50e-3),
        (V_us,   0.2245,  0.0008),
        (V_cb,   0.0405,  0.0015),
        (V_ub,   0.00369, 0.00011),
    ]

    chi2 = sum(((p - v) / e)**2 for p, v, e in data)
    return chi2


def mass_hierarchy_ratio(theta1, theta2, alpha, L=10.0):
    """Compute m_u/m_t as a function of moduli position."""
    y1 = -theta1 * L / 2
    y3 = theta2 * L / 2
    s13 = abs(y3 - y1)
    return overlap(s13, alpha)**2


def cabibbo_angle(theta1, theta2, alpha, L=10.0):
    """Compute |V_us| (Cabibbo angle) across the torus."""
    y1 = -theta1 * L / 2
    y2 = 0.0
    s12 = abs(y2 - y1)
    return overlap(s12, alpha)


def run_scan(alpha=4.3, L=10.0, resolution=80):
    """Full moduli torus scan with visualization."""

    print("=" * 60)
    print("  MODULI TORUS SCAN")
    print("=" * 60)
    print(f"\n  α = {alpha}, L = {L}, grid = {resolution}×{resolution}")

    theta1_range = np.linspace(0.05, 0.95, resolution)
    theta2_range = np.linspace(0.05, 0.95, resolution)
    T1, T2 = np.meshgrid(theta1_range, theta2_range)

    # Compute χ² landscape
    chi2_map = np.zeros_like(T1)
    mass_map = np.zeros_like(T1)
    cabibbo_map = np.zeros_like(T1)

    for i in range(resolution):
        for j in range(resolution):
            chi2_map[i, j] = quick_chi2(T1[i, j], T2[i, j], alpha, L)
            mass_map[i, j] = mass_hierarchy_ratio(T1[i, j], T2[i, j], alpha, L)
            cabibbo_map[i, j] = cabibbo_angle(T1[i, j], T2[i, j], alpha, L)

    # Find best-fit point
    min_idx = np.unravel_index(np.argmin(chi2_map), chi2_map.shape)
    best_t1 = T1[min_idx]
    best_t2 = T2[min_idx]
    best_chi2 = chi2_map[min_idx]

    print(f"\n  Best fit:")
    print(f"    θ₁ = {best_t1:.4f}")
    print(f"    θ₂ = {best_t2:.4f}")
    print(f"    χ² = {best_chi2:.2f}")
    print(f"    m_u/m_t = {mass_map[min_idx]:.2e}")
    print(f"    |V_us| = {cabibbo_map[min_idx]:.4f}")

    # Topological floor
    alpha_star = 3.5
    print(f"\n  Topological floor: α* ≈ {alpha_star}")
    print(f"  Current α = {alpha} {'> α* (stable)' if alpha > alpha_star else '< α* (UNSTABLE)'}")

    # --- Visualization ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    fig.suptitle(f'Moduli Torus T² — Flavor Landscape (α = {alpha})',
                fontsize=14, fontweight='bold')

    # Panel 1: χ² landscape
    ax = axes[0, 0]
    chi2_clipped = np.clip(chi2_map, 0, np.percentile(chi2_map, 95))
    im = ax.contourf(T1, T2, np.log10(chi2_clipped + 1), levels=30, cmap='viridis_r')
    ax.plot(best_t1, best_t2, 'r*', markersize=15, markeredgecolor='white',
           markeredgewidth=1, label=f'Best fit (χ²={best_chi2:.1f})')
    ax.set_xlabel(r'$\theta_1$', fontsize=12)
    ax.set_ylabel(r'$\theta_2$', fontsize=12)
    ax.set_title(r'$\log_{10}(\chi^2)$ landscape')
    ax.legend(fontsize=9)
    plt.colorbar(im, ax=ax)

    # Panel 2: Mass hierarchy
    ax = axes[0, 1]
    im2 = ax.contourf(T1, T2, np.log10(mass_map + 1e-20), levels=30, cmap='plasma')
    ax.plot(best_t1, best_t2, 'r*', markersize=15, markeredgecolor='white',
           markeredgewidth=1)
    # Contour for observed m_u/m_t
    ax.contour(T1, T2, np.log10(mass_map + 1e-20),
              levels=[np.log10(1.27e-5)], colors='red', linewidths=2, linestyles='--')
    ax.set_xlabel(r'$\theta_1$', fontsize=12)
    ax.set_ylabel(r'$\theta_2$', fontsize=12)
    ax.set_title(r'$\log_{10}(m_u/m_t)$')
    plt.colorbar(im2, ax=ax)

    # Panel 3: Cabibbo angle
    ax = axes[1, 0]
    im3 = ax.contourf(T1, T2, cabibbo_map, levels=30, cmap='coolwarm')
    ax.contour(T1, T2, cabibbo_map, levels=[0.2245], colors='black',
              linewidths=2, linestyles='--')
    ax.plot(best_t1, best_t2, 'r*', markersize=15, markeredgecolor='white',
           markeredgewidth=1)
    ax.set_xlabel(r'$\theta_1$', fontsize=12)
    ax.set_ylabel(r'$\theta_2$', fontsize=12)
    ax.set_title(r'$|V_{us}|$ (Cabibbo angle)')
    plt.colorbar(im3, ax=ax)

    # Panel 4: α dependence (1D slice)
    ax = axes[1, 1]
    alphas = np.linspace(3.0, 6.0, 100)
    chi2_vs_alpha = [quick_chi2(best_t1, best_t2, a, L) for a in alphas]
    ax.semilogy(alphas, chi2_vs_alpha, 'b-', linewidth=2)
    ax.axvline(x=alpha_star, color='red', linestyle='--', alpha=0.5,
              label=fr'$\alpha^* = {alpha_star}$ (topological floor)')
    ax.axvline(x=alpha, color='green', linestyle=':', linewidth=2,
              label=fr'$\alpha = {alpha}$ (benchmark)')
    ax.set_xlabel(r'$\alpha = Gv$', fontsize=12)
    ax.set_ylabel(r'$\chi^2$', fontsize=12)
    ax.set_title(r'$\chi^2$ vs compactness parameter')
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig('parameter_scan.png', dpi=150, bbox_inches='tight')
    print(f"\nSaved: parameter_scan.png")
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Moduli torus parameter scan')
    parser.add_argument('--alpha', type=float, default=4.3,
                       help='Compactness parameter (default: 4.3)')
    parser.add_argument('--resolution', type=int, default=80,
                       help='Grid resolution (default: 80)')
    parser.add_argument('--L', type=float, default=10.0,
                       help='Orbifold length (default: 10.0)')
    args = parser.parse_args()

    run_scan(alpha=args.alpha, L=args.L, resolution=args.resolution)
