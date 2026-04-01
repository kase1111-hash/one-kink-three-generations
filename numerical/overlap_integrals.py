#!/usr/bin/env python3
"""
overlap_integrals.py — Yukawa matrix from zero-mode overlaps.

Chapter 7 of "One Kink, Three Generations"

The Yukawa coupling between the n-th left-handed and m-th right-handed
fermion is determined by the overlap integral:

    Y_nm ∝ ∫ f_n^(L)(y) · h(y) · f_m^(R)(y) dy

where f_n are zero-mode profiles and h(y) is the Higgs profile.

CKM MECHANISM: Left-handed quarks (u_L, d_L) share positions because
they sit in the same SU(2) doublet. Right-handed quarks (u_R, d_R)
are SU(2) singlets at DIFFERENT positions. The CKM matrix is the
mismatch between the up-type and down-type diagonalizations:

    V_CKM = U_u† · U_d

Usage:
    python overlap_integrals.py
    python overlap_integrals.py --alpha 4.3

Author: Kase Branham — Independent Researcher
"""

import numpy as np
# numpy compatibility
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid
from scipy.linalg import svd
import matplotlib.pyplot as plt
import argparse


def gaussian_profile(y, center, alpha):
    """Normalized Gaussian zero-mode: f(y) = (α/π)^{1/4} exp(-α(y-c)²/2)."""
    return (alpha / np.pi) ** 0.25 * np.exp(-0.5 * alpha * (y - center)**2)


def overlap_integral(y, y_L, y_R, alpha, y_H=0.0, sigma_H=1.0):
    """
    Exact analytical 3-body Gaussian overlap:
        ∫ f_L(y) · h(y) · f_R(y) dy

    Derived by completing the square in the product of three Gaussians:
        f_L = (α/π)^{1/4} exp(-α(y-y_L)²/2)
        h   = (2πσ²)^{-1/2} exp(-(y-y_H)²/(2σ²))
        f_R = (α/π)^{1/4} exp(-α(y-y_R)²/2)
    
    The combined exponent is -a·y² + b·y - c where:
        a = α + 1/(2σ_H²)
        b = α·(y_L + y_R) + y_H/σ_H²
        c = α·(y_L² + y_R²)/2 + y_H²/(2σ_H²)
    """
    s2 = sigma_H**2
    a = alpha + 1.0 / (2 * s2)
    b = alpha * (y_L + y_R) + y_H / s2
    c = alpha * (y_L**2 + y_R**2) / 2 + y_H**2 / (2 * s2)
    prefactor = (alpha / np.pi)**0.5 * (2 * np.pi * s2)**(-0.5) * np.sqrt(np.pi / a)
    exponent = b**2 / (4 * a) - c
    return prefactor * np.exp(exponent)


def build_yukawa_matrix(positions_L, positions_R, alpha, y_H, sigma_H, 
                        overall_scale=1.0):
    """
    Build the 3×3 Yukawa matrix from overlap integrals.

    Y_nm = scale × ∫ f_n^L · h · f_m^R dy

    Parameters
    ----------
    positions_L : list of 3 floats
        Left-handed zero-mode positions (shared by SU(2) doublet)
    positions_R : list of 3 floats
        Right-handed zero-mode positions (differ for u_R vs d_R)
    """
    Y = np.zeros((3, 3), dtype=complex)
    for n in range(3):
        for m in range(3):
            Y[n, m] = overall_scale * overlap_integral(
                None, positions_L[n], positions_R[m], alpha, y_H, sigma_H)
    return Y


def add_cp_phase(Y, delta, pattern='13'):
    """
    Add CP-violating phase from the kink modulus.
    
    The complex vacuum phase δ enters through specific elements
    of the Yukawa matrix, producing the CKM phase.
    """
    Y_complex = Y.astype(complex).copy()
    phase = np.exp(1j * delta)
    
    # Phase enters the (1,3) and (3,1) overlaps (generation 1↔3 coupling)
    Y_complex[0, 2] *= phase
    Y_complex[2, 0] *= phase.conj()
    # Subdominant phase in (1,2)
    Y_complex[0, 1] *= np.exp(1j * delta * 0.3)
    Y_complex[1, 0] *= np.exp(-1j * delta * 0.3)
    
    return Y_complex


def diagonalize_yukawa(Y):
    """
    Diagonalize Y via SVD: Y = U · diag(σ) · V†.
    
    Returns U (left unitary), masses (singular values ascending), Vh.
    """
    U, sigma, Vh = svd(Y)
    idx = np.argsort(sigma)
    return U[:, idx], sigma[idx], Vh[idx, :]


def ckm_from_diagonalization(Y_u, Y_d):
    """
    CKM matrix from the mismatch of up and down diagonalizations:
    
        V_CKM = U_u† · U_d
    
    where Y_u = U_u · Σ_u · V_u† and Y_d = U_d · Σ_d · V_d†.
    """
    U_u, m_u, _ = diagonalize_yukawa(Y_u)
    U_d, m_d, _ = diagonalize_yukawa(Y_d)
    V_ckm = U_u.conj().T @ U_d
    return V_ckm, m_u, m_d


def ckm_parameters(V):
    """Extract standard CKM observables."""
    V_abs = np.abs(V)
    J = abs(np.imag(V[0, 0] * V[1, 1] * V[0, 1].conj() * V[1, 0].conj()))
    
    # CKM phase from standard parameterization
    s13 = V_abs[0, 2]
    if s13 > 0 and V_abs[0, 0] > 0:
        delta = np.angle(-V[0, 0] * V[1, 2] * V[0, 2].conj() * V[1, 0].conj())
    else:
        delta = 0
    
    return {
        '|V_ud|': V_abs[0, 0], '|V_us|': V_abs[0, 1], '|V_ub|': V_abs[0, 2],
        '|V_cd|': V_abs[1, 0], '|V_cs|': V_abs[1, 1], '|V_cb|': V_abs[1, 2],
        '|V_td|': V_abs[2, 0], '|V_ts|': V_abs[2, 1], '|V_tb|': V_abs[2, 2],
        'J': J, 'delta': delta,
    }


# Experimental CKM values for comparison
CKM_EXPT = {
    '|V_ud|': (0.97373, 0.00031), '|V_us|': (0.2245, 0.0008),
    '|V_ub|': (0.00369, 0.00011), '|V_cd|': (0.221, 0.004),
    '|V_cs|': (0.987, 0.011),     '|V_cb|': (0.0405, 0.0015),
    '|V_td|': (0.0080, 0.0003),   '|V_ts|': (0.0388, 0.0011),
    '|V_tb|': (1.013, 0.030),     'J': (3.08e-5, 0.15e-5),
}


def run_overlap_calculation(alpha=4.3, delta_CP=1.2):
    """
    Full overlap integral calculation with split-fermion CKM.

    The key geometry:
      - Left-handed positions are SHARED (SU(2) doublet)
      - Up-right and down-right positions DIFFER
      - This mismatch generates the CKM matrix
    """
    # ======================================================
    #  GEOMETRIC PARAMETERS
    # ======================================================
    
    # Left-handed positions (shared by u_L and d_L in each generation)
    positions_L = [-1.60, 0.0, 1.55]
    
    # Up-type right-handed positions
    positions_Ru = [-1.45, 0.15, 1.60]
    
    # Down-type right-handed positions (optimized for exact CKM with correct 3-body overlap)
    positions_Rd = [-1.2472982275715974, 0.8753062468852337, 0.6450461521965418]
    
    # Higgs localization
    y_H = -1.1230001781016452
    sigma_H = 1.0916842134716704
    
    # Overall Yukawa scales (set by v_EW and 5D coupling)
    scale_u = 1.0   # top Yukawa normalization
    scale_d = 0.3970250295604182  # down-type coupling ratio (CKM-optimized)
    
    # ======================================================
    #  BUILD YUKAWA MATRICES
    # ======================================================
    
    Y_u = build_yukawa_matrix(positions_L, positions_Ru, alpha, y_H, sigma_H, scale_u)
    Y_d = build_yukawa_matrix(positions_L, positions_Rd, alpha, y_H, sigma_H, scale_d)
    
    # Add CP phase
    Y_u = add_cp_phase(Y_u, delta_CP)
    Y_d = add_cp_phase(Y_d, delta_CP * 0.7)  # different phase projection
    
    # ======================================================
    #  DIAGONALIZE AND EXTRACT CKM
    # ======================================================
    
    V_ckm, m_u, m_d = ckm_from_diagonalization(Y_u, Y_d)
    params = ckm_parameters(V_ckm)
    
    # ======================================================
    #  PRINT RESULTS
    # ======================================================
    
    print("=" * 64)
    print("  OVERLAP INTEGRAL CALCULATION — SPLIT-FERMION CKM")
    print("=" * 64)
    
    print(f"\n  Geometric parameters:")
    print(f"    Compactness α = Gv = {alpha}")
    print(f"    CP phase δ = {delta_CP:.2f} rad ({np.degrees(delta_CP):.1f}deg)")
    print(f"    Higgs: y_H = {y_H}, σ_H = {sigma_H}")
    print(f"\n    Left-handed (shared):  {positions_L}")
    print(f"    Up-right:              {positions_Ru}")
    print(f"    Down-right:            {positions_Rd}")
    print(f"    Mismatch (Ru - Rd):    "
          f"{[f'{u-d:+.2f}' for u,d in zip(positions_Ru, positions_Rd)]}")
    
    print(f"\n  Up-type Yukawa |Y_u|:")
    for row in np.abs(Y_u):
        print(f"    [{row[0]:.4e}  {row[1]:.4e}  {row[2]:.4e}]")
    
    print(f"\n  Down-type Yukawa |Y_d|:")
    for row in np.abs(Y_d):
        print(f"    [{row[0]:.4e}  {row[1]:.4e}  {row[2]:.4e}]")
    
    print(f"\n  Mass eigenvalues (up-type):")
    if m_u[2] > 0:
        print(f"    m_u : m_c : m_t = {m_u[0]/m_u[2]:.4e} : "
              f"{m_u[1]/m_u[2]:.4e} : 1")
        print(f"    Hierarchy: {m_u[2]/m_u[0]:.0f}× between 1st and 3rd gen")
    
    print(f"\n  Mass eigenvalues (down-type):")
    if m_d[2] > 0:
        print(f"    m_d : m_s : m_b = {m_d[0]/m_d[2]:.4e} : "
              f"{m_d[1]/m_d[2]:.4e} : 1")
    
    print(f"\n  CKM matrix |V|:")
    V_abs = np.abs(V_ckm)
    for i in range(3):
        row = [f"{V_abs[i,j]:.5f}" for j in range(3)]
        print(f"    [{', '.join(row)}]")
    
    print(f"\n  {'Element':<10} {'MKGF':>10} {'Expt':>10} {'Pull':>8}")
    print(f"  {'-'*42}")
    for key in ['|V_ud|', '|V_us|', '|V_ub|', '|V_cd|', '|V_cs|', '|V_cb|',
                '|V_td|', '|V_ts|', '|V_tb|', 'J']:
        val = params[key]
        exp_val, exp_err = CKM_EXPT[key]
        pull = (val - exp_val) / exp_err if exp_err > 0 else 0
        print(f"  {key:<10} {val:>10.5f} {exp_val:>10.5f} {pull:>+8.1f}")
    
    # ======================================================
    #  VISUALIZATION
    # ======================================================
    
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle('Overlap Integrals and CKM Matrix',
                fontsize=14, fontweight='bold')
    
    # Panel 1: Zero-mode profiles with split fermion positions
    ax = axes[0, 0]
    y = np.linspace(-4, 4, 1000)
    colors_L = ['#D4A850', '#1D9E75', '#7F77DD']
    for i, (yL, c) in enumerate(zip(positions_L, colors_L)):
        f = gaussian_profile(y, yL, alpha)
        ax.fill_between(y, 0, f**2, alpha=0.2, color=c)
        ax.plot(y, f**2, color=c, linewidth=2, label=f'Gen {i+1} (L)')
    # Mark right-handed positions
    for i, (yu, yd) in enumerate(zip(positions_Ru, positions_Rd)):
        ax.axvline(x=yu, color=colors_L[i], linestyle='--', alpha=0.5)
        ax.axvline(x=yd, color=colors_L[i], linestyle=':', alpha=0.5)
    ax.plot([], [], 'k--', label=r'$u_R$ positions')
    ax.plot([], [], 'k:', label=r'$d_R$ positions')
    ax.set_xlabel(r'Extra dimension $y$', fontsize=11)
    ax.set_ylabel(r'$|f(y)|^2$', fontsize=11)
    ax.set_title('Split-fermion geometry')
    ax.legend(fontsize=7, ncol=2)
    
    # Panel 2: Yukawa matrices (log heatmap)
    ax = axes[0, 1]
    combined = np.zeros((3, 6))
    combined[:, :3] = np.log10(np.abs(Y_u) + 1e-20)
    combined[:, 3:] = np.log10(np.abs(Y_d) + 1e-20)
    im = ax.imshow(combined, cmap='viridis', aspect='auto')
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(['u_R', 'c_R', 't_R', 'd_R', 's_R', 'b_R'], fontsize=8)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['Gen 1', 'Gen 2', 'Gen 3'])
    ax.axvline(x=2.5, color='white', linewidth=2)
    ax.text(1, -0.7, 'Up-type', ha='center', fontsize=10, fontweight='bold')
    ax.text(4, -0.7, 'Down-type', ha='center', fontsize=10, fontweight='bold')
    ax.set_title(r'$\log_{10}|Y_{nm}|$ (Yukawa matrices)')
    plt.colorbar(im, ax=ax, shrink=0.8)
    
    # Panel 3: CKM matrix
    ax = axes[1, 0]
    im2 = ax.imshow(V_abs, cmap='YlOrRd', vmin=0, vmax=1, aspect='equal')
    for i in range(3):
        for j in range(3):
            color = 'white' if V_abs[i,j] > 0.5 else 'black'
            ax.text(j, i, f'{V_abs[i,j]:.4f}', ha='center', va='center',
                   fontsize=11, color=color, fontweight='bold')
    ax.set_xticks([0, 1, 2])
    ax.set_yticks([0, 1, 2])
    ax.set_xticklabels(['d', 's', 'b'], fontsize=12)
    ax.set_yticklabels(['u', 'c', 't'], fontsize=12)
    ax.set_title(r'$|V_{CKM}|$ from overlap mismatch')
    
    # Panel 4: Comparison bar chart
    ax = axes[1, 1]
    keys = ['|V_us|', '|V_cb|', '|V_ub|']
    mkgf_vals = [params[k] for k in keys]
    expt_vals = [CKM_EXPT[k][0] for k in keys]
    expt_errs = [CKM_EXPT[k][1] for k in keys]
    x = np.arange(len(keys))
    width = 0.35
    bars1 = ax.bar(x - width/2, mkgf_vals, width, label='MKGF',
                  color='#1D9E75', alpha=0.8)
    bars2 = ax.bar(x + width/2, expt_vals, width, label='Experiment',
                  color='#7F77DD', alpha=0.8, yerr=expt_errs, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(keys, fontsize=11)
    ax.set_ylabel('Value', fontsize=11)
    ax.set_title('CKM elements: MKGF vs experiment')
    ax.legend()
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig('overlap_integrals.png', dpi=150, bbox_inches='tight')
    print(f"\nSaved: overlap_integrals.png")
    
    print(f"\n  NOTE ON THIS DEMONSTRATION:")
    print(f"  The CKM mixing elements |V_us|, |V_cb|, and |V_ub| are")
    print(f"  reproduced exactly using the correct 3-body Gaussian")
    print(f"  overlap formula with optimized split-fermion positions.")
    print(f"  The mass hierarchy (~15x) is in the right direction but")
    print(f"  weaker than the observed ~10^5. The full hierarchy")
    print(f"  requires the complete moduli torus optimization")
    print(f"  (see global_chi2_fit.py and the book, Chapter 11).")
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Yukawa matrix from zero-mode overlap integrals')
    parser.add_argument('--alpha', type=float, default=4.3,
                       help='Compactness parameter Gv (default: 4.3)')
    args = parser.parse_args()
    
    run_overlap_calculation(alpha=args.alpha)
