#!/usr/bin/env python3
"""
Noise Sweep - Fixed phi_h=0, sweep phi_d only
Computes ME Noise and Analytical Noise for phi_d sweep with phi_h=0.

Parameters:
    g = 2*pi*0.2
    epsilon = 2*pi*0.04
    N = 14
    t_end = 200
    num_points = 4000
    solver: rtol=1e-6, atol=1e-8

Usage:
    python noise_sweep.py --phi_d 1.57 --output_dir results/
"""

import os
import argparse
import numpy as np
import scipy.integrate as integrate
from qutip import *
import concurrent.futures
import time
import sys

# ============================================================
# Library Compatibility
# ============================================================
if hasattr(integrate, 'cumulative_trapezoid'):
    cumtrapz = integrate.cumulative_trapezoid
else:
    cumtrapz = integrate.cumtrapz


# ============================================================
# Helper Functions for ME Noise Calculation
# ============================================================

def finalize_variance(G_N, G_raw, tlist, kappa, phi_h):
    """
    Vectorized calculation of variance.
    """
    # 1. Symmetrize G_raw - TIME-ORDERED product uses lower triangle
    lower = np.tril(G_raw)
    G_A = lower + lower.T - np.diag(np.diag(lower))
    
    # 2. Calculate Integrand Matrix
    phase_factor = np.exp(-2j * phi_h)
    integrand = np.real(G_N - phase_factor * G_A)
    
    # 3. 2D Integration using Cumulative Trapezoid
    int_1 = cumtrapz(integrand, tlist, axis=0, initial=0)
    double_integral = np.diagonal(cumtrapz(int_1, tlist, axis=1, initial=0))
    
    # 4. Final Formula
    var_list = (kappa * tlist) + (2 * kappa**2 * double_integral)
    
    return var_list


def correlation_worker(args):
    """
    Worker to calculate 2-time correlations.
    """
    task_idx, H, psi, tlist, c_ops, op1, op2, opts = args
    
    corr1 = correlation_2op_2t(H, psi, tlist, tlist, c_ops, op1, op2, 
                               reverse=False, options=opts)
    corr2 = correlation_2op_2t(H, psi, tlist, tlist, c_ops, op1, op2, 
                               reverse=True, options=opts)

    num_points = len(tlist)
    new_size = num_points * 2
    
    corr1_mapped = np.zeros((new_size, new_size), dtype=complex)
    corr2_mapped = np.zeros((new_size, new_size), dtype=complex)
    
    scale_factor = (new_size - 1) / (2 * tlist[-1])
    t_indices = (np.array(tlist) * scale_factor).astype(int)
    
    for i, ti in enumerate(t_indices):
        for j, tauj in enumerate(t_indices):
            sum_idx = ti + tauj
            if sum_idx < new_size and sum_idx < num_points:
                corr1_mapped[sum_idx, ti] = corr1[i, j]
            if sum_idx < new_size and sum_idx < num_points:
                corr2_mapped[ti, sum_idx] = corr2[i, j]

    np.fill_diagonal(corr1_mapped, 0)
    corr_sum = corr1_mapped + corr2_mapped
    
    return task_idx, corr_sum[:num_points, :num_points]


def run_simulation(phi_d, output_dir):
    """
    Run noise simulation for phi_d with phi_h=0 fixed.
    """
    phi_h = 0.0  # FIXED
    
    print("="*60)
    print("  NOISE SWEEP - Fixed phi_h=0, sweep phi_d")
    print("="*60)
    print(f"  phi_h: {phi_h:.6f} rad (FIXED)")
    print(f"  phi_d: {phi_d:.6f} rad")
    print(f"  Output: {output_dir}")
    print("="*60)
    
    start_time = time.time()
    
    # Solver options (as specified)
    solver_opts = {
        'rtol': 1e-6, 'atol': 1e-8, 'nsteps': 5000,
        'store_states': True, 'store_final_state': True
    }

    # -- Parameters (as specified) --
    wr, wa, wd = 2*np.pi*1.0, 2*np.pi*3.0, 2*np.pi*1.0
    g       = 2*np.pi*0.2      # As specified
    kappa   = 2*np.pi*0.05 
    epsilon = 2*np.pi*0.04     # As specified
    gamma   = 2*np.pi*0.001
    
    # Settings (as specified)
    N = 14           # Hilbert Space Size (as specified)
    t_end = 200
    num_points = 4000
    
    # -- Derived --
    delta   = wa - wd           
    chi     = g**2 / delta      
    phi_qb  = 2 * np.arctan(2 * chi / kappa) 
    
    print(f"\n  N={N}, t_end={t_end}, num_points={num_points}")
    print(f"  g: {g/(2*np.pi):.4f}, epsilon: {epsilon/(2*np.pi):.4f}")
    print(f"  Delta: {delta/(2*np.pi):.4f}, Chi: {chi/(2*np.pi):.4f}")
    print(f"  Phi_qb: {phi_qb:.4f} rad")

    # -- Time Settings --
    tlist = np.linspace(0, t_end, num_points)
    kappa_tau = kappa * tlist
    kt_safe = np.where(kappa_tau == 0, 1e-12, kappa_tau)

    # ============================================================
    # Build operators
    # ============================================================
    print("\n>> Building operators...")
    sys.stdout.flush()
    
    a  = tensor(destroy(N), qeye(2))
    sz = tensor(qeye(N), sigmaz())
    sm = tensor(qeye(N), sigmam())
    psi0 = tensor(basis(N, 0), basis(2, 0))
    psi1 = tensor(basis(N, 0), basis(2, 1))
    
    # Hamiltonian
    H = (
        (wr - wd) * a.dag() * a +
        (wa - wd) * sz / 2.0 +
        g * (a * sm.dag() + a.dag() * sm) +
        epsilon * (a.dag() * np.exp(1j * phi_d) + a * np.exp(-1j * phi_d))
    )
    c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma) * sm]
    
    # Output operator (Total field = cavity + input)
    alpha_in_scalar = -1j * (epsilon / np.sqrt(kappa)) * np.exp(1j * phi_d)
    a_in_op = alpha_in_scalar * tensor(qeye(N), qeye(2))
    a_out = np.sqrt(kappa) * a - a_in_op
    XX_op = 1j * (a_out.dag() * np.exp(1j * phi_h) - a_out * np.exp(-1j * phi_h))
    
    # Cavity-only Y operator (for variance subtraction)
    Y_cav_op = 1j * (a.dag() * np.exp(1j * phi_h) - a * np.exp(-1j * phi_h))

    # ============================================================
    # ME Signal Calculation
    # ============================================================
    print(">> [ME] Calculating Signal...")
    sys.stdout.flush()
    signal_start = time.time()
    
    # We compute expectation values for BOTH operators AND 'a' (for covariance)
    # index 0: XX_op (Total for Signal)
    # index 1: Y_cav_op (Cavity only for Variance)
    # index 2: a (for covariance calculation)
    res0 = mesolve(H, psi0, tlist, c_ops, e_ops=[XX_op, Y_cav_op, a], options=solver_opts)
    res1 = mesolve(H, psi1, tlist, c_ops, e_ops=[XX_op, Y_cav_op, a], options=solver_opts)

    # Extract Mean Fields <a(t)>
    mean_a_0 = res0.expect[2]
    mean_a_1 = res1.expect[2]
    
    # Signal uses total field (XX_op)
    M0_ME = np.sqrt(kappa) * cumtrapz(np.array(res0.expect[0]), tlist, initial=0)
    M1_ME = np.sqrt(kappa) * cumtrapz(np.array(res1.expect[0]), tlist, initial=0)
    Signal_ME = np.abs(M0_ME - M1_ME)
    
    # Cavity-only means (for Variance subtraction)
    M0_cav = np.sqrt(kappa) * cumtrapz(np.array(res0.expect[1]), tlist, initial=0)
    M1_cav = np.sqrt(kappa) * cumtrapz(np.array(res1.expect[1]), tlist, initial=0)
    
    print(f"   Signal completed in {time.time()-signal_start:.2f}s")
    sys.stdout.flush()

    # ============================================================
    # ME Noise Calculation (4 parallel tasks)
    # ============================================================
    print(">> [ME] Calculating Noise (4 parallel correlation tasks)...")
    sys.stdout.flush()
    noise_start = time.time()
    
    tasks = [
        (0, H, psi0, tlist, c_ops, a.dag(), a, solver_opts), 
        (1, H, psi0, tlist, c_ops, a, a, solver_opts),
        (2, H, psi1, tlist, c_ops, a.dag(), a, solver_opts), 
        (3, H, psi1, tlist, c_ops, a, a, solver_opts)
    ]
    
    unsorted_results = {}
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(correlation_worker, t) for t in tasks]
        for future in concurrent.futures.as_completed(futures):
            idx, res = future.result()
            unsorted_results[idx] = res
            elapsed = time.time() - noise_start
            print(f"   Task {idx} completed ({elapsed:.1f}s elapsed)")
            sys.stdout.flush()

    results = [unsorted_results[i] for i in range(4)]
    
    # --- COVARIANCE CALCULATION (Cov = <AB> - <A><B>) ---
    # We subtract the outer product of means BEFORE integration to avoid cancellation errors
    
    # State 0
    G_N_0_raw, G_A_0_raw = results[0], results[1]
    Cov_N_0 = G_N_0_raw - np.outer(np.conj(mean_a_0), mean_a_0)
    Cov_A_0 = G_A_0_raw - np.outer(mean_a_0, mean_a_0)

    # State 1
    G_N_1_raw, G_A_1_raw = results[2], results[3]
    Cov_N_1 = G_N_1_raw - np.outer(np.conj(mean_a_1), mean_a_1)
    Cov_A_1 = G_A_1_raw - np.outer(mean_a_1, mean_a_1)

    # Calculate Variance directly from Covariances (Note: NO subtraction of M0_cav**2 here!)
    Var0_ME = finalize_variance(Cov_N_0, Cov_A_0, tlist, kappa, phi_h)
    Var1_ME = finalize_variance(Cov_N_1, Cov_A_1, tlist, kappa, phi_h)
    
    Noise_ME = np.sqrt(np.abs(Var0_ME + Var1_ME))
    
    print(f"   Noise completed in {time.time()-noise_start:.2f}s")
    sys.stdout.flush()

    # ============================================================
    # Analytical Calculations
    # ============================================================
    print(">> [Analytic] Computing formulas...")
    sys.stdout.flush()
    
    # Analytic Signal
    angle_factor = np.sin(phi_d - phi_h)
    term_decay = 1 - (np.sin(chi * tlist + phi_qb) / np.sin(phi_qb)) * np.exp(-kappa_tau / 2)
    bracket = 1 - (4 / kt_safe) * (np.cos(phi_qb / 2)**2) * term_decay
    Signal_Analytic = np.abs(4 * epsilon * np.sin(phi_qb) * tlist * angle_factor * bracket)
    
    # Analytic Noise
    Noise_Analytic = np.sqrt(2 * kappa_tau)

    # ============================================================
    # Save Results
    # ============================================================
    os.makedirs(output_dir, exist_ok=True)
    
    filename = f"noise_phih_{phi_h:.4f}_phid_{phi_d:.4f}.npz"
    filepath = os.path.join(output_dir, filename)
    
    np.savez(filepath,
             phi_h=phi_h,
             phi_d=phi_d,
             tlist=tlist,
             kappa_tau=kappa_tau,
             Signal_ME=Signal_ME,
             Noise_ME=Noise_ME,
             Signal_Analytic=Signal_Analytic,
             Noise_Analytic=Noise_Analytic)
    
    elapsed = time.time() - start_time
    
    print(f"\n>> Results saved to: {filepath}")
    print(f">> Total elapsed time: {elapsed:.2f}s ({elapsed/3600:.2f} hours)")
    print(f">> Max ME Signal: {np.max(Signal_ME):.6f}")
    print(f">> Max ME Noise: {np.max(Noise_ME):.6f}")
    print(f">> Max Analytic Noise: {np.max(Noise_Analytic):.6f}")
    print("="*60)
    print("  SIMULATION COMPLETE")
    print("="*60)
    sys.stdout.flush()
    
    return filepath


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Noise Sweep - Fixed phi_h=0")
    parser.add_argument("--phi_d", type=float, required=True, help="Drive phase phi_d")
    parser.add_argument("--output_dir", type=str, default="results", help="Output directory")
    args = parser.parse_args()
    
    print(f"Python version: {sys.version}")
    print(f"QuTiP version: {qutip.__version__}")
    print(f"NumPy version: {np.__version__}")
    sys.stdout.flush()
    
    run_simulation(args.phi_d, args.output_dir)
