#!/usr/bin/env python3
"""
Kappa-Epsilon Grid Sweep - 2D sweep over kappa and epsilon

Sweep both kappa and epsilon each from 0.04 to 0.06 in steps of 0.01

Parameters:
    kappa   = sweep parameter
    epsilon = sweep parameter
    gamma   = 2*pi*0.001
    N = 30
    t_end = 200
    num_points = 2000
    phi_h = 0 (fixed)
    phi_d = pi/2 (fixed)
    delta = +5*2*pi (fixed)
    g = 0.341*2*pi (fixed)
    
Usage:
    python kappa_epsilon_grid_sweep.py --kappa 0.05 --epsilon 0.05 --output_dir results/
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
    lower = np.tril(G_raw)
    G_A = lower + lower.T - np.diag(np.diag(lower))
    
    phase_factor = np.exp(-2j * phi_h)
    integrand = np.real(G_N - phase_factor * G_A)
    
    int_1 = cumtrapz(integrand, tlist, axis=0, initial=0)
    double_integral = np.diagonal(cumtrapz(int_1, tlist, axis=1, initial=0))
    
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


def run_simulation(kappa_val, epsilon_val, output_dir):
    """
    Run noise simulation for a single (kappa, epsilon) pair.
    """
    print("="*60)
    print("  KAPPA-EPSILON GRID SWEEP")
    print("="*60)
    print(f"  kappa: {kappa_val:.6f} * 2pi")
    print(f"  epsilon: {epsilon_val:.6f} * 2pi")
    print(f"  Output: {output_dir}")
    print("="*60)
    
    start_time = time.time()
    
    solver_opts = {
        'rtol': 1e-6, 'atol': 1e-8, 'nsteps': 5000,
        'store_states': True, 'store_final_state': True
    }

    # -- Parameters --
    kappa   = 2*np.pi*kappa_val
    epsilon = 2*np.pi*epsilon_val
    gamma   = 2*np.pi*0.001
    
    # Fixed parameters
    g = 0.341 * 2*np.pi
    delta = 5.0 * 2*np.pi
    
    phi_h = 0.0
    phi_d = np.pi / 2.0
    
    N = 30  # Larger Hilbert space as specified
    t_end = 200
    num_points = 2000
    
    # -- Derived --
    wr = 2*np.pi*1.0
    wd = 2*np.pi*1.0
    wa = wd + delta
    
    chi = g**2 / delta
    phi_qb = 2 * np.arctan(2 * chi / kappa)
    
    print(f"\n  N={N}, t_end={t_end}, num_points={num_points}")
    print(f"  g: {g/(2*np.pi):.6f} * 2pi")
    print(f"  kappa: {kappa/(2*np.pi):.4f} * 2pi")
    print(f"  epsilon: {epsilon/(2*np.pi):.4f} * 2pi")
    print(f"  gamma: {gamma/(2*np.pi):.6f} * 2pi")
    print(f"  delta: {delta/(2*np.pi):.4f} * 2pi")
    print(f"  chi: {chi/(2*np.pi):.6f} * 2pi")
    print(f"  phi_h: {phi_h:.4f} rad, phi_d: {phi_d:.4f} rad")

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
    
    H = (
        (wr - wd) * a.dag() * a +
        (wa - wd) * sz / 2.0 +
        g * (a * sm.dag() + a.dag() * sm) +
        epsilon * (a.dag() * np.exp(1j * phi_d) + a * np.exp(-1j * phi_d))
    )
    c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma) * sm]
    
    alpha_in_scalar = -1j * (epsilon / np.sqrt(kappa)) * np.exp(1j * phi_d)
    a_in_op = alpha_in_scalar * tensor(qeye(N), qeye(2))
    a_out = np.sqrt(kappa) * a - a_in_op
    XX_op = 1j * (a_out.dag() * np.exp(1j * phi_h) - a_out * np.exp(-1j * phi_h))
    
    Y_cav_op = 1j * (a.dag() * np.exp(1j * phi_h) - a * np.exp(-1j * phi_h))

    # ============================================================
    # ME Signal Calculation
    # ============================================================
    print(">> [ME] Calculating Signal...")
    sys.stdout.flush()
    signal_start = time.time()
    
    res0 = mesolve(H, psi0, tlist, c_ops, e_ops=[XX_op, Y_cav_op, a], options=solver_opts)
    res1 = mesolve(H, psi1, tlist, c_ops, e_ops=[XX_op, Y_cav_op, a], options=solver_opts)

    mean_a_0 = res0.expect[2]
    mean_a_1 = res1.expect[2]
    
    M0_ME = np.sqrt(kappa) * cumtrapz(np.array(res0.expect[0]), tlist, initial=0)
    M1_ME = np.sqrt(kappa) * cumtrapz(np.array(res1.expect[0]), tlist, initial=0)
    Signal_ME = np.abs(M0_ME - M1_ME)
    
    M0_cav = np.sqrt(kappa) * cumtrapz(np.array(res0.expect[1]), tlist, initial=0)
    M1_cav = np.sqrt(kappa) * cumtrapz(np.array(res1.expect[1]), tlist, initial=0)
    
    print(f"   Signal completed in {time.time()-signal_start:.2f}s")
    sys.stdout.flush()

    # ============================================================
    # ME Noise Calculation
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
    
    G_N_0_raw, G_A_0_raw = results[0], results[1]
    Cov_N_0 = G_N_0_raw - np.outer(np.conj(mean_a_0), mean_a_0)
    Cov_A_0 = G_A_0_raw - np.outer(mean_a_0, mean_a_0)

    G_N_1_raw, G_A_1_raw = results[2], results[3]
    Cov_N_1 = G_N_1_raw - np.outer(np.conj(mean_a_1), mean_a_1)
    Cov_A_1 = G_A_1_raw - np.outer(mean_a_1, mean_a_1)

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
    
    angle_factor = np.sin(phi_d - phi_h)
    term_decay = 1 - (np.sin(chi * tlist + phi_qb) / np.sin(phi_qb)) * np.exp(-kappa_tau / 2)
    bracket = 1 - (4 / kt_safe) * (np.cos(phi_qb / 2)**2) * term_decay
    Signal_Analytic = np.abs(4 * epsilon * np.sin(phi_qb) * tlist * angle_factor * bracket)
    
    Noise_Analytic = np.sqrt(2 * kappa_tau)

    # ============================================================
    # Save Results
    # ============================================================
    os.makedirs(output_dir, exist_ok=True)
    
    filename = f"kappa_epsilon_grid_sweep_kappa{kappa_val:.6f}_eps{epsilon_val:.6f}.npz"
    filepath = os.path.join(output_dir, filename)
    
    np.savez(filepath,
             g=g,
             delta=delta,
             phi_h=phi_h,
             phi_d=phi_d,
             kappa=kappa,
             epsilon=epsilon,
             gamma=gamma,
             chi=chi,
             phi_qb=phi_qb,
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
    print("="*60)
    print("  SIMULATION COMPLETE")
    print("="*60)
    sys.stdout.flush()
    
    return filepath


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gamma-Kappa Sweep")
    parser.add_argument("--kappa", type=float, required=True, help="Kappa value (multiplied by 2pi)")
    parser.add_argument("--epsilon", type=float, required=True, help="Epsilon value (multiplied by 2pi)")
    parser.add_argument("--output_dir", type=str, default="results", help="Output directory")
    args = parser.parse_args()
    
    print(f"Python version: {sys.version}")
    print(f"QuTiP version: {qutip.__version__}")
    print(f"NumPy version: {np.__version__}")
    sys.stdout.flush()
    
    run_simulation(args.kappa, args.epsilon, args.output_dir)
