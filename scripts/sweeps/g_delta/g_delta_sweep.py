#!/usr/bin/env python3
"""
g-delta Grid Sweep - 101x51 grid over delta and g
Computes ME Noise and Analytical Noise for a single (g, delta) point.

Parameters:
    kappa   = 2*pi*0.05
    epsilon = 2*pi*0.05
    gamma   = 2*pi*0.001
    N = 20
    t_end = 200
    num_points = 4000
    phi_h = 0 (fixed)
    phi_d = pi/2 (fixed)
    
Sweep:
    delta: from -5*2*pi to +5*2*pi in steps of pi/5 (101 values)
    g: from 0.001*2*pi to 0.5*2*pi in steps of pi/50 (51 values)
    
Usage:
    python g_delta_sweep.py --g 0.628 --delta 6.28 --output_dir results/
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


def run_simulation(g, delta, output_dir):
    """
    Run noise simulation for a single (g, delta) pair.
    """
    print("="*60)
    print("  G-DELTA GRID SWEEP")
    print("="*60)
    print(f"  g: {g:.6f} rad/s ({g/(2*np.pi):.6f} * 2pi)")
    print(f"  delta: {delta:.6f} rad/s ({delta/(2*np.pi):.6f} * 2pi)")
    print(f"  Output: {output_dir}")
    print("="*60)
    
    start_time = time.time()
    
    # Solver options
    solver_opts = {
        'rtol': 1e-6, 'atol': 1e-8, 'nsteps': 5000,
        'store_states': True, 'store_final_state': True
    }

    # -- Parameters (as specified) --
    kappa   = 2*np.pi*0.05 
    epsilon = 2*np.pi*0.05     # As specified
    gamma   = 2*np.pi*0.001
    
    # Fixed phases
    phi_h = 0.0
    phi_d = np.pi / 2.0
    
    # Settings (as specified)
    N = 20           # Hilbert Space Size (as specified)
    t_end = 200
    num_points = 4000
    
    # -- Derived --
    # delta = wa - wd is now an input parameter
    # We set wr = wd = 2*pi*1.0 and compute wa from delta
    wr = 2*np.pi*1.0
    wd = 2*np.pi*1.0
    wa = wd + delta  # Since delta = wa - wd
    
    # Check for valid chi calculation (avoid division by zero)
    if np.abs(delta) < 1e-10:
        print("WARNING: delta ≈ 0, chi would be undefined. Using small delta.")
        delta_safe = np.sign(delta) * 1e-6 if delta != 0 else 1e-6
        chi = g**2 / delta_safe
    else:
        chi = g**2 / delta
    
    phi_qb = 2 * np.arctan(2 * chi / kappa) 
    
    print(f"\n  N={N}, t_end={t_end}, num_points={num_points}")
    print(f"  g: {g/(2*np.pi):.6f} * 2pi")
    print(f"  kappa: {kappa/(2*np.pi):.4f} * 2pi")
    print(f"  epsilon: {epsilon/(2*np.pi):.4f} * 2pi")
    print(f"  gamma: {gamma/(2*np.pi):.4f} * 2pi")
    print(f"  delta (wa-wd): {delta/(2*np.pi):.4f} * 2pi")
    print(f"  chi: {chi/(2*np.pi):.6f} * 2pi")
    print(f"  phi_h: {phi_h:.4f} rad, phi_d: {phi_d:.4f} rad")
    print(f"  phi_qb: {phi_qb:.6f} rad")

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
    
    # Hamiltonian (in rotating frame at drive frequency wd)
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
    
    # State 0
    G_N_0_raw, G_A_0_raw = results[0], results[1]
    Cov_N_0 = G_N_0_raw - np.outer(np.conj(mean_a_0), mean_a_0)
    Cov_A_0 = G_A_0_raw - np.outer(mean_a_0, mean_a_0)

    # State 1
    G_N_1_raw, G_A_1_raw = results[2], results[3]
    Cov_N_1 = G_N_1_raw - np.outer(np.conj(mean_a_1), mean_a_1)
    Cov_A_1 = G_A_1_raw - np.outer(mean_a_1, mean_a_1)

    # Calculate Variance directly from Covariances
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
    
    filename = f"g_delta_g{g/(2*np.pi):.6f}_delta{delta/(2*np.pi):.6f}.npz"
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
    print(f">> Max Analytic Noise: {np.max(Noise_Analytic):.6f}")
    print("="*60)
    print("  SIMULATION COMPLETE")
    print("="*60)
    sys.stdout.flush()
    
    return filepath


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="g-delta Grid Sweep")
    parser.add_argument("--g", type=float, required=True, help="Coupling strength g (in rad/s)")
    parser.add_argument("--delta", type=float, required=True, help="Detuning delta = wa - wd (in rad/s)")
    parser.add_argument("--output_dir", type=str, default="results", help="Output directory")
    args = parser.parse_args()
    
    print(f"Python version: {sys.version}")
    print(f"QuTiP version: {qutip.__version__}")
    print(f"NumPy version: {np.__version__}")
    sys.stdout.flush()
    
    run_simulation(args.g, args.delta, args.output_dir)
