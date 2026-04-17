import os

# ============================================================
# 0. PARALLELIZATION FLAGS (Must be set BEFORE imports)
# ============================================================
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from qutip import *
import scipy.integrate as integrate
import concurrent.futures
import time
from tqdm import tqdm 
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

# ============================================================
# 1. STYLE & COMPATIBILITY
# ============================================================
# Library Compatibility
if hasattr(integrate, 'cumulative_trapezoid'):
    cumtrapz = integrate.cumulative_trapezoid
    trapz = integrate.trapezoid
else:
    cumtrapz = integrate.cumtrapz
    trapz = integrate.trapz

# Global LaTeX Style
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 14,
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 10,
    "lines.linewidth": 2
})

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

def finalize_variance(G_N, G_raw, tlist, kappa, phi_h):
    """
    Calculates the variance from the correlation matrices.
    """
    # Extract Time Ordered part - TIME-ORDERED product uses lower triangle
    lower = np.tril(G_raw)
    G_A   = lower + lower.T - np.diag(np.diag(lower))
    
    var_list = np.zeros(len(tlist))
    phase_factor = np.exp(-2j * phi_h)
    integrand = np.real(G_N - phase_factor * G_A)
    
    for i in range(1, len(tlist)):
        sub = integrand[:i+1, :i+1]
        t_sub = tlist[:i+1]
        val = trapz(trapz(sub, t_sub, axis=0), t_sub, axis=0)
        var_list[i] = (kappa * tlist[i]) + (2 * kappa**2 * val)
    return var_list

def correlation_worker(args):
    """
    Worker to calculate 2-time correlations.
    Returns: (task_index, mapped_correlation_matrix)
    """
    task_idx, H, psi, tlist, c_ops, op1, op2, opts = args
    
    # Run QuTiP solver (Forward and Reverse)
    corr1 = correlation_2op_2t(H, psi, tlist, tlist, c_ops, op1, op2, 
                               reverse=False, options=opts)
    corr2 = correlation_2op_2t(H, psi, tlist, tlist, c_ops, op1, op2, 
                               reverse=True, options=opts)

    # Fast Numpy Mapping
    num_points = len(tlist)
    new_size = num_points * 2
    
    corr1_mapped = np.zeros((new_size, new_size), dtype=complex)
    corr2_mapped = np.zeros((new_size, new_size), dtype=complex)
    
    scale_factor = (new_size - 1) / (2 * tlist[-1])
    t_indices    = (np.array(tlist) * scale_factor).astype(int)
    
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

# ============================================================
# 3. MAIN SIMULATION
# ============================================================

def run_final_simulation():
    # Solver options
    solver_opts = {
        'rtol': 1e-4, 'atol': 1e-6, 'nsteps': 5000,
        'store_states': True, 'store_final_state': True #store_states and store_final_state must be True for the correlation_2op_2t function to work
    }

    # -- Parameters --
    wr, wa, wd = 2*np.pi*1.0, 2*np.pi*5.0, 2*np.pi*1.0
    g       = 2*np.pi*0.1      
    kappa   = 2*np.pi*0.05 
    epsilon = 2*np.pi*0.05 
    
    phi_d, phi_h = np.pi/2.0, 0.0              
    
    N = 20  # Hilbert Space Size
    
    # -- Derived --
    delta   = wa - wd           
    chi     = g**2 / delta      
    phi_qb  = 2 * np.arctan(2 * chi / kappa) 

    # -- Analytics & Heuristics --
    n_crit = delta**2 / (4 * g**2)
    n_steady_analytical = (4 * epsilon**2) / (kappa**2)
    
    # ============================================================
    # 4. SYSTEM SETUP & PRE-CHECKS
    # ============================================================
    a  = tensor(destroy(N), qeye(2))
    sz = tensor(qeye(N), sigmaz())
    psi0 = tensor(basis(N, 0), basis(2, 0)) 
    psi1 = tensor(basis(N, 0), basis(2, 1)) 
    sm = tensor(qeye(N), sigmam())

    H = (
        (wr - wd) * a.dag() * a + (wa - wd) * sz / 2
        + g * (a * sm.dag() + a.dag() * sm) 
        + epsilon * (a.dag() * np.exp(1j * phi_d) + a * np.exp(-1j * phi_d))
    )
    c_ops = [np.sqrt(kappa) * a]

    # -- Numerical Steady State Check --
    rho_ss = steadystate(H, c_ops)
    n_ss_num = expect(a.dag()*a, rho_ss)
    N_heuristic_exact = n_ss_num + 6 * np.sqrt(n_ss_num)

    # -- Dashboard Printout --
    print("\n" + "▒"*60)
    print(f"{' PARAMETERS ':^60}")
    print("▒"*60)
    print(f"  Delta (wa-wd) : {delta/(2*np.pi):<8.4f} |  g : {g/(2*np.pi):.4f}")
    print(f"  Kappa         : {kappa/(2*np.pi):<8.4f} |  Epsilon : {epsilon/(2*np.pi):.4f}")
    print(f"\n--- DERIVED QUANTITIES ---")
    print(f"  Dispersive Chi: {chi/(2*np.pi):.4f} (2π)   |  g / Delta  : {g/delta:.4f}")
    print(f"  Qubit Phase   : {phi_qb:.4f} rad")

    print(f"\n--- STABILITY & TRUNCATION CHECKS ---")
    print(f"  Analytic n_ss : {n_steady_analytical:.4f}       |  n_crit     : {n_crit:.4f}")
    print(f"  Numeric  n_ss : {n_ss_num:.4f}       |  Ratio      : {n_ss_num/n_crit:.4f}")
    print("-" * 60)
    print(f"  Heuristic N   : {N_heuristic_exact:.2f}  (<n> + 6√<n>)")
    print(f"  Chosen N      : {N}")
    
    if N < N_heuristic_exact:
        print(f"  [WARNING] N={N} is below heuristic {N_heuristic_exact:.2f}. Increase N!")
    else:
        print(f"  [OK] Hilbert space size is safe.")
        
    print("\n" + "="*60 + "\n")

    # ============================================================
    # 5. SIMULATION EXECUTION
    # ============================================================
    t_end = 50
    num_points = 50 
    tlist = np.linspace(0, t_end, num_points)

    print(">> Calculating Mean Fields (Signal)...")
    
    alpha_in_scalar = -1j * (epsilon / np.sqrt(kappa)) * np.exp(1j * phi_d)
    a_in_op = alpha_in_scalar * tensor(qeye(N), qeye(2))
    a_out = np.sqrt(kappa) * a - a_in_op
    XX_op = 1j * (a_out.dag() * np.exp(1j * phi_h) - a_out * np.exp(-1j * phi_h))
    
    res0 = mesolve(H, psi0, tlist, c_ops, e_ops=[XX_op], options=solver_opts)
    res1 = mesolve(H, psi1, tlist, c_ops, e_ops=[XX_op], options=solver_opts)
    
    M0 = np.sqrt(kappa) * cumtrapz(np.array(res0.expect[0]), tlist, initial=0)
    M1 = np.sqrt(kappa) * cumtrapz(np.array(res1.expect[0]), tlist, initial=0)
    Signal_ME = np.abs(M0 - M1)

    print(f">> Starting Parallel Noise Calculation (4 Tasks)...")
    tasks = [
        (0, H, psi0, tlist, c_ops, a.dag(), a, solver_opts), 
        (1, H, psi0, tlist, c_ops, a, a,       solver_opts),
        (2, H, psi1, tlist, c_ops, a.dag(), a, solver_opts), 
        (3, H, psi1, tlist, c_ops, a, a,       solver_opts)
    ]
    
    unsorted_results = {}
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(correlation_worker, t) for t in tasks]
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(tasks), unit="task"):
            idx, res = future.result()
            unsorted_results[idx] = res

    results = [unsorted_results[i] for i in range(4)]
    G_N_0, G_raw_0 = results[0], results[1]
    G_N_1, G_raw_1 = results[2], results[3]

    Var0 = finalize_variance(G_N_0, G_raw_0, tlist, kappa, phi_h) - M0**2
    Var1 = finalize_variance(G_N_1, G_raw_1, tlist, kappa, phi_h) - M1**2
    Noise_ME = np.sqrt(np.abs(Var0 + Var1))

    # ============================================================
    # 6. ANALYTICS & SNR
    # ============================================================
    kappa_tau = kappa * tlist
    kt_safe = np.where(kappa_tau == 0, 1e-12, kappa_tau)

    # Analytic Signal
    angle_factor = np.sin(phi_d - phi_h)
    term_decay = 1 - (np.sin(chi * tlist + phi_qb) / np.sin(phi_qb)) * np.exp(-kappa_tau / 2)
    bracket = 1 - (4 / kt_safe) * (np.cos(phi_qb / 2)**2) * term_decay
    Signal_Analytic = np.abs(4 * epsilon * np.sin(phi_qb) * tlist * angle_factor * bracket)
    
    # Analytic Noise
    Noise_Analytic = np.sqrt(2 * kappa_tau)

    # Calculate SNR (Avoiding division by zero)
    safe_noise_me = np.where(Noise_ME == 0, 1e-12, Noise_ME)
    safe_noise_an = np.where(Noise_Analytic == 0, 1e-12, Noise_Analytic)
    
    SNR_ME = Signal_ME / safe_noise_me
    SNR_Analytic = Signal_Analytic / safe_noise_an
    SNR_ratio = np.max(SNR_Analytic)/np.max(SNR_ME)

    # ============================================================
    # 7. PLOTTING (2x2 Grid)
    # ============================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_sig = axes[0, 0]
    ax_noi = axes[0, 1]
    ax_snr = axes[1, 0]
    ax_par = axes[1, 1]

    # 1. Signal (Top Left)
    ax_sig.plot(kappa_tau, Signal_ME, 'b-', label=r'\textbf{ME (Output)}', alpha=0.8)
    ax_sig.plot(kappa_tau, Signal_Analytic, 'r--', label=r'\textbf{Analytical}')
    ax_sig.set_title(r'\textbf{Signal Contrast}', pad=10)
    ax_sig.set_ylabel(r'$|M_1 - M_0|$')
    ax_sig.set_xlabel(r'$\kappa \tau$')
    ax_sig.legend(frameon=True)

    # 2. Noise (Top Right)
    ax_noi.plot(kappa_tau, Noise_ME, 'g-', label=r'\textbf{ME Noise}', alpha=0.8)
    ax_noi.plot(kappa_tau, Noise_Analytic, 'k--', label=r'\textbf{Shot Noise} $\sqrt{2\kappa\tau}$')
    ax_noi.set_title(r'\textbf{Integrated Noise}', pad=10)
    ax_noi.set_ylabel(r'$\Delta M_{\mathrm{total}}$')
    ax_noi.set_xlabel(r'$\kappa \tau$')
    ax_noi.legend(frameon=True)

    # 3. SNR (Bottom Left) - Comparison
    ax_snr.plot(kappa_tau, SNR_ME, 'purple', label=r'\textbf{ME SNR}', alpha=0.8)
    ax_snr.plot(kappa_tau, SNR_Analytic, 'k--', label=r'\textbf{Analytical SNR}')
    ax_snr.set_title(r'\textbf{Signal-to-Noise Ratio}', pad=10)
    ax_snr.set_ylabel(r'$\mathrm{SNR}(\tau)$')
    ax_snr.set_xlabel(r'$\kappa \tau$')
    ax_snr.legend(frameon=True)

    # 4. Parameter Table (Bottom Right)
    ax_par.axis('off') # Hide axis lines
    
    # Create the text block
    param_str = (
        r"\textbf{Simulation Parameters}" + "\n\n" +
        rf"$\Delta_{{ad}} / 2\pi = {delta/(2*np.pi):.1f}$" + "\n" +
        rf"$g / 2\pi = {g/(2*np.pi):.2f}$" + "\n" +
        rf"$\kappa / 2\pi = {kappa/(2*np.pi):.2f}$" + "\n" +
        rf"$\epsilon / 2\pi = {epsilon/(2*np.pi):.2f}$" + "\n" +
        rf"$\chi / 2\pi = {chi/(2*np.pi):.4f}$" + "\n" +
        rf"$\phi_{{qb}} = {phi_qb:.3f}$ rad" + "\n\n" +
        r"\textbf{Steady State Check}" + "\n" +
        rf"$n_{{ss}} (\mathrm{{Num}}) = {n_ss_num:.3f}$" + "\n" +
        rf"$N_{{\mathrm{{Hilbert}}}} = {N}$" + "\n\n" +
        r"\textbf{Final Result}" + "\n" +
        rf"$\mathrm{{max SNR Ration}} = {SNR_ratio:.3f}$"
    )
    
    # Place text in center of the quadrant
    ax_par.text(0.5, 0.5, param_str, 
                ha='center', va='center', fontsize=14, 
                bbox=dict(boxstyle="round,pad=1", fc="wheat", alpha=0.3))

    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / "Final_Quantum_Simulation.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    run_final_simulation()
