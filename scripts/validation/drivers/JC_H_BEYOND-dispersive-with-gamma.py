import os

# ============================================================
# 0. PARALLELIZATION FLAGS (Must be set BEFORE imports)
# ============================================================
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
# Note: QuTiP's "parallel" map uses joblib and will use available CPU cores.
# To limit to 4 cores, you may need to use system tools (taskset, numactl) or
# modify joblib's default n_jobs behavior programmatically if needed.

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from qutip import *
import scipy.integrate as integrate
import concurrent.futures
from tqdm import tqdm
import time
import sys
import argparse

class Tee(object):
    """
    Duplicate stdout to a log file.
    """
    def __init__(self, name, mode='w'):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        self.stderr = sys.stderr
        sys.stdout = self
        sys.stderr = self
        
    def __del__(self):
        # Handle interpreter shutdown where sys might be None
        try:
            if sys is not None:
                sys.stdout = self.stdout
                sys.stderr = self.stderr
        except (AttributeError, NameError):
            pass
            
        if hasattr(self, 'file') and not self.file.closed:
            self.file.close()
        
    def write(self, data):
        self.file.write(data)
        self.file.flush() # Ensure real-time logging
        self.stdout.write(data)
        self.stdout.flush()
        
    def flush(self):
        self.file.flush()
        self.stdout.flush()

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
# 2. HELPER FUNCTIONS (For Standard ME Noise Calculation)
# ============================================================

def finalize_variance(G_N, G_raw, tlist, kappa, phi_h):
    """
    Vectorized calculation of variance. 
    Replaces Python loops with NumPy array operations.
    """
    # 1. Symmetrize G_raw (Broadcasting) - TIME-ORDERED product uses lower triangle
    lower = np.tril(G_raw)
    G_A = lower + lower.T - np.diag(np.diag(lower))
    
    # 2. Calculate Integrand Matrix
    phase_factor = np.exp(-2j * phi_h)
    integrand = np.real(G_N - phase_factor * G_A)
    
    # 3. 2D Integration using Cumulative Trapezoid
    # We integrate along axis 0, then axis 1 to get the double integral up to time t
    # resulting in a cumulative surface integral along the diagonal.
    
    # Int(0->t) of Integrand(t', t'') dt'
    int_1 = cumtrapz(integrand, tlist, axis=0, initial=0)
    
    # Int(0->t) of [Int(0->t) ...] dt'' 
    # We only care about the diagonal where t' = t'' = t
    double_integral = np.diagonal(cumtrapz(int_1, tlist, axis=1, initial=0))
    
    # 4. Final Formula
    var_list = (kappa * tlist) + (2 * kappa**2 * double_integral)
    
    return var_list

def correlation_worker(args):
    """
    Worker to calculate 2-time correlations.
    Returns: (task_index, mapped_correlation_matrix)
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


# ============================================================
# 3. MAIN SIMULATION
# ============================================================

def run_comparison_simulation(g_val):
    # Solver options - TIGHTENED TOLERANCE
    solver_opts = {
        'rtol': 1e-6, 'atol': 1e-8, 'nsteps': 5000, 
        'store_states': True, 'store_final_state': True
    }

    # -- Parameters --
    wr, wa, wd = 2*np.pi*1.0, 2*np.pi*3.0, 2*np.pi*1.0
    g       = 2*np.pi*g_val  # g from command line
    kappa   = 2*np.pi*0.05 
    epsilon = 2*np.pi*0.05 
    gamma   = 2*np.pi*0.001 # Added for consistency 
    
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
    
    # Jaynes-Cummings Hamiltonian (BEYOND dispersive regime: wa=3.0, with gamma)
    H = (
        (wr - wd) * a.dag() * a +
        (wa - wd) * sz / 2.0 +
        g * (a * sm.dag() + a.dag() * sm) +
        epsilon * (a.dag() * np.exp(1j * phi_d) + a * np.exp(-1j * phi_d))
    )
    c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma) * sm]

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
    # 5. TIME SETTINGS
    # ============================================================
    # Time settings for Master Equation & Analytical
    t_end_me_ana = 400
    num_points_me_ana = 8000
    tlist_me_ana = np.linspace(0, t_end_me_ana, num_points_me_ana)
    
    # Time settings for Stochastic Master Equation
    t_end_sme = 100
    num_points_sme = 4000
    tlist_sme = np.linspace(0, t_end_sme, num_points_sme)
    
    # Input field for output operator (Total field = cavity + input)
    alpha_in_scalar = -1j * (epsilon / np.sqrt(kappa)) * np.exp(1j * phi_d)
    a_in_op = alpha_in_scalar * tensor(qeye(N), qeye(2))
    a_out = np.sqrt(kappa) * a - a_in_op
    XX_op = 1j * (a_out.dag() * np.exp(1j * phi_h) - a_out * np.exp(-1j * phi_h))
    
    # Cavity-only Y operator (for variance subtraction)
    Y_cav_op = 1j * (a.dag() * np.exp(1j * phi_h) - a * np.exp(-1j * phi_h))

    # ============================================================
    # 6. STANDARD MASTER EQUATION (mesolve)
    # ============================================================
    print("=" * 60)
    print(f"{' STANDARD MASTER EQUATION ':^60}")
    print("=" * 60)
    
    print(">> [ME] Calculating Mean Fields (Signal)...")
    start_time = time.time()
    
    # We compute expectation values for BOTH operators AND 'a' (for covariance)
    res0 = mesolve(H, psi0, tlist_me_ana, c_ops, e_ops=[XX_op, Y_cav_op, a], options=solver_opts)
    res1 = mesolve(H, psi1, tlist_me_ana, c_ops, e_ops=[XX_op, Y_cav_op, a], options=solver_opts)

    # Extract Mean Fields <a(t)>
    mean_a_0 = res0.expect[2]
    mean_a_1 = res1.expect[2]
    
    M0_ME = np.sqrt(kappa) * cumtrapz(np.array(res0.expect[0]), tlist_me_ana, initial=0)
    M1_ME = np.sqrt(kappa) * cumtrapz(np.array(res1.expect[0]), tlist_me_ana, initial=0)
    Signal_ME = np.abs(M0_ME - M1_ME)
    
    # Cavity-only means (for Variance subtraction)
    M0_cav = np.sqrt(kappa) * cumtrapz(np.array(res0.expect[1]), tlist_me_ana, initial=0)
    M1_cav = np.sqrt(kappa) * cumtrapz(np.array(res1.expect[1]), tlist_me_ana, initial=0)
    
    print(f"   Signal calculation completed in {time.time()-start_time:.2f}s")

    print(">> [ME] Starting Parallel Noise Calculation (4 Tasks)...")
    start_time = time.time()
    
    tasks = [
        (0, H, psi0, tlist_me_ana, c_ops, a.dag(), a, solver_opts), 
        (1, H, psi0, tlist_me_ana, c_ops, a, a, solver_opts),
        (2, H, psi1, tlist_me_ana, c_ops, a.dag(), a, solver_opts), 
        (3, H, psi1, tlist_me_ana, c_ops, a, a, solver_opts)
    ]
    
    unsorted_results = {}
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(correlation_worker, t) for t in tasks]
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(tasks), unit="task", desc="ME Noise"):
            idx, res = future.result()
            unsorted_results[idx] = res

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
    Var0_ME = finalize_variance(Cov_N_0, Cov_A_0, tlist_me_ana, kappa, phi_h)
    Var1_ME = finalize_variance(Cov_N_1, Cov_A_1, tlist_me_ana, kappa, phi_h)
    
    Noise_ME = np.sqrt(np.abs(Var0_ME + Var1_ME))
    
    print(f"   Noise calculation completed in {time.time()-start_time:.2f}s")
    print(f">> [ME] Done!\n")

    # ============================================================
    # 7. STOCHASTIC MASTER EQUATION (smesolve)
    # ============================================================
    print("=" * 60)
    print(f"{' STOCHASTIC MASTER EQUATION ':^60}")
    print("=" * 60)
    
    ntraj = 1000
    print(f">> [SME] Running {ntraj} trajectories for each qubit state...")
    
    sc_ops_sme = [-1j * np.sqrt(kappa) * a, np.sqrt(gamma) * sm]
    
    # SME solver options: Use QuTiP's built-in parallelization
    # "parallel" uses Python's standard multiprocessing (no extra packages needed)
    # "platen" method for better convergence on multiplicative noise
    sme_opts = {
        "store_measurement": True,
        "store_states": False,
        "store_final_state": False,
        "method": "platen",
        "map": "parallel",  # Standard multiprocessing, no external deps
        "num_cpus": 4,      # Respected by QuTiP
    }
    
    # Run SME for State 0
    print(f"\n>> [SME] Qubit State |0>:")
    start_time = time.time()
    res0_sme = smesolve(
        H, psi0, tlist_sme,
        c_ops=[],
        sc_ops=sc_ops_sme,
        e_ops=[],
        ntraj=ntraj,
        heterodyne=False,
        options=sme_opts
    )
    measurements_0 = np.array(res0_sme.measurement)[:, 0, :]
    print(f"   Completed in {time.time()-start_time:.2f}s")
    
    # Run SME for State 1
    print(f"\n>> [SME] Qubit State |1>:")
    start_time = time.time()
    res1_sme = smesolve(
        H, psi1, tlist_sme,
        c_ops=[],
        sc_ops=sc_ops_sme,
        e_ops=[],
        ntraj=ntraj,
        heterodyne=False,
        options=sme_opts
    )
    measurements_1 = np.array(res1_sme.measurement)[:, 0, :]
    print(f"   Completed in {time.time()-start_time:.2f}s")

    print("\n>> [SME] Processing trajectories...")
    
    offset_val = np.real(1j * (np.conj(alpha_in_scalar) - alpha_in_scalar))
    
    # Safely slice tlist_sme to match measurement output shape
    n_times_sme_actual = measurements_0.shape[1]
    tlist_sme_actual = tlist_sme[:n_times_sme_actual]

    if n_times_sme_actual != len(tlist_sme):
        print(f"   [NOTE] Slicing tlist_sme from {len(tlist_sme)} to {n_times_sme_actual} points.")

    # Process State 0
    measurements_0_corrected = measurements_0 + offset_val
    M_traj_0 = np.sqrt(kappa) * cumtrapz(measurements_0_corrected, tlist_sme_actual, axis=1, initial=0)
    
    # Process State 1
    measurements_1_corrected = measurements_1 + offset_val
    M_traj_1 = np.sqrt(kappa) * cumtrapz(measurements_1_corrected, tlist_sme_actual, axis=1, initial=0)

    mean_M_0_SME = np.mean(M_traj_0, axis=0)
    mean_M_1_SME = np.mean(M_traj_1, axis=0)
    Signal_SME = np.abs(mean_M_1_SME - mean_M_0_SME)
    
    var_M_0_SME = np.var(M_traj_0, axis=0)
    var_M_1_SME = np.var(M_traj_1, axis=0)
    Noise_SME = np.sqrt(var_M_0_SME + var_M_1_SME)
    
    # Time axis for SME
    kappa_tau_sme = kappa * tlist_sme_actual
    
    print(f">> [SME] Done!\n")

    # ============================================================
    # 8. ANALYTICAL RESULTS
    # ============================================================
    print("=" * 60)
    print(f"{' ANALYTICAL CALCULATIONS ':^60}")
    print("=" * 60)
    
    kappa_tau_ana = kappa * tlist_me_ana
    kt_safe = np.where(kappa_tau_ana == 0, 1e-12, kappa_tau_ana)

    # Analytic Signal
    angle_factor = np.sin(phi_d - phi_h)
    term_decay = 1 - (np.sin(chi * tlist_me_ana + phi_qb) / np.sin(phi_qb)) * np.exp(-kappa_tau_ana / 2)
    bracket = 1 - (4 / kt_safe) * (np.cos(phi_qb / 2)**2) * term_decay
    Signal_Analytic = np.abs(4 * epsilon * np.sin(phi_qb) * tlist_me_ana * angle_factor * bracket)
    
    # Analytic Noise
    Noise_Analytic = np.sqrt(2 * kappa_tau_ana)
    
    print(">> [Analytic] Formulas evaluated.")

    # ============================================================
    # 9. SNR CALCULATIONS
    # ============================================================
    # ME and Analytical use tlist_me_ana (4000 points, t_end=200)
    # SME uses tlist_sme_actual (up to 2000 points, t_end=100)
    # Each method has its own time column in the data file
    
    # ME/Analytical time axis
    kappa_tau_me_ana = kappa * tlist_me_ana
    
    # SME time axis  
    n_sme = len(tlist_sme_actual)
    
    # Calculate SNR for ME and Analytical (full length)
    safe_noise_me = np.where(Noise_ME == 0, 1e-12, Noise_ME)
    safe_noise_an = np.where(Noise_Analytic == 0, 1e-12, Noise_Analytic)
    
    SNR_ME = Signal_ME / safe_noise_me
    SNR_Analytic = Signal_Analytic / safe_noise_an
    
    # Calculate SNR for SME (its own length)
    safe_noise_sme = np.where(Noise_SME == 0, 1e-12, Noise_SME)
    SNR_SME = Signal_SME / safe_noise_sme
    
    SNR_ratio_ME = np.max(SNR_Analytic) / np.max(SNR_ME)
    SNR_ratio_SME = np.max(SNR_Analytic) / np.max(SNR_SME)
    
    # ============================================================
    # 9b. SAVE DATA TO FILE
    # ============================================================
    print("\n>> Saving data to file...")
    
    # ME/Analytical have num_points_me_ana points (4000)
    # SME has n_sme points (up to 2000)
    # Pad SME arrays with NaN to match ME/Analytical length
    n_me_ana = len(tlist_me_ana)
    
    def pad_with_nan(arr, target_len):
        """Pad array with NaN to reach target length."""
        if len(arr) >= target_len:
            return arr[:target_len]
        return np.concatenate([arr, np.full(target_len - len(arr), np.nan)])
    
    # Pad SME arrays
    kappa_tau_sme_padded = pad_with_nan(kappa_tau_sme, n_me_ana)
    Signal_SME_padded = pad_with_nan(Signal_SME, n_me_ana)
    Noise_SME_padded = pad_with_nan(Noise_SME, n_me_ana)
    SNR_SME_padded = pad_with_nan(SNR_SME, n_me_ana)
    
    # 12 columns: each method has its own time column
    # Analytical: Kappa*tau, Signal, Noise, SNR
    # ME: Kappa*tau, Signal, Noise, SNR
    # SME: Kappa*tau, Signal, Noise, SNR (padded with NaN beyond t_end_sme)
    
    data_matrix_12 = np.column_stack((
        kappa_tau_me_ana, Signal_Analytic, Noise_Analytic, SNR_Analytic,
        kappa_tau_me_ana, Signal_ME, Noise_ME, SNR_ME,
        kappa_tau_sme_padded, Signal_SME_padded, Noise_SME_padded, SNR_SME_padded
    ))
    
    header_12_cols = (
        f"parameters: delta={delta/(2*np.pi):.4f} g={g/(2*np.pi):.4f} "
        f"kappa={kappa/(2*np.pi):.4f} epsilon={epsilon/(2*np.pi):.4f} "
        f"gamma={gamma/(2*np.pi):.4f} N={N}\n"
        f"time_settings: ME_Ana(t_end={t_end_me_ana},n={num_points_me_ana}) SME(t_end={t_end_sme},n={num_points_sme})\n"
        "analytical-Kappa*tau analytical_signal analytical_Noise analytical_SNR "
        "ME-Kappa*tau ME_signal ME_Noise ME_SNR "
        "SME-Kappa*tau SME_signal SME_Noise SME_SNR"
    )
    
    # Save with filename based on script name and g value
    script_name = os.path.basename(__file__)
    script_base = os.path.splitext(script_name)[0]
    output_filename = f"{script_base}_g{g_val:.2f}-data.txt"
    
    np.savetxt(output_filename, data_matrix_12, header=header_12_cols, comments='', fmt='%1.6e')
    print(f"   Data saved successfully to '{output_filename}'.")
    print(f"   Note: SME columns padded with NaN beyond t={t_end_sme} (row {n_sme}).")

    # ============================================================
    # 10. PLOTTING (Skipped)
    # ============================================================
    # print("\n>> Generating plots... (Skipped)")
    
    # fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    # ax_sig = axes[0, 0]
    # ax_noi = axes[0, 1]
    # ax_snr = axes[1, 0]
    # ax_par = axes[1, 1]
    
    # ... plotting code removed ...
    
    # plt.tight_layout()
    # plt.savefig("dispersive_H_comparison.pdf", dpi=300)
    # plt.show()

    # ============================================================
    # 11. FINAL TERMINAL OUTPUT
    # ============================================================
    print("\n" + "="*60)
    print(f"{' FINAL COMPARISON ':^60}")
    print("="*60)
    print(f"  Max SNR (Analytical): {np.max(SNR_Analytic):.4f}")
    print(f"  Max SNR (ME)        : {np.max(SNR_ME):.4f}")
    print(f"  Max SNR (SME)       : {np.max(SNR_SME):.4f}")
    print("-" * 60)
    print(f"  SNR Ratio (Ana/ME) : {SNR_ratio_ME:.4f}")
    print(f"  SNR Ratio (Ana/SME): {SNR_ratio_SME:.4f}")
    print("="*60 + "\n")

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run SNR simulation with specified g value')
    parser.add_argument('--g', type=float, required=True, help='Coupling strength g (without 2*pi factor)')
    parser.add_argument('--output_dir', type=str, default='.', help='Output directory for results')
    args = parser.parse_args()
    
    # Change to output directory if specified
    if args.output_dir != '.':
        os.makedirs(args.output_dir, exist_ok=True)
        os.chdir(args.output_dir)
    
    # Setup automatic logging
    script_name = os.path.basename(__file__)
    script_base = os.path.splitext(script_name)[0]
    log_filename = f"{script_base}_g{args.g:.2f}.log"
    
    # Redirect stdout and stderr to both console and log file
    logger = Tee(log_filename)
    sys.stderr = logger
    
    print(f">> Logging output to: {log_filename}")
    print(f">> Running with g = {args.g}")
    run_comparison_simulation(args.g)
