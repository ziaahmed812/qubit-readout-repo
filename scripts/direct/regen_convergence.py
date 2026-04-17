#!/usr/bin/env python3
"""
Regenerate dynamic_convergence_plot.pdf with 1.5x larger fonts, no title.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qutip import *
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = REPO_ROOT / "figures" / "appendix"

matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 21,
    "axes.labelsize": 22,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
    "lines.linewidth": 2,
    "axes.linewidth": 1.2,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.spines.left": True,
    "axes.spines.bottom": True,
})

wc      = 0.0
wa      = 2.0
g       = 0.2
kappa   = 0.05
gamma   = 0.0
epsilon = 0.05

tlist = np.linspace(0, 50 / kappa, 1500)
tolerance = 1e-5
N_test = 21

def compute_n_and_Ptop(N: int):
    a = destroy(N); sm = sigmam(); I_atom = qeye(2); I_cav = qeye(N)
    H = (wc * tensor(a.dag() * a, I_atom) + 0.5 * wa * tensor(I_cav, sigmaz())
         + g * (tensor(a.dag(), sm) + tensor(a, sm.dag())) + epsilon * tensor(a + a.dag(), I_atom))
    c_ops = [np.sqrt(kappa) * tensor(a, I_atom)]
    if gamma > 0.0:
        c_ops.append(np.sqrt(gamma) * tensor(I_cav, sm))
    psi0 = tensor(basis(N, 0), basis(2, 0))
    n_op = tensor(num(N), I_atom)
    top_ket = basis(N, N - 1)
    P_top_op = tensor(top_ket * top_ket.dag(), I_atom)
    result = mesolve(H, psi0, tlist, c_ops, e_ops=[n_op, P_top_op])
    return np.array(result.expect[0], dtype=float), np.array(result.expect[1], dtype=float)

Delta_r = wc
n_ss = epsilon**2 / ((kappa/2)**2 + Delta_r**2)
N_heur = int(np.ceil(n_ss + 6.0 * np.sqrt(n_ss)))

print(f"Running N={N_test}...")
n_t_N, Ptop_N = compute_n_and_Ptop(N_test)
print(f"Running N={N_test-1}...")
n_t_Nm1, Ptop_Nm1 = compute_n_and_Ptop(N_test - 1)
print(f"Running N={N_heur}...")
n_t_heur, Ptop_heur = compute_n_and_Ptop(N_heur)

delta_N_vs_Nm1 = np.abs(n_t_N - n_t_Nm1)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

mstep = 50
midx = np.arange(0, len(tlist), mstep)

# (a) Photon number dynamics
ax[0].plot(kappa * tlist, n_t_N, color='C0', linewidth=2.0, label=rf"$N={N_test}$")
ax[0].plot((kappa * tlist)[midx], n_t_Nm1[midx], 'x', linestyle='None', color='C1', markersize=8, label=rf"$N={N_test-1}$")
ax[0].plot((kappa * tlist)[midx], n_t_heur[midx], 'o', linestyle='None', color='C2', markersize=5, label=rf"$N={N_heur}$ (heuristic)")
ax[0].text(0.03, 0.95, "(a)", transform=ax[0].transAxes, fontsize=21, fontweight="bold", va="top")
ax[0].set_xlabel(r"Time ($\kappa t$)", fontsize=22)
ax[0].set_ylabel(r"$\langle n(t)\rangle$", fontsize=22)
ax[0].grid(True, linestyle=":", alpha=0.7)
ax[0].legend(loc="lower right", fontsize=14)

# (b) Convergence diagnostics
ax[1].semilogy(kappa * tlist, delta_N_vs_Nm1, color='C0', linewidth=2.0,
               label=rf"$|\Delta n(t)|$: $N={N_test}$ vs $N={N_test-1}$")
Ptop_Nm1_env = np.maximum.accumulate(Ptop_Nm1)
ax[1].semilogy(kappa * tlist, Ptop_Nm1_env, color='C1', linewidth=2.0,
               label=rf"$\max_{{t'\leq t}} P_{{{N_test-2}}}(t')$, $N={N_test-1}$")
Ptop_N_env = np.maximum.accumulate(Ptop_N)
ax[1].semilogy(kappa * tlist, Ptop_N_env, color='C2', linewidth=2.0,
               label=rf"$\max_{{t'\leq t}} P_{{{N_test-1}}}(t')$, $N={N_test}$")
Ptop_heur_env = np.maximum.accumulate(Ptop_heur)
ax[1].semilogy(kappa * tlist, Ptop_heur_env, color='C3', linewidth=2.0,
               label=rf"$\max_{{t'\leq t}} P_{{{N_heur-1}}}(t')$, $N={N_heur}$ (heuristic)")
ax[1].axhline(tolerance, color="gray", linestyle=":", linewidth=1.5,
              label=r"Tolerance = $1.0 \times 10^{-5}$")
ax[1].text(0.03, 0.95, "(b)", transform=ax[1].transAxes, fontsize=21, fontweight="bold", va="top")
ax[1].set_xlabel(r"Time ($\kappa t$)", fontsize=22)
ax[1].set_ylabel(r"Error / Population", fontsize=22)
ax[1].grid(True, which="both", linestyle=":", alpha=0.7)
ax[1].legend(loc="lower right", fontsize=14)

ax[0].tick_params(axis='both', which='major', labelsize=18)
ax[1].tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "dynamic_convergence_plot.pdf", dpi=300, bbox_inches="tight")
plt.close()
print("Saved dynamic_convergence_plot.pdf")
