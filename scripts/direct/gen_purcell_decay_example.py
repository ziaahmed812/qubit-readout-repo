#!/usr/bin/env python3
"""
Generate purcell_decay_example.pdf — Time-domain Purcell decay at Δ = 4g.
Based on user-provided code. Key: sm = tensor(qeye(N), destroy(2)), proj_e = sm.dag()*sm.
"""
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from qutip import *
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"

mpl.rcParams.update({
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

# ================================================================
# PARAMETERS (from solve_purcell.py / user code)
# ================================================================
g     = 0.05 * 2 * np.pi
kappa = 0.05 * 2 * np.pi
N     = 30

# ================================================================
# OPERATORS
# ================================================================
a  = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))
sp = sm.dag()
sz = tensor(qeye(N), sigmaz())
proj_e = sm.dag() * sm

c_ops = [np.sqrt(kappa) * a]

# ================================================================
# Solve at Δ = 4g
# ================================================================
Delta_example = 4.0 * g
gamma_p_ex = (g**2 * kappa) / (Delta_example**2 + (kappa / 2.0)**2)
print(f"Δ = {Delta_example:.4f}, γ_P = {gamma_p_ex:.6f}")

H_ex = (Delta_example / 2.0) * sz + g * (a * sp + a.dag() * sm)
psi0 = tensor(basis(N, 0), basis(2, 1))
tmax_ex = 8.0 / gamma_p_ex
tlist_ex = np.linspace(0.0, tmax_ex, 1000)

print("Solving master equation...")
res_ex = mesolve(H_ex, psi0, tlist_ex, c_ops, [proj_e])
Pe_ex = np.array(res_ex.expect[0], dtype=float)

# Analytic Purcell decay
Pe_purcell = np.exp(-gamma_p_ex * tlist_ex)

# Dimensionless x-axis: κ t
kt = kappa * tlist_ex

# ================================================================
# PLOT (1.5× fonts, no title)
# ================================================================
plt.figure(figsize=(6.0, 5.0))

plt.plot(kt, Pe_ex, color='black', linewidth=2, label=r"Numerical $P_e(t)$")
plt.plot(kt, Pe_purcell, color='green', linestyle='--', linewidth=2.0,
         label=r"Analytic $e^{-\Gamma_{\rm P} t}$")

plt.xlabel(r"$\kappa t$", fontsize=22)
plt.ylabel(r"$P_e(t)$", fontsize=22)
plt.grid(True, linestyle=':')
plt.legend(loc="best", fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "purcell_decay_example.pdf", dpi=300, bbox_inches="tight")
plt.close()
print("Saved purcell_decay_example.pdf")
