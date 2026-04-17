#!/usr/bin/env python3
"""
Generate purcell_rates_loglog.pdf — Purcell decay rate vs Δ/g on log-log scale.
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

Delta_list = np.linspace(0.1, 15.0, 30) * g

# ================================================================
# OPERATORS
# ================================================================
a  = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))
sp = sm.dag()
sz = tensor(qeye(N), sigmaz())
proj_e = sm.dag() * sm

c_ops = [np.sqrt(kappa) * a]

gamma_num_list = []
gamma_purcell_list = []

# ================================================================
# SWEEP DETUNINGS
# ================================================================
print("Computing Purcell rates...")
for idx, Delta in enumerate(Delta_list):
    gamma_p = (g**2 * kappa) / (Delta**2 + (kappa / 2)**2)
    gamma_purcell_list.append(gamma_p)

    tmax = 8.0 / gamma_p
    tlist = np.linspace(0, tmax, 1200)

    psi0 = tensor(basis(N, 0), basis(2, 1))
    res = mesolve(H=(Delta / 2) * sz + g * (a * sp + a.dag() * sm),
                  rho0=psi0, tlist=tlist, c_ops=c_ops, e_ops=[proj_e])
    Pe = np.array(res.expect[0], dtype=float)

    # Restrict to valid region and use middle 60% for fitting
    mask = Pe > 1e-5
    t_fit = tlist[mask]
    P_fit = Pe[mask]
    n = len(t_fit)
    t_fit = t_fit[n // 5: 4 * n // 5]
    P_fit = P_fit[n // 5: 4 * n // 5]

    logP = np.log(P_fit)
    slope, intercept = np.polyfit(t_fit, logP, 1)
    gamma_num = -slope
    gamma_num_list.append(gamma_num)

    if idx % 5 == 0:
        print(f"  [{idx+1}/{len(Delta_list)}] Δ/g={Delta/g:.1f}, "
              f"γ_P(analytic)={gamma_p:.2e}, γ(numeric)={gamma_num:.2e}")

gamma_num_list = np.array(gamma_num_list)
gamma_purcell_list = np.array(gamma_purcell_list)
print("Done computing rates.")

# ================================================================
# PLOT: log-log with LaTeX style (1.5× fonts, no title)
# ================================================================
plt.figure(figsize=(6.0, 5.0))

x = Delta_list / g
y_analytic = gamma_purcell_list
y_numeric  = gamma_num_list

mask = x >= 2.0
x_m     = x[mask]
y_m     = y_analytic[mask]
y_num_m = y_numeric[mask]

plt.loglog(x_m, y_m, color='black', marker='o', markersize=7,
           linewidth=3, label=r"Analytic Purcell")
plt.loglog(x_m, y_num_m, color='green', linestyle='-', marker='+',
           markersize=10, linewidth=1.8, label=r"Numerical (JC + Lindblad ME)")

# Reference slope −2
mid = len(x_m) // 2
x_ref0, y_ref0 = x_m[mid], y_m[mid]
x_ref = np.array([x_ref0 / 2, x_ref0 * 2])
y_ref = y_ref0 * (x_ref / x_ref0)**(-2)
plt.loglog(x_ref, y_ref, color='red', linestyle='--', linewidth=2.2,
           label=r"Slope $-2$")

plt.xlabel(r"Detuning $\Delta / g$", fontsize=22)
plt.ylabel(r"Decay rate $\Gamma$", fontsize=22)
plt.grid(True, which='both', linestyle=':')
plt.ylim(top=1e-1)
plt.legend(loc="upper right", fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "purcell_rates_loglog.pdf", dpi=300, bbox_inches="tight")
plt.close()
print("Saved purcell_rates_loglog.pdf")
