#!/usr/bin/env python3
"""
Regenerate JC_H_dispersive_g0.05_SNR.pdf with 1.5x larger fonts, no title.
"""
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, cumulative_trapezoid
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"
DATA_FILE = REPO_ROOT / "data" / "selected" / "validation" / "JC_H_dispersive_g0.05-data.txt"

# Global LaTeX Style — 1.5× font sizes
plt.rcParams.update({
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
})

# Load ME/SME data
data = np.loadtxt(DATA_FILE, skiprows=3)

kappa_tau_me  = data[:, 4]
snr_me        = data[:, 7]
kappa_tau_sme = data[:, 8]
snr_sme       = data[:, 11]
valid_sme     = ~np.isnan(snr_sme)

# ODE-based SNR with γ_κ
epsilon    = 2 * np.pi * 0.05
kappa      = 2 * np.pi * 0.05
g          = 2 * np.pi * 0.05
Delta      = 2 * np.pi * 4.0

chi        = g**2 / Delta
gamma_k    = (g / Delta)**2 * kappa
sqrt_kappa = np.sqrt(kappa)

phi_d = np.pi / 2.0
drive = -1j * epsilon * np.exp(1j * phi_d)
drive_re = np.real(drive)
drive_im = np.imag(drive)
alpha_in = -1j * (epsilon / sqrt_kappa) * np.exp(1j * phi_d)

def solve_and_get_M(sz0, t_eval):
    def sigma_z(t):
        return (1.0 + sz0) * np.exp(-gamma_k * t) - 1.0
    def ode(t, y):
        x, yv = y
        sz = sigma_z(t)
        return [chi*sz*yv - 0.5*kappa*x + drive_re,
               -chi*sz*x  - 0.5*kappa*yv + drive_im]
    sol = solve_ivp(ode, (t_eval[0], t_eval[-1]), [0.0, 0.0],
                    t_eval=t_eval, method='RK45', rtol=1e-10, atol=1e-12)
    a_out_im = sqrt_kappa * sol.y[1] - np.imag(alpha_in)
    XX = 2.0 * a_out_im
    return sqrt_kappa * cumulative_trapezoid(XX, sol.t, initial=0.0)

t_end = 400
t = np.linspace(0, t_end, 20000)
kappa_tau_ode = kappa * t

M_p = solve_and_get_M(+1.0, t)
M_m = solve_and_get_M(-1.0, t)
S_ode = np.abs(M_p - M_m)
N_ode = np.sqrt(2.0 * kappa * t)
N_safe = np.where(N_ode == 0, 1e-12, N_ode)
snr_ode = S_ode / N_safe

c_ode = '#1f77b4'
c_me  = '#ff7f0e'
c_sme = '#d62728'

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(kappa_tau_ode, snr_ode, '-', color=c_ode, lw=4,
        label='dispersive Hamiltonian')
ax.plot(kappa_tau_me, snr_me, '--', color=c_me, lw=3,
        label='driven JC model (ME)')
ax.plot(kappa_tau_sme[valid_sme], snr_sme[valid_sme], '-', color=c_sme, lw=3, alpha=0.8,
        label='driven JC model (SME)')

if np.any(valid_sme):
    ax.scatter(kappa_tau_sme[valid_sme][-1], snr_sme[valid_sme][-1],
               color=c_sme, s=150, zorder=5, edgecolors='none')

ax.set_xlabel(r'$\kappa \tau$', fontsize=22)
ax.set_ylabel(r'SNR', fontsize=22)
ax.legend(loc='best', fontsize=18)
ax.grid(True, alpha=0.3)
ax.tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()

plt.savefig(OUTPUT_DIR / "JC_H_dispersive_g0.05_SNR.pdf", format='pdf', bbox_inches='tight', dpi=300)
plt.close()
print("Saved JC_H_dispersive_g0.05_SNR.pdf")
