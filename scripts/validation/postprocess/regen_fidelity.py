#!/usr/bin/env python3
"""
Regenerate fidelity_JC_BEYOND_dispersive_g0.20.pdf with 1.5x larger fonts, no title.
"""
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import erf
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"
DATA_FILE = REPO_ROOT / "data" / "selected" / "validation" / "JC_H_BEYOND-dispersive-with-gamma_g0.20-data.txt"

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

def fidelity(signal, noise):
    with np.errstate(divide='ignore', invalid='ignore'):
        snr = np.where(noise > 0, signal / noise, 0.0)
    return 0.5 * (1.0 + erf(snr / (2.0 * np.sqrt(2.0))))

print(f"Loading: {DATA_FILE}")
data = np.loadtxt(str(DATA_FILE), skiprows=3)

kappa_tau_ana = data[:, 0]; signal_ana = data[:, 1]; noise_ana = data[:, 2]
kappa_tau_me = data[:, 4]; signal_me = data[:, 5]; noise_me = data[:, 6]
kappa_tau_sme = data[:, 8]; signal_sme = data[:, 9]; noise_sme = data[:, 10]

F_ana = fidelity(signal_ana, noise_ana); F_me = fidelity(signal_me, noise_me)
F_sme = fidelity(signal_sme, noise_sme)
valid_sme = ~np.isnan(signal_sme) & ~np.isnan(noise_sme)

max_F_me = np.max(F_me); idx_max_me = np.argmax(F_me)

c_ana = '#1f77b4'; c_me = '#ff7f0e'; c_sme = '#d62728'

fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
ax.plot(kappa_tau_ana, F_ana, '-', color=c_ana, lw=2.5, label='Dispersive (no decay)')
ax.plot(kappa_tau_me, F_me, '--', color=c_me, lw=2.0, label=r'Driven JC (ME, with $\gamma$)')
ax.plot(kappa_tau_me[idx_max_me], max_F_me, 'ko', ms=10, zorder=5)
ax.annotate(rf"$F_{{\max}} = {max_F_me:.4f}$" + "\n" + rf"$\kappa\tau = {kappa_tau_me[idx_max_me]:.1f}$",
            xy=(kappa_tau_me[idx_max_me], max_F_me),
            xytext=(kappa_tau_me[idx_max_me] + 25.0, max_F_me - 0.15),
            fontsize=18,
            arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.4', alpha=0.95))
ax.axhline(0.5, color='#6b7280', linestyle='--', lw=1.8, alpha=0.9, label='Random guess', zorder=1)
ax.set_xlabel(r'$\kappa \tau$', fontsize=22)
ax.set_ylabel(r'Fidelity $F$', fontsize=22)
ax.set_xlim(0, kappa_tau_me[-1]); ax.set_ylim(0.48, 1.02)
ax.grid(True, linestyle=':', alpha=0.5)
ax.legend(loc='lower right', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=18)

plt.savefig(OUTPUT_DIR / "fidelity_JC_BEYOND_dispersive_g0.20.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved fidelity_JC_BEYOND_dispersive_g0.20.png")
