#!/usr/bin/env python3
"""
Regenerate dispersive_H_g0.20_SNR.pdf with 1.5x larger fonts, no title.
"""
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"
DATA_FILE = REPO_ROOT / "data" / "selected" / "validation" / "dispersive_H_g0.20-data.txt"

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

# Load data
data = np.loadtxt(DATA_FILE, skiprows=3)

kappa_tau_ana = data[:, 0]
snr_ana = data[:, 3]
kappa_tau_me = data[:, 4]
snr_me = data[:, 7]
kappa_tau_sme = data[:, 8]
snr_sme = data[:, 11]

valid_sme = ~np.isnan(snr_sme)

c_ana = '#1f77b4'
c_me = '#ff7f0e'
c_sme = '#d62728'

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(kappa_tau_ana, snr_ana, '-', color=c_ana, lw=4, label='Analytical')
ax.plot(kappa_tau_me, snr_me, '--', color=c_me, lw=3, label='ME')
ax.plot(kappa_tau_sme[valid_sme], snr_sme[valid_sme], '-', color=c_sme, lw=3, alpha=0.8, label='SME')

if np.any(valid_sme):
    ax.scatter(kappa_tau_sme[valid_sme][-1], snr_sme[valid_sme][-1], 
               color=c_sme, s=150, zorder=5, edgecolors='none')

ax.set_xlabel(r'$\kappa \tau$', fontsize=22)
ax.set_ylabel(r'SNR', fontsize=22)
ax.legend(loc='best', fontsize=18)
ax.grid(True, alpha=0.3)

ax.tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "dispersive_H_g0.20_SNR.pdf", format='pdf', bbox_inches='tight', dpi=300)
plt.close()
print("Saved dispersive_H_g0.20_SNR.pdf")
