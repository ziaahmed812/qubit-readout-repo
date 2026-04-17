#!/usr/bin/env python3
"""
Regenerate cavity_shift.pdf and stark_shift.pdf with 1.5x larger fonts, no titles.
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"

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

# --- Parameters ---
chi = -1.0
omega_c = 0.0
omega_q = 0.0
kappa = 0.5
gamma = 0.2

w_cav = np.linspace(-4, 4, 1000)
w_qub = np.linspace(-8, 2, 1000)

def lorentzian(w, w0, width):
    return (width**2 / 4) / ((w - w0)**2 + width**2 / 4)

# --- Plot 1: Cavity Frequency Shift ---
freq_g = omega_c - chi
freq_e = omega_c + chi

L_g = lorentzian(w_cav, freq_g, kappa)
L_e = lorentzian(w_cav, freq_e, kappa)

fig1, ax1 = plt.subplots(figsize=(7, 4.5))
ax1.plot(w_cav, L_g, label=r'Qubit in $|g\rangle$', color='#1f77b4', lw=2.5)
ax1.plot(w_cav, L_e, label=r'Qubit in $|e\rangle$', color='#d62728', linestyle='--', lw=2.5)
ax1.fill_between(w_cav, L_g, alpha=0.15, color='#1f77b4')
ax1.fill_between(w_cav, L_e, alpha=0.15, color='#d62728')

ax1.set_xlabel(r'Frequency ($\omega$)', fontsize=22)
ax1.set_ylabel(r'Transmission (a.u.)', fontsize=22)
ax1.set_yticks([])
ax1.set_ylim(0, 1.35)

ax1.set_xticks([freq_e, omega_c, freq_g])
ax1.set_xticklabels([r'$\omega_c + \chi$', r'$\omega_c$', r'$\omega_c - \chi$'], fontsize=21)

peak_height = 1.0
arrow_y = 1.15
ax1.annotate('', xy=(freq_e - 0.05, arrow_y), xytext=(freq_g + 0.05, arrow_y),
             arrowprops=dict(arrowstyle='<|-|>', color='black', lw=1.2, shrinkA=0, shrinkB=0, mutation_scale=15))
ax1.text((freq_e + freq_g)/2, arrow_y + 0.08, r'$2\chi$', ha='center', va='bottom', fontsize=21)

ax1.axvline(freq_g, ymax=0.74, color='#1f77b4', linestyle=':', alpha=0.7, lw=1.2)
ax1.axvline(freq_e, ymax=0.74, color='#d62728', linestyle=':', alpha=0.7, lw=1.2)

ax1.legend(loc='upper right', frameon=False, fontsize=18)
ax1.grid(True, alpha=0.3)
ax1.tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'cavity_shift.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved cavity_shift.pdf")

# --- Plot 2: ac Stark Shift ---
ns = [0, 1, 2, 3]
colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#d62728']

fig2, ax2 = plt.subplots(figsize=(7, 4.5))

peak_positions = []
for i, n in enumerate(ns):
    shift = 2 * chi * n
    peak_pos = omega_q + shift
    peak_positions.append(peak_pos)
    L_q = lorentzian(w_qub, peak_pos, gamma)
    ax2.plot(w_qub, L_q, color=colors[i], lw=2.5)
    ax2.fill_between(w_qub, L_q, alpha=0.2, color=colors[i])
    ax2.text(peak_pos, 1.08, rf'$n={n}$', ha='center', va='bottom', fontsize=21,
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.8))
    ax2.axvline(peak_pos, ymax=0.74, color=colors[i], linestyle=':', alpha=0.6, lw=1.2)

ax2.set_xlabel(r'Frequency ($\omega$)', fontsize=22)
ax2.set_ylabel(r'Absorption (a.u.)', fontsize=22)
ax2.set_yticks([])
ax2.set_ylim(0, 1.35)
ax2.set_xlim(-8, 2)

ax2.set_xticks(peak_positions)
ax2.set_xticklabels([r'$\omega_q$', r'$\omega_q{+}2\chi$', r'$\omega_q{+}4\chi$', r'$\omega_q{+}6\chi$'], fontsize=21)

arrow_y = 0.25
ax2.annotate('', xy=(peak_positions[-1] - 0.8, arrow_y), xytext=(peak_positions[0] + 0.5, arrow_y),
             arrowprops=dict(arrowstyle='-|> ', color='black', lw=1.2, shrinkA=0, shrinkB=0, mutation_scale=15))
ax2.text((peak_positions[0] + peak_positions[-1])/2, arrow_y + 0.04, 
         r'Increasing $\bar{n}$', ha='center', va='bottom', fontsize=18, style='italic')

ax2.grid(True, alpha=0.3)
ax2.tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'stark_shift.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved stark_shift.pdf")
