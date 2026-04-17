#!/usr/bin/env python3
"""
Generate the Figure 3.10 homodyne-phase sensitivity artifact from externally regenerated sweep outputs.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import glob
from qutip import *
import scipy.integrate as integrate
import sys

if hasattr(integrate, 'cumulative_trapezoid'):
    cumtrapz = integrate.cumulative_trapezoid
else:
    cumtrapz = integrate.cumtrapz

REPO_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"
RESULTS_DIR = REPO_ROOT / "external_results" / "delta_phi" / "results"

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

def load_snr_data(results_dir):
    files = sorted(glob.glob(str(results_dir / "noise_phih_0.0000_phid_*.npz")))
    phi_d_list, max_snr_list = [], []
    for f in files:
        data = np.load(f)
        phi_d = float(data['phi_d'])
        signal = data['Signal_ME']; noise = data['Noise_ME']
        with np.errstate(divide='ignore', invalid='ignore'):
            snr = np.where(noise > 1e-10, signal / noise, 0)
        max_snr_list.append(np.nanmax(snr))
        phi_d_list.append(phi_d)
    return np.array(phi_d_list), np.array(max_snr_list)

def main():
    if not RESULTS_DIR.exists():
        raise FileNotFoundError(f"Expected externally regenerated sweep outputs in {RESULTS_DIR}. See external_results/README.md.")
    print(">> Loading Full JC HPC results...")
    phi_d_jc, max_snr_jc = load_snr_data(RESULTS_DIR)
    sort_idx = np.argsort(phi_d_jc)
    phi_d_jc = phi_d_jc[sort_idx]; max_snr_jc = max_snr_jc[sort_idx]
    delta_phi_pi = phi_d_jc / np.pi
    max_idx = np.argmax(max_snr_jc)
    max_delta_phi_pi = delta_phi_pi[max_idx]
    max_snr_val = max_snr_jc[max_idx]

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    line_color = '#2E86AB'; fill_color = '#A6E1FA'
    ax.fill_between(delta_phi_pi, 0, max_snr_jc, alpha=0.3, color=fill_color)
    ax.plot(delta_phi_pi, max_snr_jc, '-', color=line_color, linewidth=2.5)
    ax.plot(max_delta_phi_pi, max_snr_val, 'ko', ms=12, zorder=5)
    ax.annotate(f'Max = {max_snr_val:.1f}\n$\\Delta\\phi = {max_delta_phi_pi:.2f}\\pi$',
                 xy=(max_delta_phi_pi, max_snr_val),
                 xytext=(max_delta_phi_pi + 0.3, max_snr_val * 0.75),
                 fontsize=18,
                 arrowprops=dict(arrowstyle='->', color='black', lw=2),
                 bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.4', alpha=0.95))
    ax.set_xlabel(r'$\Delta\phi = \phi_d - \phi_h$ (units of $\pi$)', fontsize=22)
    ax.set_ylabel(r'Max SNR', fontsize=22)
    ax.set_xlim(0, 2); ax.set_ylim(0, max_snr_val * 1.15)
    ax.set_xticks([0, 0.5, 1, 1.5, 2])
    ax.set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.savefig(OUTPUT_DIR / "snr_vs_phi_d.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved snr_vs_phi_d.pdf")

if __name__ == "__main__":
    main()
