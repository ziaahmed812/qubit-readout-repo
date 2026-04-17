#!/usr/bin/env python3
"""
Generate the Figure 3.9 gamma-sweep artifact from externally regenerated sweep outputs.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"
RESULTS_DIR = REPO_ROOT / "external_results" / "gamma_sweep" / "results"

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

def main():
    if not RESULTS_DIR.exists():
        raise FileNotFoundError(f"Expected externally regenerated sweep outputs in {RESULTS_DIR}. See external_results/README.md.")
    fig, ax = plt.subplots(figsize=(8, 6))
    files = sorted(RESULTS_DIR.glob("*.npz"))
    print(f"Found {len(files)} result files")
    gamma_values = []
    for f in files:
        data = np.load(f)
        gamma_values.append(float(data['gamma']) / (2 * np.pi))
    cmap = plt.cm.viridis
    norm = Normalize(vmin=min(gamma_values), vmax=max(gamma_values))
    for f in sorted(files):
        data = np.load(f)
        gamma = float(data['gamma']) / (2 * np.pi)
        kappa_tau = data['kappa_tau']; signal = data['Signal_ME']; noise = data['Noise_ME']
        with np.errstate(divide='ignore', invalid='ignore'):
            snr = np.where(noise > 0, signal / noise, 0)
        ax.plot(kappa_tau, snr, color=cmap(norm(gamma)), linewidth=1.5, alpha=0.8)
    ax.set_xlabel(r'$\kappa \tau$', fontsize=22)
    ax.set_ylabel(r'SNR', fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    sm = ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label(r'$\gamma/(2\pi)$', fontsize=22)
    cbar.ax.tick_params(labelsize=18)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "gamma_sweep_snr_vs_kappa_tau.pdf", format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved gamma_sweep_snr_vs_kappa_tau.pdf")

if __name__ == "__main__":
    main()
