#!/usr/bin/env python3
"""
Generate the Figure 3.8b (kappa, epsilon) heatmap from externally regenerated grid outputs.

The thesis artifact filename `gamma_kappa_max_snr_heatmap.pdf` is preserved for consistency with
the original thesis source, but the physical content is the Figure 3.8b heatmap in the
(kappa, epsilon) plane.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"
RESULTS_DIR = REPO_ROOT / "external_results" / "kappa_epsilon_grid" / "results"

matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 21,
    "axes.labelsize": 22,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
    "lines.linewidth": 2,
    "axes.linewidth": 1.4,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.spines.left": True,
    "axes.spines.bottom": True,
})

def main():
    if not RESULTS_DIR.exists():
        raise FileNotFoundError(f"Expected externally regenerated sweep outputs in {RESULTS_DIR}. See external_results/README.md.")
    kappa_vals = np.array([0.04, 0.045, 0.05, 0.055, 0.06])
    epsilon_vals = np.array([0.04, 0.045, 0.05, 0.055, 0.06])
    max_snr_grid = np.zeros((len(epsilon_vals), len(kappa_vals)))
    files = list(RESULTS_DIR.glob("*.npz"))
    print(f"Found {len(files)} result files")
    for f in files:
        data = np.load(f)
        kappa = float(data['kappa']) / (2 * np.pi)
        epsilon = float(data['epsilon']) / (2 * np.pi)
        signal = data['Signal_ME']; noise = data['Noise_ME']
        with np.errstate(divide='ignore', invalid='ignore'):
            snr = np.where(noise > 0, signal / noise, 0)
        max_snr = np.nanmax(snr)
        i_kappa = np.argmin(np.abs(kappa_vals - kappa))
        i_epsilon = np.argmin(np.abs(epsilon_vals - epsilon))
        max_snr_grid[i_epsilon, i_kappa] = max_snr

    print(f"Max SNR range: {max_snr_grid.min():.4f} to {max_snr_grid.max():.4f}")

    fig, ax = plt.subplots(figsize=(7, 5.5))
    im = ax.imshow(max_snr_grid, origin='lower', aspect='auto', cmap='viridis',
                   extent=[kappa_vals[0] - 0.0025, kappa_vals[-1] + 0.0025,
                           epsilon_vals[0] - 0.0025, epsilon_vals[-1] + 0.0025])

    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r'max(SNR)', fontsize=22)
    cbar.ax.tick_params(labelsize=18)

    ax.set_xlabel(r'$\kappa/(2\pi)$', fontsize=22)
    ax.set_ylabel(r'$\epsilon/(2\pi)$', fontsize=22)

    ax.set_xticks(kappa_vals)
    ax.set_yticks(epsilon_vals)
    ax.tick_params(axis='both', which='major', labelsize=18)

    plt.tight_layout()
    legacy_output = OUTPUT_DIR / "gamma_kappa_max_snr_heatmap.pdf"
    reviewer_alias = OUTPUT_DIR / "figure_3_8b_kappa_epsilon_heatmap.pdf"
    plt.savefig(legacy_output, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(reviewer_alias, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {legacy_output.name} and {reviewer_alias.name}")

if __name__ == "__main__":
    main()
