#!/usr/bin/env python3
"""
Generate the Figure 3.6a max-SNR heatmap from externally regenerated (g, Delta) sweep outputs.
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"
RESULTS_DIR = REPO_ROOT / "external_results" / "g_delta_grid" / "results"

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

KAPPA = 0.05 * 2 * np.pi

def load_results(results_dir):
    results_path = Path(results_dir)
    npz_files = sorted(results_path.glob("g_delta_*.npz"))
    print(f"Found {len(npz_files)} result files")
    g_values, delta_values, max_snr_values, is_monotonic = [], [], [], []
    for npz_file in npz_files:
        try:
            data = np.load(npz_file)
            g = float(data['g'])
            delta = float(data['delta'])
            signal = data['Signal_ME']
            noise = data['Noise_ME']
            noise_safe = np.where(noise > 0, noise, np.inf)
            snr = signal / noise_safe
            max_idx = np.nanargmax(snr)
            max_snr = snr[max_idx]
            mono = max_idx >= len(snr) - int(0.01 * len(snr))
            g_values.append(g); delta_values.append(delta)
            max_snr_values.append(max_snr); is_monotonic.append(mono)
        except Exception as e:
            print(f"Error loading {npz_file.name}: {e}")
    return np.array(g_values), np.array(delta_values), np.array(max_snr_values), np.array(is_monotonic)

def create_heatmap(g_vals, delta_vals, max_snr_vals, is_monotonic, output_dir):
    unique_g = np.sort(np.unique(g_vals))
    unique_delta = np.sort(np.unique(delta_vals))
    snr_grid = np.full((len(unique_g), len(unique_delta)), np.nan)
    g_to_idx = {g: i for i, g in enumerate(unique_g)}
    delta_to_idx = {d: i for i, d in enumerate(unique_delta)}
    for g, delta, snr, mono in zip(g_vals, delta_vals, max_snr_vals, is_monotonic):
        i, j = g_to_idx.get(g), delta_to_idx.get(delta)
        if i is not None and j is not None:
            snr_grid[i, j] = np.nan if mono else snr
    delta_2pi = unique_delta / (2 * np.pi)
    g_2pi = unique_g / (2 * np.pi)

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(snr_grid, aspect='auto', origin='lower',
                   extent=[delta_2pi.min(), delta_2pi.max(), g_2pi.min(), g_2pi.max()],
                   cmap='viridis', interpolation='nearest')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r'$\mathrm{max}(\mathrm{SNR})$', fontsize=22)
    cbar.ax.tick_params(labelsize=18)

    # Dispersive regime boundary
    delta_pos = np.linspace(0.01, delta_2pi.max(), 200)
    g_boundary_pos = np.sqrt(KAPPA / 2 * delta_pos * 2 * np.pi) / (2 * np.pi)
    delta_neg = np.linspace(delta_2pi.min(), -0.01, 200)
    g_boundary_neg = np.sqrt(KAPPA / 2 * np.abs(delta_neg * 2 * np.pi)) / (2 * np.pi)
    mask_pos = (g_boundary_pos >= g_2pi.min()) & (g_boundary_pos <= g_2pi.max())
    ax.plot(delta_pos[mask_pos], g_boundary_pos[mask_pos], 'k--', linewidth=2,
            label=r'$g^2/|\Delta| = \kappa/2$')
    mask_neg = (g_boundary_neg >= g_2pi.min()) & (g_boundary_neg <= g_2pi.max())
    ax.plot(delta_neg[mask_neg], g_boundary_neg[mask_neg], 'k--', linewidth=2)

    # Max SNR contours
    delta_zero_idx = np.searchsorted(unique_delta, 0)
    max_g_for_each_delta, delta_for_max_g = [], []
    for j, delta in enumerate(unique_delta):
        col = snr_grid[:, j]
        if not np.all(np.isnan(col)):
            max_g_for_each_delta.append(unique_g[np.nanargmax(col)] / (2 * np.pi))
            delta_for_max_g.append(delta / (2 * np.pi))

    max_delta_pos_for_each_g, g_for_max_delta_pos = [], []
    for i, g in enumerate(unique_g):
        row_pos = snr_grid[i, delta_zero_idx:]
        if not np.all(np.isnan(row_pos)) and len(row_pos) > 0:
            max_delta_pos_for_each_g.append(unique_delta[delta_zero_idx + np.nanargmax(row_pos)] / (2 * np.pi))
            g_for_max_delta_pos.append(g / (2 * np.pi))

    max_delta_neg_for_each_g, g_for_max_delta_neg = [], []
    for i, g in enumerate(unique_g):
        row_neg = snr_grid[i, :delta_zero_idx]
        if not np.all(np.isnan(row_neg)) and len(row_neg) > 0:
            max_delta_neg_for_each_g.append(unique_delta[np.nanargmax(row_neg)] / (2 * np.pi))
            g_for_max_delta_neg.append(g / (2 * np.pi))

    ax.plot(delta_for_max_g, max_g_for_each_delta, '-', color='#d62728', linewidth=2.5,
            label=r'max SNR along $g$')
    ax.plot(max_delta_pos_for_each_g, g_for_max_delta_pos, '--', color='#d62728', linewidth=2.5,
            label=r'max SNR along $|\Delta|$')
    ax.plot(max_delta_neg_for_each_g, g_for_max_delta_neg, '--', color='#d62728', linewidth=2.5)

    ax.set_xlabel(r'$\Delta/(2\pi)$', fontsize=22)
    ax.set_ylabel(r'$g/(2\pi)$', fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=18)
    margin = 0.3
    ax.set_xlim(delta_2pi.min() + margin, delta_2pi.max() - margin)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=18)
    plt.tight_layout()
    plt.savefig(Path(output_dir) / "max_snr_heatmap.pdf", format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved max_snr_heatmap.pdf")

def main():
    if not RESULTS_DIR.exists():
        raise FileNotFoundError(f"Expected externally regenerated sweep outputs in {RESULTS_DIR}. See external_results/README.md.")
    g_vals, delta_vals, max_snr_vals, is_monotonic = load_results(RESULTS_DIR)
    if len(g_vals) == 0:
        print("No results found!"); return
    create_heatmap(g_vals, delta_vals, max_snr_vals, is_monotonic, OUTPUT_DIR)

if __name__ == "__main__":
    main()
