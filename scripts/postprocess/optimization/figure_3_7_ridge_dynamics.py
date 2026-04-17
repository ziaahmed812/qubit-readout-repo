#!/usr/bin/env python3
"""
Generate the Figure 3.7 ridge-dynamics panels from externally regenerated (g, Delta) sweep outputs.
Optimized: two-pass loading to avoid keeping the full grid in RAM.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 21,
    "axes.labelsize": 22,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
    "lines.linewidth": 1,
    "axes.linewidth": 1.2,
    "axes.spines.top": True,
    "axes.spines.right": True,
})

RESULTS_DIR = REPO_ROOT / "external_results" / "g_delta_grid" / "results"

def main():
    if not RESULTS_DIR.exists():
        raise FileNotFoundError(f"Expected externally regenerated sweep outputs in {RESULTS_DIR}. See external_results/README.md.")
    npz_files = sorted(RESULTS_DIR.glob("g_delta_*.npz"))
    print(f"Found {len(npz_files)} result files")

    # ── Pass 1: only extract (g, delta, max_snr, filename) ──
    print("Pass 1: scanning max SNR values...")
    records = []
    for i, npz_file in enumerate(npz_files):
        try:
            data = np.load(npz_file)
            g = float(data['g']); delta = float(data['delta'])
            signal = data['Signal_ME']; noise = data['Noise_ME']
            noise_safe = np.where(noise > 0, noise, np.inf)
            max_snr = np.nanmax(signal / noise_safe)
            records.append((g, delta, max_snr, npz_file))
            data.close()
        except Exception as e:
            print(f"Error: {npz_file.name}: {e}")
        if (i + 1) % 500 == 0:
            print(f"  scanned {i+1}/{len(npz_files)}...", flush=True)

    print(f"  scanned {len(records)} files total")

    all_g = sorted(set(r[0] for r in records))
    all_delta = sorted(set(r[1] for r in records))
    print(f"  unique g: {len(all_g)}, unique delta: {len(all_delta)}")

    # Build lookup: (g, delta) -> (max_snr, filepath)
    lookup = {}
    for g, delta, max_snr, fpath in records:
        key = (g, delta)
        if key not in lookup or max_snr > lookup[key][0]:
            lookup[key] = (max_snr, fpath)

    # ── Determine which files to load for each plot ──
    # Plot 1: for each delta, find the g with max SNR
    best_for_delta = {}
    for delta in all_delta:
        best_g, best_snr, best_file = None, -np.inf, None
        for g in all_g:
            if (g, delta) in lookup:
                ms, fp = lookup[(g, delta)]
                if ms > best_snr:
                    best_snr = ms; best_g = g; best_file = fp
        if best_file is not None:
            best_for_delta[delta] = best_file

    # Plot 2: for each g, find the delta with max SNR
    best_for_g = {}
    for g in all_g:
        best_delta, best_snr, best_file = None, -np.inf, None
        for delta in all_delta:
            if (g, delta) in lookup:
                ms, fp = lookup[(g, delta)]
                if ms > best_snr:
                    best_snr = ms; best_delta = delta; best_file = fp
        if best_file is not None:
            best_for_g[g] = (best_delta, best_file)

    files_to_load = set(best_for_delta.values()) | set(v[1] for v in best_for_g.values())
    print(f"Pass 2: loading {len(files_to_load)} selected files for plotting...")

    # ── Pass 2: load only selected files ──
    file_data = {}
    for fp in files_to_load:
        data = np.load(fp)
        kappa_tau = data['kappa_tau']; signal = data['Signal_ME']; noise = data['Noise_ME']
        noise_safe = np.where(noise > 0, noise, np.inf)
        snr = signal / noise_safe
        file_data[fp] = (kappa_tau, snr)
        data.close()

    cmap = plt.cm.viridis

    # ── Plot 1: optimal g curves colored by |delta| ──
    abs_delta_2pi = [abs(d) / (2 * np.pi) for d in all_delta]
    norm_delta = Normalize(vmin=min(abs_delta_2pi), vmax=max(abs_delta_2pi))
    fig, ax = plt.subplots(figsize=(8, 6))
    for delta in all_delta:
        if delta in best_for_delta:
            fp = best_for_delta[delta]
            kappa_tau, snr = file_data[fp]
            ax.plot(kappa_tau, snr, color=cmap(norm_delta(abs(delta) / (2 * np.pi))),
                    alpha=0.8, linewidth=1.2)
    ax.set_xlabel(r'$\kappa\tau$', fontsize=22); ax.set_ylabel(r'SNR', fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=18)
    sm = ScalarMappable(cmap=cmap, norm=norm_delta); sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax); cbar.set_label(r'$|\Delta|/(2\pi)$', fontsize=22)
    cbar.ax.tick_params(labelsize=18)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "snr_curves_along_g.pdf", format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved snr_curves_along_g.pdf")

    # ── Plot 2: optimal delta curves colored by g ──
    g_2pi = [g / (2 * np.pi) for g in all_g]
    norm_g = Normalize(vmin=min(g_2pi), vmax=max(g_2pi))
    fig, ax = plt.subplots(figsize=(8, 6))
    for g in all_g:
        if g in best_for_g:
            _, fp = best_for_g[g]
            kappa_tau, snr = file_data[fp]
            ax.plot(kappa_tau, snr, color=cmap(norm_g(g / (2 * np.pi))),
                    alpha=0.8, linewidth=1.2)
    ax.set_xlabel(r'$\kappa\tau$', fontsize=22); ax.set_ylabel(r'SNR', fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=18)
    sm = ScalarMappable(cmap=cmap, norm=norm_g); sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax); cbar.set_label(r'$g/(2\pi)$', fontsize=22)
    cbar.ax.tick_params(labelsize=18)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "snr_curves_along_delta.pdf", format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved snr_curves_along_delta.pdf")

if __name__ == "__main__":
    main()
