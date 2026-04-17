#!/usr/bin/env python3
"""
Regenerate thesis_wigner_snapshots.pdf with 1.5x larger fonts, no title.
SPECIAL: Keep x-axis and y-ticks only for first column of each row.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qutip import *
from qutip.wigner import wigner
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = REPO_ROOT / "figures" / "chapter3"

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

def run_phase_space_sim(dispersive=False):
    h_type = "dispersive" if dispersive else "full JC"
    print(f"\n>> Running simulation ({h_type})...")
    sys.stdout.flush()
    solver_opts = {
        'rtol': 1e-6, 'atol': 1e-8, 'nsteps': 5000,
        'store_states': True, 'store_final_state': True
    }
    wr = 2*np.pi * 1.0; wa = 2*np.pi * 3.0; wd = 2*np.pi * 1.0
    g = 2*np.pi * 0.2; kappa = 2*np.pi * 0.05; epsilon = 2*np.pi * 0.04; gamma = 2*np.pi * 0.001
    phi_h = 0.0; phi_d = np.pi / 2.0
    N = 14; t_end = 200; num_points = 4000
    delta = wa - wd; chi = g**2 / delta
    tlist = np.linspace(0, t_end, num_points)
    snapshot_times = [0, 10, 28.4, 50, 100, 200]
    snapshot_indices = [np.argmin(np.abs(tlist - t)) for t in snapshot_times]
    actual_times = [tlist[idx] for idx in snapshot_indices]
    print(f"   Snapshots at t = {[f'{t:.1f}' for t in actual_times]}")
    a = tensor(destroy(N), qeye(2)); sz = tensor(qeye(N), sigmaz()); sm = tensor(qeye(N), sigmam())
    psi0 = tensor(basis(N, 0), basis(2, 0)); psi1 = tensor(basis(N, 0), basis(2, 1))
    if dispersive:
        H = chi * a.dag() * a * sz + epsilon * (a.dag() * np.exp(1j * phi_d) + a * np.exp(-1j * phi_d))
        c_ops = [np.sqrt(kappa) * a]
    else:
        H = ((wr - wd) * a.dag() * a + (wa - wd) * sz / 2.0 +
             g * (a * sm.dag() + a.dag() * sm) +
             epsilon * (a.dag() * np.exp(1j * phi_d) + a * np.exp(-1j * phi_d)))
        c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma) * sm]
    print("   Solving for |0⟩..."); res0 = mesolve(H, psi0, tlist, c_ops, e_ops=[a], options=solver_opts)
    print("   Solving for |1⟩..."); res1 = mesolve(H, psi1, tlist, c_ops, e_ops=[a], options=solver_opts)
    states_0 = [res0.states[idx].ptrace(0) for idx in snapshot_indices]
    states_1 = [res1.states[idx].ptrace(0) for idx in snapshot_indices]
    return states_0, states_1, actual_times

def main():
    states0_disp, states1_disp, t_snapshots = run_phase_space_sim(dispersive=True)
    states0_full, states1_full, _ = run_phase_space_sim(dispersive=False)
    n_snapshots = len(t_snapshots)
    print(f"\n>> Generating Thesis Snapshot Plots...")
    a_op = destroy(states0_disp[0].dims[0][0])

    def plot_wigner_snapshot(ax, state0, state1, t_val, show_ylabel=True, show_xlabel=True, show_legend=False):
        c0 = expect(a_op, state0); c1 = expect(a_op, state1)
        mid_x = np.real(c0 + c1) / 2.0; mid_y = np.imag(c0 + c1) / 2.0
        R = 2.2
        grid_center_x = mid_x * np.sqrt(2); grid_center_y = mid_y * np.sqrt(2)
        grid_span = R * np.sqrt(2)
        w_xvec = np.linspace(grid_center_x - grid_span, grid_center_x + grid_span, 200)
        w_yvec = np.linspace(grid_center_y - grid_span, grid_center_y + grid_span, 200)
        W0 = wigner(state0, w_xvec, w_yvec); W1 = wigner(state1, w_xvec, w_yvec)
        scale = 1.0 / np.sqrt(2)
        x_cav = w_xvec * scale; y_cav = w_yvec * scale
        ax.contourf(x_cav, y_cav, W0, levels=15, cmap='Blues', alpha=0.6)
        ax.contourf(x_cav, y_cav, W1, levels=15, cmap='Reds', alpha=0.6)
        contour_frac = 0.1
        ax.contour(x_cav, y_cav, W0, levels=[contour_frac * np.max(W0)], colors='blue', linewidths=1.0, zorder=3)
        ax.contour(x_cav, y_cav, W1, levels=[contour_frac * np.max(W1)], colors='red', linewidths=1.0, zorder=3)
        ax.axhline(y=0, color='green', linestyle='--', linewidth=1.0, alpha=0.5, zorder=2)
        ax.axvline(x=0, color='green', linestyle='--', linewidth=1.0, alpha=0.5, zorder=2)
        ax.plot(np.real(c0), np.imag(c0), 'o', color='blue', ms=5, label=r'$|g\rangle$', zorder=5)
        ax.plot(np.real(c1), np.imag(c1), 'o', color='red', ms=5, label=r'$|e\rangle$', zorder=5)
        ax.plot([np.real(c0), np.real(c1)], [np.imag(c0), np.imag(c1)], 'k-', lw=1.0, zorder=4)
        ax.set_xlim(mid_x - R, mid_x + R); ax.set_ylim(mid_y - R, mid_y + R)
        ax.set_aspect('equal'); ax.set_facecolor('white')
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.text(0.05, 0.95, f'$t = {t_val:.1f}$', transform=ax.transAxes,
                fontsize=18, verticalalignment='top',
                bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.3', alpha=0.9))
        # SPECIAL: Only show y-label and y-ticks on first column
        if show_ylabel:
            ax.set_ylabel(r'$\mathrm{Im}[\langle a \rangle]$', fontsize=22)
        else:
            ax.set_ylabel('')
            ax.set_yticklabels([])
        # SPECIAL: Only show x-label and x-ticks on first column of each row
        if show_xlabel:
            ax.set_xlabel(r'$\mathrm{Re}[\langle a \rangle]$', fontsize=22)
        else:
            ax.set_xlabel('')
            ax.set_xticklabels([])
        if show_legend:
            ax.legend(loc='upper right', fontsize=18)

    fig, axes = plt.subplots(2, n_snapshots, figsize=(16, 6), constrained_layout=True)
    # Top row: Dispersive — x-axis labels only on first column
    for i, (s0, s1, t) in enumerate(zip(states0_disp, states1_disp, t_snapshots)):
        plot_wigner_snapshot(axes[0, i], s0, s1, t,
                            show_ylabel=(i==0), show_xlabel=False, show_legend=(i==0))
        if i > 0:
            axes[0, i].set_xticklabels([])
    # Bottom row: Full JC — x-axis and y-ticks only on first column
    for i, (s0, s1, t) in enumerate(zip(states0_full, states1_full, t_snapshots)):
        plot_wigner_snapshot(axes[1, i], s0, s1, t,
                            show_ylabel=(i==0), show_xlabel=(i==0), show_legend=False)

    plt.savefig(OUTPUT_DIR / "thesis_wigner_snapshots.pdf", dpi=300)
    plt.close()
    print("Saved thesis_wigner_snapshots.pdf")

if __name__ == "__main__":
    main()
