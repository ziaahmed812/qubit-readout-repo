#!/usr/bin/env python3
"""
Regenerate thesis_wigner_final.pdf with 1.5x larger fonts, no title.
SPECIAL: Change |0⟩ → |g⟩ and |1⟩ → |e⟩.
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
    N = 14; t_max_snr = 28.41; t_end = 35; num_points = 1000
    delta = wa - wd; chi = g**2 / delta
    tlist = np.linspace(0, t_end, num_points)
    t_idx = np.argmin(np.abs(tlist - t_max_snr))
    print(f"   Using t = {tlist[t_idx]:.2f} (index {t_idx})")
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
    print("   Solving for |g⟩..."); res0 = mesolve(H, psi0, tlist, c_ops, e_ops=[a], options=solver_opts)
    print("   Solving for |e⟩..."); res1 = mesolve(H, psi1, tlist, c_ops, e_ops=[a], options=solver_opts)
    state_0 = res0.states[t_idx].ptrace(0); state_1 = res1.states[t_idx].ptrace(0)
    return state_0, state_1, kappa, tlist[t_idx]

def main():
    state0_disp, state1_disp, kappa, t_snr = run_phase_space_sim(dispersive=True)
    state0_full, state1_full, _, _ = run_phase_space_sim(dispersive=False)
    print(f"\n>> Generating Final Thesis Plot at t = {t_snr:.2f}...")
    a_op = destroy(state0_disp.dims[0][0])
    cent0_disp = expect(a_op, state0_disp); cent1_disp = expect(a_op, state1_disp)
    cent0_full = expect(a_op, state0_full); cent1_full = expect(a_op, state1_full)

    def plot_wigner_centered(ax, state0, state1, c0, c1, show_legend=False):
        mid_x = np.real(c0 + c1) / 2.0; mid_y = np.imag(c0 + c1) / 2.0
        R = 2.2
        grid_center_x = mid_x * np.sqrt(2); grid_center_y = mid_y * np.sqrt(2)
        grid_span = R * np.sqrt(2)
        w_xvec = np.linspace(grid_center_x - grid_span, grid_center_x + grid_span, 500)
        w_yvec = np.linspace(grid_center_y - grid_span, grid_center_y + grid_span, 500)
        W0 = wigner(state0, w_xvec, w_yvec); W1 = wigner(state1, w_xvec, w_yvec)
        scale = 1.0 / np.sqrt(2)
        x_cav = w_xvec * scale; y_cav = w_yvec * scale
        ax.contourf(x_cav, y_cav, W0, levels=20, cmap='Blues', alpha=0.6)
        ax.contourf(x_cav, y_cav, W1, levels=20, cmap='Reds', alpha=0.6)
        contour_frac = 0.1
        ax.contour(x_cav, y_cav, W0, levels=[contour_frac * np.max(W0)], colors='blue', linewidths=1.5, zorder=3)
        ax.contour(x_cav, y_cav, W1, levels=[contour_frac * np.max(W1)], colors='red', linewidths=1.5, zorder=3)
        ax.axhline(y=0, color='green', linestyle='--', linewidth=1.5, alpha=0.6, zorder=2)
        ax.axvline(x=0, color='green', linestyle='--', linewidth=1.5, alpha=0.6, zorder=2)
        # CHANGED: |0⟩ → |g⟩ and |1⟩ → |e⟩
        ax.plot(np.real(c0), np.imag(c0), 'o', color='blue', ms=8, label=r'$|g\rangle$', zorder=5)
        ax.plot(np.real(c1), np.imag(c1), 'o', color='red', ms=8, label=r'$|e\rangle$', zorder=5)
        ax.plot([np.real(c0), np.real(c1)], [np.imag(c0), np.imag(c1)], 'k-', lw=1.5, zorder=4)
        diff = c1 - c0
        angle_from_pos_x = np.degrees(np.angle(diff))
        angle = 180 - angle_from_pos_x
        if angle < 0: angle += 360
        if angle >= 360: angle -= 360
        ax.text(0.05, 0.95, f'$\\theta = {angle:.1f}^\\circ$', transform=ax.transAxes,
                fontsize=18, verticalalignment='top',
                bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.5', alpha=0.9))
        ax.set_xlabel(r'$\mathrm{Re}[\langle a \rangle]$', fontsize=22)
        ax.set_ylabel(r'$\mathrm{Im}[\langle a \rangle]$', fontsize=22)
        ax.set_xlim(mid_x - R, mid_x + R); ax.set_ylim(mid_y - R, mid_y + R)
        ax.set_aspect('equal'); ax.set_facecolor('white')
        if show_legend:
            ax.legend(loc='upper right', fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=18)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.5), constrained_layout=True)
    plot_wigner_centered(ax1, state0_disp, state1_disp, cent0_disp, cent1_disp, show_legend=True)
    plot_wigner_centered(ax2, state0_full, state1_full, cent0_full, cent1_full, show_legend=False)
    plt.savefig(OUTPUT_DIR / "thesis_wigner_final.pdf", dpi=300)
    plt.close()
    print("Saved thesis_wigner_final.pdf")

if __name__ == "__main__":
    main()
