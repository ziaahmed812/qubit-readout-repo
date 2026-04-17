#!/usr/bin/env python3
"""
Create animations of cavity field evolution for Dispersive H and Full JC H.

Parameters (matching noise_sweep.py):
  - t_end = 200, num_points = 4000
  - g = 2π × 0.2, kappa = 2π × 0.05, epsilon = 2π × 0.04
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from qutip import *
from qutip.wigner import wigner
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MEDIA_DIR = REPO_ROOT / "media" / "animations"

# ============================================================
# Global LaTeX Style
# ============================================================
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "lines.linewidth": 2,
    "axes.linewidth": 1.2,
})

def run_simulation(dispersive=False):
    """Run simulation and return all states."""
    h_type = "dispersive" if dispersive else "full JC"
    print(f"\n>> Running simulation ({h_type})...")
    sys.stdout.flush()
    
    solver_opts = {
        'rtol': 1e-6, 'atol': 1e-8, 'nsteps': 5000,
        'store_states': True, 'store_final_state': True
    }

    # Parameters (matching noise_sweep.py exactly)
    wr = 2*np.pi * 1.0
    wa = 2*np.pi * 3.0
    wd = 2*np.pi * 1.0
    g       = 2*np.pi * 0.2
    kappa   = 2*np.pi * 0.05
    epsilon = 2*np.pi * 0.04
    gamma   = 2*np.pi * 0.001
    
    phi_h = 0.0
    phi_d = np.pi / 2.0
    
    N = 14
    t_end = 500
    num_points = 5000
    
    delta = wa - wd
    chi = g**2 / delta
    
    print(f"   Parameters: g={g/(2*np.pi):.2f}, kappa={kappa/(2*np.pi):.2f}, epsilon={epsilon/(2*np.pi):.2f}")
    print(f"   t_end={t_end}, num_points={num_points}")

    tlist = np.linspace(0, t_end, num_points)

    # Build operators
    a  = tensor(destroy(N), qeye(2))
    sz = tensor(qeye(N), sigmaz())
    sm = tensor(qeye(N), sigmam())
    psi0 = tensor(basis(N, 0), basis(2, 0))
    psi1 = tensor(basis(N, 0), basis(2, 1))
    
    if dispersive:
        H = (
            chi * a.dag() * a * sz +
            epsilon * (a.dag() * np.exp(1j * phi_d) + a * np.exp(-1j * phi_d))
        )
        c_ops = [np.sqrt(kappa) * a]
    else:
        H = (
            (wr - wd) * a.dag() * a +
            (wa - wd) * sz / 2.0 +
            g * (a * sm.dag() + a.dag() * sm) +
            epsilon * (a.dag() * np.exp(1j * phi_d) + a * np.exp(-1j * phi_d))
        )
        c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma) * sm]
    
    print("   Solving for |0⟩...")
    res0 = mesolve(H, psi0, tlist, c_ops, e_ops=[a], options=solver_opts)
    print("   Solving for |1⟩...")
    res1 = mesolve(H, psi1, tlist, c_ops, e_ops=[a], options=solver_opts)
    
    # Extract all cavity states
    states_0 = [res0.states[i].ptrace(0) for i in range(num_points)]
    states_1 = [res1.states[i].ptrace(0) for i in range(num_points)]
    
    # Expectation values for trajectory
    a_expect_0 = res0.expect[0]
    a_expect_1 = res1.expect[0]
    
    print(f"   Done!")
    
    return states_0, states_1, a_expect_0, a_expect_1, tlist


def create_animation(states_0, states_1, a_expect_0, a_expect_1, tlist, h_type, output_file, show_max_snr=False):
    """Create animation of cavity field evolution."""
    print(f"\n>> Creating animation for {h_type}...")
    
    num_frames = len(tlist)
    a_op = destroy(states_0[0].dims[0][0])
    
    # Precompute all centers
    print("   Precomputing centers...")
    centers_0 = np.array([expect(a_op, s) for s in states_0])
    centers_1 = np.array([expect(a_op, s) for s in states_1])
    
    # Find global bounds for consistent view
    all_real = np.concatenate([np.real(centers_0), np.real(centers_1)])
    all_imag = np.concatenate([np.imag(centers_0), np.imag(centers_1)])
    
    # Use steady-state region for view (after t=50)
    steady_idx = np.argmin(np.abs(tlist - 50))
    mid_x = np.mean([np.real(centers_0[steady_idx:]).mean(), np.real(centers_1[steady_idx:]).mean()])
    mid_y = np.mean([np.imag(centers_0[steady_idx:]).mean(), np.imag(centers_1[steady_idx:]).mean()])
    
    R = 2.5  # View radius
    
    # Precompute Wigner grid (fixed for all frames)
    grid_center_x = mid_x * np.sqrt(2)
    grid_center_y = mid_y * np.sqrt(2)
    grid_span = R * np.sqrt(2)
    
    w_xvec = np.linspace(grid_center_x - grid_span, grid_center_x + grid_span, 150)
    w_yvec = np.linspace(grid_center_y - grid_span, grid_center_y + grid_span, 150)
    
    scale = 1.0 / np.sqrt(2)
    x_cav = w_xvec * scale
    y_cav = w_yvec * scale
    
    # Setup figure
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Initial empty plots
    contourf0 = None
    contourf1 = None
    contour0 = None
    contour1 = None
    
    # Trajectory lines (will grow) - more prominent
    traj0_line, = ax.plot([], [], 'b-', alpha=0.9, lw=2)
    traj1_line, = ax.plot([], [], 'r-', alpha=0.9, lw=2)
    
    # Center dots
    dot0, = ax.plot([], [], 'o', color='blue', ms=10, label=r'$|0\rangle$', zorder=10)
    dot1, = ax.plot([], [], 'o', color='red', ms=10, label=r'$|1\rangle$', zorder=10)
    
    # Separation vector
    sep_line, = ax.plot([], [], 'k-', lw=2, zorder=9)
    
    # Coordinate axes
    ax.axhline(y=0, color='green', linestyle='--', linewidth=1.5, alpha=0.5, zorder=1)
    ax.axvline(x=0, color='green', linestyle='--', linewidth=1.5, alpha=0.5, zorder=1)
    
    # Time text
    time_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, fontsize=14,
                        verticalalignment='top', 
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    # Max SNR marker
    t_max_snr = 28.41
    snr_text = ax.text(0.02, 0.88, '', transform=ax.transAxes, fontsize=12,
                       verticalalignment='top', color='green',
                       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
    
    ax.set_xlim(mid_x - R, mid_x + R)
    ax.set_ylim(mid_y - R, mid_y + R)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\mathrm{Re}[\langle a \rangle]$')
    ax.set_ylabel(r'$\mathrm{Im}[\langle a \rangle]$')
    ax.set_title(f'{h_type} Hamiltonian - Cavity Field Evolution', fontsize=14)
    ax.legend(loc='upper right')
    ax.set_facecolor('white')
    
    # Frame skip for reasonable file size (every 10th frame = 400 frames total)
    frame_skip = 10
    frame_indices = list(range(0, num_frames, frame_skip))
    
    # Find frame closest to max SNR and add pause frames there
    t_max_snr_idx = np.argmin(np.abs(tlist - t_max_snr))
    max_snr_frame = t_max_snr_idx // frame_skip
    
    # Insert extra copies of max SNR frame for pause effect (90 extra frames = 3 second pause)
    # Only for full JC, not dispersive
    pause_frames = 90
    if show_max_snr:
        insert_pos = frame_indices.index(max_snr_frame * frame_skip) if (max_snr_frame * frame_skip) in frame_indices else None
        if insert_pos is not None:
            for _ in range(pause_frames):
                frame_indices.insert(insert_pos + 1, max_snr_frame * frame_skip)
    
    print(f"   Total frames: {len(frame_indices)} (from {num_frames} time points)")
    
    def init():
        traj0_line.set_data([], [])
        traj1_line.set_data([], [])
        dot0.set_data([], [])
        dot1.set_data([], [])
        sep_line.set_data([], [])
        time_text.set_text('')
        snr_text.set_text('')
        return traj0_line, traj1_line, dot0, dot1, sep_line, time_text, snr_text
    
    def animate(frame_num):
        idx = frame_indices[frame_num]
        t = tlist[idx]
        
        # Get current states
        state0 = states_0[idx]
        state1 = states_1[idx]
        c0 = centers_0[idx]
        c1 = centers_1[idx]
        
        # Clear previous contours
        for coll in ax.collections[:]:
            coll.remove()
        
        # Recalculate Wigner
        W0 = wigner(state0, w_xvec, w_yvec)
        W1 = wigner(state1, w_xvec, w_yvec)
        
        # Plot contours
        ax.contourf(x_cav, y_cav, W0, levels=15, cmap='Blues', alpha=0.6, zorder=2)
        ax.contourf(x_cav, y_cav, W1, levels=15, cmap='Reds', alpha=0.6, zorder=2)
        
        # Contour borders
        if np.max(W0) > 0:
            ax.contour(x_cav, y_cav, W0, levels=[0.1 * np.max(W0)], 
                      colors='blue', linewidths=1.5, zorder=3)
        if np.max(W1) > 0:
            ax.contour(x_cav, y_cav, W1, levels=[0.1 * np.max(W1)], 
                      colors='red', linewidths=1.5, zorder=3)
        
        # Update trajectory (up to current time)
        traj0_line.set_data(np.real(centers_0[:idx+1]), np.imag(centers_0[:idx+1]))
        traj1_line.set_data(np.real(centers_1[:idx+1]), np.imag(centers_1[:idx+1]))
        
        # Update dots
        dot0.set_data([np.real(c0)], [np.imag(c0)])
        dot1.set_data([np.real(c1)], [np.imag(c1)])
        
        # Update separation vector
        sep_line.set_data([np.real(c0), np.real(c1)], [np.imag(c0), np.imag(c1)])
        
        # Update time text
        time_text.set_text(f'$t = {t:.1f}$')
        
        # Max SNR indicator - only show during pause at max SNR (full JC only)
        if show_max_snr and abs(t - t_max_snr) < 0.5:
            snr_text.set_text(f'MAX SNR ($t = {t_max_snr:.1f}$)')
        else:
            snr_text.set_text('')
        
        if frame_num % 50 == 0:
            print(f"   Frame {frame_num}/{len(frame_indices)} (t={t:.1f})")
        
        return traj0_line, traj1_line, dot0, dot1, sep_line, time_text, snr_text
    
    print("   Rendering animation...")
    anim = FuncAnimation(fig, animate, init_func=init, frames=len(frame_indices),
                        interval=50, blit=False)
    
    # Save as MP4 (15 fps for slower playback - 2x slower than before)
    writer = FFMpegWriter(fps=15, metadata=dict(artist='QuTiP'), bitrate=3000)
    anim.save(output_file, writer=writer)
    
    plt.close()
    print(f"   Saved: {output_file}")


def main():
    from pathlib import Path
    SCRIPT_DIR = Path(__file__).parent
    
    # Run dispersive simulation (no max SNR pause)
    states0_disp, states1_disp, a0_disp, a1_disp, tlist = run_simulation(dispersive=True)
    create_animation(states0_disp, states1_disp, a0_disp, a1_disp, tlist,
                    "Dispersive", MEDIA_DIR / "dispersive_evolution.mp4", show_max_snr=False)
    
    # Run full JC simulation (with max SNR pause)
    states0_full, states1_full, a0_full, a1_full, tlist = run_simulation(dispersive=False)
    create_animation(states0_full, states1_full, a0_full, a1_full, tlist,
                    "Full JC", MEDIA_DIR / "full_jc_evolution.mp4", show_max_snr=True)
    
    print("\n" + "="*60)
    print("  ANIMATIONS COMPLETE")
    print("="*60)


if __name__ == "__main__":
    main()
