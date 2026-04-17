#!/bin/bash
#SBATCH --job-name=g-delta
#SBATCH --partition=small-short,large-short,small-long,large-long
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --array=0-5150
#SBATCH --output=logs/g_delta_%A_%a.out
#SBATCH --error=logs/g_delta_%A_%a.err

# ============================================================
# g-delta Grid Sweep: 101 x 51 = 5151 jobs
#
# delta: 101 values from -5*2*pi to +5*2*pi (step pi/5)
# g: 51 values from 0.001*2*pi to 0.5*2*pi (step pi/50)
#
# Fixed parameters:
#   phi_h = 0, phi_d = pi/2
#   kappa = 0.05*2*pi, epsilon = 0.05*2*pi, gamma = 0.001*2*pi
#   N = 20, t_end = 200, num_points = 4000
#
# Preserved as an example launcher derived from the thesis runs.
# Adjust partition, wall time, and VENV_ACTIVATE for the target cluster.
#
# Expected runtime: ~6-10 hours per job
# Memory: ~8-10 GB per job
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [[ -n "${VENV_ACTIVATE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${VENV_ACTIVATE}"
fi

# Enforce threading invariant (4 cores per job)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Create directories
mkdir -p logs
mkdir -p results

# Map array index -> (g, delta)
# Grid: 101 delta values x 51 g values = 5151 total
# k -> (i_delta, j_g) where i_delta = k // 51, j_g = k % 51
read g_val delta_val < <(python3 - <<EOF
import numpy as np
k = int(${SLURM_ARRAY_TASK_ID})

# Delta: from -5*2*pi to +5*2*pi in steps of pi/5 (101 values)
delta_start = -5 * 2 * np.pi
delta_end = 5 * 2 * np.pi
delta_step = np.pi / 5
delta_values = np.arange(delta_start, delta_end + delta_step/2, delta_step)

# g: from 0.001*2*pi to 0.5*2*pi in steps of pi/50 (51 values)
g_start = 0.001 * 2 * np.pi
g_end = 0.5 * 2 * np.pi
g_step = np.pi / 50
g_values = np.arange(g_start, g_end + g_step/2, g_step)

N_g = len(g_values)  # 51

i_delta = k // N_g
j_g = k % N_g

g_val = g_values[j_g]
delta_val = delta_values[i_delta]

print(f"{g_val:.10f} {delta_val:.10f}")
EOF
)

echo "=============================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "g: ${g_val}"
echo "delta: ${delta_val}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}"
echo "Host: $(hostname)"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo "Started: $(date)"
echo "=============================================="

# Run simulation
python g_delta_sweep.py \
    --g "${g_val}" \
    --delta "${delta_val}" \
    --output_dir "results"

echo "Finished: $(date)"
