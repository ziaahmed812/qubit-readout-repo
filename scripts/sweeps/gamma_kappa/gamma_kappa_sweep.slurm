#!/bin/bash
#SBATCH --job-name=gam-kap
#SBATCH --partition=large-long
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --array=0-24
#SBATCH --output=logs/gamma_kappa_%A_%a.out
#SBATCH --error=logs/gamma_kappa_%A_%a.err

# ============================================================
# Legacy filename preserved: this launcher actually sweeps kappa and epsilon.
# Preserved as an example launcher derived from the thesis runs.
# Adjust partition, wall time, and VENV_ACTIVATE for the target cluster.
#
# kappa: from 0.04 to 0.06 in steps of 0.005 (5 values)
# epsilon: from 0.04 to 0.06 in steps of 0.005 (5 values)
#
# Grid: k -> (i_kappa, j_epsilon) where i_kappa = k // 5, j_epsilon = k % 5
#
# Fixed parameters:
#   g = 0.341*2*pi, delta = +5*2*pi
#   phi_h = 0, phi_d = pi/2
#   gamma = 0.001*2*pi
#   N = 30, t_end = 200, num_points = 2000
#
# Expected runtime: longer due to N=30
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [[ -n "${VENV_ACTIVATE:-}" ]]; then
  # shellcheck disable=SC1090
  source "${VENV_ACTIVATE}"
fi

# Enforce threading invariant
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs
mkdir -p results

# Map array index -> (kappa, epsilon) values
read kappa_val epsilon_val < <(python3 - <<EOF
import numpy as np
k = int(${SLURM_ARRAY_TASK_ID})

# Kappa: from 0.04 to 0.06 in steps of 0.005 (5 values)
kappa_values = np.array([0.04, 0.045, 0.05, 0.055, 0.06])

# Epsilon: from 0.04 to 0.06 in steps of 0.005 (5 values)
epsilon_values = np.array([0.04, 0.045, 0.05, 0.055, 0.06])

N_eps = len(epsilon_values)  # 5

i_kappa = k // N_eps
j_epsilon = k % N_eps

kappa_val = kappa_values[i_kappa]
epsilon_val = epsilon_values[j_epsilon]

print(f"{kappa_val:.6f} {epsilon_val:.6f}")
EOF
)

echo "=============================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "kappa: ${kappa_val}"
echo "epsilon: ${epsilon_val}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}"
echo "Host: $(hostname)"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo "Started: $(date)"
echo "=============================================="

python gamma_kappa_sweep.py \
    --kappa "${kappa_val}" \
    --epsilon "${epsilon_val}" \
    --output_dir "results"

echo "Finished: $(date)"
