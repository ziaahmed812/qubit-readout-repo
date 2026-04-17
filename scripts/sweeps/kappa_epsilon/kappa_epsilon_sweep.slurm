#!/bin/bash
#SBATCH --job-name=kappa-eps
#SBATCH --partition=small-long
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --array=0-10
#SBATCH --output=logs/kappa_eps_%A_%a.out
#SBATCH --error=logs/kappa_eps_%A_%a.err

# ============================================================
# Kappa=Epsilon Sweep: 11 jobs
#
# epsilon (=kappa): from 0.01 to 0.51 in steps of 0.05 (11 values)
#        0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.51
#
# Fixed parameters:
#   g = 0.341*2*pi, delta = +5*2*pi
#   phi_h = 0, phi_d = pi/2
#   gamma = 0.001*2*pi
#   N = 20, t_end = 200, num_points = 2000
#
# Preserved as an example launcher derived from the thesis runs.
# Adjust partition, wall time, and VENV_ACTIVATE for the target cluster.
#
# Expected runtime: ~6-12 hours per job
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

# Map array index -> epsilon value
read epsilon_val < <(python3 - <<EOF
import numpy as np
k = int(${SLURM_ARRAY_TASK_ID})

# Epsilon: from 0.01 to 0.51 in steps of 0.05 (11 values)
epsilon_start = 0.01
epsilon_step = 0.05
epsilon_val = epsilon_start + k * epsilon_step

print(f"{epsilon_val:.6f}")
EOF
)

echo "=============================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "epsilon (=kappa): ${epsilon_val}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}"
echo "Host: $(hostname)"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo "Started: $(date)"
echo "=============================================="

python kappa_epsilon_sweep.py \
    --epsilon "${epsilon_val}" \
    --output_dir "results"

echo "Finished: $(date)"
