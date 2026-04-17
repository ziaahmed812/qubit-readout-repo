#!/bin/bash
#SBATCH --job-name=gamma-sweep
#SBATCH --partition=small-long
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --array=0-20
#SBATCH --output=logs/gamma_sweep_%A_%a.out
#SBATCH --error=logs/gamma_sweep_%A_%a.err

# ============================================================
# Gamma Sweep: 21 jobs
#
# gamma: from 0.001 to 0.101 in steps of 0.005 (21 values)
#        0.001, 0.006, 0.011, ..., 0.096, 0.101
#
# Fixed parameters:
#   g = 0.341*2*pi, delta = +5*2*pi
#   phi_h = 0, phi_d = pi/2
#   kappa = 0.05*2*pi, epsilon = 0.05*2*pi
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

# Map array index -> gamma value
read gamma_val < <(python3 - <<EOF
import numpy as np
k = int(${SLURM_ARRAY_TASK_ID})

# Gamma: from 0.001 to 0.101 in steps of 0.005 (21 values)
gamma_start = 0.001
gamma_step = 0.005
gamma_val = gamma_start + k * gamma_step

print(f"{gamma_val:.6f}")
EOF
)

echo "=============================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "gamma: ${gamma_val}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}"
echo "Host: $(hostname)"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo "Started: $(date)"
echo "=============================================="

python gamma_sweep.py \
    --gamma "${gamma_val}" \
    --output_dir "results"

echo "Finished: $(date)"
