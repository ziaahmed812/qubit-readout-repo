#!/bin/bash
#SBATCH --job-name=fixed-phih
#SBATCH --partition=small-short
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-360
#SBATCH --output=logs/noise_%A_%a.out
#SBATCH --error=logs/noise_%A_%a.err

# ============================================================
# Fixed phi_h=0, sweep phi_d: 361 jobs
# phi_d: 361 values from 0 to 2*pi
# phi_h: 0 (fixed)
#
# Preserved as an example launcher derived from the thesis runs.
# Adjust partition, wall time, and VENV_ACTIVATE for the target cluster.
#
# Expected runtime: ~6 hours per job
# Memory: ~6-8 GB per job
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

# Map array index -> phi_d (361 steps from 0 to 2*pi)
read phi_d < <(python3 - <<EOF
import numpy as np
k = int(${SLURM_ARRAY_TASK_ID})

N_phi = 361
phi_d_values = np.linspace(0, 2*np.pi, N_phi)
phi_d = phi_d_values[k]

print(f"{phi_d:.6f}")
EOF
)

echo "=============================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "phi_h: 0.000000 (FIXED)"
echo "phi_d: ${phi_d}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "=============================================="

# Run simulation
python noise_sweep.py \
    --phi_d "${phi_d}" \
    --output_dir "results"

echo "Finished: $(date)"
