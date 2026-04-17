# Reproducibility Guide

This guide explains how much of the thesis material can be rerun directly from this repository, and where the limits are. The repository is meant for reviewer inspection and technical transparency, not as a full archival dump of every heavy parameter sweep.

## Three Reproducibility Levels

The material in this repository falls into three groups.

### 1. Directly reproducible here

These figures and media can be regenerated directly from the scripts in the repository:

- Figures 3.1a-3.2b
- Figures 3.11-3.12
- Figure B.1
- the two supplementary animations

### 2. Reproducible from included data

These figures can be regenerated because the required compact validation datasets are included under `data/selected/validation/`:

- Figures 3.3a-3.5

### 3. Workflow preserved; full rerun needs regenerated sweep outputs

For these figures, the repository includes:

- the final thesis figure
- the original sweep script
- the script that turns sweep outputs into the final figure

What it does **not** include is the full saved output of the large sweeps themselves:

- Figures 3.6a-3.10

That omission is deliberate. The goal here is to make the workflow inspectable and understandable without shipping bulky sweep archives that are not needed for review.

## Recorded Thesis Compute Stack

The larger QuTiP runs used in the thesis were recorded on a CPU-only environment with:

- **Python** 3.12
- **QuTiP** 5.2.2
- **NumPy** 2.3.x
- **SciPy** 1.16.x
- **Matplotlib**
- **tqdm**

Additional system dependencies:

- a LaTeX distribution, because several scripts use `matplotlib` with `text.usetex=True`
- `ffmpeg`, for the MP4 animations

This is the recorded thesis compute stack, not a promise that every preserved heavy-sweep path has been fully revalidated end-to-end inside this review repository.

## Suggested Local Environment

The following is a reasonable reconstruction of the recorded thesis environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "qutip==5.2.2" "numpy>=2.3,<2.4" "scipy>=1.16,<1.17" matplotlib tqdm
```

If `python3.12` is not available locally, another recent Python 3 version may still be adequate for code inspection and smaller reruns, but the configuration above is the closest documented match.

## Figures You Can Regenerate Directly

Run the following from the repository root.

### Dispersive interpretation and Purcell benchmark

```bash
python scripts/direct/regen_dispersive_H_interpretations.py
python scripts/direct/gen_purcell_decay_example.py
python scripts/direct/gen_purcell_rates_loglog.py
```

Expected outputs:

- `figures/chapter3/cavity_shift.pdf`
- `figures/chapter3/stark_shift.pdf`
- `figures/chapter3/purcell_decay_example.pdf`
- `figures/chapter3/purcell_rates_loglog.pdf`

### Phase-space figures

```bash
python scripts/direct/regen_thesis_wigner_final.py
python scripts/direct/regen_thesis_wigner_snapshots.py
```

Expected outputs:

- `figures/chapter3/thesis_wigner_final.pdf`
- `figures/chapter3/thesis_wigner_snapshots.pdf`

### Appendix convergence figure

```bash
python scripts/direct/regen_convergence.py
```

Expected output:

- `figures/appendix/dynamic_convergence_plot.pdf`

### Supplementary animations

```bash
python scripts/animations/make_animations.py
```

Expected outputs:

- `media/animations/dispersive_evolution.mp4`
- `media/animations/full_jc_evolution.mp4`

## Validation Figures Reproduced From Included Data

These scripts regenerate the thesis validation figures from compact data files that are included in the repository:

```bash
python scripts/validation/postprocess/regen_dispersive_H_SNR.py
python scripts/validation/postprocess/regen_JC_dispersive_g005.py
python scripts/validation/postprocess/regen_JC_beyond_dispersive_g020.py
python scripts/validation/postprocess/regen_fidelity.py
```

Expected outputs:

- `figures/chapter3/dispersive_H_g0.20_SNR.pdf`
- `figures/chapter3/JC_H_dispersive_g0.05_SNR.pdf`
- `figures/chapter3/JC_H_BEYOND-dispersive-with-gamma_g0.20_SNR.pdf`
- `figures/chapter3/fidelity_JC_BEYOND_dispersive_g0.20.png`

These correspond to the validation parameter sets summarized in **Tables 3.1 and 3.2** of the accompanying thesis PDF.

## Larger Optimization Figures

The following figures depend on larger parameter sweeps:

- Figures 3.6a-3.7b: `scripts/sweeps/g_delta/g_delta_sweep.py`
- Figure 3.8a: `scripts/sweeps/kappa_epsilon/kappa_epsilon_sweep.py`
- Figure 3.8b: `scripts/sweeps/kappa_epsilon_grid/kappa_epsilon_grid_sweep.py`
- Figure 3.9: `scripts/sweeps/gamma/gamma_sweep.py`
- Figure 3.10: `scripts/sweeps/delta_phi/noise_sweep.py`

The corresponding figure-generation scripts are preserved in:

- `scripts/postprocess/optimization/`

If these figures need to be rerun in full, regenerated sweep outputs should be placed under:

- `external_results/`

The expected directory layout is documented in [../external_results/README.md](../external_results/README.md).

## Notes On The Heavier Runs

The recorded thesis runner notes for the larger scans describe a CPU-only setup with:

- 4 CPU cores
- 15 GiB RAM
- Python 3.12
- QuTiP 5.2.2

That information matters mainly for the larger sweep scripts rather than the smaller self-contained figure scripts.
