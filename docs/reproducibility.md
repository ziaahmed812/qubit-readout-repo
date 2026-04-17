# Reproducibility Guide

This repository supports external review of the readout results in the accompanying thesis PDF. It therefore distinguishes carefully between:

- results that can be regenerated directly from the code included here
- results that can be regenerated here from included compact data files
- results for which the repository preserves the original computational workflow and final artifact, but not the full heavy sweep outputs

## Recorded Thesis Compute Stack

The larger QuTiP runs used in the thesis were recorded on a CPU-only compute environment with:

- **Python**: 3.12
- **QuTiP**: 5.2.2
- **NumPy**: 2.3.x
- **SciPy**: 1.16.x
- **Matplotlib**
- **tqdm**

Additional system dependencies:

- a LaTeX distribution, because many scripts use `matplotlib` with `text.usetex=True`
- `ffmpeg`, for the MP4 animations

This is the **recorded thesis compute stack**, not a fully revalidated end-to-end environment guarantee for every preserved heavy-sweep path in this submission repository.

## Suggested Local Environment

The following is a best-effort reconstruction of the recorded thesis stack:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "qutip==5.2.2" "numpy>=2.3,<2.4" "scipy>=1.16,<1.17" matplotlib tqdm
```

If `python3.12` is not available locally, another recent Python 3 version may still be adequate for code inspection and some smaller reruns, but the configuration above is the closest documented match to the thesis compute environment.

## Reproducibility Classes

### 1. Directly reproducible here

These figures and media can be regenerated directly from the scripts included in this repository:

- Figures 3.1a-3.2b
- Figures 3.11-3.12
- Figure B.1
- the two supplementary animations

### 2. Reproducible here from included selected data

These figures can be regenerated here because the required compact validation data files are included in `data/selected/validation/`:

- Figures 3.3a-3.5

### 3. Provenance-preserved; full rerun requires regenerating omitted sweep outputs

These figures are represented here by:

- the final thesis artifact
- the original scientific sweep driver
- the postprocess script that turns sweep outputs into the final artifact

but **not** by the large archived sweep outputs themselves:

- Figures 3.6a-3.10

This is a deliberate design choice for the submission repository. The heavy output bundles are not required for reviewers to inspect the methodology or the final thesis artifacts, but a complete rerun of those figures would require regenerating the omitted sweep results externally.

## Directly Reproducible Figures

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

## Selected-Data Validation Figures

These scripts use compact included data files rather than omitted sweep bundles.

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

These figures correspond to the validation parameter sets summarized in **Tables 3.1 and 3.2** of the accompanying thesis PDF.

## Heavy Sweep Figures

The following thesis figures depend on larger parameter sweeps:

- Figures 3.6a-3.7b: `scripts/sweeps/g_delta/g_delta_sweep.py`
- Figure 3.8a: `scripts/sweeps/kappa_epsilon/kappa_epsilon_sweep.py`
- Figure 3.8b: `scripts/sweeps/kappa_epsilon_grid/kappa_epsilon_grid_sweep.py`
- Figure 3.9: `scripts/sweeps/gamma/gamma_sweep.py`
- Figure 3.10: `scripts/sweeps/delta_phi/noise_sweep.py`

The corresponding artifact-generation scripts are preserved in:

- `scripts/postprocess/optimization/`

If a full rerun of these figures is desired, the omitted sweep outputs should be regenerated and placed under:

- `external_results/`

The expected directory structure is documented in `external_results/README.md`.

## Compute Notes For Heavy Sweeps

The recorded thesis runner notes for the larger scans describe a CPU-only environment with:

- 4 CPU cores
- 15 GiB RAM
- Python 3.12
- QuTiP 5.2.2

That information is mainly relevant to the heavy sweep drivers rather than the smaller self-contained figure scripts.
