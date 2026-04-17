# Thesis Figure Map

This page is meant to be read alongside the accompanying thesis PDF. It answers a simple question: for each thesis figure, what file in this repository corresponds to it, and what script or workflow produced it?

The figure numbering follows the LaTeX build metadata of the thesis source, so it matches the numbering seen in the compiled PDF.

Validation figures use the parameter sets summarized in **Tables 3.1 and 3.2** of the accompanying PDF. The main optimization outcomes are summarized in **Table 3.3**.

## Figure Overview

| Thesis figure | What it shows | What is included here | Reproducibility |
| --- | --- | --- | --- |
| Figure 3.1 | two viewpoints on the dispersive interaction | Figure 3.1a and 3.1b output files plus the direct plotting script | directly reproducible here |
| Figure 3.2 | Purcell benchmark for the driven Jaynes-Cummings model | Figure 3.2a and 3.2b output files plus the direct plotting scripts | directly reproducible here |
| Figure 3.3 | validation of SNR extraction in the dispersive limit | output files, validation drivers, selected data, and figure scripts | reproducible from included data |
| Figure 3.4 | breakdown of the dispersive approximation | output file, validation driver, selected data, and figure script | reproducible from included data |
| Figure 3.5 | fidelity beyond the dispersive regime | output file, selected data, and figure script | reproducible from included data |
| Figure 3.6 | optimization across detuning and coupling | output files, sweep script, and figure scripts | workflow preserved; full rerun needs regenerated sweep outputs |
| Figure 3.7 | time dependence along the optimal ridge | output files, sweep script, and figure script | workflow preserved; full rerun needs regenerated sweep outputs |
| Figure 3.8 | optimization across linewidth and drive strength | output files, sweep scripts, and figure scripts | workflow preserved; full rerun needs regenerated sweep outputs |
| Figure 3.9 | effect of intrinsic qubit decay | output file, sweep script, and figure script | workflow preserved; full rerun needs regenerated sweep outputs |
| Figure 3.10 | dependence on the homodyne phase | output file, sweep script, and figure script | workflow preserved; full rerun needs regenerated sweep outputs |
| Figure 3.11 | phase-space comparison of the readout states | output file plus direct plotting script | directly reproducible here |
| Figure 3.12 | time-resolved phase-space snapshots | output file plus direct plotting script | directly reproducible here |
| Figure B.1 | convergence diagnostics for the driven Jaynes-Cummings master equation | output file plus direct plotting script | directly reproducible here |

## Figure-By-Figure Map

| Thesis item | Scientific calculation | Figure-generation step | Output file in this repo | Reproducibility | Note |
| --- | --- | --- | --- | --- | --- |
| Figure 3.1a | `scripts/direct/regen_dispersive_H_interpretations.py` | same script | `figures/chapter3/cavity_shift.pdf` | directly reproducible here | cavity frequency shift panel |
| Figure 3.1b | `scripts/direct/regen_dispersive_H_interpretations.py` | same script | `figures/chapter3/stark_shift.pdf` | directly reproducible here | ac Stark shift panel |
| Figure 3.2a | `scripts/direct/gen_purcell_decay_example.py` | same script | `figures/chapter3/purcell_decay_example.pdf` | directly reproducible here | Purcell decay example |
| Figure 3.2b | `scripts/direct/gen_purcell_rates_loglog.py` | same script | `figures/chapter3/purcell_rates_loglog.pdf` | directly reproducible here | Purcell rate benchmark |
| Figure 3.3a | `scripts/validation/drivers/dispersive_H.py` | `scripts/validation/postprocess/regen_dispersive_H_SNR.py` using `data/selected/validation/dispersive_H_g0.20-data.txt` | `figures/chapter3/dispersive_H_g0.20_SNR.pdf` | reproducible from included data | dispersive SNR validation |
| Figure 3.3b | `scripts/validation/drivers/JC_H_dispersive.py` | `scripts/validation/postprocess/regen_JC_dispersive_g005.py` using `data/selected/validation/JC_H_dispersive_g0.05-data.txt` | `figures/chapter3/JC_H_dispersive_g0.05_SNR.pdf` | reproducible from included data | driven Jaynes-Cummings validation near the dispersive limit |
| Figure 3.4 | `scripts/validation/drivers/JC_H_BEYOND-dispersive-with-gamma.py` | `scripts/validation/postprocess/regen_JC_beyond_dispersive_g020.py` using `data/selected/validation/JC_H_BEYOND-dispersive-with-gamma_g0.20-data.txt` | `figures/chapter3/JC_H_BEYOND-dispersive-with-gamma_g0.20_SNR.pdf` | reproducible from included data | beyond-dispersive benchmark |
| Figure 3.5 | `scripts/validation/drivers/JC_H_BEYOND-dispersive-with-gamma.py` with `data/selected/validation/JC_H_BEYOND-dispersive-with-gamma_g0.20-data.txt` | `scripts/validation/postprocess/regen_fidelity.py` using `data/selected/validation/JC_H_BEYOND-dispersive-with-gamma_g0.20-data.txt` | `figures/chapter3/fidelity_JC_BEYOND_dispersive_g0.20.png` | reproducible from included data | fidelity extracted from the same beyond-dispersive benchmark dataset used for Figure 3.4 |
| Figure 3.6a | `scripts/sweeps/g_delta/g_delta_sweep.py` | `scripts/postprocess/optimization/figure_3_6a_max_snr_heatmap.py` | `figures/chapter3/max_snr_heatmap.pdf` | workflow preserved; full rerun needs regenerated sweep outputs | maximum-SNR heatmap |
| Figure 3.6b | `scripts/sweeps/g_delta/g_delta_sweep.py` | `scripts/postprocess/optimization/figure_3_6b_optimal_time_heatmap.py` | `figures/chapter3/optimal_time_heatmap.pdf` | workflow preserved; full rerun needs regenerated sweep outputs | optimal-time heatmap |
| Figure 3.7a | `scripts/sweeps/g_delta/g_delta_sweep.py` | `scripts/postprocess/optimization/figure_3_7_ridge_dynamics.py` | `figures/chapter3/snr_curves_along_g.pdf` | workflow preserved; full rerun needs regenerated sweep outputs | time traces along the coupling optimum at fixed detuning |
| Figure 3.7b | `scripts/sweeps/g_delta/g_delta_sweep.py` | `scripts/postprocess/optimization/figure_3_7_ridge_dynamics.py` | `figures/chapter3/snr_curves_along_delta.pdf` | workflow preserved; full rerun needs regenerated sweep outputs | time traces along the detuning optimum at fixed coupling |
| Figure 3.8a | `scripts/sweeps/kappa_epsilon/kappa_epsilon_sweep.py` | `scripts/postprocess/optimization/figure_3_8a_kappa_equals_epsilon_cut.py` | `figures/chapter3/kappa_epsilon_sweep_snr_vs_kappa_tau.pdf` | workflow preserved; full rerun needs regenerated sweep outputs | cut taken along equal drive strength and linewidth |
| Figure 3.8b | `scripts/sweeps/kappa_epsilon_grid/kappa_epsilon_grid_sweep.py` | `scripts/postprocess/optimization/figure_3_8b_kappa_epsilon_heatmap.py` | `figures/chapter3/gamma_kappa_max_snr_heatmap.pdf` and reviewer-facing alias `figures/chapter3/figure_3_8b_kappa_epsilon_heatmap.pdf` | workflow preserved; full rerun needs regenerated sweep outputs | the thesis filename is a legacy leftover; the actual figure is the linewidth-drive heatmap used for Figure 3.8b |
| Figure 3.9 | `scripts/sweeps/gamma/gamma_sweep.py` | `scripts/postprocess/optimization/figure_3_9_gamma_sweep.py` | `figures/chapter3/gamma_sweep_snr_vs_kappa_tau.pdf` | workflow preserved; full rerun needs regenerated sweep outputs | intrinsic-decay sweep |
| Figure 3.10 | `scripts/sweeps/delta_phi/noise_sweep.py` | `scripts/postprocess/optimization/figure_3_10_homodyne_phase_sensitivity.py` | `figures/chapter3/snr_vs_phi_d.pdf` | workflow preserved; full rerun needs regenerated sweep outputs | homodyne-phase sensitivity |
| Figure 3.11 | `scripts/direct/regen_thesis_wigner_final.py` | same script | `figures/chapter3/thesis_wigner_final.pdf` | directly reproducible here | phase-space comparison |
| Figure 3.12 | `scripts/direct/regen_thesis_wigner_snapshots.py` | same script | `figures/chapter3/thesis_wigner_snapshots.pdf` | directly reproducible here | time-resolved phase-space snapshots |
| Figure B.1 | `scripts/direct/regen_convergence.py` | same script | `figures/appendix/dynamic_convergence_plot.pdf` | directly reproducible here | convergence diagnostics |

## Supplementary Animations

The accompanying thesis PDF does not include the two MP4 files below, but they are part of this repository because they help visualize the same readout dynamics in a more direct way.

| Supplementary item | Scientific generation | Output file in this repo | Note |
| --- | --- | --- | --- |
| Supplementary Animation 1 | `scripts/animations/make_animations.py` | `media/animations/dispersive_evolution.mp4` | cavity-field evolution in the dispersive model |
| Supplementary Animation 2 | `scripts/animations/make_animations.py` | `media/animations/full_jc_evolution.mp4` | cavity-field evolution in the full Jaynes-Cummings model |
