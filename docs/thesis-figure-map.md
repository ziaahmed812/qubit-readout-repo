# Thesis Figure Map

This document maps the repository outputs to the figure numbering used in the accompanying thesis PDF. The numbering below follows the LaTeX build metadata of the thesis source, so it matches the printed figure numbers seen by reviewers.

Validation figures draw on the parameter sets summarized in **Table 3.1** and **Table 3.2** of the accompanying PDF. The optimization outcomes are summarized in **Table 3.3**.

## Parent Figure Overview

| Thesis figure | Thesis title / scientific focus | Repository support | Reproducibility class |
| --- | --- | --- | --- |
| Figure 3.1 | interpretations of the dispersive interaction | Figure 3.1a and Figure 3.1b artifacts plus direct plotting script | directly reproducible here |
| Figure 3.2 | Purcell benchmark for the driven Jaynes-Cummings model | Figure 3.2a and Figure 3.2b artifacts plus direct plotting scripts | directly reproducible here |
| Figure 3.3 | validation of SNR extraction in the dispersive limit | Figure 3.3a and Figure 3.3b artifacts, validation drivers, selected data, and postprocess scripts | reproducible here from included selected data |
| Figure 3.4 | deviation from the dispersive approximation | final artifact, validation driver, selected data, and postprocess script | reproducible here from included selected data |
| Figure 3.5 | measurement fidelity beyond the dispersive regime | final artifact, selected data, and postprocess script | reproducible here from included selected data |
| Figure 3.6 | optimisation map in the detuning-coupling plane | final artifacts, original sweep driver, and optimization postprocess scripts | provenance-preserved; full rerun requires regenerating omitted sweep outputs |
| Figure 3.7 | SNR dynamics along the optimal ridge | final artifacts, original sweep driver, and optimization postprocess script | provenance-preserved; full rerun requires regenerating omitted sweep outputs |
| Figure 3.8 | drive and linewidth optimisation in the linewidth-drive plane | final artifacts, original sweep drivers, and optimization postprocess scripts | provenance-preserved; full rerun requires regenerating omitted sweep outputs |
| Figure 3.9 | impact of intrinsic qubit decay | final artifact, original sweep driver, and optimization postprocess script | provenance-preserved; full rerun requires regenerating omitted sweep outputs |
| Figure 3.10 | sensitivity of the maximum SNR to the homodyne phase | final artifact, original sweep driver, and optimization postprocess script | provenance-preserved; full rerun requires regenerating omitted sweep outputs |
| Figure 3.11 | phase-space view of the readout states | final artifact plus direct plotting script | directly reproducible here |
| Figure 3.12 | temporal phase-space dynamics of the readout states | final artifact plus direct plotting script | directly reproducible here |
| Figure B.1 | convergence diagnostics for the driven Jaynes-Cummings master equation | final artifact plus direct plotting script | directly reproducible here |

## Panel And Artifact Map

| Thesis item | Scientific data generation | Artifact generation / postprocessing | Final repository artifact | Reproducibility class | Notes |
| --- | --- | --- | --- | --- | --- |
| Figure 3.1a | `scripts/direct/regen_dispersive_H_interpretations.py` | same script | `figures/chapter3/cavity_shift.pdf` | directly reproducible here | cavity frequency shift panel |
| Figure 3.1b | `scripts/direct/regen_dispersive_H_interpretations.py` | same script | `figures/chapter3/stark_shift.pdf` | directly reproducible here | ac Stark shift panel |
| Figure 3.2a | `scripts/direct/gen_purcell_decay_example.py` | same script | `figures/chapter3/purcell_decay_example.pdf` | directly reproducible here | Purcell decay example |
| Figure 3.2b | `scripts/direct/gen_purcell_rates_loglog.py` | same script | `figures/chapter3/purcell_rates_loglog.pdf` | directly reproducible here | Purcell rate benchmark |
| Figure 3.3a | `scripts/validation/drivers/dispersive_H.py` | `scripts/validation/postprocess/regen_dispersive_H_SNR.py` using `data/selected/validation/dispersive_H_g0.20-data.txt` | `figures/chapter3/dispersive_H_g0.20_SNR.pdf` | reproducible here from included selected data | dispersive SNR validation |
| Figure 3.3b | `scripts/validation/drivers/JC_H_dispersive.py` | `scripts/validation/postprocess/regen_JC_dispersive_g005.py` using `data/selected/validation/JC_H_dispersive_g0.05-data.txt` | `figures/chapter3/JC_H_dispersive_g0.05_SNR.pdf` | reproducible here from included selected data | driven JC dispersive-limit validation |
| Figure 3.4 | `scripts/validation/drivers/JC_H_BEYOND-dispersive-with-gamma.py` | `scripts/validation/postprocess/regen_JC_beyond_dispersive_g020.py` using `data/selected/validation/JC_H_BEYOND-dispersive-with-gamma_g0.20-data.txt` | `figures/chapter3/JC_H_BEYOND-dispersive-with-gamma_g0.20_SNR.pdf` | reproducible here from included selected data | beyond-dispersive benchmark |
| Figure 3.5 | `scripts/validation/drivers/JC_H_BEYOND-dispersive-with-gamma.py` with `data/selected/validation/JC_H_BEYOND-dispersive-with-gamma_g0.20-data.txt` | `scripts/validation/postprocess/regen_fidelity.py` using `data/selected/validation/JC_H_BEYOND-dispersive-with-gamma_g0.20-data.txt` | `figures/chapter3/fidelity_JC_BEYOND_dispersive_g0.20.png` | reproducible here from included selected data | fidelity derived from the beyond-dispersive benchmark data |
| Figure 3.6a | `scripts/sweeps/g_delta/g_delta_sweep.py` | `scripts/postprocess/optimization/figure_3_6a_max_snr_heatmap.py` | `figures/chapter3/max_snr_heatmap.pdf` | provenance-preserved; full rerun requires regenerating omitted sweep outputs | maximum-SNR heatmap |
| Figure 3.6b | `scripts/sweeps/g_delta/g_delta_sweep.py` | `scripts/postprocess/optimization/figure_3_6b_optimal_time_heatmap.py` | `figures/chapter3/optimal_time_heatmap.pdf` | provenance-preserved; full rerun requires regenerating omitted sweep outputs | optimal-time heatmap |
| Figure 3.7a | `scripts/sweeps/g_delta/g_delta_sweep.py` | `scripts/postprocess/optimization/figure_3_7_ridge_dynamics.py` | `figures/chapter3/snr_curves_along_g.pdf` | provenance-preserved; full rerun requires regenerating omitted sweep outputs | ridge dynamics along the coupling optimum at fixed detuning |
| Figure 3.7b | `scripts/sweeps/g_delta/g_delta_sweep.py` | `scripts/postprocess/optimization/figure_3_7_ridge_dynamics.py` | `figures/chapter3/snr_curves_along_delta.pdf` | provenance-preserved; full rerun requires regenerating omitted sweep outputs | ridge dynamics along the detuning optimum at fixed coupling |
| Figure 3.8a | `scripts/sweeps/kappa_epsilon/kappa_epsilon_sweep.py` | `scripts/postprocess/optimization/figure_3_8a_kappa_equals_epsilon_cut.py` | `figures/chapter3/kappa_epsilon_sweep_snr_vs_kappa_tau.pdf` | provenance-preserved; full rerun requires regenerating omitted sweep outputs | cut along equal drive strength and linewidth |
| Figure 3.8b | `scripts/sweeps/kappa_epsilon_grid/kappa_epsilon_grid_sweep.py` | `scripts/postprocess/optimization/figure_3_8b_kappa_epsilon_heatmap.py` | `figures/chapter3/gamma_kappa_max_snr_heatmap.pdf` and reviewer-facing alias `figures/chapter3/figure_3_8b_kappa_epsilon_heatmap.pdf` | provenance-preserved; full rerun requires regenerating omitted sweep outputs | the thesis artifact filename is a legacy name; the actual content is the linewidth-drive heatmap described in Figure 3.8b |
| Figure 3.9 | `scripts/sweeps/gamma/gamma_sweep.py` | `scripts/postprocess/optimization/figure_3_9_gamma_sweep.py` | `figures/chapter3/gamma_sweep_snr_vs_kappa_tau.pdf` | provenance-preserved; full rerun requires regenerating omitted sweep outputs | intrinsic-decay sweep |
| Figure 3.10 | `scripts/sweeps/delta_phi/noise_sweep.py` | `scripts/postprocess/optimization/figure_3_10_homodyne_phase_sensitivity.py` | `figures/chapter3/snr_vs_phi_d.pdf` | provenance-preserved; full rerun requires regenerating omitted sweep outputs | homodyne-phase sensitivity |
| Figure 3.11 | `scripts/direct/regen_thesis_wigner_final.py` | same script | `figures/chapter3/thesis_wigner_final.pdf` | directly reproducible here | phase-space comparison |
| Figure 3.12 | `scripts/direct/regen_thesis_wigner_snapshots.py` | same script | `figures/chapter3/thesis_wigner_snapshots.pdf` | directly reproducible here | time-resolved phase-space snapshots |
| Figure B.1 | `scripts/direct/regen_convergence.py` | same script | `figures/appendix/dynamic_convergence_plot.pdf` | directly reproducible here | convergence diagnostics |

## Supplementary Reviewer Media

The accompanying thesis PDF does not embed the two MP4 files below, but they are included here as supplementary reviewer material.

| Supplementary item | Scientific generation | Final repository artifact | Notes |
| --- | --- | --- | --- |
| Supplementary Animation 1 | `scripts/animations/make_animations.py` | `media/animations/dispersive_evolution.mp4` | dispersive model cavity-field evolution |
| Supplementary Animation 2 | `scripts/animations/make_animations.py` | `media/animations/full_jc_evolution.mp4` | full Jaynes-Cummings cavity-field evolution |
