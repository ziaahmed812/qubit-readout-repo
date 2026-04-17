# External Sweep Outputs

This directory is intentionally empty in the public review repository.

The optimization figures in the thesis come from larger parameter sweeps. This repository includes:

- the final thesis figures
- the original sweep scripts
- the scripts that turn sweep outputs into figures

What it does **not** include is the bulky saved output of those sweeps.

If those figures ever need to be rerun in full, regenerated sweep outputs should be placed here using the following layout:

- `g_delta_grid/results/`
- `kappa_epsilon/results/`
- `kappa_epsilon_grid/results/`
- `gamma_sweep/results/`
- `delta_phi/results/`

So this directory is simply a placeholder for those external sweep results. The workflow is preserved in the repository, but the heavy saved outputs themselves are not bundled here.
