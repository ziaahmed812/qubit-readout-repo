# Supplementary Media

This repository includes two MP4 animations as supplementary reviewer material. They are not embedded in the thesis figures, but they provide a compact visual summary of the cavity-field dynamics in the dispersive and full Jaynes-Cummings models.

## Included Files

| Asset | Output path | Generator | Summary |
| --- | --- | --- | --- |
| Supplementary Animation 1 | `media/animations/dispersive_evolution.mp4` | `scripts/animations/make_animations.py` | cavity-field evolution in the dispersive model |
| Supplementary Animation 2 | `media/animations/full_jc_evolution.mp4` | `scripts/animations/make_animations.py` | cavity-field evolution in the full Jaynes-Cummings model |

## How To Regenerate

From the repository root:

```bash
python scripts/animations/make_animations.py
```

This writes:

- `media/animations/dispersive_evolution.mp4`
- `media/animations/full_jc_evolution.mp4`

The script uses `matplotlib.animation` together with `ffmpeg`, so a working `ffmpeg` installation is required.

## Review Context

These animations should be read as supplementary visualizations of the same readout physics discussed in Chapter 3:

- the dispersive animation provides a reference evolution in the simpler regime
- the full Jaynes-Cummings animation highlights the stronger non-dispersive deformation of the cavity-field dynamics

They are included to make the readout dynamics easier to inspect during external review, not as a substitute for the thesis figures.
