# Supplementary Media

This repository includes two MP4 animations that are not part of the thesis figure set, but are useful for inspecting the same readout physics in a more direct visual form.

## Included Files

| Asset | Output file | Generator | What it shows |
| --- | --- | --- | --- |
| Supplementary Animation 1 | `media/animations/dispersive_evolution.mp4` | `scripts/animations/make_animations.py` | cavity-field evolution in the dispersive model |
| Supplementary Animation 2 | `media/animations/full_jc_evolution.mp4` | `scripts/animations/make_animations.py` | cavity-field evolution in the full Jaynes-Cummings model |

## How To Regenerate Them

From the repository root:

```bash
python scripts/animations/make_animations.py
```

This produces:

- `media/animations/dispersive_evolution.mp4`
- `media/animations/full_jc_evolution.mp4`

The script uses `matplotlib.animation` together with `ffmpeg`, so `ffmpeg` needs to be available locally.

## How To Read Them

These animations should be read as supplementary visualizations of the same readout dynamics discussed in Chapter 3:

- the dispersive animation gives the simpler reference picture
- the full Jaynes-Cummings animation shows the stronger deformation and non-dispersive behavior

They are included to help a reviewer see the dynamics more directly, not to replace the thesis figures.
