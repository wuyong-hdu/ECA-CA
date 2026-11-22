# Efficient Passive Radar Clutter Suppression using Correlation-Accumulation-based Extensive Cancellation Algorithm

ECA-CA accumulates blockwise correlations to estimate one global clutter
coefficient vector, preserving full-CPI ECA processing with lower computational
and memory requirements. This repository contains the 14 experimental
configurations, measured data and results used for Figures 3-8.

## Reproduce figures and display saved results

In MATLAB R2024a, open this directory and run:

```matlab
results = reproduce_results();
```

This prints PTR, mean execution time, peak working memory and relative residual
error for all 14 configurations, and displays Figures 3-8 using the saved results.
It does not rerun the clutter suppression algorithms. Figures (PDF/PNG/FIG) and
the numerical summary are saved under `output/`. Use `reproduce_results(false)`
to export without keeping figure windows open.

The published figures and numerical results can also be viewed directly in
`results/figures/`, `results/summary.txt` and `results/metrics.csv`.

## Run the experiments

```matlab
results = run_experiments();
```

The included `data/analysis_frame.mat` contains the reference and surveillance
channels for the selected 200 ms frame, plus the boundary samples required by
ECA-FS (278 ms in total). New results are saved under `runs/`.

## Configuration

The three parameter series contain 14 unique configurations:

| Varied parameter | Values | Fixed parameters |
|---|---|---|
| B: number of processing blocks | 25, 50, 100, 200, 400, 800, 1600, 3200 | N=1600000, K=400 |
| K: clutter suppression order | 100, 200, 400, 800 | N=1600000, B=200 |
| N: signal length | 200000, 400000, 800000, 1600000 | K=400, B=200 |

Sampling rate: 8 MHz. Target evaluation: range bin 9, Doppler -55 Hz.
PTR threshold factor: 4.5. Timing: one warm-up and 50 repetitions, using eight
computational-library threads. Settings are in `src/section5_config.m`.

Figures 3, 4, 7 and 8 show relative error, PTR, execution time and working memory.
Figure 5 uses N=1600000, K=400, B=25; Figure 6 compares B=25, 100, 400 and 1600
at N=1600000 and K=400.

Noncommercial use only under the [PolyForm Noncommercial License 1.0.0](LICENSE).
