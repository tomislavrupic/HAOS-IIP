# Transverse Backreaction Maps

## Purpose

Visualize how the lowest restricted transverse mode responds in space to a localized scalar/background deformation.

## Setup

- branch variants: `['baseline', 'puncture']`
- deformation family: `bump`
- amplitudes: `[0.05, 0.1]`
- lattice size: `12`

## Sensitivity summary

| variant | eta | mean sensitivity | correlation with scalar profile |
| --- | ---: | ---: | ---: |
| baseline | 0.05 | 0.008449 | -0.0799 |
| baseline | 0.1 | 0.008907 | -0.0934 |
| puncture | 0.05 | 0.004770 | 0.1507 |
| puncture | 0.1 | 0.005312 | 0.2862 |

## Direct result

- observation: localized scalar/background bumps produce nonuniform changes in the lowest restricted transverse mode rather than a flat redistribution across the lattice
- conclusion: for the baseline branch at eta=0.10, the mean sensitivity is 0.008907 with profile correlation -0.0934, so the response is spatially structured but not simply proportional to the scalar bump profile

## Artifacts

- results: `data/20260310_164414_transverse_backreaction_maps.json`
- plots: `plots/20260310_164414_transverse_backreaction_maps.png`, `plots/20260310_164414_transverse_sensitivity_vs_scalar_profile.png`
- timestamp: `20260310_164414`
