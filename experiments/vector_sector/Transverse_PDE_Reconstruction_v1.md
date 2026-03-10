# Transverse PDE Reconstruction

## Purpose

Test whether the lowest restricted modes admit a local coarse-grained PDE of the form `curl curl A ≈ lambda A` with `div A ≈ 0`.

## Setup

- branch: baseline periodic torus
- sizes: `[12, 16]`
- restricted modes reconstructed: 3

## Residual summary

| `n` | mode | eigenvalue | relative curl-curl residual | relative divergence norm |
| --- | ---: | ---: | ---: | ---: |
| 12 | 1 | 0.263337 | 0.050648 | 0.000000 |
| 12 | 2 | 0.263337 | 0.050648 | 0.000000 |
| 12 | 3 | 0.263337 | 0.050648 | 0.000000 |
| 16 | 1 | 0.150761 | 0.028620 | 0.000000 |
| 16 | 2 | 0.150761 | 0.028620 | 0.000000 |
| 16 | 3 | 0.301523 | 0.028620 | 0.000000 |

## Direct result

- observation: after coarse-graining to a periodic vector field, the lowest restricted modes remain nearly divergence-free and the curl-curl residual decreases with lattice size
- conclusion: for the lowest mode, the relative curl-curl residual decreases from 0.051 at n=12 to 0.029 at n=16, consistent with an emergent local coarse-grained PDE in the tested range

## Artifacts

- results: `data/20260310_154831_transverse_pde_reconstruction.json`
- plots: `plots/20260310_154831_pde_reconstruction_residuals.png`, `plots/20260310_154831_transverse_field_residual_maps.png`, `plots/20260310_154831_divergence_curl_phase_pde_reconstruction.png`
- timestamp: `20260310_154831`
