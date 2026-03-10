# Scalar-Transverse Coupling Matrix

## Purpose

Measure the low-dimensional response matrix `C_ij(eta) = <a_i, (T_eta - T_0) a_j>` for the low restricted transverse basis.

## Setup

- branch: `baseline`
- deformation family: `anisotropic`
- lattice size: `12`
- matrix modes: `10`

## Matrix summary

| eta | diagonal fraction | off-diagonal fraction | nearest-neighbor mean |
| ---: | ---: | ---: | ---: |
| 0.02 | 0.0738 | 0.9262 | 0.000015 |
| 0.05 | 0.1297 | 0.8703 | 0.000027 |
| 0.1 | 0.3365 | 0.6635 | 0.000038 |

## Direct result

- observation: the deformation-induced coupling matrix remains structured, but its weight is broadly distributed across the low-mode block rather than sharply diagonal-dominant
- conclusion: at eta=0.10, the diagonal fraction is 0.3365 and the off-diagonal fraction is 0.6635, indicating structured but non-diagonal mode mixing in the low transverse sector

## Artifacts

- results: `data/20260310_164216_scalar_transverse_coupling_matrix.json`
- plots: `plots/20260310_164216_transverse_coupling_matrix_heatmap.png`, `plots/20260310_164216_transverse_mode_mixing_vs_eta.png`
- timestamp: `20260310_164216`
