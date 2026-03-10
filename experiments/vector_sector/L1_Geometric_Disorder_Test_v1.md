# L1 Geometric Disorder Test

## Purpose

Test whether the restricted transverse band survives loss of exact lattice symmetry when only the node embedding is perturbed.

## Setup

- sizes: `[12, 16, 20]`
- baseline periodic topology preserved
- perturbation: `x_i -> x_i + sigma xi`
- `sigma = 0.05 h`
- deterministic seed base: `20260310`

## Lowest restricted mode comparison

| `n` | ordered `lambda_1` | disordered `lambda_1` | ordered IPR | disordered IPR | first-10 spectral shift |
| --- | ---: | ---: | ---: | ---: | ---: |
| 12 | 0.263337 | 0.263251 | 0.000454 | 0.000434 | 0.0503 |
| 16 | 0.150761 | 0.150738 | 0.000163 | 0.000185 | 0.0501 |
| 20 | 0.097277 | 0.097268 | 0.000119 | 0.000100 | 0.0501 |

## Direct result

- observation: the restricted transverse spectrum persists under weak geometric disorder, with the ordered and disordered scaled bands remaining close and the lowest restricted mode staying divergence-free
- conclusion: in the tested weak-disorder regime, loss of exact lattice symmetry does not destroy the low restricted transverse band

## Artifacts

- results: `data/20260310_144621_L1_geometric_disorder_test.json`
- plots: `plots/20260310_144621_transverse_band_disorder.png`, `plots/20260310_144621_divergence_curl_phase_disorder.png`
- timestamp: `20260310_144621`
