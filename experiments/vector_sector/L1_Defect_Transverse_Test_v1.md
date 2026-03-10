# L1 Defect Transverse Test

## Purpose

Test whether punctures or line defects weaken harmonic dominance enough for the restricted transverse sector

`d1* d1` on `ker(d0*) intersect (H1)^perp`

to produce a descending low spectral band.

## Setup

- substrate: periodic cubic lattice with oriented `x`, `y`, `z` edges
- sizes: `[6, 8, 10]`
- variants:
  - baseline periodic torus
  - single cubic puncture: centered `2x2x2` node removal
  - line defect: removed `z`-edge column at central `(x, y)` index
- kernel parameter: `epsilon = 0.2`
- operator definitions unchanged:
  - `d0 = W1^(1/2) B0`
  - `d1 = W2^(1/2) C`
  - `L1 = d0 d0* + d1* d1`

## Restricted transverse floor

### Baseline periodic torus

| `n` | nodes | edges | harmonic dim | lowest restricted eigenvalue |
| --- | ---: | ---: | ---: | ---: |
| 6 | 216 | 648 | 3 | 0.932912 |
| 8 | 512 | 1536 | 3 | 0.563345 |
| 10 | 1000 | 3000 | 3 | 0.372535 |

### Single cubic puncture

| `n` | nodes | edges | harmonic dim | lowest restricted eigenvalue |
| --- | ---: | ---: | ---: | ---: |
| 6 | 208 | 612 | 3 | 0.684887 |
| 8 | 504 | 1500 | 3 | 0.441271 |
| 10 | 992 | 2964 | 3 | 0.317111 |

### Line defect

| `n` | nodes | edges | harmonic dim | lowest restricted eigenvalue |
| --- | ---: | ---: | ---: | ---: |
| 6 | 216 | 642 | 3 | 0.766501 |
| 8 | 512 | 1528 | 3 | 0.514298 |
| 10 | 1000 | 2990 | 3 | 0.350303 |

## Diagnostics

- divergence-curl phase plot for the low `L1` window
- restricted transverse spectrum for each `(n, variant)`
- lowest restricted eigenvalue versus lattice size

## Direct result

- observation: the lowest restricted transverse eigenvalue decreases with lattice size in all three substrate branches, and both defect branches lie below the baseline torus at every tested n
- conclusion: the restricted transverse sector develops a descending low band under the tested puncture and line-defect substrate modifications

## Artifacts

- results: `data/20260310_130326_L1_defect_scan.json`
- plots: `plots/20260310_130326_divergence_curl_phase.png`, `plots/20260310_130326_restricted_transverse_spectrum.png`, `plots/20260310_130326_restricted_eigenvalue_vs_n.png`
- timestamp: `20260310_130326`
