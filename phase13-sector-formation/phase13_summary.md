# Phase XIII - Multi-Mode Statistical Sector Formation

## Objective

Test whether dilute populations of the frozen persistent localized candidate preserve approximate identity, develop reproducible spacing structure, and remain distinguishable from the deterministic control hierarchy.

## Frozen Inputs

- Candidate: `low_mode_localized_wavepacket` from Phase XI and Phase XII
- Refinement levels: `[48, 60, 72, 84]` with `h = 1 / n_side`
- Ensemble sizes: `[3, 5, 7]`
- Layout seeds: `[1301, 1302, 1303]` with minimum physical separation `0.12`
- Evaluation tau: `0.8`
- Persistence tau grid: `[0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 3.6, 4.2]`
- Phase VIII short-time window retained by contract: `[0.02, 0.03, 0.05, 0.075, 0.1]`

## Key Results

- Mean branch evaluation survival fraction: `1.0`.
- Mean control evaluation survival fraction: `0.990476190476`.
- Mean branch evaluation cluster frequency: `0.238095238095`.
- Mean control evaluation cluster frequency: `0.228571428571`.
- Mean branch evaluation spectral distortion: `-0.2213869077`.
- Mean control evaluation spectral distortion: `-0.232292030413`.
- Mean branch nearest-neighbor scale: `0.263125976359`.
- Mean control nearest-neighbor scale: `0.2636216034`.
- Branch seed std (survival, spacing, cluster): `0.0`, `0.100627125854`, `0.249443825785`.
- Branch successive-refinement distances (survival, spacing): `0.0`, `0.0`.
- Branch refinement drift within fixed seeded layouts (spacing, cluster, exclusion, pair-hist L1): `0.0`, `0.0`, `0.0`, `0.0`.
- Maximum branch-control gaps across the tau grid (survival, spectral distortion): `0.071164021164`, `0.046202260651`.
- Terminal tau `4.2` survival fractions (branch, control): `0.957936507937`, `0.886772486773`.

## Bounded Interpretation

The branch is treated as feasible only if each deterministic layout family remains compact across refinement, and if the same seeded layouts on the control hierarchy separate under time-resolved survival or spectral occupancy diagnostics. Seed-to-seed spacing variation is treated as deterministic geometry variation, not as refinement instability. These results do not assert thermodynamics, particle ensembles, or kinetic laws.

Phase XIII establishes multi-mode statistical sector formation feasibility for the frozen operator hierarchy.
