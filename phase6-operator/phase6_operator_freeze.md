# Phase VI Operator Freeze

Freeze timestamp: `2026-03-21T14:41:26Z`.

## Objective

Freeze one branch-valid operator family across a deterministic refinement hierarchy and test whether later spectral work is well-posed on the frozen branch.

## Selected Operator Class

- Selected class: `cochain_laplacian`.
- Selection rationale: it is the simplest branch-local family already exposed by the frozen DK2D complex, with direct symmetry structure, direct semiboundedness diagnostics, and no additional auxiliary machinery.
- Not selected: weighted graph Laplacian, because it would discard the higher-degree cochain structure already frozen in the branch.
- Not selected: discrete Dirac / Hodge-Dirac as the primary family, because the squared cochain Laplacian provides the clearer semibounded feasibility test at this stage.

## Frozen Inputs

- Phase IV bundle: `phase4-sector-freeze/runs/phase4_sector_freeze_bundle_latest.json`.
- Phase IV config: `phase4-sector-freeze/configs/stage24_4_stable_smeared_law_and_ordering_extension_boundary_runs.json`.
- Phase V authority ledger: `phase5-readout/runs/phase5_runs_ledger_latest.json`.
- Frozen sector identifier: `phaseIV_stable_smeared_sector`.
- Branch operator sector: `dirac_kahler_minimal_0_1_form` with `periodic` `regular_lattice` and `static_gaussian` construction discipline.

## Refinement Hierarchy

- Refinement parameter: `h = 1 / n_side`.
- Interpretation of `h`: lattice-spacing surrogate on the unit-periodic branch-local square cell complex.
- Base resolution: `12`.
- Refinement multipliers: `[1, 2, 3, 4]`.
- Deterministic rule: for each multiplier `m`, set `n_side = base_resolution * m`, build the periodic DK2D complex with frozen `base_epsilon`, and define `O_h = delta_h`.

## Admissibility Diagnostics

| level | n_side | h | dim | nnz | symmetry defect | min eigen probe | first positive | spectral radius | nullspace | reorder defect | admissible |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| R1 | 12 | 0.083333333333 | 576 | 2880 | 0.000000000000 | -0.000000000000 | 0.263337445092 | 7.862309796963 | 4 | 0.000000000000 | True |
| R2 | 24 | 0.041666666667 | 2304 | 11520 | 0.000000000000 | -0.000000000000 | 0.067853205626 | 7.965353020924 | 4 | 0.000000000000 | True |
| R3 | 36 | 0.027777777778 | 5184 | 25920 | 0.000000000000 | -0.000000000000 | 0.030325938407 | 7.984582776023 | 4 | 0.000000000000 | True |
| R4 | 48 | 0.020833333333 | 9216 | 46080 | 0.000000000000 | -0.000000000000 | 0.017091721482 | 7.991324152244 | 4 | 0.000000000000 | True |

All diagnostics are computed, not asserted. Hermiticity defect is the relative Frobenius defect of `O_h - O_h^T`; semiboundedness is probed by sampled smallest eigenvalues and random Rayleigh quotients; node reordering consistency is tested by a deterministic torus translation of the full cochain basis followed by restoration to the original order.

## Basic Spectral Sanity Probes

- All levels admissible: `True`.
- Nullspace estimate values: `[4, 4, 4, 4]`.
- Spectral radius values: `[7.862309796963, 7.965353020924, 7.984582776023, 7.991324152244]`.
- First positive eigenvalue values: `[0.263337445092, 0.067853205626, 0.030325938407, 0.017091721482]`.
- Spectral radius non-decreasing: `True`.
- Node reordering consistent: `True`.

## Bounded Interpretation

- The frozen operator family is branch-local only.
- The diagnostics support self-adjointness, semiboundedness, stable kernel estimation, and deterministic refinement construction on the tested hierarchy.
- This is an operator-feasibility result, not a continuum or geometry claim.

## Explicit Non-Claims

- No continuum Laplacian is claimed.
- No geometric coefficients are computed.
- No heat-kernel asymptotics are computed.
- No spectral invariants are claimed beyond the finite diagnostic samples listed here.
- No continuum or GR correspondence is asserted.

Phase VI now has a frozen operator family suitable for spectral feasibility testing.
