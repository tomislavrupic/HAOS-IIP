# HAOS_IIP_GitHub_Main_47_2_Comparison_and_Local_Release_50_1_v1

Status:

- release-comparison note
- purpose: compare the live GitHub `main` state against the current local HAOS-IIP worktree before freezing release `50.1`

## 1. Live GitHub baseline

As checked against the live remote:

- repository: `tomislavrupic/HAOS-IIP`
- branch: `main`
- live commit: `cf5d48fb2178b92879a98dcdf0d52a4ec90ce60e`
- live release headline: `Freeze geometry_emergence release 47.2`

Interpretation:

- the public GitHub branch is still frozen at the geometry-emergence tranche
- none of the continuum-scaling, HAOS-to-harmonic derivation, or active-sector transport work is live on GitHub yet

## 2. Local delta relative to live GitHub

The current local worktree contains one coherent release-sized tranche above `47.2`.

### 2.1 Continuum-scaling and operator updates

Main additions:

- strengthened scalar operator-convergence runner with perturbed point-set and boundary-condition controls
- updated scalar continuum note and latest scalar result payload
- expanded transverse comparison with explicit harmonic / coexact sector readout
- sector-window and direct-gap scan runners plus corresponding notes and result bundles

Primary artifacts:

- `numerics/simulations/kernel_operator_convergence.py`
- `numerics/simulations/transverse_continuum_comparison.py`
- `numerics/simulations/transverse_sector_window_scan.py`
- `numerics/simulations/transverse_sector_gap_scaling.py`
- `experiments/operators/Kernel_Laplacian_Convergence_v1.md`
- `experiments/vector_sector/Transverse_Continuum_Comparison_v1.md`
- `experiments/vector_sector/Transverse_Sector_Window_Scan_v1.md`

### 2.2 HAOS-to-harmonic foundations stack

Main additions:

- bounded derivation-program spine from `H1-H8` through `T1-T8`
- standalone executed notes for `T1`, `T2`, `T4`, `T5`, `T6`, `T7`, and `T8`
- bounded scalar-locality cleanup `F1`
- scalar toy bridge note linking the abstract quadratic step to the first Dirichlet / Hodge lift

Primary artifacts:

- `docs/notes/foundations/HAOS_to_Harmonic_Operator_Derivation_Program_v1.md`
- `docs/notes/foundations/HAOS_Quadratic_Defect_Theorem_T1_v1.md`
- `docs/notes/foundations/HAOS_Relational_Dirichlet_Theorem_T2_v1.md`
- `docs/notes/foundations/HAOS_Incidence_Theorem_T4_v1.md`
- `docs/notes/foundations/HAOS_Hodge_Defect_Theorem_T5_v1.md`
- `docs/notes/foundations/HAOS_Sector_Decomposition_Theorem_T6_v1.md`
- `docs/notes/foundations/HAOS_Physical_Sector_Theorem_T7_v1.md`
- `docs/notes/foundations/HAOS_Continuum_Stability_Theorem_T8_v1.md`
- `docs/notes/foundations/HAOS_Strict_Scalar_Locality_From_Local_Failure_Assembly_F1_v1.md`

### 2.3 Physical-sector transport execution tranche

Main additions:

- explicit active-sector comparison map `J_n`
- restricted coexact transport runner and first bounded transport note
- mild-disorder extension through `n = 28`
- execution-program note translating the future-work section into bounded work packages
- disorder-threshold sweep runner added but not yet frozen as a completed result bundle

Primary artifacts:

- `numerics/simulations/transverse_active_sector_transport.py`
- `numerics/simulations/transverse_active_sector_disorder_sweep.py`
- `experiments/vector_sector/Transverse_Active_Sector_Transport_v1.md`
- `docs/notes/foundations/HAOS_Continuum_Physics_Execution_Program_v1.md`

Current bounded read:

- the clean torus remains the weak branch-identity case under refinement
- mild disorder stabilizes the same transported coexact family strongly under the fixed `J_n` map
- the threshold sweep has begun, but only the baseline endpoint was confirmed before the run was intentionally interrupted

## 3. Frozen paper delta relative to live GitHub

Two numbered papers exist locally but are not live on GitHub:

1. `48.1` Scalar Continuum Control, Harmonic-Coexact Sector Separation, and a Bounded Roadmap Toward Continuum Physics on HAOS-IIP
2. `49.1` From Recoverable Coherence to Harmonic Operator Structure: A Bounded Derivation Program, Physical-Sector Restriction, and Continuum-Stability Roadmap for HAOS-IIP

This release adds one more paper to package the repository state itself:

3. `50.1` GitHub-State Comparison and Continuum-Foundations Synchronization Release for HAOS-IIP

## 4. Release scope for `50.1`

Included in the intended release scope:

- HAOS-IIP numerics, notes, experiment writeups, data bundles, plots, paper sources, compiled PDFs, and release-index updates associated with the continuum / derivation / transport tranche

Explicitly excluded from the intended release scope:

- unrelated `Images/` assets
- any unfinished disorder-threshold result bundle that was not cleanly stamped to `data/`, `plots/`, and `EXPERIMENT_LOG.md`

## 5. Current authority boundary

After comparison against live GitHub, the honest repository state is:

- scalar continuum control is materially stronger than the live `47.2` branch suggests
- the HAOS-to-harmonic derivation ladder now exists as a bounded theorem stack instead of only a conceptual bridge
- active-sector transport is real and executable, but vector-sector continuum closure is still not complete
- the disorder-threshold question is active work, not a frozen claim

## 6. Recommended publish action

Publish the HAOS-IIP tranche as a single synchronization release that:

1. lands papers `48.1`, `49.1`, and `50.1`
2. lands the foundations notes and execution-program note
3. lands the active-sector transport numerics and bounded result bundles
4. excludes unrelated image assets

That action would move GitHub from the frozen `47.2` geometry-emergence state to the present continuum / harmonic-derivation / active-sector-transport state without overstating unfinished disorder-threshold work.
