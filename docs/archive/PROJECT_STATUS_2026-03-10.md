# PROJECT_STATUS

Date: March 10, 2026

## Status summary

The repository has been restructured into a research layout with:

- canonical theory files
- classified legacy source documents
- initialized Git repository
- reusable starter numerics
- reproducible orchestration layer
- one explicit 3D eigenmode experiment preserved under `experiments/`
- a periodic/twisted `L1` gauge-sector experiment with coexact projection

The theory architecture has now been tightened further:

- same substrate
- same kernel
- three operator sectors

This replaces the earlier compressed reading that one node Laplacian might explain all sectors directly.

## Current gauge-sector verdict

The present HAOS/IIP edge operator supports a genuine non-scalar flux-responsive branch, but explicit harmonic-vs-coexact separation shows that the lowest modes remain harmonic or mixed topological structure rather than a low coexact vector band. In the current periodic/twisted range, this is a negative Maxwell test.

- `n = 5`, `m = 1`: full lowest `0.013057`, harmonic-projected `0.509808`, coexact-projected `0.458457`
- `n = 4`, `m = 1`: full lowest `0.014110`, harmonic-projected `0.611897`, coexact-projected `0.679582`

## FILES_MOVED

### Foundations -> `docs/notes/foundations/`

- `HAOS_Frozen_Specification_v1.docx`
- `HAOS_IIP_Minimal_Theory_Statement_v2.docx`
- `Interaction_Invariant_Physics_Full_v1.docx`

### Kernel -> `theory/kernel/`

- `HAOS_IIP_Kernel_v1.docx`
- `HAOS_Spectral_Recoverability_Principle_v1.docx`
- `Discrete_Scale_Invariance_and_Coherence_Stability_v1.docx`

### Geometry -> `theory/geometry/`

- `Interaction_Kernels_to_Emergent_Geometry_v1.docx`
- `HAOS_IIP_Spectral_Geometry_Note_v1.docx`
- `HAOS_IIP_Spectral_Gravity_Note_v1.docx`
- `HAOS_IIP_Singularities_v1.docx`

### Gauge -> `theory/gauge/`

- `HAOS_IIP_Emergent_Gauge_Sector_Note_v1.md`

### Particles -> `theory/particles/`

- `Rebuilding_Dudas_Picture_in_HAOS_IIP_v1.docx`

### Applications -> `docs/notes/applications/`

- `HAOS_IIP_Microtubules_v1.docx`
- `HAOS_IIP_Fluid_Derivation_v1.docx`
- `The_Meissner_Effect_v1.docx`

### Archive -> `docs/archive/`

- `Interaction_Invariant_Foundations_of_Consciousness_v1.docx`
- `Mathematical_Foundations_of_QATC_v1.docx`
- `QATC_2026_Creative_Dynamics_v1.docx`
- `The_Physics_of_Imagination_v1.docx`
- `ALi_22_v1.docx`
- `HAOPS_v1.docx`
- `STATEFLOW_v1.docx`
- `The_Hierarchical_Autonomy_of_Reality_v1.docx`

### Experiments -> `experiments/eigenmodes/haos_iip_3d_low_mode_study/`

- `HAOS_IIP_3D_Low_Mode_Study_v1.md`
- `haos_iip_3d_low_modes_v1.py`
- `haos_iip_results_v1.json`
- `haos_iip_mode_plots/`

## FILES_VERSIONED

- existing active notes were given `_v1` or `_v2` suffixes
- the minimal theory statement retained `_v2` because the source title already marked it as a later-stage document
- canonical files were created separately under `docs/canon/` rather than renaming source drafts into canon

## TERMINOLOGY_CONFLICTS

Primary conflicts found in legacy material:

1. `epsilon`

- used as kernel width in geometry/kernel notes
- also used as coherence threshold in recoverability inequalities
- canonical standard:
  - `epsilon_k`: kernel width
  - `epsilon_c`: coherence threshold

2. `L`

- used as graph Laplacian
- also loosely used in some drafts for Lagrangian-like objects
- canonical standard:
  - `L0`: node Laplacian
  - `L1`: edge/Hodge Laplacian
  - `L`: generic Laplacian only when sector is obvious
  - `mathcal{L}`: continuum operator

3. charge language

- some drafts treat charge as winding, some as coupling asymmetry, some as placeholder gauge quantity
- canonical status:
  - winding/defect charge is `[P]`
  - no physical charge derivation is `[E]`

## CANON_CREATED

Created under `docs/canon/`:

- `HAOS_IIP_CORE_THEORY.md`
- `KERNEL_DEFINITION.md`
- `OPERATOR_STRUCTURE.md`
- `EMERGENT_GEOMETRY.md`
- `GAUGE_PROGRAM.md`
- `OPEN_PROBLEMS.md`

## NUMERICS_INITIALIZED

Created under `numerics/simulations/`:

- `laplacian_modes.py`
- `gauge_modes.py`
- `hodge_modes.py`
- `parameter_sweep.py`
- `periodic_twisted_l1.py`

These are small reusable scripts, separate from experiment-specific code.

Added orchestration layer:

- `scripts/run_experiments.py`
- `config.json`
- `experiments/EXPERIMENT_LOG.md`
- `Makefile`
- `numerics/simulations/hodge_modes.py`
- `numerics/simulations/parameter_sweep.py`

The repository now supports a single command workflow:

```bash
python3 scripts/run_experiments.py
```

or

```bash
make run
```

This loads shared parameters, runs the starter simulations, stores structured JSON results in `data/`, writes plots into `plots/`, and appends a dated run record to `experiments/EXPERIMENT_LOG.md`.

## NEXT_RESEARCH_TARGET

The next technical target should be:

1. introduce punctures or controlled defects to separate torus-cycle modes from local transverse circulation
2. extend the periodic scan to larger lattices and check whether the coexact floor moves down toward a cleaner low band
3. only after that, build the minimal Dirac-type branch `D_H`
