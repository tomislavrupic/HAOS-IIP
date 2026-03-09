# HAOS-IIP

HAOS-IIP is a research repository for a weighted interaction-fabric program.

The mature architecture is:

```text
one substrate
one kernel
three operator sectors
```

Those three sectors are:

- `L0`: node Laplacian for scalar / geometry sector
- `L1`: edge / Hodge Laplacian for gauge / vector sector
- `D_H`: Dirac-type lift for fermion / spinor sector

## Canonical entry points

- [docs/canon/HAOS_IIP_CORE_THEORY.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/HAOS_IIP_CORE_THEORY.md)
- [docs/canon/KERNEL_DEFINITION.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/KERNEL_DEFINITION.md)
- [docs/canon/OPERATOR_STRUCTURE.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/OPERATOR_STRUCTURE.md)
- [docs/canon/EMERGENT_GEOMETRY.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/EMERGENT_GEOMETRY.md)
- [docs/canon/GAUGE_PROGRAM.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/GAUGE_PROGRAM.md)
- [docs/canon/OPEN_PROBLEMS.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/OPEN_PROBLEMS.md)

## Repository layout

- `docs/canon/`: stabilized technical statements
- `docs/notes/`: retained working theory and applications
- `docs/archive/`: meta, philosophical, and legacy material
- `theory/`: topic-grouped source material
- `numerics/`: reusable simulation starters and future notebooks
- `experiments/`: experiment-specific code, plots, and result notes
- `data/`: structured JSON results from orchestrated runs
- `plots/`: current and timestamped plot outputs
- `scripts/`: thin orchestration layer for reproducible runs
- `papers/`: future paper drafts and figures
- `ai/`: Codex prompts and agent definitions

## Current active experiment

- [experiments/eigenmodes/haos_iip_3d_low_mode_study/HAOS_IIP_3D_Low_Mode_Study_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/experiments/eigenmodes/haos_iip_3d_low_mode_study/HAOS_IIP_3D_Low_Mode_Study_v1.md)

It currently gives the strongest explicit numerical result in the repo:

- the bare 3D node Laplacian `L0` shows scalar geometry-like low modes
- a phase-dressed scalar branch shows connection sensitivity
- no autonomous gauge or fermion sector has been derived yet

## Current architecture

The repository no longer assumes that one node Laplacian explains every sector.

Instead:

```text
same substrate
same kernel
different operator lifts
```

## Working rules

- Do not delete source material.
- Preserve old versions with `_vN` suffixes.
- Use `_CANON` only for stabilized canonical files.
- Keep speculative claims out of canonical files unless they are labeled `[P]` or `[O]`.
- Put experiment-specific outputs under `experiments/`, not into theory folders.

## Quick start

Run the repo as a small experiment lab:

```bash
cd /Volumes/Samsung\ T5/2026/HAOS/HAOS\ DOCS/HAOS-IIP
python3 scripts/run_experiments.py
```

Or use:

```bash
cd /Volumes/Samsung\ T5/2026/HAOS/HAOS\ DOCS/HAOS-IIP
make run
```

This does the following:

- loads shared parameters from `config.json`
- runs the reusable scalar and gauge starter experiments
- runs the weighted edge/Hodge `L1` experiment
- runs a parameter sweep over `epsilon`, `N`, and substrate choice
- writes JSON outputs into `data/`
- writes plots into `plots/`
- appends a record to `experiments/EXPERIMENT_LOG.md`

Adjust experiment parameters in `config.json` rather than hard-coding them in the simulation modules.

## Direct module runs

Run the current 3D low-mode experiment:

```bash
cd /Volumes/Samsung\ T5/2026/HAOS/HAOS\ DOCS/HAOS-IIP/experiments/eigenmodes/haos_iip_3d_low_mode_study
python3 haos_iip_3d_low_modes_v1.py
```

Run the reusable starter simulations:

```bash
cd /Volumes/Samsung\ T5/2026/HAOS/HAOS\ DOCS/HAOS-IIP
python3 numerics/simulations/laplacian_modes.py
python3 numerics/simulations/gauge_modes.py
python3 numerics/simulations/hodge_modes.py
python3 numerics/simulations/parameter_sweep.py
```

## Next explicit target

Build the minimal Dirac-type branch `D_H` after the `L1` edge/Hodge branch has been mapped more cleanly across flux backgrounds and defect perturbations.
