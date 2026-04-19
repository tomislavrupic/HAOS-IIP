# HAOS-IIP

Frozen branch-local cochain-Laplacian hierarchy: bounded emergence diagnostics and HAOS-to-harmonic derivation program (`T1-T8` closed).

Success criterion: running the public reproduction script should produce identical table values and an identical plot to the frozen baseline.

<p align="center">
  <img src="Images/HAOS-IIP%20Continuum-Bridge%20Feasibility%20Path.png" alt="HAOS-IIP continuum-bridge feasibility path" width="100%" />
</p>

## Overview

HAOS-IIP is a structured numerical research program investigating whether a reproducible emergence ladder can be identified inside a frozen operator hierarchy on discrete graph domains.

The project focuses on:

- bounded feasibility diagnostics
- controlled refinement hierarchies
- frozen telemetry definitions
- reproducible phase bundles

This repository does not assert geometric structure, physical correspondence, or continuum validity. It provides a computational environment for testing whether such structures could stabilize under disciplined constraints.

Public-scope boundary: this reproduces a bounded mesoscopic-to-proto-geometric feasibility arc. It does not claim a continuum limit or physical ontology.

## Visual Entry Points

The landing page uses only a small curated media set so a cold reader can orient quickly without getting buried in assets.

- [Deterministic validation hierarchy, Phases V-X](Images/HAOS-IIP.png)
- [Emergent feasibility hierarchy, Phases X-XVIII](Images/Hierarchy%20of%20Emergent%20Feasibility%20Phases.png)
- [Phase V diagnostic data readout](Images/HAOS-IIP%20Phase%20V%20Diagnostic%20Data%20Readout.png)
- [Morphological stability network design](Images/Morphological%20Stability%20Network%20Design.png)
- [Proto-particle feasibility study diagram](Images/Proto-Particle%20Feasibility%20Study%20Diagram.png)
- [HAOS-IIP system ledger](Images/HAOS-IIP_System_Ledger.pdf)
- [HAOS-IIP emergence telemetry](Images/HAOS-IIP_EMERGENCE_TELEMETRY.pdf)
- [System Check: The Cautious Bridge](Images/System_Check__The_Cautious_Bridge.mp4)
- [Media index](docs/media/MEDIA_INDEX.md)
- [Image vault](Images/README.md)

## Scope Of The Program

The program studies emergence inside a frozen branch regime defined by:

- a branch-local cochain-Laplacian operator family
- deterministic initialization routines
- matched altered-connectivity controls
- shared telemetry metrics across phases

Primary research questions include:

- Does a stable propagation band form under refinement?
- Can temporal ordering emerge from branch-internal dynamics?
- Does a consistent causal influence structure close?
- Can proto-geometric distance surrogates stabilize?

All conclusions are framed as feasibility statements, not ontological claims.

## Repository Structure

```text
phase3-stability/ ... phase18-distance-surrogate/  frozen phase bundles and diagnostics
continuum-sketch/                                  low-cost post-processing scaling test
examples/                                          one-command public reproduction spine
telemetry/                                         frozen emergence metrics
haos_core/                                         shared core primitives
papers/                                            numbered PDF releases and drafts
ai/                                                saved Codex prompts and workflow support
```

Each phase bundle follows a uniform contract:

```text
build_*.py
phase*_manifest.json
phase*_summary.md
runs/
diagnostics/
(optional) plots/
```

## Quick Start

For a cold external read, start with the narrow public reproduction path:

```bash
python3 examples/quick_reproduce.py
```

This one command:

- checks the frozen Phase XV-XVIII bundles
- rebuilds the artifact-only continuum sketch
- writes a compact public table and figure under `examples/output/`
- verifies both against `examples/expected_output.csv` and `examples/expected_plot.svg`

If you only want to validate a single frozen bundle, run an integrity check:

```bash
python3 run_phase.py 18 --check
```

Example:

```bash
python3 run_phase.py 10 --check
```

No new simulations are required to validate frozen results.

## Telemetry Demo

HAOS-IIP also exposes a small telemetry-based structural-stability diagnostic. This demo illustrates how HAOS-IIP telemetry can be used as a structural-stability diagnostic for evolving graph-like systems. It is not a claim about physical spacetime emergence.

Run the default demo bundle:

```bash
python3 -m haos_iip.demo stability
```

Request a machine-readable payload for one scenario:

```bash
python3 -m haos_iip.demo stability baseline --json
```

Run a deterministic scan grid:

```bash
python3 -m haos_iip.demo stability --scan noise=0.00:0.10:0.05 connectivity=0.00:0.20:0.10
```

Run a deterministic generated variant:

```bash
python3 -m haos_iip.demo stability baseline --noise 0.05 --connectivity-drop 0.2 --cluster-split
```

Public metric aliases:

- `structural_retention` = persistence
- `temporal_consistency` = ordering
- `causal_deformation` = depth drift
- `geometric_integrity` = distance coherence

![Trajectory stability demo](stability_demo_overview.svg)

Short note:

- [A Minimal Structural-Stability Oracle Based on Frozen HAOS-IIP Telemetry](docs/notes/applications/A_Minimal_Structural_Stability_Oracle_Based_on_Frozen_HAOS_IIP_Telemetry_v1.md)

## Program Status

**Current milestone (50.2 Foundations Cleanup):**  
The bounded HAOS-to-harmonic derivation program is fully executed.  
`T1 -> T2 -> T4 -> T5 -> T6 -> T7 -> T8` ladder plus `F2-F4` closures are now complete and frozen in the public repository.

The frozen branch-local cochain-Laplacian hierarchy now supports a reproducible ladder from recoverable coherence to harmonic operator structure, including Laplacian form, Hodge-type defect structure, sector decomposition, physical-sector restriction, and refinement-stability criteria under explicit HAOS-aligned constraints.

Milestone anchors:

- [50.2 Foundations-Cleanup Synchronization Release](papers/pdf_releases/50.2%20Foundations-Cleanup%20Synchronization%20Release%20for%20HAOS-IIP.pdf)
- [49.1 From Recoverable Coherence to Harmonic Operator Structure](papers/pdf_releases/49.1%20From%20Recoverable%20Coherence%20to%20Harmonic%20Operator%20Structure%3A%20A%20Bounded%20Derivation%20Program%2C%20Physical-Sector%20Restriction%2C%20and%20Continuum-Stability%20Roadmap%20for%20HAOS-IIP.pdf)
- [PROGRAM_STATE_MILESTONE_18.md](PROGRAM_STATE_MILESTONE_18.md) (original geometry-emergence baseline)

## Reproducibility Contract

Core experimental layers are frozen and defined in:

- [API_CONTRACT.md](API_CONTRACT.md)
- [telemetry/frozen_metrics.py](telemetry/frozen_metrics.py)

This guarantees:

- stable operator definitions
- consistent initialization rules
- invariant telemetry formulas
- matched control construction

Emergence diagnostics rely only on these frozen interfaces.

## Continuum-Sketch Layer

The bounded continuum-facing tranche, including scalar convergence controls, harmonic/coexact sector separation, and active-sector transport machinery, is now frozen as part of the `50.1/50.2` synchronization releases.

A minimal post-processing protocol for low-cost scaling inspection is provided in:

- [continuum-sketch/](continuum-sketch/)

This stage performs:

- refinement trend checks
- propagation-band stability inspection
- proto-distance scaling diagnostics

It does not introduce new dynamics.

## Papers

Latest releases:

- [50.2 Foundations-Cleanup Synchronization Release](papers/pdf_releases/50.2%20Foundations-Cleanup%20Synchronization%20Release%20for%20HAOS-IIP.pdf)
- [50.1 GitHub-State Comparison and Continuum-Foundations Synchronization Release](papers/pdf_releases/50.1%20GitHub-State%20Comparison%20and%20Continuum-Foundations%20Synchronization%20Release%20for%20HAOS-IIP.pdf)
- [49.1 From Recoverable Coherence to Harmonic Operator Structure](papers/pdf_releases/49.1%20From%20Recoverable%20Coherence%20to%20Harmonic%20Operator%20Structure%3A%20A%20Bounded%20Derivation%20Program%2C%20Physical-Sector%20Restriction%2C%20and%20Continuum-Stability%20Roadmap%20for%20HAOS-IIP.pdf)

Numbered synthesis papers are released in:

- [papers/pdf_releases/](papers/pdf_releases/)
- [papers/pdf_releases/INDEX.md](papers/pdf_releases/INDEX.md)

Supporting repository visuals and media are indexed in:

- [docs/media/MEDIA_INDEX.md](docs/media/MEDIA_INDEX.md)
- [Images/README.md](Images/README.md)

Earlier milestone example:

- [43.1 Ordering, Causal Closure, and Proto-Geometric Distance-Surrogate Feasibility on a Frozen Branch-Local Cochain-Laplacian Hierarchy.pdf](papers/pdf_releases/43.1%20Ordering,%20Causal%20Closure,%20and%20Proto-Geometric%20Distance-Surrogate%20Feasibility%20on%20a%20Frozen%20Branch-Local%20Cochain-Laplacian%20Hierarchy.pdf)

## Limitations

This repository intentionally avoids:

- physical interpretation claims
- geometric reconstruction claims
- continuum field-theory derivations
- cosmological or ontological assertions

The program is strictly a numerical emergence feasibility study.

## Citation

If you use results or methods from this repository, please cite this work using:

- [CITATION.cff](CITATION.cff)

## License

See:

- [LICENSE](LICENSE)
- [COPYRIGHT.md](COPYRIGHT.md)
- [THEORY_NOTICE.md](THEORY_NOTICE.md)

## Author

Tomislav Rupić  
Independent multimedia researcher and computational emergence practitioner.
Website: [tomislav-rupic.com](https://tomislav-rupic.com)  
Email: [tom.d.vox@gmail.com](mailto:tom.d.vox@gmail.com)
