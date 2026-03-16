# Stage 23.10 -- DK Braid Topology Mechanism Isolation Probe

## Purpose

Classify the braid-like exchange window seen in the Dirac-Kahler branch as either:

- an intrinsic exchange mechanism of the frozen operator-kernel architecture
- or a clustered seed texture that does not survive controlled symmetry destruction

This stage is diagnostic, not exploratory.
It is meant to decide the status of the Phase III braid line, not to open a broader atlas.

## Scientific question

On the frozen HAOS-IIP Dirac-Kahler packet architecture, does braid-like exchange persist once clustered symmetry is deliberately degraded while the operator class remains fixed?

Equivalently:

- if the braid is intrinsic, it should reappear outside the tight clustered control and remain stable under refinement and mild phase reseeding
- if the braid is seed-dependent, it should stay confined to the clustered control or collapse under small geometric damage

## Frozen architecture

Do not modify:

- DK grading rules
- Gaussian kernel class
- interaction-invariance update
- timestep policy
- global normalization

Allowed variation only in:

- initial packet geometry
- mild support-skew / degree-spectrum surrogate
- motif structure
- weak local anisotropy

## Minimal decisive lattice

Run exactly 9 primary families:

1. clustered baseline control
2. clustered degraded
3. clustered topology-broken
4. corridor baseline
5. corridor motif-injected
6. corridor anisotropic
7. triad baseline
8. triad degree-skewed
9. distributed random seed

Use one common operator parameter set for all runs.

## Required observables

Per primary run compute:

- `topology_class`
- `braid_survival_time`
- `flow_concentration_index`
- `grade_transfer_asymmetry`
- `refinement_topology_stability`
- `topology_repeatability_score`

Allowed topology labels:

- `braid_like_exchange`
- `transfer_smeared`
- `dispersive_pass`
- `localized_capture`

## Stability checks

For every primary family:

- rerun at doubled resolution `n -> 2n`
- rerun the base geometry under small phase jitter

This stage does not claim persistence law.
It only classifies whether the braid window is architectural or clustered.

## Decision rule

Classify as `intrinsic_exchange_mechanism` only if:

- braid topology appears in at least 2 non-clustered families
- those braid hits satisfy `refinement_topology_stability <= 1.1`
- and `topology_repeatability_score >= 0.7`

Classify as `clustered_texture_artefact` if:

- braid remains confined to clustered controls
- or the clustered braid collapses under mild geometry degradation

Otherwise return `mechanism_boundary_inconclusive`.

## Required outputs

Per-run:

- stamped JSON
- stamped CSV
- trajectory visualization

Summary:

- topology survival matrix
- family-comparison persistence panel
- mechanism classification note

## Interpretation rule

Do not claim persistence, binding, exclusion law, or emergent particles.

The only allowed end statement is one of:

- `intrinsic_exchange_mechanism`
- `clustered_texture_artefact`
- `mechanism_boundary_inconclusive`
