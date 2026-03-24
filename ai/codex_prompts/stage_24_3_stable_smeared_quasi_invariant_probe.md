# Stage 24.3 / IV-A3 - Stable-Smeared Quasi-Invariant Probe

## Purpose
Probe whether the law-defined stable smeared sector on the frozen clustered Dirac-Kahler branch carries any branch-stable quasi-invariant, bounded composite, or preserved ordering relation.

This is a reduction-only stage. It must not reopen braid, enlarge the observable core, or introduce a new operator branch.

## Authoritative inputs
- `23.12 / III-C3`: one-sided smeared-dominant phase closure
- `23.14 / III-C5`: frozen Phase III claim ledger
- `24.1 / IV-A1`: accepted smeared-sector observable core
- `24.2 / IV-A2`: compact effective law

## Positive and contrast sets
- positive set: rows on the Stage 23.12 closure lattice that satisfy both the Stage 24.2 flow-threshold law and the frozen `stable_closed_smeared` truth label
- contrast set: near-threshold failures and transient mixed controls above the threshold boundary

## Allowed observables
Primary/core only:
- `topology_survival_time`
- `refinement_stability_flag`
- `weak_coupling_stability_flag`
- `reverse_stability_flag`
- `flow_concentration_index`

Context only:
- `phase_corridor_position`
- `phase_corridor_width`

Forbidden:
- `transfer_asymmetry`
- `coherence_decay_rate`
- `topology_return_error`
- `braid_like_exchange_indicator`
- `mechanism_sensitivity_score`

## Candidate families
At minimum compare:
- exact closure-support signatures built from the three accepted stability flags
- bounded single-observable envelopes on accepted core quantities
- branch-local composite quantities built only from accepted core quantities
- preserved ordering relations built only from accepted core quantities

Direct restatements of the Stage 24.2 threshold law may be reported, but they may not be selected as the final quasi-invariant read.

## Success criterion
The stage succeeds only if at least one candidate quantity or relation is significantly more stable inside the law-defined stable smeared sector than in the contrast controls, while using only accepted core observables and adding structure beyond the Stage 24.2 threshold rule.

## Required outputs
- `Stage_IV_A3_Stable_Smeared_Quasi_Invariant_Probe_v1.md`
- `stage24_3_stable_smeared_quasi_invariant_probe.py`
- `stage24_3_stable_smeared_quasi_invariant_probe_runs.json`
- stamped JSON summary
- stamped CSV candidate ledger
- summary visuals for candidate comparison, the selected relation, and stable-vs-contrast separation

## Interpretation boundary
Do not claim a universal invariant, conservation law, or continuum correspondence.

At most claim:
the law-defined stable smeared sector carries a branch-local quasi-invariant, bounded composite, or preserved ordering relation on the frozen clustered DK branch.
