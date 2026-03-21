# Phase V Authoritative Summary

Authority freeze timestamp: `2026-03-21T14:23:11Z`.

## Objective

Freeze Phase V as a deterministic recovery-histogram readout over the frozen Phase IV stable sector bundle, with bounded reporting of reproducibility, robustness, and control discrimination.

## Frozen Input Schema

```json
{
  "sector_identifier": "...",
  "graph_reference_or_seed": {},
  "kernel_configuration_snapshot": {},
  "selector_configuration": {},
  "invariant_baseline_vector": {},
  "perturbation_policy_descriptor": {}
}
```

The authoritative frozen run uses the current stable Phase IV bundle and rejects oversized upstream payloads before compression into this schema.

## Recovery Score Definition

The Phase V recovery score is the arithmetic mean of five normalized invariant-recovery components: `law_fit_score`, `stable_subset_size`, `stable_rho`, `ordering_gap`, and `monotonic_consistency_score`. Each component is clamped to `[0, 1]`. No alternate score is introduced in the authority freeze.

## Run Modes

- `frozen_sector`: authoritative mode over the frozen Phase IV stable sector bundle.
- `dummy_sector`: synthetic control mode retained for scaffolding and smoke tests.
- Control classes: `stable_frozen_sector`, `degraded_sector_control`, `shuffled_null_control`.

## Reproducibility Results

- Stable frozen sector was swept across `5` deterministic schedules.
- Dominant stable class: `linear_like`.
- Class consistency fraction: `1.000000`.
- Mean recovery score average/range: `0.710931` / `0.006224`.
- Pairwise histogram L1 average/max: `0.059765` / `0.087890`.

## Robustness Results

- Baseline stable class: `linear_like`.
- Overall classifier stability fraction across narrow policy variation: `0.888889`.
- Max histogram L1 from baseline: `0.150390`.
- Max mean recovery score span across the tested band: `0.025673`.
- Tight-cluster test: `True`.
- Recorded edge-case classifier flip:
- `amplitude_jitter = 0.14` produced `quadratic_like` with mean recovery `0.704820` and histogram L1 from baseline `0.095705`.

## Control Separation Results

- Stable aggregate: mean `0.710931`, dominant class `linear_like`.
- Degraded aggregate: mean `0.375590`, dominant class `unstable`, observed classes `mixed, unstable`.
- Null aggregate: mean `0.061472`, dominant class `unstable`.
- Stable vs degraded: histogram L1 `1.553125`, mean-score gap `0.335341`.
- Stable vs null: histogram L1 `2.000000`, mean-score gap `0.649459`.
- Degraded vs null: histogram L1 `1.861328`, mean-score gap `0.314118`.

## Bounded Interpretation

- Within the tested configuration band, the stable frozen sector remains a high-recovery regime with a dominant `linear_like` label.
- Degraded and null controls remain clearly separated in both recovery level and distributional distance.
- The recorded amplitude-jitter edge flip is treated as a boundary marker for the descriptive classifier, not as evidence for a new mechanism or theory layer.

## Explicit Non-Claims

- No claim is made that the descriptive scaling classifier is invariant outside the tested policy band.
- No claim is made that the observed labels establish a universal recovery law.
- No claim is made that degraded or null controls encode physical mechanisms beyond their role as structural controls.
- No claim is made beyond the frozen Phase IV bundle, the configured perturbation schedules, and the current runner implementation.

## Authority Artifacts

- Baseline frozen run: `linear_like` with mean recovery `0.713270`.
- Manifest: `phase5-readout/phase5_authoritative_manifest.json`.
- Runs ledger: `phase5-readout/runs/phase5_runs_ledger_latest.json`.
- Figures: `phase5-readout/phase5_summary_figures/phase5_histogram_controls.svg`, `phase5-readout/phase5_summary_figures/phase5_reproducibility_sweep_summary.svg`, `phase5-readout/phase5_summary_figures/phase5_robustness_variation_summary.svg`.

Claim boundary: Phase V authority is limited to deterministic recovery-histogram behavior on the frozen Phase IV stable sector bundle under the configured perturbation schedules, policy bands, and control classes. It does not claim universal scaling, theory expansion, or classifier invariance outside this tested regime.
