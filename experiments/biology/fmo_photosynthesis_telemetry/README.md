# FMO Photosynthesis Telemetry

Phase 56 sidecar experiment for HAOS-IIP.

This folder tests whether the spectral-dynamics telemetry instrument developed in
Phase 55 can preserve recoverable site/pathway identity on a small FMO-inspired
weighted transfer graph.

The model uses a standard 7-site FMO-like Hamiltonian only as a deterministic
interaction topology. It is not a molecular simulation, exciton-transport model,
quantum-biology result, or claim about photosynthesis. The point is narrower:
test whether HAOS spectral dynamics and the level-5 null ladder remain specific
on a compact transfer network that is very different from the microtubule lattice.

Run from the repository root:

```bash
uv run --with numpy --with matplotlib python experiments/biology/fmo_photosynthesis_telemetry/run_fmo_spectral_telemetry.py
```

Outputs are written to `experiments/biology/fmo_photosynthesis_telemetry/outputs/`.

For a seed/noise robustness pass:

```bash
uv run --with numpy --with matplotlib python experiments/biology/fmo_photosynthesis_telemetry/run_fmo_robustness.py
```

The robustness runner writes CSV, JSON, Markdown, and a compact summary figure
under `outputs/robustness_56/`.

For the Phase 56.2 pathway-strengthening sweep:

```bash
uv run --with numpy --with matplotlib python experiments/biology/fmo_photosynthesis_telemetry/run_fmo_variant_sweep.py
```

This compares spectral, local, hybrid, sink-aware, and environment-assisted
address variants. The tracked summary is `fmo_pathway_strengthening_56_2.md`,
with a compact figure at `figures/fmo_variant_summary_56_2.png`.

For the Phase 56.3 pathway/flux modeling sweep:

```bash
uv run --with numpy --with matplotlib python experiments/biology/fmo_photosynthesis_telemetry/run_fmo_pathway_sweep.py
```

This adds explicit signed-flux restoration along the declared FMO pathway edges.
The tracked summary is `fmo_pathway_flux_56_3.md`, with a compact figure at
`figures/fmo_pathway_flux_summary_56_3.png`.

For the Phase 56.4 intrinsic pathway dynamics sweep:

```bash
uv run --with numpy --with matplotlib python experiments/biology/fmo_photosynthesis_telemetry/run_fmo_emergent_sweep.py
```

This uses weak directed/temporal pathway bias plus a level-6 directed trajectory
null. The tracked summary is `fmo_intrinsic_pathway_56_4.md`, with a compact
figure at `figures/fmo_intrinsic_pathway_summary_56_4.png`.

## Initial Result

The first 50-run robustness pass is intentionally kept as an honest falsifier:

- runs: `50`
- pass_rate: `0.040000`
- recoverability_mean: `0.803724`
- site_identity_mean: `0.984934`
- pathway_identity_mean: `0.546597`
- active_null_z_mean: `2.472024`
- control_pass_max: `2`

Interpretation: the compact FMO-like graph preserves site identity strongly, but
pathway identity is not yet specific enough under the level-5 null ladder.
Controls still pass often enough that Phase 56 should be treated as an initial
MARGINAL/FAIL diagnostic, not an external validation result like Phase 55.2.

See `fmo_robustness_56_initial.md` for the tracked summary.

## Phase 56.2 Result

The first pathway-strengthening sweep improved the best variant but did not clear
the target bar:

- variants: `6`
- total runs: `300`
- best variant: `sink_spectral`
- best pass_rate: `0.080000`
- best pathway_identity_mean: `0.622057`
- best recoverability_mean: `0.827528`
- best active_null_z_mean: `2.659739`
- best control_pass_max: `2`

Interpretation: adding a reaction-center sink improves pathway retention over
the spectral baseline (`0.622` vs `0.547`) but does not produce a robust PASS.
Local and hybrid address terms degrade the compact FMO transfer graph. The next
real improvement likely needs a better pathway/flux model, not more address
smoothing.

## Phase 56.3 Result

The explicit pathway/flux model solves the pathway-retention metric but still
does not solve strict control specificity:

- specs: `6`
- total runs: `300`
- best pathway spec: `flux_strong`
- best pathway_identity_mean: `0.906155`
- best recoverability_mean: `0.930908`
- best active_null_z_mean: `3.227554`
- best pass_rate: `0.060000` (`flux_light`)
- control_pass_max: `2`

Interpretation: explicit flux restoration can preserve the declared transfer
pathway, but the controls can also exploit that engineered pathway term. The FMO
bottleneck has therefore moved from pathway retention to branch-specific
control discrimination. This remains a diagnostic result, not a validation pass.

## Phase 56.4 Result

The intrinsic pathway dynamics plus level-6 directed trajectory null reduce the
control leak somewhat, but they do not preserve enough pathway identity:

- specs: `6`
- total runs: `300`
- best spec: `modes5_directed`
- best pass_rate: `0.160000`
- best pathway_identity_mean: `0.671540`
- best temporal_pathway_mean: `0.694755`
- best recoverability_mean: `0.848259`
- best active_null_z_mean: `2.612894`
- best control_pass_mean: `1.320000`

Interpretation: compared with explicit flux restoration, intrinsic directed
dynamics are less over-engineered and improve the pass rate slightly, but the
pathway-retention target is not met. The FMO sidecar remains a useful diagnostic
failure: it separates site identity, pathway retention, and strict specificity
instead of collapsing them into one optimistic score.

PASS requires:

- observed recoverability above threshold,
- site and pathway identity retention above threshold,
- the active level-5 null to pass,
- zero strict controls passing,
- observed recoverability to beat the best control by a margin.
