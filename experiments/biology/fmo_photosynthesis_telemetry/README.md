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

PASS requires:

- observed recoverability above threshold,
- site and pathway identity retention above threshold,
- the active level-5 null to pass,
- zero strict controls passing,
- observed recoverability to beat the best control by a margin.
