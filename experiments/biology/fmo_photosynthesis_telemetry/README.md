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

PASS requires:

- observed recoverability above threshold,
- site and pathway identity retention above threshold,
- the active level-5 null to pass,
- zero strict controls passing,
- observed recoverability to beat the best control by a margin.
