# HAOS-IIP Physics Bridge 53.0

This directory is an external post-processing layer over frozen scalar-carrier artifacts.

It does not modify HAOS core, change scalar operators, re-run the scalar simulations, or claim physical ontology. Its job is narrower: translate already-frozen HAOS-IIP scalar measurements into physics-facing proxy observables with explicit `PASS`, `OPEN`, and `WATCH` boundaries.

## Run

```bash
python3 experiments/physics_bridge/physics_observables.py
```

or through the public reproduction entry point:

```bash
python3 examples/quick_reproduce.py --physics-observables-only
```

## Outputs

- `experiments/physics_bridge/results/physics_observables.csv`
- `experiments/physics_bridge/results/physics_observables.json`
- `experiments/physics_bridge/results/physics_observables_summary.md`

## Discipline

- `PASS` means the proxy passes its existing frozen scalar threshold.
- `OPEN` means the proxy fails a threshold and must not be promoted.
- `WATCH` means the aggregate read passes, but a narrower margin remains unresolved.

The bridge dictionary is a test protocol, not a theory feature.

Any root-note or harmonic-partial reading must remain external to this directory and cannot change the bridge rows.

The inverse-square row is intentionally split:

- raw recoverability-gradient inverse-square closure remains `OPEN`;
- shell-native law-aware inverse-square closure is `PASS` under its stated reconstruction thresholds.

The power-law scaling row is intentionally split:

- raw local-gradient power scaling remains `OPEN`;
- shell-native law-aware power scaling is `PASS` as a bounded continuum-like proxy.

The localized-bump row is also intentionally split:

- weak localized bump response is `PASS` after excluding the earliest source-core transient layer;
- stronger localized bump response remains `OPEN` under the same transport thresholds.

The disorder-native row is also intentionally split:

- aggregated mild-disorder closure is `PASS`;
- seed-universal disorder closure is separately tracked and is currently `PASS` on the tested delayed fit window `[0.006, 0.024]`.

## Optional Sidecar Audit

An additional external sidecar audit is available for the raw local-gradient boundary:

```bash
python3 experiments/physics_bridge/raw_gradient_shape_audit.py
```

or through the public reproduction entry point:

```bash
python3 examples/quick_reproduce.py --raw-gradient-audit-only
```

This sidecar audit does not change the canonical bridge rows. Its strongest bounded result is narrower: the full raw local-gradient observable remains `OPEN`, while a source-core-excluded, amplitude-normalized shape-only read can pass the same raw thresholds as a sidecar proxy.

## Celestial Boundary Audit

The Phase 57 celestial-holography boundary audit is available as an external
sidecar:

```bash
python3 experiments/physics_bridge/celestial_boundary_audit/run_celestial_boundary_audit.py
```

Generated outputs:

- `experiments/physics_bridge/celestial_boundary_audit/haos_vs_celestial_requirements.csv`
- `experiments/physics_bridge/celestial_boundary_audit/bridge_status.json`
- `experiments/physics_bridge/celestial_boundary_audit/celestial_boundary_audit.md`

The saved roadmap remains at:

`docs/notes/foundations/HAOS_IIP_Celestial_Holography_Audit_Roadmap_v1.md`

The audit treats celestial holography as a high-standard boundary check for
HAOS-IIP physics-facing language. Its expected status is `OPEN`: it is not a
core change and not a claim of BMS, Virasoro, S-matrix, soft-theorem,
collinear-limit, or gravitational-memory recovery.

## Spherical Harmonic Control Probe

The Phase 58/58.1 spherical-harmonic control probe is available as a known-target
boundary-geometry sanity check:

```bash
python3 experiments/physics_bridge/spherical_harmonic_control_probe/run_spherical_harmonic_probe.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/spherical_harmonic_control_probe/run_spherical_harmonic_probe.py
```

This probe asks whether HAOS-style spectral telemetry can distinguish low-order
S2 mode organization from generic graph smoothness under controls. The 58.1
runner hardens the score with geodesic edge-weight signatures and l(l+1)-like
band-spacing checks. It is not a claim of celestial holography.

## Soft-Structure Proxy Test

The Phase 59 soft-structure proxy test is available as a toy amplitude-like
sidecar:

```bash
python3 experiments/physics_bridge/soft_structure_proxy_test/run_soft_structure_proxy_test.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/soft_structure_proxy_test/run_soft_structure_proxy_test.py
```

This probe asks whether soft-specific telemetry can detect known toy
pole/residue/factorization breakage better than generic smoothness metrics. It
is not a claim of gravitational soft-theorem recovery, celestial amplitudes, or
S-matrix reconstruction.

## Claim-Gated Physics Bridge Update

The Phase 60 claim-gated bridge update consolidates the celestial-facing sidecars:

```bash
python3 experiments/physics_bridge/claim_gated_bridge_update/run_claim_gated_bridge_update.py
```

Generated outputs:

- `experiments/physics_bridge/claim_gated_bridge_update/claim_gate_table.csv`
- `experiments/physics_bridge/claim_gated_bridge_update/bridge_status.json`
- `experiments/physics_bridge/claim_gated_bridge_update/claim_gated_physics_bridge.md`

This update locks the language boundary: the spherical and toy-soft sidecars may
be cited as bounded proxy PASS rows, while celestial holography, BMS charge
recovery, Virasoro/CCFT recovery, S-matrix reconstruction, real soft-theorem
recovery, and gravitational-memory observables remain `OPEN`.

## Bell/CHSH Computational Reference Probe

The Bell/CHSH reference probe is available as a bounded computational sidecar:

```bash
uv run python experiments/physics_bridge/bell_chsh_reference/run_bell_chsh_reference.py
```

Generated outputs:

- `experiments/physics_bridge/bell_chsh_reference/precommitment_contract.json`
- `experiments/physics_bridge/bell_chsh_reference/bell_chsh_diagnostics.csv`
- `experiments/physics_bridge/bell_chsh_reference/bridge_status.json`
- `experiments/physics_bridge/bell_chsh_reference/bell_chsh_reference_report.md`

This probe only checks a standard CHSH reference harness against classical
controls. It is not a laboratory Bell experiment, not a loophole-free Bell test,
and not a CST or HAOS-IIP derivation of Bell correlations.

## Bell HAOS-IIP Candidate Bridge

The first isolated HAOS-IIP-to-Bell candidate bridge is available as a bounded
candidate bundle:

```bash
uv run python experiments/physics_bridge/bell_haos_candidate/run_bell_haos_candidate.py
```

Generated outputs:

- `experiments/physics_bridge/bell_haos_candidate/precommitment_contract.json`
- `experiments/physics_bridge/bell_haos_candidate/assumption_ledger.json`
- `experiments/physics_bridge/bell_haos_candidate/candidate_trials.csv`
- `experiments/physics_bridge/bell_haos_candidate/candidate_correlations.json`
- `experiments/physics_bridge/bell_haos_candidate/no_signalling_diagnostics.csv`
- `experiments/physics_bridge/bell_haos_candidate/setting_independence_diagnostics.csv`
- `experiments/physics_bridge/bell_haos_candidate/control_results.csv`
- `experiments/physics_bridge/bell_haos_candidate/bell_candidate_report.md`
- `experiments/physics_bridge/bell_haos_candidate/bell_candidate_result.json`

This bridge implements B0/B1/B2 only. The frozen Bell reference sidecar remains
a scoreboard convention; this bundle does not import the quantum reference
curve, does not implement B3 joint closure, and does not establish HAOS-IIP or
CST Bell recovery.

The candidate branch has since been frozen as Bell Bridge v0.3. B3.1, B3.2,
B3.2.1, and B3.2.2 are retained as negative/provenance artifacts. The terminal
classification is `GENERIC_RELATIONAL_GEOMETRY`: the current `G / J_lambda`
construction produces genuine relational structure, but not HAOS-specific
structure sufficient to support a Bell bridge. CHSH scoring from this branch is
not authorized.

## Synthetic Relational Geometry Calibration

The synthetic relational-geometry calibration suite is available as an
instrument-only proving ground:

```bash
uv run python experiments/physics_bridge/synthetic_relational_geometry_calibration/run_synthetic_calibration.py
```

Generated outputs:

- `experiments/physics_bridge/synthetic_relational_geometry_calibration/precommitment_contract.json`
- `experiments/physics_bridge/synthetic_relational_geometry_calibration/synthetic_source_manifest.json`

## HBP Continuity

The HBP registry is the synthetic bridge program used by the 67.1 companion
release. It inherits the hardened 66.5 boundary discipline and keeps the row
labels terminal:

- PB-01: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`, `CONTROL_INVALID`,
  `HOLDOUT_TRANSFER_PASS`, `MIXED_OPEN`, `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
- PB-03: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`, `CONTROL_INVALID`,
  `MIXED_OPEN`
- PB-04: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`, `HAOS_BELL_STATUS_OPEN`,
  `MIXED_OPEN`

For the compact continuity note, see
`docs/notes/foundations/post_67_1_bridge_audit.md`.
- `experiments/physics_bridge/synthetic_relational_geometry_calibration/semantic_relation_matrix.csv`
- `experiments/physics_bridge/synthetic_relational_geometry_calibration/calibration_geometry_matrix.csv`
- `experiments/physics_bridge/synthetic_relational_geometry_calibration/calibration_control_results.csv`
- `experiments/physics_bridge/synthetic_relational_geometry_calibration/calibration_refinement_persistence.csv`
- `experiments/physics_bridge/synthetic_relational_geometry_calibration/synthetic_calibration_report.md`
- `experiments/physics_bridge/synthetic_relational_geometry_calibration/synthetic_calibration_result.json`

This suite tests whether the semantic/refinement instrument can recover known
synthetic relational structure, reject destructive controls, and leave ambiguous
partial-preservation cases open. It is not a physics bridge and does not
authorize Bell scoring.

## Literature Component Bridge

The literature-derived component bridge is available as a bounded method
calibration sidecar:

```bash
uv run python experiments/physics_bridge/literature_component_bridge/run_literature_component_bridge.py
```

Generated outputs:

- `experiments/physics_bridge/literature_component_bridge/precommitment_contract.json`
- `experiments/physics_bridge/literature_component_bridge/source_manifest.json`
- `experiments/physics_bridge/literature_component_bridge/component_scores.csv`
- `experiments/physics_bridge/literature_component_bridge/control_results.csv`
- `experiments/physics_bridge/literature_component_bridge/literature_component_bridge_report.md`
- `experiments/physics_bridge/literature_component_bridge/literature_component_bridge_result.json`

This sidecar combines spectral identity, Hodge sector, graph curvature,
transport, and kernel-distance components extracted from the bridge-literature
corpus. It is an operational method calibration only. A component-level `PASS`
does not establish a physical bridge, continuum limit, field theory, quantum
result, or gravity result.

## Hydrogen Spectra Computational Reference Probe

The hydrogen spectra reference probe is available as a bounded gross-spectrum
sidecar:

```bash
uv run python experiments/physics_bridge/hydrogen_spectra_reference/run_hydrogen_spectra_reference.py
```

Generated outputs:

- `experiments/physics_bridge/hydrogen_spectra_reference/precommitment_contract.json`
- `experiments/physics_bridge/hydrogen_spectra_reference/transition_catalog.json`
- `experiments/physics_bridge/hydrogen_spectra_reference/hydrogen_spectra_diagnostics.csv`
- `experiments/physics_bridge/hydrogen_spectra_reference/bridge_status.json`
- `experiments/physics_bridge/hydrogen_spectra_reference/hydrogen_spectra_reference_report.md`

This probe only checks a declared Rydberg-series reference harness and component
controls. It is not a laboratory hydrogen spectrum, not a fine-structure
calculation, and not a CST or HAOS-IIP derivation of hydrogen spectra.

## Semiconductor Band-Structure Computational Reference Probe

The semiconductor band-structure reference probe is available as a bounded toy
band sidecar:

```bash
uv run python experiments/physics_bridge/semiconductor_band_reference/run_semiconductor_band_reference.py
```

Generated outputs:

- `experiments/physics_bridge/semiconductor_band_reference/precommitment_contract.json`
- `experiments/physics_bridge/semiconductor_band_reference/band_catalog.json`
- `experiments/physics_bridge/semiconductor_band_reference/semiconductor_band_diagnostics.csv`
- `experiments/physics_bridge/semiconductor_band_reference/bridge_status.json`
- `experiments/physics_bridge/semiconductor_band_reference/semiconductor_band_reference_report.md`

This probe only checks a declared toy direct-gap band reference and component
controls. It is not a physical semiconductor calculation, not an ab-initio
band-structure calculation, and not a CST or HAOS-IIP derivation of
semiconductor physics.

## Gravitational Memory Toy Probe

The Phase 61 gravitational-memory toy probe is available as a synthetic
permanent-displacement benchmark on a sampled `S2` boundary:

```bash
python3 experiments/physics_bridge/gravitational_memory_toy_probe/run_gravitational_memory_toy_probe.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/gravitational_memory_toy_probe/run_gravitational_memory_toy_probe.py
```

This probe asks whether HAOS-style spectral telemetry can detect a known
memory-like step deformation better than controls that preserve smoothness,
low-mode structure, or time-series scale. It remains inside the Phase 60 claim
gate: it is not a claim of real gravitational-memory recovery, BMS soft hair,
celestial holography, or S-matrix reconstruction.

## Multi-Pole Supertranslation Toy Probe

The Phase 62 multi-pole supertranslation toy probe extends the Phase 61 memory
benchmark from one permanent deformation to an ordered synthetic `l=2,3,4`
mode-shift sequence:

```bash
python3 experiments/physics_bridge/multipole_supertranslation_probe/run_multipole_supertranslation_probe.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/multipole_supertranslation_probe/run_multipole_supertranslation_probe.py
```

This probe asks whether HAOS-style spectral telemetry can recover cumulative
field identity, multi-pole harmonic address, event ordering, and permanence
better than smooth or partial-statistic controls. It remains a toy benchmark:
it is not a claim of BMS supertranslation recovery, BMS charge recovery, real
gravitational memory, celestial holography, or S-matrix reconstruction.

## Supertranslation + Memory Composition Toy Probe

The Phase 63 composition probe combines the Phase 61 memory toy and Phase 62
multi-pole supertranslation toy into a two-stage synthetic benchmark:

```bash
python3 experiments/physics_bridge/supertranslation_memory_composition_probe/run_supertranslation_memory_composition_probe.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/supertranslation_memory_composition_probe/run_supertranslation_memory_composition_probe.py
```

This probe asks whether HAOS-style spectral telemetry can recover a synthetic
supertranslation-like shift, a later induced memory-like response, their
temporal order, and the fixed toy response relation better than controls that
preserve only pieces of the construction. It remains inside the Phase 60 claim
gate and is not a claim of real BMS supertranslation recovery, gravitational
memory, celestial holography, or S-matrix reconstruction.

## GW Memory Entry Gate

The Phase 64 GW memory entry gate moves from pure toy boundary fields toward
strain-like time-series data:

```bash
python3 experiments/physics_bridge/gw_memory_entry_gate/run_gw_memory_proxy.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/gw_memory_entry_gate/run_gw_memory_proxy.py
```

The default run uses a deterministic GW150914-like surrogate so the probe is
offline and reproducible. Optional local GWOSC-style HDF5, CSV/TXT, or NPY
strain files can be supplied with `--input-file`.

This probe asks whether HAOS-style memory/composition telemetry can distinguish
a structured strain-derived proxy event from time-shuffle, phase-scramble,
amplitude-preserving, event-window, off-event, and noise controls. It remains
inside the claim gate: it is not a claim of real gravitational-memory detection,
BMS charge recovery, soft theorem recovery, celestial holography, or S-matrix
reconstruction.

## GW Event-Window Hardening

The Phase 65 event-window hardening runner attacks the main Phase 64 leakage
path with stricter event-window controls and multiple deterministic surrogate
events:

```bash
python3 experiments/physics_bridge/gw_memory_entry_gate/run_gw_event_window_hardening.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/gw_memory_entry_gate/run_gw_event_window_hardening.py
```

This runner adds sliding-window, partial-overlap, chirp-reversal,
envelope-locked phase, timing-randomization, and micro-chunk controls. The
current status is `MARGINAL`: the target remains stable across replicates, but
event-window controls still compete strongly enough that the leak is not closed.
This is a hardening diagnostic, not a gravitational-memory claim.

## GW Event-Window Specificity

The Phase 66 event-window specificity runner adds event-internal ridge-coherence
and analytic-phase continuity metrics:

```bash
python3 experiments/physics_bridge/gw_memory_entry_gate/run_gw_event_window_specificity.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/gw_memory_entry_gate/run_gw_event_window_specificity.py
```

The current status is `MARGINAL`: the new metrics improve separation over Phase
65, but partial-overlap, micro-window, and sliding-window controls still prevent
a bounded PASS. This remains a strain-derived proxy diagnostic only.
