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
