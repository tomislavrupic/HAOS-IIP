# Phase XI -- Localized Mode Protection and Failure-Mechanism Program

This phase tests whether the frozen Phase X localized candidate
`low_mode_localized_wavepacket` is protected by reproducible structural features
or is only an accidental configuration.

Frozen inputs:
- `phase6-operator/phase6_operator_manifest.json`
- `phase7-spectral/phase7_spectral_manifest.json`
- `phase8-trace/phase8_trace_manifest.json`
- `phase9-invariants/phase9_manifest.json`
- `phase10-bridge/phase10_manifest.json`
- `phaseX-proto-particle/phaseX_integrated_manifest.json`
- `phaseX-proto-particle/runs/proto_particle_candidates.json`

Primary outputs:
- `runs/phase11_perturbation_survival_ledger.csv`
- `runs/phase11_failure_classification.csv`
- `runs/phase11_persistence_scaling.csv`
- `phase11_summary.md`
- `phase11_manifest.json`
- `runs/phase11_runs.json`

The builder is deterministic and self-contained:

```bash
python3 phase11-protection/build_localized_mode_protection.py
```

Validation:

```bash
python3 phase11-protection/diagnostics/check_phase11_bundle.py
```
