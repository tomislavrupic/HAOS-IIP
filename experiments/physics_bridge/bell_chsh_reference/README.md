# Bell/CHSH Computational Reference Probe

This sidecar runs a bounded CHSH reference calculation:

- deterministic local response functions must obey the local CHSH bound;
- a standard singlet-reference correlation table should violate that bound;
- seeded finite-sample telemetry should reproduce the reference behavior;
- local and uncorrelated controls should remain non-violating or near zero.

It is not a laboratory Bell experiment, not a loophole-free Bell test, and not a
CST or HAOS-IIP derivation of Bell correlations.

## Run

```bash
uv run python experiments/physics_bridge/bell_chsh_reference/run_bell_chsh_reference.py
```

## Outputs

- `precommitment_contract.json`
- `bell_chsh_diagnostics.csv`
- `bridge_status.json`
- `bell_chsh_reference_report.md`

## Status Boundary

`BELL_REFERENCE_SANITY_PASS` means only that the computational CHSH reference
harness reproduces its precommitted analytic expectations and controls. CST Bell
status remains `OPEN`.
