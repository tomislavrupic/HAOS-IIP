# Semiconductor Band-Structure Computational Reference Probe

This sidecar runs a bounded toy band-structure reference calculation:

- a declared two-band direct-gap catalog should reproduce its own band table;
- controls preserve some confounds while breaking gap, directness, curvature,
  dispersion, band identity, or k-ordering;
- all scores remain component-visible.

It is not a physical semiconductor calculation, not an ab-initio band-structure
calculation, and not a CST or HAOS-IIP derivation of semiconductor physics.

## Run

```bash
uv run python experiments/physics_bridge/semiconductor_band_reference/run_semiconductor_band_reference.py
```

## Outputs

- `precommitment_contract.json`
- `band_catalog.json`
- `semiconductor_band_diagnostics.csv`
- `bridge_status.json`
- `semiconductor_band_reference_report.md`

## Status Boundary

`SEMICONDUCTOR_REFERENCE_SANITY_PASS` means only that the toy band-structure
reference harness reproduces its declared reference and rejects its declared
controls. CST semiconductor status remains `OPEN`.
