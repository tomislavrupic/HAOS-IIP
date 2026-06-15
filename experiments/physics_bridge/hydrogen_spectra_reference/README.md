# Hydrogen Spectra Computational Reference Probe

This sidecar runs a bounded gross-spectrum reference calculation:

- a declared hydrogen Rydberg-series catalog should reproduce its own line table;
- controls preserve line counts, spectral ranges, ordering, exact line sets, or
  transition labels while breaking the tested structure;
- a minimal `delta_l = 1` selection-rule check is included as a reference
  bookkeeping gate.

It is not a laboratory hydrogen spectrum, not a fine-structure calculation, and
not a CST or HAOS-IIP derivation of hydrogen spectra.

## Run

```bash
uv run python experiments/physics_bridge/hydrogen_spectra_reference/run_hydrogen_spectra_reference.py
```

## Outputs

- `precommitment_contract.json`
- `transition_catalog.json`
- `hydrogen_spectra_diagnostics.csv`
- `bridge_status.json`
- `hydrogen_spectra_reference_report.md`

## Status Boundary

`HYDROGEN_REFERENCE_SANITY_PASS` means only that the computational gross-spectrum
reference harness reproduces its precommitted line law and rejects its declared
controls. CST hydrogen status remains `OPEN`.
