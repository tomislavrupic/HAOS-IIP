# Phase 59 Soft-Structure Proxy Test

This physics-bridge sidecar is a toy soft-pole and factorization benchmark.

It does not use real scattering amplitudes and does not claim recovery of
gravitational soft theorems, celestial amplitudes, BMS charges, Virasoro
structure, S-matrix data, or gravitational memory.

## Run

```bash
python3 experiments/physics_bridge/soft_structure_proxy_test/run_soft_structure_proxy_test.py
```

If your local `python3` does not have NumPy/Matplotlib:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/soft_structure_proxy_test/run_soft_structure_proxy_test.py
```

## Outputs

- `soft_diagnostics.csv`
- `bridge_status.json`
- `soft_structure_report.md`
- `figures/control_comparison.png`
- `figures/metric_breakdown.png`
- `figures/residue_recovery.png`
- `figures/smoothness_vs_soft_specificity.png`

## Test Logic

The target is synthetic amplitude-like data with:

- a known leading soft pole at `omega = 0`;
- known residues;
- rank-one residue factorization across six toy channels.

Controls preserve some broad properties while breaking the soft structure:

- per-frequency spectrum shuffle;
- wrong residues with the same pole;
- shifted pole locations;
- non-factorizing residues;
- smooth-but-wrong amplitudes.

`PASS` means the toy target is separable from controls on residue recovery, pole
location, and factorization. Even a `PASS` is only a toy soft-structure sanity
check, not a physics claim.
