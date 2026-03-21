# Phase VII Spectral Feasibility

Primary authoritative track: Track A (`spectral_feasibility`).

This phase consumes the frozen Phase VI operator manifest and the frozen Phase V authority manifest without modifying either freeze. The operator under test is the branch-local periodic DK2D block cochain Laplacian `delta_h` on the frozen refinement hierarchy `h = 1 / n_side`.

Run:

```bash
python3 phase7-spectral/build_spectral_feasibility.py
```

Artifacts:

- `runs/phase7_spectral_scaling_ledger.csv`
- `phase7_spectral_summary.md`
- `phase7_spectral_manifest.json`
- `plots/phase7_eigenvalue_vs_h.svg`
- `plots/phase7_spectral_radius_vs_h.svg`
- `plots/phase7_trace_proxy_vs_h.svg`

Track B is not authoritative in this phase bundle.
