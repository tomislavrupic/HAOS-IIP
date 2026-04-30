# Phase 64 - GW Memory Entry Gate

This sidecar moves the claim-gated celestial bridge from pure toy fields toward
strain-like time-series data.

It asks one bounded question:

Can HAOS-style memory/composition telemetry distinguish a structured
strain-derived proxy event from controls?

It does not claim real gravitational memory, BMS charge recovery, soft theorem
recovery, celestial holography, or S-matrix reconstruction.

## Run

Default offline surrogate:

```bash
python3 experiments/physics_bridge/gw_memory_entry_gate/run_gw_memory_proxy.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/gw_memory_entry_gate/run_gw_memory_proxy.py
```

Optional local strain file:

```bash
python3 experiments/physics_bridge/gw_memory_entry_gate/run_gw_memory_proxy.py \
  --input-file path/to/local/GWOSC_file.hdf5 \
  --event-time 1126259462.4 \
  --detector H1
```

Supported local formats:

- GWOSC-style HDF5/H5 with a strain dataset such as `strain/Strain`
- CSV/TXT with either one strain column or two columns `time,strain`
- NPY with shape `(N,)` or `(N, 2)`

## Outputs

- `gw_memory_diagnostics.csv`
- `bridge_status.json`
- `gw_memory_entry_report.md`
- `prepared_strain_preview.csv`
- `prepared_strain_metadata.json`
- `figures/strain_proxy_trace.png`
- `figures/control_comparison.png`
- `figures/metric_breakdown.png`
- `figures/time_frequency_proxy.png`

## Controls

- global time shuffle
- FFT phase scramble
- amplitude-preserving spectrum surrogate
- event-window chunk scramble
- structured event shifted away from the claimed event time
- noise-only same-scale baseline

## Claim Gate

`PASS` means the strain-derived proxy separates a structured strain-like event
from the listed controls under the current telemetry.

`PASS` does not mean:

- real gravitational memory was detected
- BMS charges or soft hair were recovered
- celestial holography was recovered
- real soft theorems were tested
- an S-matrix was reconstructed

The default run uses a deterministic GW150914-like surrogate so the phase is
reproducible without network access or external data.

## Phase 65 Event-Window Hardening

Phase 65 hardens the main Phase 64 leak by running multiple deterministic
surrogate events against stricter event-window controls:

```bash
python3 experiments/physics_bridge/gw_memory_entry_gate/run_gw_event_window_hardening.py
```

If NumPy/Matplotlib are not installed in the default Python runtime:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/gw_memory_entry_gate/run_gw_event_window_hardening.py
```

Additional outputs:

- `event_window_hardening_diagnostics.csv`
- `event_window_hardening_summary.csv`
- `event_window_hardening_status.json`
- `gw_event_window_hardening_report.md`
- `figures/phase65_event_scores.png`
- `figures/phase65_margin_distribution.png`
- `figures/phase65_control_ladder.png`

The current Phase 65 status is `MARGINAL`: target scores remain high across
replicates, but envelope-locked and sliding/micro-window controls still compete.
