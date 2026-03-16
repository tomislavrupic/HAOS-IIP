# Stage C0.12 Harmonic Detuning Continuum Scan v1

Timestamped JSON: `data/20260315_214453_stage_c0_12_harmonic_detuning_continuum_scan.json`
Timestamped CSV: `data/20260315_214453_stage_c0_12_harmonic_detuning_continuum_scan.csv`

Continuity assessment: `sharp_transition_threshold`

Per-run summary:

- `C012_detune_000`
  - delta: `0.000`
  - topology class: `braid_like_exchange`
  - flow concentration index: `0.8989`
  - address selectivity index: `1.0000`
  - address sorting score: `1.0000`
  - topology survival time: `0.3759`
  - transfer asymmetry: `0.3130`
- `C012_detune_0125`
  - delta: `0.125`
  - topology class: `braid_like_exchange`
  - flow concentration index: `0.9006`
  - address selectivity index: `0.9790`
  - address sorting score: `0.5378`
  - topology survival time: `0.5011`
  - transfer asymmetry: `0.3090`
- `C012_detune_0250`
  - delta: `0.250`
  - topology class: `braid_like_exchange`
  - flow concentration index: `0.9022`
  - address selectivity index: `0.9595`
  - address sorting score: `0.5542`
  - topology survival time: `0.5011`
  - transfer asymmetry: `0.3052`
- `C012_detune_0375`
  - delta: `0.375`
  - topology class: `braid_like_exchange`
  - flow concentration index: `0.9037`
  - address selectivity index: `0.9411`
  - address sorting score: `0.5678`
  - topology survival time: `0.5011`
  - transfer asymmetry: `0.3014`
- `C012_detune_0500`
  - delta: `0.500`
  - topology class: `braid_like_exchange`
  - flow concentration index: `0.9051`
  - address selectivity index: `0.9237`
  - address sorting score: `0.5805`
  - topology survival time: `0.5011`
  - transfer asymmetry: `0.2977`
- `C012_detune_0625`
  - delta: `0.625`
  - topology class: `braid_like_exchange`
  - flow concentration index: `0.9073`
  - address selectivity index: `0.8959`
  - address sorting score: `0.6016`
  - topology survival time: `0.5011`
  - transfer asymmetry: `0.2915`
- `C012_detune_0750`
  - delta: `0.750`
  - topology class: `transfer_smeared`
  - flow concentration index: `0.9091`
  - address selectivity index: `0.8706`
  - address sorting score: `0.6206`
  - topology survival time: `0.5011`
  - transfer asymmetry: `0.2853`
- `C012_detune_0875`
  - delta: `0.875`
  - topology class: `transfer_smeared`
  - flow concentration index: `0.9103`
  - address selectivity index: `0.8471`
  - address sorting score: `0.6363`
  - topology survival time: `0.0000`
  - transfer asymmetry: `0.2791`
- `C012_detune_1000`
  - delta: `1.000`
  - topology class: `transfer_smeared`
  - flow concentration index: `0.9108`
  - address selectivity index: `0.8248`
  - address sorting score: `0.6496`
  - topology survival time: `0.0000`
  - transfer asymmetry: `0.2727`

Continuity note:
- Topology stays in one protected class and then switches once near delta=0.750. Selectivity decays monotonically, while the sorting curve shows secondary texture inside the protected branch rather than a second topology transition.

## Plots
- `plots/20260315_214453_stage_c0_12_C012_detune_000_topology_vs_time_trace.png`
- `plots/20260315_214453_stage_c0_12_C012_detune_0125_topology_vs_time_trace.png`
- `plots/20260315_214453_stage_c0_12_C012_detune_0250_topology_vs_time_trace.png`
- `plots/20260315_214453_stage_c0_12_C012_detune_0375_topology_vs_time_trace.png`
- `plots/20260315_214453_stage_c0_12_C012_detune_0500_topology_vs_time_trace.png`
- `plots/20260315_214453_stage_c0_12_C012_detune_0625_topology_vs_time_trace.png`
- `plots/20260315_214453_stage_c0_12_C012_detune_0750_topology_vs_time_trace.png`
- `plots/20260315_214453_stage_c0_12_C012_detune_0875_topology_vs_time_trace.png`
- `plots/20260315_214453_stage_c0_12_C012_detune_1000_topology_vs_time_trace.png`
- `plots/20260315_214453_stage_c0_12_topology_phase_diagram.png`
- `plots/20260315_214453_stage_c0_12_selectivity_decay_curve.png`
- `plots/20260315_214453_stage_c0_12_sorting_degradation_curve.png`
- `plots/20260315_214453_stage_c0_12_braid_survival_panel.png`
