# EXPERIMENT_LOG

This file records orchestrated experiment runs.

## Laplacian geometry test
- Date: 2026-03-09T22:31:38
- Config: substrate=random_geometric, nodes=500, epsilon=0.2, seed=42
- Results: `data/20260309_223137_laplacian_modes.json`
- Plots: `plots/20260309_223137_laplacian_modes_laplacian_modes_spectrum.png`, `plots/20260309_223137_laplacian_modes_laplacian_modes_mode1.png`
- Observation: low spectrum contains localized modes
- Conclusion: geometry sector mixed with substrate or boundary artifacts

## Gauge sector test
- Date: 2026-03-09T22:31:38
- Config: substrate=cubic_lattice, nodes=216, epsilon=0.2, seed=42
- Results: `data/20260309_223137_gauge_modes.json`
- Plots: `plots/20260309_223137_gauge_modes_gauge_modes_spectrum.png`, `plots/20260309_223137_gauge_modes_gauge_modes_phase_mode.png`
- Observation: phase dressing produces circulating low modes
- Conclusion: connection-sensitive gauge branch present, but still scalar in background transport

## Laplacian geometry test
- Date: 2026-03-09T22:33:52
- Config: substrate=random_geometric, nodes=500, epsilon=0.2, seed=42
- Results: `data/20260309_223352_laplacian_modes.json`
- Plots: `plots/20260309_223352_laplacian_modes_laplacian_modes_spectrum.png`, `plots/20260309_223352_laplacian_modes_laplacian_modes_mode1.png`
- Observation: low spectrum contains localized modes
- Conclusion: geometry sector mixed with substrate or boundary artifacts

## Gauge sector test
- Date: 2026-03-09T22:33:52
- Config: substrate=cubic_lattice, nodes=216, epsilon=0.2, seed=42
- Results: `data/20260309_223352_gauge_modes.json`
- Plots: `plots/20260309_223352_gauge_modes_gauge_modes_spectrum.png`, `plots/20260309_223352_gauge_modes_gauge_modes_phase_mode.png`
- Observation: phase dressing produces circulating low modes
- Conclusion: connection-sensitive gauge branch present, but still scalar in background transport
