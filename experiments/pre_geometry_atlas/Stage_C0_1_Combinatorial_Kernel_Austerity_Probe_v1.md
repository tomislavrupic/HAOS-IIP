# Stage C0.1 Combinatorial Kernel Austerity Probe v1

Timestamped summary: `data/20260315_185133_stage_c0_1_combinatorial_kernel_probe.json`
Timestamped run table: `data/20260315_185133_stage_c0_1_combinatorial_kernel_probe.csv`

Architecture notice: The frozen regular-lattice Hodge branch already uses uniform nearest-neighbor edge weights, so this probe instantiates an explicit edge-graph kernel operator to make the metric-removal question non-vacuous.

Kernel comparison table:

- `clustered_composite_anchor` / `Gaussian baseline`
  - interaction label: `metastable composite regime`
  - coarse label: `diffuse coarse field regime`
  - composite lifetime: `1.4764`
  - binding persistence: `1.0000`
  - coarse persistence sigma=4: `1.0000`
  - delta vs Gaussian composite lifetime: `0.0000`
  - delta vs Gaussian binding persistence: `0.0000`
- `clustered_composite_anchor` / `Pure adjacency`
  - interaction label: `metastable composite regime`
  - coarse label: `diffuse coarse field regime`
  - composite lifetime: `1.4764`
  - binding persistence: `1.0000`
  - coarse persistence sigma=4: `1.0000`
  - delta vs Gaussian composite lifetime: `0.0000`
  - delta vs Gaussian binding persistence: `0.0000`
- `clustered_composite_anchor` / `Graph-shell`
  - interaction label: `metastable composite regime`
  - coarse label: `diffuse coarse field regime`
  - composite lifetime: `1.4764`
  - binding persistence: `1.0000`
  - coarse persistence sigma=4: `1.0000`
  - delta vs Gaussian composite lifetime: `0.0000`
  - delta vs Gaussian binding persistence: `0.0000`
- `counter_propagating_corridor` / `Gaussian baseline`
  - interaction label: `transient binding regime`
  - coarse label: `diffuse coarse field regime`
  - composite lifetime: `0.0000`
  - binding persistence: `0.0000`
  - coarse persistence sigma=4: `1.0000`
  - delta vs Gaussian composite lifetime: `0.0000`
  - delta vs Gaussian binding persistence: `0.0000`
- `counter_propagating_corridor` / `Pure adjacency`
  - interaction label: `transient binding regime`
  - coarse label: `diffuse coarse field regime`
  - composite lifetime: `0.0000`
  - binding persistence: `0.0000`
  - coarse persistence sigma=4: `1.0000`
  - delta vs Gaussian composite lifetime: `0.0000`
  - delta vs Gaussian binding persistence: `0.0000`
- `counter_propagating_corridor` / `Graph-shell`
  - interaction label: `transient binding regime`
  - coarse label: `diffuse coarse field regime`
  - composite lifetime: `0.0000`
  - binding persistence: `0.0000`
  - coarse persistence sigma=4: `1.0000`
  - delta vs Gaussian composite lifetime: `0.0000`
  - delta vs Gaussian binding persistence: `0.0000`
- `phase_ordered_symmetric_triad` / `Gaussian baseline`
  - interaction label: `transient binding regime`
  - coarse label: `diffuse coarse field regime`
  - composite lifetime: `0.0000`
  - binding persistence: `0.0000`
  - coarse persistence sigma=4: `1.0000`
  - delta vs Gaussian composite lifetime: `0.0000`
  - delta vs Gaussian binding persistence: `0.0000`
- `phase_ordered_symmetric_triad` / `Pure adjacency`
  - interaction label: `transient binding regime`
  - coarse label: `diffuse coarse field regime`
  - composite lifetime: `0.0000`
  - binding persistence: `0.0000`
  - coarse persistence sigma=4: `1.0000`
  - delta vs Gaussian composite lifetime: `0.0000`
  - delta vs Gaussian binding persistence: `0.0000`
- `phase_ordered_symmetric_triad` / `Graph-shell`
  - interaction label: `transient binding regime`
  - coarse label: `diffuse coarse field regime`
  - composite lifetime: `0.0000`
  - binding persistence: `0.0000`
  - coarse persistence sigma=4: `1.0000`
  - delta vs Gaussian composite lifetime: `0.0000`
  - delta vs Gaussian binding persistence: `0.0000`

Interpretation boundary:
- this is an explicit edge-graph control branch, not a replacement of the frozen Hodge operator stack
- on the frozen regular-lattice branch, nearest-neighbor edge weights are already uniform, so adjacency-only is not a meaningful falsification unless a separate edge-graph operator is instantiated
- positive outcomes here mean graph-structural survival inside this control branch, not a blanket claim about every existing stage
