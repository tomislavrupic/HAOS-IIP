# Stage C0.5 - Emergent Compression Preorder on the No-Distance Branch

Purpose:

Stage C0.5 keeps the explicit no-distance edge-graph branch fixed, but replaces packet evolution with a relational compression rule. Initial packet support is used only to seed the active edge set. Every later step is generated combinatorially from survivor-class refinement.

Frozen ingredients:

- same graph size and projected-transverse setup as C0.1-C0.4
- same representative seeds
- same adjacency and graph-shell kernel families
- no Euclidean distance
- no Gaussian embedding distance
- no adaptive reinforcement or external clock

Core diagnostics:

1. Behavioral-congruence survivor count

- Build edge congruence classes on the line graph from adjacency or graph-shell signatures.
- Survivor support = edges that remain in non-singleton congruence classes.
- Plot `W(step)` = survivor count vs compression step.
- Define emergent `tau` as the first stabilization depth of the survivor support.

2. Defect ledger

- At each transition record survivor edges whose local motif-incidence signature changed while their behavioral congruence class did not.
- Report defect density `delta = |ledger| / |edges|`.
- Record the first quiet step where `delta` drops below threshold.

3. Non-functoriality probe

- Compare adjacency and graph-shell congruence partitions on the same survivor graph.
- Count pairwise congruence mismatches as non-functoriality defects.
- In the alternating branch, choose the next refinement that reduces this mismatch count most.

4. PASS separation test

- In the alternating branch, forbid any refinement that would merge previously distinct survivor classes.
- The first blocked merge defines the emergent `Z2` arrow step.

Minimal run design:

- 3 representatives
- 3 compression modes
  - adjacency preorder
  - graph-shell preorder
  - alternating PASS preorder

Total: 9 runs

Required outputs:

- stamped JSON
- stamped CSV
- per-run compression trace plot with `W`, defect density, mismatch count, and obstruction rank
- per-run survivor overlay
- summary class matrix and depth panels
- atlas note

Interpretation boundary:

- This stage tests whether a purely combinatorial compression preorder can sequence the no-distance branch.
- It does not claim physical time, geometry emergence, particles, or universal irreversibility.
