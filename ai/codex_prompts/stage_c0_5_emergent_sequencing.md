# Stage C0.5 - Emergent Sequencing from Incidence-Mismatch Accumulation and Protected Topology Stability

Purpose:

Keep the Paper 24.1 no-distance edge-graph control branch fixed, but replace the assumed packet-evolution clock with a combinatorial sequencing rule derived from incidence-mismatch accumulation on adjacency and graph-shell support refinements.

Frozen ingredients:

- same graph size as C0.1-C0.4
- same adjacency and graph-shell kernel families
- same representative seeds
- same projected-transverse state sector
- same shell-weight measurement vocabulary
- no Euclidean distance and no Gaussian embedding distance

New observables:

1. Persistence-class stabilization depth

- After each support update, compare the current topology label with the adjacency-next and graph-shell-next labels.
- Count the number of distinct active persistence classes at that step.
- The first step where this count stabilizes is the emergent sequencing depth.

2. Incidence-mismatch ledger

- Record local motif or degree changes when the chosen next support does not preserve the protected topology label of the representative seed.
- Report mismatch density as ledger entries divided by total edge count.
- The first quiet step is the first step where mismatch density falls below threshold.

3. Kernel-update consistency probe

- Compare adjacency-next and graph-shell-next topology labels on the same current support.
- Any disagreement is a kernel-update consistency failure.

4. Protected topology-stability constraint

- In the protected mode, use the PASS-style blocked-merge rule from the survivor-preorder branch.
- The first blocked merge is the emergent protected ordering step.

Minimal run design:

- 3 representatives
- 3 sequencing modes
  - adjacency sequencing
  - graph-shell sequencing
  - protected alternating sequencing

Total: 9 runs

Required outputs:

- stamped JSON
- stamped CSV
- per-run sequencing trace plot
- per-run support overlay
- summary sequencing-class matrix
- depth and mismatch panels
- atlas note

Interpretation boundary:

- This stage asks whether sequencing can be derived from incidence-mismatch accumulation and protected topology stability on the no-distance branch.
- It does not claim physical time, geometry emergence, particles, or universal irreversibility.
