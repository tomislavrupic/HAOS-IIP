# Stage C0.3 Topology-Exchange Robustness Under Incidence Noise v1

Timestamped summary: `data/20260315_192917_stage_c0_3_incidence_noise_robustness.json`
Timestamped run table: `data/20260315_192917_stage_c0_3_incidence_noise_robustness.csv`

Architecture notice: This branch varies only the combinatorial incidence structure of the explicit edge-graph control operator. Packet initialization, projected-transverse state sector, graph size, and shell-weight kernel family remain frozen.
Baseline notice: The current C0.1-C0.2 combinatorial branch does not yet contain a clean braid_like_exchange baseline. C0.3 therefore anchors to the nearest exchange-like protected surrogate, counter_propagating_corridor, and treats braid detection as an allowed outcome rather than an assumed starting point.

Per-run summary:

- `Baseline no-noise`
  - topology class: `smeared_transfer`
  - flow concentration: `46.1431`
  - exchange coherence: `0.2158`
  - channel count: `10`
  - loop count: `0`
  - degree spectrum shift: `0.0000`
  - clustering coefficient shift: `0.0000`
- `Ultra-weak rewiring`
  - topology class: `smeared_transfer`
  - flow concentration: `42.9971`
  - exchange coherence: `0.2330`
  - channel count: `10`
  - loop count: `0`
  - degree spectrum shift: `0.0786`
  - clustering coefficient shift: `-0.0008`
- `Weak rewiring`
  - topology class: `smeared_transfer`
  - flow concentration: `40.6109`
  - exchange coherence: `0.2560`
  - channel count: `10`
  - loop count: `0`
  - degree spectrum shift: `0.1361`
  - clustering coefficient shift: `-0.0024`
- `Moderate rewiring`
  - topology class: `smeared_transfer`
  - flow concentration: `40.1626`
  - exchange coherence: `0.2528`
  - channel count: `10`
  - loop count: `0`
  - degree spectrum shift: `0.2357`
  - clustering coefficient shift: `-0.0071`
- `Moderate deletion-insertion`
  - topology class: `smeared_transfer`
  - flow concentration: `40.7380`
  - exchange coherence: `0.2449`
  - channel count: `10`
  - loop count: `0`
  - degree spectrum shift: `0.2722`
  - clustering coefficient shift: `-0.0046`
- `Clustered rewiring patch`
  - topology class: `smeared_transfer`
  - flow concentration: `39.2678`
  - exchange coherence: `0.2529`
  - channel count: `10`
  - loop count: `0`
  - degree spectrum shift: `0.1811`
  - clustering coefficient shift: `-0.0027`
- `Random localized defect pair`
  - topology class: `smeared_transfer`
  - flow concentration: `37.1418`
  - exchange coherence: `0.2443`
  - channel count: `11`
  - loop count: `0`
  - degree spectrum shift: `0.2097`
  - clustering coefficient shift: `-0.0010`
- `Degree-biased perturbation`
  - topology class: `smeared_transfer`
  - flow concentration: `40.1548`
  - exchange coherence: `0.2450`
  - channel count: `10`
  - loop count: `0`
  - degree spectrum shift: `0.2178`
  - clustering coefficient shift: `-0.0059`
- `Connectivity-edge stress test`
  - topology class: `smeared_transfer`
  - flow concentration: `42.8851`
  - exchange coherence: `0.2350`
  - channel count: `10`
  - loop count: `0`
  - degree spectrum shift: `0.0278`
  - clustering coefficient shift: `-0.0000`

Interpretation boundary:
- this scan asks whether no-distance combinatorial exchange-like transport is robust to controlled incidence perturbations
- it does not claim geometry emergence, particles, or bound-state physics
