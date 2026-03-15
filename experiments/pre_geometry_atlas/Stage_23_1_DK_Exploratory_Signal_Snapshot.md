# Stage 23.1 DK Exploratory Signal Snapshot

## Scope
This snapshot freezes the corrected Stage 23.1 read after three passes:
- corrected `n = 12` Dirac-Kaehler collision scan using combined-field peak tracking
- matched `n = 12` scalar bosonic control matrix on the same geometry-phase grid
- `n = 24` Dirac-Kaehler clustered-pair refinement only

## Corrected DK base scan
Timestamped JSON: `data/20260315_140539_stage23_1_dk_collision.json`

Persistence-first result:
- candidate bounded collision regimes: `0 / 10`
- weak persistence gain: `2 / 10`
- no persistence gain: `8 / 10`

Collision texture result:
- `pass-through with grade exchange`: `9 / 10`
- `unresolved / mixed`: `1 / 10`

Clustered-family signal:
- `S23_1_clustered_in_phase`
  - collision label: `unresolved / mixed`
  - persistence label: `weak persistence gain`
  - composite lifetime: `0.7517`
  - binding persistence: `1.0000`
  - grade-transfer amplitude: `0.1303`
- `S23_1_clustered_out_of_phase`
  - collision label: `pass-through with grade exchange`
  - persistence label: `weak persistence gain`
  - composite lifetime: `0.1253`
  - binding persistence: `0.1667`
  - grade-transfer amplitude: `0.2303`

## Matched scalar control
Timestamped JSON: `data/20260315_140932_stage23_1_scalar_collision_control.json`

Persistence-first result:
- candidate bounded collision regimes: `1 / 10`
- no persistence gain: `9 / 10`

Key control outcome:
- `S23_1_scalar_clustered_in_phase`
  - collision label: `metastable composite`
  - persistence label: `candidate bounded collision regime`
  - composite lifetime: `0.6058`
  - binding persistence: `1.0000`
- `S23_1_scalar_clustered_out_of_phase`
  - collision label: `pass-through dispersive`
  - persistence label: `no persistence gain`

Control consequence:
- the strongest clustered persistence signature is not uniquely Dirac-Kaehler
- the DK sector does not outperform the matched scalar control on persistence at `n = 12`

## Clustered-pair refinement at n = 24
Timestamped JSON: `data/20260315_141011_stage23_1_dk_clustered_refinement.json`

Result:
- candidate bounded collision regimes: `0 / 2`
- `S23_1_clustered_in_phase_n24`
  - collision label: `unresolved / mixed`
  - persistence label: `weak persistence gain`
  - composite lifetime: `0.7489`
  - binding persistence: `1.0000`
  - grade-transfer amplitude: `0.0414`
- `S23_1_clustered_out_of_phase_n24`
  - collision label: `pass-through dispersive`
  - persistence label: `no persistence gain`
  - composite lifetime: `0.0000`
  - binding persistence: `0.0000`
  - grade-transfer amplitude: `0.0819`

Refinement consequence:
- the in-phase clustered anomaly persists only as an unresolved close-packed encounter
- the out-of-phase clustered anomaly does not survive refinement as a persistence signal
- grade-transfer amplitude weakens under refinement

## Freeze line
Stage 23.1 is exploratory-positive only in the following limited sense:
- Dirac-Kaehler propagation changes collision texture materially through widespread grade exchange and altered clustered-pair encounter structure

Stage 23.1 is not positive in the stronger persistence sense:
- no DK run earns a bounded collision regime under the corrected read
- matched scalar control is at least as strong in clustered-pair persistence at `n = 12`
- clustered-pair refinement does not promote the DK signal into a clean bounded regime

## Operational conclusion
Freeze Stage 23.1 as:
- a real operator-sector collision-phenomenology signal
- not yet a validated persistence-landscape shift
- sufficient to justify Phase III continuation only if the next branch targets DK-specific collision structure rather than generic persistence claims
