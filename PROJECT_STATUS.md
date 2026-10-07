# HAOS-IIP Current Research Status

Reconciled: **7 October 2026**.

This is the current navigation index. It summarizes the authorities linked below;
it does not replace their contracts, change scientific classifications, or
unblock implementation. Release 66.5 remains a frozen baseline and 66.4 remains
the recommended entry paper, rather than dates for the latest local research.

## Current evidence and lifecycle

| Work | Recorded outcome | Development reading | Evidence / authority |
| --- | --- | --- | --- |
| EL-R3-THEORY-REVISION-01 | `FROZEN_THEORY_REVISION` | Complete; do not schedule the revision again | [Theory](docs/roadmaps/EL_R3_THEORY_REVISION_01.md) |
| EL-R3-MESOSCOPIC-FIELD-01 | `PARTIAL_RECOVERY_ONLY` | `TERMINAL_NEGATIVE`; Rung 3 unsupported | [Frozen result](experiments/emergence_ladder/rung3_mesoscopic_field/final/summary.md), [manifest](experiments/emergence_ladder/rung3_mesoscopic_field/final/result_manifest.json) |
| HAOS-AF-AR-01 | `BOUNDED_TOY_ANTIFRAGILITY_SUPPORTED` | Preserve the bounded sidecar; no higher-rung promotion | [Result](experiments/antifragility/adaptive_recoverability_v1/final/summary.md), [hostile audit](experiments/antifragility/adaptive_recoverability_v1/audit/hostile_audit.md) |
| PIX-7 finite-capacity exchange | `INTERNAL_FINITE_GRAPH_PROBE`; diagnostics `PASS_SCOPED` | Finite conservative model; physical gravity not demonstrated | [Stored result](experiments/physics_bridge/pix7_finite_capacity_exchange/outputs/results.json), [derivation and boundaries](docs/notes/foundations/PIX_7_Constructed_Finite_Capacity_Exchange_v1.md) |
| EL-R4-OC-02 / EL-R5-CSR-01 | Registry candidates; implementation blocked | Missing lower-rung support and required frozen precommitments | [Lifecycle registry](docs/branch_governance/branch_lifecycle_summary.md) |
| HBP PB-01 through PB-04 | `QUARANTINED_INVALID` | Historical artifacts only; no in-place rehabilitation | [HBP status](experiments/hbp/hbp_status_snapshot.md) |
| HBP-IR-01 | `INSTRUMENT_VALID` | Supporting synthetic calibration; no external prediction claim | [Integrity result](experiments/hbp/integrity_repair_v2/results/hbp_ir_01_report.md) |
| Current scale-bridge 66.5 mechanism | CP2 / CP3 / comparative / CP5 gates `OPEN` | `TERMINAL_NEGATIVE` for development; lower-level evidence retained | [Lifecycle registry](docs/branch_governance/branch_lifecycle_summary.md), [claim contract](docs/claim_governance/scale_bridge_66_5.claim.json) |
| HAOS-GEN synthetic generative line | Verified negative; paused pending external task | No further synthetic versions without the external-task gate | [Status](HAOS_GEN_STATUS.json) |

The mesoscopic field restores two perturbation families at each tested level,
but its within-cell nullspace family remains unrecovered. The resulting `2/3`
rate is below the frozen `0.75` gate. A valid checker result does not promote
this partial scientific outcome.

The antifragility sidecar's positive result remains limited by its symmetric
ring, designer-specified update, structural capacity metric, and static-hardening
reference. PIX-7 numerical consistency is not a gravity gate pass. Neither
sidecar changes the emergence ladder or creates an authorized successor.

## Next gate

The highest supported emergence rung remains Rung 2. The theory revision and
its mesoscopic successor are finished; no new Rung 3 candidate is currently
authorized. A successor requires independent motivation, a new identifier, and
a frozen precommitment under the registry's reopening policy. Existing failures,
thresholds, controls, and result identifiers remain preserved.

The two candidates in the registry's active queue are dependency-blocked.
Their presence in that queue is not permission to implement operational closure
or cross-scale recovery. See the [reconciled emergence ladder](docs/roadmaps/HAOS_IIP_EMERGENCE_LADDER_2026-07-11.md).

## Preservation and verification

The [reviewed local snapshot](docs/snapshots/RESEARCH_SNAPSHOT_2026-10-07.md)
records the Git anchor, exact file inventory, verification scope, and restoration
instructions. It preserves 63 formerly untracked files and seven existing
modified files, including raw rows, figures, arrays, contracts, tests, and audits.
It is a local snapshot, not a release or independent replication.

From the repository root:

```bash
uv run python scripts/check_branch_lifecycle.py
uv run python scripts/check_claim_contracts.py
uv run python experiments/emergence_ladder/rung3_mesoscopic_field/check_bundle.py
uv run python experiments/antifragility/adaptive_recoverability_v1/check_bundle.py
uv run python examples/quick_reproduce.py
```

These commands check saved artifacts and governance. Do not rerun single-use
final partitions as a substitute for their bundle checkers.

## Authority and historical notes

- Scientific outcomes: each frozen contract, result, and manifest.
- Development decisions: [lifecycle JSON](docs/branch_governance/branch_lifecycle_registry.json);
  its Markdown summary is generated from that authority.
- Public claim ceilings: [claim governance](docs/claim_governance/README.md).
- Current navigation: this page, with the README and ladder pointing to it.
- [March 10 project status](docs/archive/PROJECT_STATUS_2026-03-10.md) is preserved
  byte-for-byte as history, including its old proposed next steps.
- The [post-67.1 snapshot](docs/notes/foundations/post_67_1_status_snapshot.md)
  is historical; its earlier HBP readiness descriptions are superseded.
- The [one-bit scale falsifier M0 record](experiments/emergence_ladder/one_bit_scale_falsifier_01/README.md)
  preserves the missing-theory gate as it stood at that time. The later completed
  revision does not retrospectively execute that falsifier or turn its draft
  into a frozen contract.
