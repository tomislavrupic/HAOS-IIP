# Geometry Bridge Recoverability Dashboard

## Status

- Result: `PASS_WITH_FRAGILITY`
- Classification: `COMPOSITE_DIAGNOSTIC_DASHBOARD_ONLY`
- Result hash: `geometry_bridge_dashboard_dc14cbaf7f66529e66a57f49`

## Labels

- `GEOMETRY_BRIDGE_RECOVERABILITY_DASHBOARD_BUILT`
- `MOONSHINE_ARITHMETIC_CHANNEL_INCLUDED`
- `GAUSSIAN_PRIME_SUPPORT_CHANNEL_PARTIAL`
- `BETTI_VECTOR_CHANNEL_INCLUDED`
- `THRESHOLD_ROBUSTNESS_INCLUDED`
- `NULL_RARITY_INCLUDED`
- `CROSS_CHANNEL_COVARIANCE_INCLUDED`
- `FORMAL_TARGET_STATUS_INCLUDED`
- `CLAIM_BOUNDARY_CHECKS_PASS`
- `GLOBAL_CLAIM_OPEN`
- `PHYSICAL_BRIDGE_NOT_ESTABLISHED`
- `PASS_WITH_FRAGILITY`

## Channels

| Channel | Status | Classification |
| --- | --- | --- |
| Moonshine arithmetic diagnostic | `PASS` | `ARITHMETIC_SUPPORT_DIAGNOSTIC` |
| Gaussian-prime norm-lift support | `PARTIAL` | `REPRESENTATIVE_INPUT_TO_BETTI_GRAPH` |
| Betti_0 / Betti_1 diagnostic | `PASS` | `TOPOLOGICAL_DIAGNOSTIC_SIDECAR` |
| Threshold robustness | `PASS_WITH_FRAGILITY` | `LOCAL_ROBUSTNESS_BAND` |
| Null ensemble rarity | `OPEN` | `FALSE_POSITIVE_RATE_ESTIMATE` |
| Cross-channel perturbation covariance | `PASS` | `SHARED_SUPPORT_DIAGNOSTIC_COUPLING` |
| Pinned source validation | `PASS` | `PINNED_SOURCE_VALIDATION_ONLY` |
| Formal Lean targets | `OPEN` | `FORMAL_TARGET_SCAFFOLD_ONLY` |
| Failure ledger | `PASS` | `RECOVERABILITY_FAILURE_CONTRACT` |
| Claim-boundary checks | `PASS` | `NO_FORBIDDEN_UPGRADE_LABELS` |

## Interpretation

`PASS_WITH_FRAGILITY` means the local executable diagnostics recover declared structure, but null rarity, standalone norm-lift recovery, and local Lean certification remain incomplete or open.

## Boundary

This dashboard aggregates bounded diagnostics only. It does not prove Monstrous Moonshine, recover a Monster action, establish a physical bridge, or derive continuum, quantum, or gravity structure.
