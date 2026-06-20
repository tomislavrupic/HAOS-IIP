# Moonshine-Betti Bridge Report

## Status

- Result: `PASS`
- Classification: `SHARED_SUPPORT_DIAGNOSTIC_COUPLING`
- Result hash: `moonshine_betti_bridge_a7d2b62c0c8be6d8c6b757e9`
- Moonshine input hash: `moonshine_diag_3ea77659b5d076fc6c03e4c6`
- Betti vector input hash: `betti_vector_1d5373fa671fd513fe2d5ac0`

## Bridge Type

`shared-support diagnostic coupling`

Useful phrase: `shared-support perturbation covariance`.

Dangerous phrase: `Moonshine-Betti mechanism`.

## Labels

- `MOONSHINE_BETTI_BRIDGE_MANIFEST_BUILT`
- `SHARED_SUPPORT_DIAGNOSTIC_COUPLING_BUILT`
- `SHARED_SUPPORT_PERTURBATION_COVARIANCE_REPORTED`
- `MOONSHINE_BETTI_MECHANISM_NOT_ESTABLISHED`
- `MOONSHINE_PROOF_NOT_ESTABLISHED`
- `PHYSICAL_BRIDGE_NOT_ESTABLISHED`
- `CLAIM_GATED_SHARED_SUPPORT_DIAGNOSTIC`
- `SHARED_SUPPORT_COVARIANCE_PASS`

## Matched Perturbation Response

| Case | Moonshine condition | Betti condition | Moonshine total | Betti vector | Betti edge | Status |
| --- | --- | --- | ---: | ---: | ---: | --- |
| self_recovery | known_positive | known_positive | 0 | 0 | 0 | `PASS` |
| shared_support_replacement | nonsupersingular_prime_control | support_replacement_control | 0.522389629805 | 7 | 0.175 | `PASS` |

## Boundary

This bridge records whether a declared shared-support perturbation moves both bounded diagnostic channels.
It does not derive the Betti graph from Moonshine data, prove Monstrous Moonshine, construct a Monster action, or establish a physical bridge.

Unmatched channel-specific controls remain visible in `moonshine_betti_bridge_result.json` and must not be compressed into a mechanism claim.
