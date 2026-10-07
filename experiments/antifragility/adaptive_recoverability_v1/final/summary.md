# HAOS-AF-AR-01 Final Summary

Classification: `BOUNDED_TOY_ANTIFRAGILITY_SUPPORTED`

This finite toy experiment tests an operational antifragility criterion. It does not derive general antifragility from HAOS-IIP.

## Demonstrated

- `budget_invariant`
- `collapse_never_earlier`
- `collapse_shift_present`
- `dose_response_nondecreasing`
- `family_generalization`
- `held_out_locations_disjoint`
- `identity_preserved`
- `no_target_regression`
- `primary_capacity_preserved`
- `sham_inert`
- `tail_gain`
- `target_beats_passive`
- `target_beats_random`
- `target_beats_reversed`

## Failed

- None under the frozen finite schedule.

## Core Result

- median held-out recovery-capacity gain: `0.340433542952`
- tenth-percentile gain: `0.169267214738`
- paired target-minus-passive CI: `[0.35377550397, 0.35377550397]`
- paired target-minus-random CI: `[0.157985497913, 0.177898836799]`
- minimum operator-identity cosine: `0.954839748452`
- minimum primary-capacity retention: `0.768421052632`

## Static Reference Boundary

A statically hardened graph can reach the same endpoint without exposure; it is a noncausal ceiling reference, not evidence for the adaptive mechanism.
The experiment therefore attributes a bounded stress-triggered transition, not a unique or optimal design principle.

## Open

- Other graph families, adaptive laws, resource budgets, functions, and stress distributions.
- Durability beyond the declared single-update held-out test schedule.
- Independent implementation and replication.
- Whether the criterion remains discriminative in empirical systems.

## Not Claimed

- No general antifragility, physical law, biological adaptation, economic advantage, agency, universality, continuum limit, or new physics.
- `NO_HIGHER_RUNG_PROMOTION`: this sidecar does not alter the emergence ladder or any frozen parent result.
