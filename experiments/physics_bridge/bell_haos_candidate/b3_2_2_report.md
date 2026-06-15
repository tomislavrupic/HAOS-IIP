# B3.2.2 Semantic and Refinement Specificity Audit

Implemented fact: this harness keeps B3.2 `G`, `J_lambda`, settings, and thresholds frozen.
Design choice: S1 and S2 are independent provenance gates.
Heuristic: semantic affinity is token-overlap plus inverse canonical token distance.
Unverified hypothesis: HAOS-specific Bell geometry remains unestablished unless both gates pass.

## Labels
- SETTING_SEMANTICS_NOT_ESTABLISHED
- REFINEMENT_SPECIFICITY_NOT_ESTABLISHED
- GENERIC_RELATIONAL_GEOMETRY
- CHSH_SCORING_NOT_AUTHORIZED
- HAOS_BELL_DERIVATION_NOT_ESTABLISHED

## S1 Semantic Ordering
- target composite: `0.4222758711105997`
- label permutation composite: `0.08224299065420561`
- margin: `0.3400328804563941`
- status: `SETTING_SEMANTICS_NOT_ESTABLISHED`

## S2 Refinement Persistence
- target pair-order stability: `0.9206652376715403`
- control pair-order stability: `0.8071536560643975`
- target aligned distance: `0.12238126985783403`
- control aligned distance: `0.30334483749399194`
- degradation points: `3`
- status: `REFINEMENT_SPECIFICITY_NOT_ESTABLISHED`

## Boundary
- No CHSH score is computed.
- B3.2 geometry is not modified.
- A pass here would only permit consideration of the next gate, not a Bell derivation claim.
