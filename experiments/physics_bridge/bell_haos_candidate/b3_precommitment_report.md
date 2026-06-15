# B3 Precommitment Contract

Implemented fact: this report is generated before B3 scoring from the same contract used to authorize execution.
Design choice: B3 tests one explicit joint-closure-cost candidate rather than generic global closure language.
Heuristic: the cost uses frozen Phase 18 closure features and chain-order token displacements.
Unverified hypothesis: this may or may not produce CHSH violation; no HAOS Bell derivation is established.

## Authorization
- status: `AUTHORIZED_PRECOMMITMENT`
- authorization: `IMPLEMENTATION_AUTHORIZED_WITH_EXPLICIT_INVOICE`
- implementation authorized: `True`
- contract hash: `b3_precommitment_contract_74db225269ebe1705396f518`

## Joint Probability Model
- product distribution: `P(Q=q | a,b,lambda) = exp(-beta*C(q,a,b,lambda)) / sum_{r in {-1,+1}} exp(-beta*C(r,a,b,lambda))`
- local registration: `P(A=x,B=y | a,b,lambda) = 0.5 * P(Q=x*y | a,b,lambda)`
- cost: `C(q,a,b,lambda) = base_instability(lambda) + 0.5*(1 - q*preferred_product(a,b,lambda)*closure_strength(lambda))`
- preferred product: `sign of chain-order displacement between setting token pairs`

## Assumption Invoice
- selected invoice: `FACTORIZATION_RELAXED; OUTCOME_INDEPENDENCE_RELAXED`
- factorization: `RELAXED_BY_EXPLICIT_JOINT_OUTCOME_PRODUCT_DISTRIBUTION`
- measurement independence: `PRESERVED_AND_TESTED_BY_SOURCE_SETTING_MI`
- parameter independence: `NOT_RELAXED; TESTED_BY_NO_SIGNALLING_MARGINALS`
- outcome independence: `RELAXED_BY_PRODUCT_LEVEL_JOINT_DISTRIBUTION`
- forward-only causality: `NO_SPACETIME_CAUSAL_CLAIM; COMPUTATIONAL_BATCH_CANDIDATE_AFTER_SETTINGS_ASSIGNED`

## Falsification Conditions
- NO_SIGNALLING_FAIL invalidates a standard Bell-nonlocal interpretation
- SETTING_INDEPENDENCE_FAIL means measurement independence is not preserved
- POST_SELECTION_DETECTED disqualifies the run
- TARGET_LEAKAGE_DETECTED disqualifies the run
- S<=2 or CI overlapping 2 means CHSH violation is not detected
- matched local controls exceeding the bound robustly indicates implementation leakage or scoreboard failure

## Boundary
- This is not a physical Bell experiment.
- This is not a loophole-free Bell test.
- This is not a derivation of quantum mechanics or the Born rule.
- A clean `S <= 2` remains a valid result.
