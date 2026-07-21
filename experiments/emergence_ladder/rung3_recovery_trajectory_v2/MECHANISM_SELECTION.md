# EL-R3-RT-02 Mechanism Selection

Status: frozen before implementation and final evaluation.

Selection rule: rank mechanism families by independence from the recovery
metric, clarity, leakage resistance, falsifiability, locality, infrastructure
compatibility, parameter economy, distinction from passive relaxation, ability
to form a recovery basin, and value for later closure/cross-scale tests. No
candidate was executed before this ranking was frozen.

## Candidate A - Local Constraint-Error Feedback

- Motivation: a cochain-compatible state satisfies local edge relations. A
  disruption creates measurable local compatibility errors.
- Available information: current endpoint values, locally stored signs of edge
  differences, edge incidence, and current inequality violations.
- Update: `x <- x + eta D^T M[s max(0,m-sDx)] / degree`, where `D` is the local
  incidence map, `s` is one-bit relational orientation memory, `m` is a global
  non-state-specific separation margin, and `M` masks unavailable or already
  satisfied constraints.
- Privileged reference: no full state or future trajectory; local relational
  memory is stored before intervention and fixes state only up to a global
  offset.
- Scale: local updates with a mesoscopic redundant constraint graph.
- Free parameters: correction gain and inequality margin.
- Circularity risk: medium-low. The mechanism minimizes compatibility error;
  trajectory fidelity is evaluation-only.
- Expected failure: relation memory corruption, wrong topology, insufficient
  redundancy, overly large perturbation, or unstable gain.
- Destructive control: relation permutation, topology rewiring, signal-blind
  direction, and memory ablation.
- HAOS-IIP relation: cochain compatibility, local closure inconsistency, and
  perturbation/recovery logic.
- Cost: low.
- Negative-result value: high; it tests whether local structural memory is
  sufficient for endogenous repair.

## Candidate B - Distributed Parity Reconstruction

- Motivation: replicated parity checks can reconstruct discrete local damage.
- Available information: local binary checks and neighboring votes.
- Update: asynchronous bit/phase flips that reduce violated parity checks.
- Privileged reference: no full state, but a fixed codebook is required.
- Scale: local/distributed.
- Free parameters: parity layout, update schedule, and vote threshold.
- Circularity risk: low.
- Expected failure: correlated errors beyond code distance.
- Destructive control: parity-check deletion and check-label permutation.
- HAOS-IIP relation: discrete invariant redundancy and branch identity.
- Cost: medium.
- Negative-result value: high, but the discrete codebook is less continuous
  than the current cochain transport regime.

## Candidate C - Adaptive Operator Repair

- Motivation: disruption-induced inconsistency could temporarily modify local
  edge weights and then relax them when consistency returns.
- Available information: edge residuals and local weight state.
- Update: coupled state/weight dynamics with residual-driven edge adaptation.
- Privileged reference: nominal local weights must be stored.
- Scale: local to mesoscopic.
- Free parameters: state gain, weight gain, weight decay, clipping, and timescale.
- Circularity risk: medium.
- Expected failure: adaptation chases noise, locks in damaged topology, or
  creates oscillation.
- Destructive control: freeze, randomize, or reverse weight adaptation.
- HAOS-IIP relation: operator hierarchy and perturbation-sensitive transport.
- Cost: high.
- Negative-result value: high, but too many parameters make first attribution
  difficult.

## Candidate D - Slow-Memory Homeostasis

- Motivation: slow variables may retain identity while fast variables absorb
  shocks.
- Available information: current fast variables and locally coupled slow traces.
- Update: fast relaxation plus slow-variable feedback.
- Privileged reference: no exact checkpoint, but slow traces approximate a
  pre-intervention manifold.
- Scale: local with two timescales.
- Free parameters: two gains, two timescales, coupling, and saturation.
- Circularity risk: medium-high because the slow manifold can become a hidden
  reference if not tightly constrained.
- Expected failure: frozen slow layer, memory drift, or feedback oscillation.
- Destructive control: slow-layer shuffle, disconnect, or timescale inversion.
- HAOS-IIP relation: persistence and prospective cross-scale organization.
- Cost: medium-high.
- Negative-result value: high, but leakage interpretation is harder than A.

## Frozen Ranking

Scores are `1` (weak) to `5` (strong). Parameter economy is higher when fewer
free parameters are required.

| Mechanism | Metric independence | Clarity | Leakage resistance | Falsifiability | Locality | Infrastructure fit | Parameter economy | Beyond passive | Recovery basin | Later-rung value | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| A Constraint-error feedback | 5 | 5 | 4 | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 49 |
| B Distributed parity | 5 | 5 | 4 | 5 | 5 | 3 | 3 | 5 | 5 | 4 | 44 |
| C Adaptive operator repair | 4 | 4 | 4 | 5 | 5 | 5 | 2 | 5 | 4 | 5 | 43 |
| D Slow-memory homeostasis | 4 | 4 | 3 | 5 | 4 | 4 | 2 | 5 | 4 | 5 | 40 |

Selected: **Candidate A - local constraint-error feedback with redundant edge
memory**.

The selection is based on mechanism clarity and leakage resistance, not a
preliminary recovery score. Candidate B remains the next genuinely distinct
family if A fails.
