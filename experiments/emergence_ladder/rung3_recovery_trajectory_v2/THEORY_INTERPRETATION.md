# EL-R3-RT-02 Theory Interpretation

Final classification: `PARTIAL_RECOVERY_ONLY`.

## Supported Mechanistic Signal

One-bit local relational memory plus compatibility-error feedback performs real
constraint repair:

- median state-region gain: `0.997665`;
- median trajectory-corridor gain: `0.506329`;
- median edge-sign identity: `1.000000`;
- median node-rank identity: `0.978777`;
- median operator-restoration score: `0.919147`;
- paired target-minus-passive gain CI: `[0.997453, 0.998072]`.

This is not passive relaxation. Passive gain is zero. Operator-only filtering
repairs some local smoothness but loses trajectory and identity. Signal-blind
correction is unstable and consumes energy in the wrong direction. The oracle
demonstrates that the evaluation can register full recovery but is excluded
from evidence.

## Why Rung 3 Did Not Pass

No final target run satisfied every recovery dimension simultaneously.

- functional recovery median: `0.000000`;
- recovery rate: `0.000000`;
- all `96` target runs classify as `identity-preserving adaptation`, not
  organizational-regime recovery;
- qualifying recovery basin: empty;
- first failure threshold: `0.25`, the smallest evaluated magnitude;
- validation gate: failed;
- operator-only and signal-blind superiority gates fail at the categorical
  recovery-rate level because target and admissible controls all have zero fully
  recovered runs;
- memory and redundancy ablations reduce raw relational gains, but cannot
  degrade a target recovery rate that is already zero;
- the trivial-attractor rejection gate cannot promote a non-recovering target,
  even though the trivial control visibly collapses variance.

The mechanism restores relational ordering but does not reconstruct the
predeclared function or reliably re-enter the full trajectory corridor. Local
orientation bits retain identity information, but not enough amplitude and
functional information to qualify as organizational-regime recovery.

The reported median latency of `20` steps is the first state-error plus identity
re-entry time. Because no run passes all dimensions, full recovery latency is
unreached rather than `20`.

## Theory Narrowing

- Memory was present, compressed, distributed, and non-oracular.
- Constraint signals were sufficient for relational repair but insufficient
  for functional restoration.
- Restoration did not conflict with identity; identity was the strongest
  surviving dimension.
- Dynamics were not too dissipative in the target, but the one-bit code was too
  coarse to specify the required functional manifold.
- The system therefore has persistence plus endogenous relational repair, not
  demonstrated recoverability.

## Next Distinct Mechanism

The next legitimate Rung 3 family is distributed parity reconstruction with a
finite error-correction radius. It must encode function-relevant information
without storing a continuous state checkpoint, and it must retain the same
passive, filtering, topology, trivial-attractor, signal-blind, parameter-budget,
oracle, and ablation boundaries.

Do not retune RT-02 margins, gains, functions, or thresholds. Its partial result
is frozen.
