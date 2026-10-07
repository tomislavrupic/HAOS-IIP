# HAOS-AF-AR-01 Post-Final Hostile Audit

Status:

- written after the single completed final run;
- does not modify the frozen contract, implementation, final artifacts, result
  hash, or internal classification;
- narrows interpretation only.

## Verified Internal Result

The frozen bundle checker recomputed all `3888` rows and returned:

```text
classification: BOUNDED_TOY_ANTIFRAGILITY_SUPPORTED
result_hash: haos_af_ar_01_4458994589914f9df321a879
checker_status: PASS
```

All fourteen preregistered gates passed. The target's median held-out
recovery-capacity gain over passive was `0.340433542952`; its tenth-percentile
gain was `0.169267214738`; and its minimum row-level gain was
`0.116263749814`. The smallest family median was `0.246176832723` for
backup-channel damage. Minimum operator-identity cosine was `0.954839748452`,
and minimum retained primary capacity was `0.768421052632` under an invariant
total interaction-weight budget.

## Strongest Legitimate Reading

This is a constructive existence witness for the frozen operational criterion:
on this finite graph family, a stress-triggered persistent update can improve
later recovery capacity, preserve the declared identity and resource bounds,
avoid tail regression and earlier collapse, and separate from passive,
equal-budget random, reversed, and sham controls.

It demonstrates that interaction invariance and beneficial adaptation are not
logically incompatible. The invariant object is the audit contract and
recognizable graph identity; the adaptive object is the internal allocation of
interaction weight.

## Structural Limitations

### 1. Rotational Symmetry Weakens Location Holdout

The ring is rotationally symmetric. Displaced damage locations are disjoint in
the schedule, but target and passive scores are invariant under rotation. This
is why their seed-block bootstrap interval is effectively degenerate. The eight
final seeds do not supply eight independent target-versus-passive outcomes.

Family holdout is more informative than location holdout here: distributed
primary damage and backup-channel damage differ structurally from the exposure
family and both retain positive margins. A future claim should use irregular
graphs or heterogeneous local tasks so location generalization becomes
nontrivial.

### 2. The Adaptive Direction Is Designer-Specified

Exposure does not discover a new repair law. The policy is frozen in advance:
greater observed primary damage transfers more weight to uniformly distributed
distance-two backups. The result tests stress-triggered beneficial
reconfiguration, not autonomous learning, mechanism discovery, or optimality.

### 3. Static Hardening Reaches the Same Endpoint

The noncausal static-hardened reference can instantiate the same maximum backup
fraction without exposure. The experiment therefore does not show that damage
is necessary to find the endpoint or that adaptation beats the best static
design. It shows that the declared stress signal causally triggers movement
from the baseline to a more recoverable bounded configuration.

### 4. Metric Alignment Is Built Into the Mechanism Class

Local backup weight is expected to improve algebraic connectivity and effective
conductance. Random and reversed controls show that the allocation pattern and
direction matter, but the target was selected with these graph metrics in mind.
A stronger test requires an independently frozen function that can conflict
with connectivity improvement, plus development and final graph families that
are structurally separated.

### 5. Durability Is One Update Window

The experiment tests one persistent update followed by the declared held-out
stress ladder. It does not test repeated adaptation, forgetting, hysteresis,
resource exhaustion, adversarial stress sequences, or recovery after the backup
structure itself changes.

## Hostile Verdict

Retain the internal classification exactly as frozen, but phrase the public
conclusion as:

```text
HAOS-AF-AR-01 supplies bounded finite-toy support for an operational
antifragility criterion. It does not show that HAOS-IIP generally derives
antifragility or that the tested update is learned, necessary, optimal, or
empirically valid.
```

The next discriminating experiment should use heterogeneous graphs, an
independently frozen task function, repeated exposure/test cycles, and a static
robust-design frontier. Until that passes, the present result remains a
methodological proof of compatibility and testability.
