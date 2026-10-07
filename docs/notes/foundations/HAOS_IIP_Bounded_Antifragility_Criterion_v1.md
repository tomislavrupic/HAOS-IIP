# HAOS-IIP Bounded Antifragility Criterion v1

Status:

- downstream operational derivation;
- no modification of HAOS core or frozen telemetry;
- no new physics or universal antifragility claim;
- paired with the finite toy sidecar `HAOS-AF-AR-01`.

## 1. Result of the Derivation

Recoverability does not imply antifragility.

A system can return exactly to the same viable state after every bounded
perturbation while gaining no future capacity. Such a system is recoverable,
but its post-exposure improvement is zero. A system can also remain nearly
unchanged under stress; that is robustness, not antifragility.

HAOS-IIP can nevertheless define a bounded antifragility criterion by adding a
cross-exposure update while keeping the interaction and audit contract frozen.
The invariant contract is the ruler that makes improvement distinguishable
from a silent change of rules.

## 2. Frozen Objects

Let:

- `C` be a frozen evaluation contract;
- `S_0` be the pre-exposure interaction system;
- `P_e` be a bounded exposure perturbation;
- `O_e` be the information legitimately observable during exposure;
- `A(S_0, O_e)` be a persistent adaptive update producing `S_1`;
- `H` be a predeclared held-out perturbation distribution;
- `J_C(S, P)` be recovery capacity under the unchanged contract `C`.

The update may change allowed state, memory, or interaction weights. It may not
change the metric, thresholds, held-out schedule, resource accounting, or claim
language after outcomes are known.

## 3. Antifragility Functional

Define the bounded held-out gain:

```text
AF_C(S_0, P_e, A; H)
  = E_{P_h ~ H}[J_C(A(S_0, O_e), P_h) - J_C(S_0, P_h)].
```

A scalar positive mean is insufficient. Promotion additionally requires:

1. **exposure causality**: sham exposure does not produce the gain;
2. **held-out generalization**: `P_h` is unavailable to the update and differs
   by declared location, seed, severity, or family;
3. **identity preservation**: the declared system identity remains above its
   frozen threshold;
4. **functional preservation**: the target function is preserved or improved;
5. **resource invariance**: improvement is not purchased by an undeclared
   increase in capacity, energy, payload, or intervention budget;
6. **control separation**: the target beats passive, signal-blind or shuffled,
   and reversed-update controls;
7. **tail safety**: mean improvement does not hide a worse declared tail or an
   earlier collapse threshold;
8. **durability**: the improvement persists for the declared evaluation window;
9. **reproducibility**: the complete schedule and decision are recomputable
   from committed artifacts.

When these gates pass, the bounded classification is:

```text
identity-preserving improvement of held-out recoverability caused by prior
bounded stress under an unchanged interaction and audit contract.
```

## 4. Strict Hierarchy

The operational hierarchy is strict:

```text
robustness
  -> recoverability
  -> adaptive recovery
  -> bounded antifragility.
```

- Robustness permits little immediate degradation.
- Recoverability permits degradation followed by re-entry into a viable class.
- Adaptive recovery permits a stress-triggered internal update.
- Bounded antifragility requires that update to improve later held-out recovery
  without violating identity, function, resource, control, or tail gates.

No lower rung logically entails a higher one.

## 5. Why Interaction Invariance Matters

Antifragility requires change, so it cannot mean that every system variable is
invariant. The invariant object is the comparison contract and the recognizable
system identity. The adaptive object is the permitted state, memory, or
interaction structure.

If the score, stress distribution, resource budget, or identity definition
changes after exposure, apparent improvement is not measurable antifragility;
it may be target movement. If identity or function is lost, the result may be
replacement, mutation, or drift rather than improvement of the same system.

## 6. HAOS-AF-AR-01 Instantiation

The first sidecar uses a weighted ring with local backup interactions. Exposure
to primary-edge damage transfers a bounded fraction of the same total
interaction-weight budget from primary edges to local distance-two backups.
The update stores only one scalar backup-budget fraction and does not store the
damage location or future test information.

Held-out evaluation includes displaced primary damage, distributed primary
damage, and backup-channel damage. The recovery-capacity score combines
algebraic connectivity and global effective conductance. Separate gates audit
operator identity, primary-capacity retention, total budget, tail gain,
collapse ordering, dose response, sham response, and matched controls.

The sidecar's maximum claim is finite toy support for this operational
criterion. A passing result would not show that HAOS-IIP generally generates
antifragility, that the update is optimal, or that any empirical system is
antifragile.

## 7. Falsifiers

The bounded claim fails if any of the following occurs:

- held-out gain is absent or confidence intervals include the frozen null;
- gain disappears outside the exposure location;
- shuffled or reversed updates perform equivalently;
- sham exposure produces the same update;
- identity or function falls below threshold;
- the resource budget increases;
- the lower tail regresses or collapse occurs earlier;
- missing rows, post-hoc thresholds, or test leakage enter the bundle.

## 8. Claim Boundary

This note derives an operational test from HAOS-IIP's recoverability and
claim-gating discipline. It does not derive antifragility as a necessary law of
interaction, prove a universal theorem, or authorize physical, biological,
economic, autonomous, or ontological claims.
