# Geometry Bridge Chain

This folder contains a synthetic, frozen bridge chain that separates four
questions:

1. intrinsic geometry recovery
2. transformation semantics
3. probability-rule recovery
4. observable prediction
5. chain-level audit
6. hidden synthetic geometry recovery
7. hidden synthetic transformation recovery
8. hidden synthetic probability-rule recovery
9. Gaussian-prime norm-lift arithmetic geometry

The chain is intentionally bounded. It is not a Bell experiment and it does
not claim a physical mechanism.

## Rungs

### GEO-01 - Synthetic intrinsic geometry recovery

- path: `synthetic_intrinsic_geometry_recovery/`
- purpose: test whether operator-only transport features recover latent
  intrinsic geometry better than conventional graph/spectral baselines
- current state: open / non-dominant on holdout
- useful signal: pairwise geometry exists, but the coarse holdout target does
  not beat the best baseline
- failure localization: [failure_localization_note.md](./synthetic_intrinsic_geometry_recovery/failure_localization_note.md)
- spiral localization: [spiral_failure_localization.md](./synthetic_intrinsic_geometry_recovery/spiral_failure_localization.md)
- spiral feature specificity study: [spiral_feature_specificity_study.md](./spiral_feature_specificity_study.md)
- spiral feature specificity result: [spiral_feature_specificity_report.md](./spiral_feature_specificity_report.md)

### GEO-02 - Synthetic transformation semantics recovery

- path: `synthetic_transformation_semantics_recovery/`
- purpose: test whether frozen transport/holonomy observables distinguish
  identity, pure gauge, flux, and orientation reversal
- current state: open
- useful signal: transport/holonomy observables are distinct enough to support
  a transformation-semantics layer

### GEO-03 - Synthetic probability-rule recovery

- path: `synthetic_probability_rule_recovery/`
- purpose: test whether a frozen probability rule over transport/holonomy
  features predicts held-out transformation classes better than a null baseline
- current state: open
- useful signal: the frozen rule beats the null on probabilistic scores, but it
  is not yet strong enough for a closed claim

### GEO-04 - Synthetic observable prediction

- path: `synthetic_observable_prediction/`
- purpose: turn the frozen probability rule into an observable forecast on a
  held-out synthetic target
- current state: valid
- useful signal: the observable prediction now clears holdout class accuracy
  and pairwise ordering against null on the synthetic target

### GEO-HIDDEN-01 - Synthetic hidden geometry recovery

- path: `synthetic_hidden_geometry_recovery/`
- purpose: test whether frozen operator and diffusion features recover
  intrinsic distance, orientation, transformation class, and held-out
  relations from a hidden synthetic geometry
- current state: open
- useful signal: distance, orientation, and held-out relations recover well;
  transformation recovery remains the weak point and is the declared open
  criterion
- spectral hardening summary: [spectral_diagnostics_summary.md](./spectral_diagnostics_summary.md)
- terminal reading: `BENCHMARK_OPEN`, `TRANSFORMATION_RECOVERY_BOUNDARY_OPEN`,
  `TRANSFORMATION_CLASS_NOT_ROBUST_ON_HOLDOUT`

#### Spectral Hardening Pass (GEO-HIDDEN-01)

- normalized Laplacian variants now run as the primary spectral path:
  `L_sym` and `L_rw`
- low-mode diagnostics now report Fiedler sign stability separately from holdout
  transformation transfer
- Cheeger sweep conductance is recorded as a bottleneck diagnostic
- current frozen numbers:
  - `transform_accuracy`: `0.500000`
  - `fiedler_transform_accuracy`: `0.250000`
  - `fiedler_sym_accuracy`: `0.250000`
  - `fiedler_rw_accuracy`: `0.312500`
  - `fiedler_sign_stability`: `0.513889`
  - `cheeger_conductance`: `0.304962`
- interpretation: the diagnostic layer improved, but
  `TRANSFORMATION_RECOVERY_BOUNDARY_OPEN` remains the correct terminal reading
- summary note: [spectral_diagnostics_summary.md](./spectral_diagnostics_summary.md)

### GEO-GP-01 - Gaussian-prime norm-lift geometry

- path: `gaussian_prime_norm_lift/`
- purpose: test whether Gaussian-prime ramified / inert / split classes and
  the norm lift `N(a + bi) = a^2 + b^2` produce stable shell, spectral,
  norm-lift, and cochain-flow telemetry under matched controls
- current state: valid synthetic arithmetic calibration
- useful signal: rotation/self controls remain invariant; class shuffling
  degrades prime-shell semantics; norm shuffling moves norm-lift and weighted
  spectral telemetry; topology destruction moves spectral and cochain telemetry
- report: [gaussian_prime_norm_lift_report.md](./gaussian_prime_norm_lift/gaussian_prime_norm_lift_report.md)
- terminal reading: `GAUSSIAN_PRIME_NORM_LIFT_BUILT`,
  `COMPONENT_CONTROLS_PASS`, `PHYSICAL_BRIDGE_NOT_ESTABLISHED`,
  `MONSTER_MOONSHINE_NOT_TESTED`

### GEO-T1-01 - Synthetic hidden transformation recovery

- path: `synthetic_hidden_transformation_recovery/`
- purpose: test whether a frozen observer can recover identity, inverse,
  composition closure, orbits, stabilizers, equivalence classes, holdout
  compositions, and refinement persistence from a hidden finite
  transformation system
- current state: valid
- useful signal: identity, inverse, composition closure, orbit structure,
  stabilizer size, equivalence classes, held-out compositions, and refinement
  persistence recover from the hidden transformation algebra

### GEO-P1-01 - Synthetic hidden probability rule recovery

- path: `synthetic_hidden_probability_rule_recovery/`
- purpose: test whether a frozen observer can recover a hidden Bernoulli law
  from transformation-conditioned observations with frozen holdout, null
  controls, and calibration
- current state: valid
- useful signal: the benchmark recovers a synthetic probability rule with
  frozen holdout transfer, calibration, and null rejection

## Current chain reading

Implemented:

`geometry -> transformation semantics -> probability rule -> observable prediction`

Verified only in a synthetic calibration sense. The chain is present, and the
observable layer now clears its frozen synthetic holdout. That still does not
justify a stronger Bell claim.

The hidden probability-rule benchmark is a separate synthetic calibration task.
It checks whether the probability layer itself can be recovered before any
observable-prediction or Bell-adjacent questions are revisited.

The hidden-geometry benchmark is a separate synthetic calibration task. It is
intended to check the pregeometric layer before any Bell-adjacent questions are
revisited.

The hidden-transformation benchmark is a separate synthetic calibration task.
It checks whether the algebra of change can be recovered from concealed
structure before any Bell-adjacent questions are revisited.

The chain-level audit at `run_geometry_chain_audit.py` keeps the boundary
explicit. It treats the chain as pregeometric / proto-geometric, not as a Bell
mechanism.

## Chain Snapshot

| Rung | Current state | Terminal reading | Short reading |
| --- | --- | --- | --- |
| GEO-01 | `open / non-dominant on holdout` | `GEOMETRY_NOT_DEMONSTRATED`, `MIXED_OPEN` | Pairwise geometry exists, but the coarse holdout target does not beat the best baseline. |
| GEO-02 | `open` | `TRANSFORMATION_SEMANTICS_OPEN`, `MIXED_OPEN` | Transport / holonomy observables separate classes, but the transformation layer is still open. |
| GEO-03 | `open` | `PROBABILITY_RULE_OPEN`, `MIXED_OPEN` | The frozen rule beats the null probabilistically, but not enough for closure. |
| GEO-04 | `valid` | `OBSERVABLE_PREDICTION_VALID` | Observable prediction clears holdout class accuracy and pairwise ordering against null. |
| GEO-HIDDEN-01 | `open` | `BENCHMARK_OPEN`, `TRANSFORMATION_RECOVERY_BOUNDARY_OPEN`, `TRANSFORMATION_CLASS_NOT_ROBUST_ON_HOLDOUT` | Distance, orientation, and held-out relations recover well; transformation recovery remains the weak point. |
| GEO-GP-01 | `valid synthetic arithmetic calibration` | `GAUSSIAN_PRIME_NORM_LIFT_BUILT`, `COMPONENT_CONTROLS_PASS`, `PHYSICAL_BRIDGE_NOT_ESTABLISHED` | Gaussian-prime shell semantics, norm-lift spectral telemetry, and cochain-flow diagnostics respond to matched arithmetic controls. |
| GEO-T1-01 | `valid` | `BENCHMARK_VALID` | Hidden transformation algebra recovers identity, inverse, composition closure, orbits, stabilizers, equivalence classes, holdout compositions, and refinement persistence. |
| GEO-P1-01 | `valid` | `PROBABILITY_RULE_VALID` | Hidden Bernoulli-law recovery passes calibration, holdout transfer, and null rejection. |

## Boundary

This folder is an instrumentation scaffold. It preserves the separation between
synthetic calibration and unsupported physics claims.
