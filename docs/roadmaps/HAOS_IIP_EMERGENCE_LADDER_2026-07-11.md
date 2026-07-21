# HAOS-IIP Emergence Ladder - 2026-07-11

Status: evidence-derived project roadmap and claim boundary.

This ladder separates increasingly difficult claims. A result at one rung does
not promote a higher rung. Classifications describe repository evidence, not an
ontology.

## One-Page Ladder

| Rung | Operational question | Classification | Confidence | Decisive reading |
| ---: | --- | --- | --- | --- |
| 0 | Can the instruments measure their declared mechanisms reliably? | `PARTIALLY_SUPPORTED` | High | Frozen phase checks and HBP-IR-01 are valid; historical HBP PB-01 to PB-04 remain quarantined. |
| 1 | Is organization distinguishable from matched nulls? | `SUPPORTED` | High | Stable Phase V recovery distributions separate strongly from degraded and shuffled controls. |
| 2 | Does organization persist under bounded perturbation? | `SUPPORTED` | High | Phase XI maps survival basins, perturbation thresholds, and repeatable failure channels across refinements. |
| 3 | Does the system recover toward a constrained organizational region after disruption? | `NEGATIVE_RESULT` | High for tested mechanisms | RT-01 is negative; RT-02 repairs relational identity only; RP-01 adds finite-radius parity repair but fails the functional validation gate and does not consume final seeds. |
| 4 | Do internal relations predict future evolution beyond external descriptions? | `NEGATIVE_RESULT` | High for current mechanism | EL-R4-OC-01 shows aggregate advantage but its run-level confidence interval includes zero. |
| 5 | Are local, mesoscopic, and global structures coupled nontrivially? | `PARTIALLY_SUPPORTED` | Medium | Refinement, propagation, temporal, and distance-surrogate structure exist; explicit cross-scale dependency transfer is missing. |
| 6 | Do macrovariables provide a stable, efficient future description? | `PARTIALLY_SUPPORTED` | Medium-low | Effective propagation and EOT-01 synthetic term recovery exist; equal-information out-of-sample macro prediction is absent. |
| 7 | Does structure predict across unseen contexts? | `INSTRUMENT_INVALID` | High | Historical HBP predictive artifacts have invalid controls or model paths; HBP-IR-01 repairs only the instrument scaffold. |
| 8 | Do validated emergent units compose into new stable units? | `NOT_TESTED` | High | Interaction and collective phases do not perform decomposition/recomposition against aggregation controls. |
| 9 | Do distinct microscopic families converge on shared macro behavior? | `INCONCLUSIVE` | Medium-high | CP5 and universality-facing diagnostics remain open and family coverage is insufficient. |
| 10 | Do effective structures converge under size and resolution change? | `INCONCLUSIVE` | High | Phase X supports one bounded extension; CP2 and later same-surrogate convergence gates remain open. |
| 11 | Does HAOS-IIP make an empirically discriminating physical prediction? | `NOT_TESTED` | High | Reference sidecars and bounded external probes do not contain an independent HAOS-derived physical prediction. |

## Emergence Frontier

The highest genuinely supported rung is **Rung 2 - Persistence under
perturbation**.

The first unsupported rung is **Rung 3 - Recoverability**. Phase V and Phase XI
contain useful recovery-score and persistence evidence, but they do not compare
full post-intervention trajectories against passive relaxation, ordinary
operator filtering, and trivial attractors. EL-R3-RT-01 supplied that missing
comparison and returned `NEGATIVE_RESULT`. EL-R3-RT-02 then tested independently
motivated local relation feedback and returned `PARTIAL_RECOVERY_ONLY`.
EL-R3-RP-01 tested distributed local parity and stopped at
`VALIDATION_GATE_FAILED`: relation and operator repair succeed inside the code
radius, but the independent function does not recover. None promotes the rung.

The minimum blocker is primarily **theory and information architecture**.
One-bit orientations and local parity over order relations store enough
information for identity-preserving constraint repair, but not the mesoscopic
amplitude/phase information required by the independent function. More gain,
decoder iterations, or threshold tuning cannot supply the missing information.

## Rung Details

### Rung 0 - Instrument Integrity

- Operational definition: deterministic, mechanism-valid measurement with
  canonical features, explicit missingness, no leakage, same-path controls, and
  semantic checkers.
- Repository evidence: phase checkers, artifact hashes, public reproduction,
  CST failure localization, and HBP-IR-01.
- Relevant areas: frozen phases 3-19, `run_phase.py`, CST, HBP.
- Strongest supporting artifact:
  `experiments/hbp/integrity_repair_v2/results/hbp_ir_01_result.json`.
- Strongest counterevidence: PB-01 through PB-04 contain invalid controls,
  holdout selection, baseline aliasing, or incomplete mechanism routing.
- Classification: `PARTIALLY_SUPPORTED`.
- Confidence: high.
- Missing validation: migrate a real successor benchmark through the repaired
  instrument without reusing historical result identifiers.
- Next decisive experiment: use HBP-IR-01 on one new, precommitted operational
  dataset with untouched holdout.
- Non-claim: instrument validity is not predictive validity.

### Rung 1 - Non-Random Organization

- Operational definition: organization separates from matched nulls with
  effect sizes across repeated deterministic schedules.
- Repository evidence: Phase V stable, degraded, and shuffled recovery
  distributions; spectral and invariant control comparisons in later phases.
- Relevant phases: V, VII-IX.
- Strongest supporting artifact:
  `phase5-readout/phase5_authoritative_summary.md`.
- Strongest counterevidence: descriptive class identity flips at the recorded
  amplitude-jitter boundary.
- Classification: `SUPPORTED`.
- Confidence: high inside the frozen Phase IV/V regime.
- Missing validation: independent operator families and non-frozen replication.
- Next decisive experiment: repeat the same null separation on a distinct
  operator family without changing the recovery score.
- Non-claim: structured behavior alone is not emergence or physics.

### Rung 2 - Persistence Under Perturbation

- Operational definition: bounded perturbation-response curves retain identity
  until measured failure thresholds are crossed.
- Repository evidence: Phase XI survival basins across stiffness,
  connectivity, phase, and boundary perturbations with matched controls.
- Relevant phases: XI and the Phase V robustness audit.
- Strongest supporting artifact: `phase11-protection/phase11_summary.md`.
- Strongest counterevidence: phase-randomization survival is only `0.444444`
  on the branch, and thresholds are branch-local.
- Classification: `SUPPORTED`.
- Confidence: high for the tested candidate and perturbation grid.
- Missing validation: broader candidates and adversarial perturbations.
- Next decisive experiment: preregister a new candidate family and reuse the
  same perturbation channels unchanged.
- Non-claim: persistence does not establish recovery, closure, or ontology.

### Rung 3 - Recoverability

- Operational definition: after disruption, trajectories return toward a
  constrained invariant region and outperform passive, degraded, and null
  alternatives.
- Repository evidence: Phase V recovery histograms, Phase XI persistence,
  EL-R3-RT-01, EL-R3-RT-02, and EL-R3-RP-01.
- Relevant phases/experiments: V, XI, and `experiments/emergence_ladder/`.
- Strongest supporting artifact: `phase5-readout/phase5_authoritative_summary.md`.
- Strongest new partial signal:
  `experiments/emergence_ladder/rung3_recovery_trajectory_v2/final/aggregate_result.json`.
- Strongest counterevidence: RP-01 correctly decodes all 16 within-radius
  validation targets but has functional restoration median `0.0`, full recovery
  rate `0.0`, and no separation from RT-02 or matched controls. Its validation
  gate fails and final seeds remain untouched.
- Classification: `NEGATIVE_RESULT` at the ladder level; RT-02 is
  `PARTIAL_RECOVERY_ONLY` and cannot promote the rung.
- Confidence: high for the three tested mechanism families.
- Missing validation: an independently justified mesoscopic information
  representation that is function-relevant but cannot reconstruct the full
  continuous state.
- Next decisive experiment: none authorized until EL-R3-THEORY-REVISION-01
  explains what restorative variable stores the missing function-relevant
  information and freezes new destructive controls.
- Non-claim: recoverability does not imply operational closure or a law of
  nature.

### Rung 4 - Operational Closure

- Operational definition: internal relational state predicts the next system
  event better than external metadata and naive sequence baselines, and the
  advantage fails under transition destruction.
- Repository evidence: Phase XVII causal-closure feasibility and EL-R4-OC-01.
- Relevant phases/experiments: XVI, XVII, and
  `experiments/emergence_ladder/rung4_operational_closure/`.
- Strongest supporting artifact:
  `experiments/emergence_ladder/rung4_operational_closure/results/operational_closure_result.json`.
- Strongest counterevidence: holdout aggregate advantage is `0.067901`, but the
  run-level 95% interval is `[-0.037037, 0.046296]`.
- Classification: `NEGATIVE_RESULT` for the current event-transition mechanism;
  this cannot promote the ladder while Rung 3 remains unsupported.
- Confidence: high.
- Missing validation: raw state trajectories, explicit external forcing, and
  intervention-aware internal/external predictor budgets.
- Next decisive experiment: EL-R4-OC-02 over intervention-native state
  trajectories rather than compressed chain labels.
- Non-claim: Phase XVII's name does not by itself establish operational closure.

### Rung 5 - Cross-Scale Organization

- Operational definition: predefined local, mesoscopic, and global variables
  exchange predictive information under a fixed coarse-graining protocol.
- Repository evidence: Phase X refinement descriptors, Phase XV propagation,
  Phase XVI order, and Phase XVIII distance-surrogate scaling.
- Relevant phases: X and XV-XVIII.
- Strongest supporting artifact: `phase18-distance-surrogate/phase18_summary.md`.
- Strongest counterevidence: no scale-shuffled or cross-scale-disrupted control
  tests coordinated recovery; current same-surrogate gates remain open.
- Classification: `PARTIALLY_SUPPORTED`.
- Confidence: medium.
- Missing validation: explicit scale variables, lag/information transfer, and
  controls preserving within-scale statistics while destroying coupling.
- Next decisive experiment: EL-R5-CSR-01 cross-scale recovery coordination,
  gated on a valid Rung 4 mechanism.
- Non-claim: refinement stability is not cross-scale causal organization.

### Rung 6 - Effective Dynamics

- Operational definition: independently defined macrovariables predict future
  behavior under equal information budgets and yield stable effective
  operators on holdout regimes.
- Repository evidence: Phase XV effective propagation descriptors and EOT-01's
  synthetic Laplacian plus squared-Laplacian recovery.
- Relevant areas: Phase XV and `experiments/effective_operator_theory/`.
- Strongest supporting artifact:
  `experiments/effective_operator_theory/results/effective_operator_result.json`.
- Strongest counterevidence: EOT-01 recovers terms used to generate its own
  synthetic data; it is calibration, not artifact-derived discovery.
- Classification: `PARTIALLY_SUPPORTED`.
- Confidence: medium-low.
- Missing validation: macrovariables defined independently of targets,
  equal-dimensional random/compressed baselines, and unseen refinements.
- Next decisive experiment: EOT-02 artifact-derived effective expansion with a
  power gate for the currently sparse refinement levels.
- Non-claim: an EFT-like fitting discipline is not a physical EFT.

### Rung 7 - Cross-Context Predictive Transfer

- Operational definition: calibration-selected structure predicts untouched
  systems, seeds, topologies, or regimes through mechanism-valid controls.
- Repository evidence: PB-01 to PB-04 attempts and HBP-IR-01.
- Relevant area: `experiments/hbp/`.
- Strongest supporting artifact:
  `experiments/hbp/integrity_repair_v2/results/hbp_ir_01_result.json`.
- Strongest counterevidence: historical HBP controls, model selection, and
  candidate/baseline paths are invalid.
- Classification: `INSTRUMENT_INVALID` for historical predictive results.
- Confidence: high.
- Missing validation: one new dataset routed through HBP-IR-01 with an untouched
  holdout and no historical identifier reuse.
- Next decisive experiment: a fresh network-recovery contract after lower-rung
  macrovariables are independently validated.
- Non-claim: HBP-IR-01's synthetic score is not external transfer evidence.

### Rung 8 - Compositional Emergence

- Operational definition: independently validated units combine into a new,
  non-additive macrostate that can be decomposed and recovered hierarchically.
- Repository evidence: Phase XII interactions and Phases XIII-XIV collective
  sectors are precursors only.
- Relevant phases: XII-XIV.
- Strongest supporting artifact: `phase12-interactions/phase12_summary.md`.
- Strongest counterevidence: there is no uncoupled/random/static aggregation
  comparison with decomposition and recomposition.
- Classification: `NOT_TESTED`.
- Confidence: high.
- Missing validation: validated lower units, coupling sweep, non-additivity,
  decomposition, and hierarchical perturbation recovery.
- Next decisive experiment: two-unit composition benchmark after operational
  closure is established.
- Non-claim: interaction or collective behavior alone is not compositional
  emergence.

### Rung 9 - Universality Class Evidence

- Operational definition: distinct microscopic families converge under one
  predefined coarse-graining to shared dimensionless behavior.
- Repository evidence: CP5, scale-bridge comparisons, and multiple synthetic
  domains.
- Relevant areas: continuum-sketch and bridge diagnostics.
- Strongest supporting artifact:
  `continuum-sketch/post_66_5_roadmap_run_summary.md`.
- Strongest counterevidence: family coverage and hardened controls remain
  insufficient; the current 66.5 mechanism is terminal-negative.
- Classification: `INCONCLUSIVE`.
- Confidence: medium-high.
- Missing validation: genuinely distinct microscopic families, competing flow
  hypotheses, and finite-size uncertainty.
- Next decisive experiment: only after Rungs 4-6 provide shared macrovariables.
- Non-claim: visual or spectral similarity is not a universality class.

### Rung 10 - Continuum-Limit Evidence

- Operational definition: effective observables converge under independent
  size, resolution, and discretization sweeps with uncertainty and artifact
  controls.
- Repository evidence: Phase X's one-level extension and later refinement
  summaries.
- Relevant phases: X, XV-XVIII, continuum-sketch CP2.
- Strongest supporting artifact: `phase10-bridge/phase10_summary.md`.
- Strongest counterevidence: same-surrogate recovery and universality gates
  remain open; only a narrow family and few refinements are available.
- Classification: `INCONCLUSIVE`.
- Confidence: high.
- Missing validation: more independent refinements, discretization families,
  uncertainty, and convergence-model comparison.
- Next decisive experiment: finite-size scaling after a stable macro-dynamic
  observable exists.
- Non-claim: compact scaling over one extension is not a continuum limit.

### Rung 11 - Empirically Discriminable Physical Bridge

- Operational definition: a preregistered HAOS-derived observable distinguishes
  itself from established alternatives on independent physical data.
- Repository evidence: Bell, hydrogen, semiconductor reference harnesses and
  bounded external/materials sidecars.
- Relevant areas: physics bridge, materials bridge, external validation.
- Strongest supporting artifact:
  `experiments/physics_bridge/bell_chsh_reference/bridge_status.json` as a
  reference scoreboard only.
- Strongest counterevidence: Bell, hydrogen, and semiconductor mechanism slots
  remain empty; no independent discriminating HAOS prediction exists.
- Classification: `NOT_TESTED`.
- Confidence: high.
- Missing validation: an independently generated observable, competing model
  predictions, preregistration, and external data.
- Next decisive experiment: none until lower rungs establish a transferable
  effective variable and prediction rule.
- Non-claim: reference reproduction is not mechanism derivation or empirical
  confirmation.

## Frontier Candidate Ranking

Scores use `1` (weak) to `5` (strong). Circularity is scored as resistance to
circularity, so higher is better.

| Candidate | Claim importance | Falsifiability | Metric independence | Control quality | Cost | Circularity resistance | Infrastructure reuse | Negative-result value | Unlock value | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| EL-R3-RT-01 trajectory recovery vs passive/trivial controls | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 45 |
| EL-R4-OC-01 internal-vs-external operational closure | 5 | 5 | 4 | 5 | 5 | 4 | 5 | 5 | 5 | 43 |
| EL-R5-CSR-01 cross-scale recovery coordination | 5 | 5 | 5 | 4 | 3 | 4 | 4 | 5 | 5 | 40 |

EL-R3-RT-01 was selected after the hostile ladder read showed that earlier
recovery evidence lacked passive and trivial-attractor trajectory controls. It
returned a valid negative. EL-R4-OC-01 was also executed and remains useful
failure localization, but it cannot be used to skip Rung 3.

## Following-Rung Priority

1. **EL-R3-THEORY-REVISION-01.** Explain and precommit a new mesoscopic
   restorative representation; do not retune RT-02 or RP-01.
2. **EL-R4-OC-02 intervention-native closure.** Execute only after Rung 3 is
   supported; use raw internal state and external forcing under equal budgets.
3. **EL-R5-CSR-01 cross-scale recovery.** Execute only after Rung 4 survives
   intervention; retain scale-shuffled and passive controls.

The next cycle should not retune EL-R3-RT-01, EL-R3-RT-02, EL-R3-RP-01, or
EL-R4-OC-01. Their negative/partial results are frozen.
