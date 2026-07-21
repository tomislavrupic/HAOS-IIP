# HAOS-IIP Project Audit - 2026-07-11

Status: repository-state and failure-localization audit

Scope: current checkout at `8d8a49c`, including frozen phases, CST, geometry
and physics bridges, HBP, effective-operator experiments, tests, generated
artifacts, and claim boundaries.

## 1. Executive Summary

The project is coherent at its frozen phase and public reproduction layers.
The phase authority chain, artifact-only reproduction spine, CST negative
result, reference sidecars, synthetic calibration bundles, and their bounded
non-claims reproduce in the inspected environment.

The project is not coherent as one uniformly validated experimental system.
The principal break is in HBP PB-02 through PB-04: result hashes reproduce, but
several controls do not execute their declared perturbations, PB-02 selects and
reconstructs its best baseline through inconsistent paths, draft contracts are
executed as if authorized, and PB-04 compares identically implemented baseline
and HAOS paths. These artifacts cannot support predictive-bridge promotion.

Answers to the required state questions:

- Is the project currently coherent? **Partially.** Frozen phases and bounded
  calibration bundles are coherent; HBP execution is not methodologically
  coherent with its contracts.
- Is it reproducible? **Yes for the checked artifact and test paths.** This is
  computational reproducibility, not scientific validity. Several HBP defects
  reproduce exactly.
- Are its strongest claims supported? **The bounded claims are.** The evidence
  supports deterministic feasibility and negative/open bridge statements. It
  does not support physical correspondence, Bell derivation, hydrogen or
  semiconductor derivation, quantum gravity, continuum ontology, or an
  empirical HAOS mechanism.
- Highest-risk weakness: **HBP hash-valid artifacts can pass checkers even when
  their experimental comparisons and controls are internally invalid.**
- Single most valuable next improvement: **a versioned HBP integrity repair
  that makes contracts executable, selects models without holdout access, and
  routes target and controls through the same model path.**
- What must not be touched: `haos_core` public APIs,
  `telemetry/frozen_metrics.py`, frozen phase manifests, Phase V readout
  definitions, the Phase VI operator hierarchy, the Phase VIII short-time
  window, historical result artifacts, and archived pre-refactor code.

Overall classification:

- frozen phase architecture: `SIGNAL STABLE`
- public reproduction: `RS`
- CST v0.2.2: `OPEN / TARGET NOT DISTINCT / UNDERPOWERED`
- geometry chain: `PARTIAL / OPEN`
- Bell mechanism: `NOT ESTABLISHED`
- HBP predictive interpretation: `CONTROL INVALID / REPAIR REQUIRED`
- effective-operator scaffold: `SYNTHETIC CALIBRATION PASS`
- physical mechanism: `NOT ESTABLISHED`

## 2. Repository and Git State

State before this audit:

- branch: `review/phase18-cst-physicsbridge-hbp`
- HEAD: `8d8a49c feat(effective-operator): add EFT-like expansion scaffold`
- worktree: clean
- staged files: none
- unstaged files: none
- untracked files: none
- merge conflicts: none found
- configured remote: `origin git@github.com:tomislavrupic/HAOS-IIP.git`
- relation to tracked review branch: ahead 12, behind 0
- relation to `origin/main`: ahead 0, behind 0
- local `main` relation to `origin/main`: behind 35, ahead 0

Remote relations above use the configured local remote refs. No fetch, merge,
rebase, commit, or push was performed.

Repository scale:

- tracked files: 5,068
- Python files outside `.venv` and `.git`: 379
- JSON files outside `.venv` and `.git`: 1,608
- checkout size: approximately 2.6 GB
- `Images/`: approximately 770 MB
- `docs/`: approximately 561 MB
- `experiments/`: approximately 223 MB
- `papers/`: approximately 137 MB

Large tracked examples include a 42.7 MB video, a 42.4 MB literature PDF,
multiple 12-20 MB PDFs, and the 14.4 MB Bell `candidate_trials.csv`. Two metric
field JSON files are byte-identical 9.2 MB blobs under different names. This is
not a correctness failure, but it is a long-term reviewability and clone-cost
risk.

The recurring Git fsmonitor warning
`fsmonitor_ipc__send_query: unspecified error` did not alter status results. It
appears environment-local and was not repaired in repository files.

Recent active development is concentrated in:

- `experiments/effective_operator_theory/`
- `experiments/geometry_bridge/monstrous_moonshine_diagnostic/`
- `papers/bridge_literature/`
- `experiments/physics_bridge/literature_component_bridge/`

Apparently unfinished or terminal-open areas include HBP PB-02 through PB-04,
the geometry transformation bottleneck, Bell B3.2.2, CST power sufficiency,
the formal Lean target scaffold, and draft papers containing explicit
placeholder citations or figures.

## 3. Architecture Map

| Layer | Source of truth | Role | Current state |
| --- | --- | --- | --- |
| Repository rules | `AGENTS.md`, `ROOT_EXECUTION_GUIDE.md`, `MIGRATION_MAP.md` | Freeze and import boundaries | Active |
| Frozen core | `haos_core/core.py`, `haos_core/io.py` | Graph, transport, selector, invariant primitives | Frozen |
| Frozen telemetry | `telemetry/frozen_metrics.py` | Recovery, persistence, temporal, causal metrics | Frozen |
| Phase orchestration | `run_phase.py` | Builder/checker dispatch | Active; guard repaired |
| Frozen phase chain | `phase3-stability/` through `phase19-spectral-address/` | Measurement rigs and authority artifacts | Reproducing |
| Public reproduction | `examples/quick_reproduce.py` | Artifact-only Phase XV-XVIII spine | Reproducing |
| Prediction ledger | `predictions/prediction_ledger.json` | Bounded prediction records | Active, external tests mostly open |
| CST | `experiments/cst/` | Closure stability instrument | Frozen experimental v0.2.2 |
| Geometry bridge | `experiments/geometry_bridge/` | Synthetic geometry and arithmetic calibration | Mixed/open |
| Physics bridge | `experiments/physics_bridge/` | Reference harnesses and claim-gated proxies | Reference-valid; derivations open |
| HBP | `experiments/hbp/` | Operational and predictive bridge registry | Registry valid; PB execution flawed |
| Effective operator | `experiments/effective_operator_theory/` | EFT-like synthetic hierarchy extraction | Synthetic pass only |
| External/materials/biology | `external_validation/`, `experiments/materials_bridge/`, `experiments/biology/` | Sidecar validation and demos | Bounded; partially exercised |
| Speculative material | `speculative_bridge/` | Explicit analogy playground | Non-authoritative |
| Historical runtime | `archive-pre-refactor/` | Provenance only | Archived; no active imports found |

The phase import rule is respected in the inspected phase builders: active phase
code imports standard libraries, declared numerical dependencies, and
`haos_core`; it consumes earlier phases through artifacts rather than importing
their implementation modules.

Source-of-truth ambiguity remains in `API_CONTRACT.md`: it describes future
module paths such as `operators/cochain_laplacian.py`, while the implemented and
repository-frozen callable surface is `haos_core`. The document says it is a
contract for later implementation, but `ROOT_EXECUTION_GUIDE.md` calls it the
current callable surface. This should be resolved in a documentation-only
versioned contract clarification.

## 4. Active Experiment and Phase Map

| Subsystem | Inputs | Controls/baselines | Output status | Claim ceiling |
| --- | --- | --- | --- | --- |
| Phases VI-XIX | Frozen manifests and phase-local artifacts | Phase-specific matched controls | Checkers pass where registered | Bounded feasibility only |
| CST v0.2.2 | Frozen Phase XV-XVIII ledgers, two seeds | Component-specific controls | Partially discriminative; target not distinct; underpowered | Experimental closure diagnostic |
| Bell reference | Declared local and quantum reference laws | Local/random controls | Reference sanity pass | Scoreboard only |
| Bell candidate B0-B3.2.2 | Frozen Phase XVIII features and precommitted synthetic mechanisms | Leakage, post-selection, geometry nulls | Bell-local then terminal generic geometry | No Bell derivation |
| Hydrogen reference | Declared Rydberg relation | Malformed formula controls | Reference sanity pass | No HAOS derivation |
| Semiconductor reference | Declared toy direct-gap model | Malformed band controls | Reference sanity pass | Toy reference only |
| Geometry chain | Synthetic hidden geometry and transformations | Permutation, topology and refinement nulls | Mixed/open | Synthetic calibration only |
| Moonshine/Betti | Pinned arithmetic constants and declared relation graph | D4/isomorphism and destructive nulls | Pass with fragility; Lean target open | Arithmetic/topological diagnostic |
| HBP PB-01 | Synthetic network recovery | Graph and label controls | Not distinct; control invalid | Synthetic calibration |
| HBP PB-02 | PowerGraph IEEE24 | Eight declared controls | Stored `CONTROL_INVALID`; quarantined | No predictive claim |
| HBP PB-03 | PowerGraph cascade sample | Four declared controls | Reproducible non-distinct result; implementation incomplete | No predictive claim |
| HBP PB-04 | Synthetic gene network to PowerGraph | Four declared controls | Reproducible non-distinct result; candidate path non-discriminative | No transfer claim |
| Effective operator | Synthetic ring Laplacian | Identity and unstable-sign controls | Synthetic pass | No EFT-QG or physical field theory |

Negative and open outcomes are present in the generated results and generally
remain visible in reports. CST and Bell have the strongest explicit failure
localization. HBP has explicit negative labels, but its checker and tests do not
currently distinguish deterministic reproduction from experimental validity.

## 5. Commands Executed

Repository and architecture inspection included Git state/history/branches,
tracked and ignored files, directory maps, large-file inventory, frozen API and
manifest reads, import scans, conflict-marker scans, placeholder scans,
absolute-path scans, broad-exception scans, seed scans, and result-status reads.

Primary executable checks:

```bash
uv run python examples/quick_reproduce.py
uv run python examples/quick_reproduce.py --physics-observables-only
uv run python run_phase.py 3 --check
uv run python phase4-sector-freeze/diagnostics/check_phase4_bundle.py
uv run python phase5-readout/diagnostics/check_phase5_bundle.py
uv run python run_phase.py 10 --check
uv run python run_phase.py X --check
uv run python run_phase.py 11 --check
uv run python run_phase.py 12 --check
uv run python run_phase.py 13 --check
uv run python run_phase.py 14 --check
uv run python run_phase.py 15 --check
uv run python run_phase.py 16 --check
uv run python run_phase.py 17 --check
uv run python run_phase.py 18 --check
uv run python run_phase.py 19 --check
uv run python experiments/cst/diagnostics/check_cst_bundle.py
uv run python experiments/hbp/check_hbp_bundle.py
uv run python experiments/biology/external_dataset_bridge/quick_reproduce.py --check
uv run python -m haos_iip.demo stability baseline --json
uv run python experiments/physics_bridge/bell_chsh_reference/run_bell_chsh_reference.py
uv run python experiments/physics_bridge/hydrogen_spectra_reference/run_hydrogen_spectra_reference.py
uv run python experiments/physics_bridge/semiconductor_band_reference/run_semiconductor_band_reference.py
uv run python experiments/effective_operator_theory/run_effective_operator_expansion.py
uv run python -m compileall -q haos_core telemetry experiments examples run_phase.py
```

Unit-test families:

```bash
uv run python -m unittest discover tests -v
uv run python -m unittest discover experiments/cst/tests -v
uv run python -m unittest discover experiments/hbp/tests -v
uv run python -m unittest discover experiments/effective_operator_theory/tests -v
uv run python -m unittest discover experiments/geometry_bridge/tests -v
uv run python -m unittest discover experiments/geometry_bridge/gaussian_prime_norm_lift/tests -v
uv run python -m unittest discover experiments/geometry_bridge/monstrous_moonshine_diagnostic/tests -v
uv run python -m unittest discover experiments/geometry_bridge/synthetic_hidden_geometry_recovery/tests -v
uv run python -m unittest discover experiments/geometry_bridge/synthetic_hidden_probability_rule_recovery/tests -v
uv run python -m unittest discover experiments/geometry_bridge/synthetic_hidden_transformation_recovery/tests -v
uv run python -m unittest discover experiments/geometry_bridge/synthetic_intrinsic_geometry_recovery/tests -v
uv run python -m unittest discover experiments/geometry_bridge/synthetic_observable_prediction/tests -v
uv run python -m unittest discover experiments/geometry_bridge/synthetic_probability_rule_recovery/tests -v
uv run python -m unittest discover experiments/geometry_bridge/synthetic_transformation_semantics_recovery/tests -v
uv run python -m unittest discover experiments/physics_bridge/bell_haos_candidate/tests -v
uv run python -m unittest discover experiments/physics_bridge/synthetic_relational_geometry_calibration/tests -v
uv run python -m unittest discover experiments/physics_bridge/literature_component_bridge/tests -v
```

Dependency and structured-data checks:

```bash
find . -name '*.json' -not -path './.git/*' -not -path './.venv/*' -exec jq empty {} \;
UV_PROJECT_ENVIRONMENT=/tmp/haos-iip-audit-venv uv sync --locked
/tmp/haos-iip-audit-venv/bin/python -c 'from experiments.hbp.pb02_external_network_recovery import pb02_holdout'
```

No repository-wide linter, type checker, coverage gate, Markdown link checker,
or documentation build is configured in `pyproject.toml` or the root execution
guide. Those checks were therefore recorded as missing rather than improvised
as authority.

## 6. Test and Reproducibility Results

- 230 discovered unit tests passed after the safe fixes.
- 41 HBP test executions include ten PB-02 tests twice because
  `test_pb02.py` wildcard-imports `test_pb02_holdout.py`.
- All 1,608 JSON files parsed successfully.
- Python compilation completed without syntax errors for the active inspected
  trees.
- Registered phase/check bundles 3, 4, 5, 10, X, and 11-19 passed.
- Phases 6-9 have builders but no registered checkers. They were not rebuilt.
- Public reproduction passed without changing tracked outputs.
- Biology external-dataset stored-output check passed 41/41 internal checks.

Representative stable hashes:

- CST: `cst_result_7ff2cde2ba1264027997629a` with verdict `FAIL`
- HBP registry: `hbp_result_1fdc653880468c32a182cc95`
- PB-02 stored artifact: `pb02_result_c84f151554dc98485a16bee4`
  with status `CONTROL_INVALID`
- Effective operator: `effective_operator_a740ce933dd1dbd479931847`
- Bell reference: `bell_chsh_result_0639c6444f710c890c8df5fa`
- Hydrogen reference: `hydrogen_spectra_result_6a7fe5cf85a776c6b674315e`
- Semiconductor reference: `semiconductor_band_result_75dbf20a24707f7da1114240`

The Phase 2.2 CST result references a different CST benchmark hash because it
deliberately runs the benchmark with scalar compression disabled. This was
checked in source and is not stale provenance.

## 7. Confirmed Bugs and Defects

| ID | Priority / severity | Category | Affected files | Evidence | Consequence | Confidence | Action | Frozen core affected? |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AUD-001 | P1 / High | Orchestration | `run_phase.py` | `--check` fell back to builder when checker missing | A check command could rewrite frozen artifacts | Confirmed | Fixed fail-closed; tests added | No |
| AUD-002 | P1 / High | Control validity | `pb02_holdout.py:565-590` | Shuffled target is created as tuple element 3, discarded as `_`, while element 2 is original `y` | Shuffled-label control is guaranteed to score 1.0 | Confirmed | Quarantine result; repair only in versioned rerun | No |
| AUD-003 | P1 / High | Artifact consistency | `pb02_holdout.py:483-495`, `687-706` | Baseline feature slices differ between scoring and stored prediction reconstruction | Same result reports baseline Spearman 0.991454 and 0.990555 | Confirmed | Centralize frozen feature registry in PB-02 v2 | No |
| AUD-004 | P1 / High | Holdout leakage | `pb02_holdout.py:618`, `678` | Best baseline is selected by maximum holdout Spearman | Holdout is used for model selection, violating the contract | Confirmed | Select on calibration only in new precommitment | No |
| AUD-005 | P1 / High | Contract enforcement | PB-02 contracts and runner | Contracts remain `DRAFT`/`PRECOMMITMENT_DRAFT`; required output names differ; runner checks only IDs and execution-mode string | Draft protocol executed as frozen benchmark | Confirmed | Add executable contract validator and version gate | No |
| AUD-006 | P1 / High | HBP discrimination | PB-03/PB-04 runners and manifests | PB-03 implements one coarse baseline instead of manifest set; PB-04 baseline and HAOS calls are identical | Comparison cannot establish incremental HAOS value | Confirmed | HBP v2 repair; no result promotion | No |
| AUD-007 | P1 / High | Control validity | PB-03/PB-04 runners | Seed repeat compares arrays to themselves; topology controls replace predictions with mean/sorted predictions rather than perturbing source topology through same path | Declared controls do not test their contracts | Confirmed | Rebuild controls through target execution path | No |
| AUD-008 | P1 / High | Dependency reproducibility | `pyproject.toml`, PB-02/PB-03 loaders | Clean locked environment failed with `ModuleNotFoundError: h5py` | HBP tests/runners depended on machine residue | Confirmed | Fixed dependency and lockfile | No |
| AUD-009 | P2 / Medium | Portability | HBP source and tests | Six Python files hard-coded an author-local absolute volume path | Cold clone fails outside author machine | Confirmed | Fixed with shared path resolver and env override | No |
| AUD-010 | P2 / Medium | Numerical semantics | `benchmark_utils.py`, PB-02/PB-03 loaders | NaN and infinity are converted to observed zero | Missingness is aliased with signal zero and can affect rankings | Confirmed | Preserve mask/counts in versioned HBP repair | No |
| AUD-011 | P2 / Medium | Uncertainty labeling | `pb02_holdout.py:744-755` | `seed_variance` is variance between two model scores, not across seeds | Reported uncertainty is mislabeled | Confirmed | Remove or compute actual independent-seed variance in v2 | No |
| AUD-012 | P2 / Medium | Frozen metric edge cases | `telemetry/frozen_metrics.py` | Empty persistence grid indexes `tau_grid[0]`; self-loop is scored acyclic because singleton SCC loops are not counted | External adapter can fail or overstate acyclicity on edge cases | Confirmed | Do not edit frozen metric; guard and test at adapter boundary | Yes, if changed |
| AUD-013 | P2 / Medium | Checker weakness | PB-02 checker/tests | Checker validates existence and result hash only; tests assert current invalid status and duplicate suite | Deterministic invalid artifacts appear healthy | Confirmed | Add semantic integrity checker in PB-02 v2 | No |
| AUD-014 | P2 / Medium | Documentation drift | root status docs and HBP readiness/status docs | Release 66.5 text coexists with post-67.1 and later experiments; readiness files say runners are missing | Readers cannot infer actual project state reliably | Confirmed | HBP status narrowed now; broader release-state consolidation remains | No |
| AUD-015 | P2 / Medium | Missing validation | phases 6-9 | Builders exist, no independent checkers | Frozen authority cannot be checked through root contract | Confirmed | Add read-only manifest/hash checkers, no rebuilds | No |
| AUD-016 | P2 / Medium | Source-of-truth ambiguity | `API_CONTRACT.md`, `haos_core/`, root guide | Contract names unimplemented future modules while rules freeze implemented `haos_core` APIs | New contributors may implement against wrong surface | Confirmed | Versioned documentation clarification | No |

## 8. Suspected Bugs Requiring Deeper Testing

- Eigenvector sign and degenerate-subspace ambiguity are handled locally in
  several geometry probes but lack one repository-wide invariant policy.
  Sign-stable tests do not necessarily establish eigenspace stability.
- Several reference and toy probes use fixed `1e-12` denominator floors. Most
  are numerically reasonable, but sensitivity is not consistently reported.
- Some generated results contain status under different keys or nested objects.
  There is no common result schema, making automated registry-wide auditing
  fragile.
- The prediction ledger contains field-facing thresholds seeded from one
  materials sidecar. The ledger states this limitation, but independent
  replication and multiple-comparison discipline remain open.
- The large media and artifact corpus may contain stale or duplicate
  presentation files not connected to current source hashes. This was not
  exhaustively regenerated.

## 9. Numerical Risks

- HBP missing values are imputed as zero without a missingness channel.
- PB-02's incremental prediction interval estimates mean raw prediction
  differences, not uncertainty in the Spearman performance margin.
- PB-02 reports only two model scores as `seed_variance`; no independent seed
  variance is measured.
- Frozen telemetry uses an epsilon floor that can map degenerate zero inputs to
  finite values; the behavior is not fully tested for empty, singleton, or
  self-loop cases.
- Constant-vector Spearman returns zero in HBP utilities. This is useful for
  computation but must remain distinguishable from evidence of no association.
- Synthetic geometry tests often use one fixed seed family and modest holdout
  sizes. Passing deterministic controls is calibration evidence, not broad
  robustness.

## 10. Experimental and Control Weaknesses

HBP is the dominant weakness:

- PB-02 control implementation is wrong and its result path is internally
  inconsistent.
- PB-02 chooses the best conventional baseline on holdout.
- PB-03 does not implement the baseline families declared in its manifest.
- PB-03 and PB-04 seed-repeat controls are tautologies.
- PB-03/PB-04 topology controls do not transform topology through the same
  feature and model path as the target.
- PB-04 names identical source/target feature paths as baseline and HAOS.
- PB-02 tests are duplicated and mostly validate stored artifact shape/hash,
  not scientific contract compliance.

Other boundaries:

- CST is underpowered at two independent seeds. Bootstrap uncertainty is not
  independent evidence; the bundle correctly says so.
- CST controls are component-specific; some controls preserve the signal seen
  by other components. Phase 2.2 correctly keeps the instrument partial.
- Bell B3.2.2 terminates at generic relational geometry. No CHSH promotion is
  authorized.
- Geometry recovery remains weak on the transformation/spiral holdout path.
  Strong synthetic subcomponents do not close external semantics.
- Effective-operator recovery detects the terms used to generate its synthetic
  dynamics. This is instrument calibration, not independent law discovery.
- Moonshine/Betti diagnostics depend on a declared relation graph and remain an
  arithmetic/topological sidecar. The local Lean artifact is an open target,
  not a proof.

## 11. Claim-Boundary Review

The strongest repository-level boundary language is sound: the project does
not claim a continuum limit, physical correspondence, spacetime, Bell recovery,
hydrogen derivation, semiconductor derivation, quantum gravity, consciousness,
or ontology.

Supported statements:

- frozen discrete operator and artifact chains reproduce;
- multiple synthetic instruments detect predeclared structures and controls;
- CST currently fails to distinguish its target and is underpowered;
- Bell reference and local-candidate scoreboards work, while HAOS Bell
  derivation remains absent;
- geometry and effective-operator bundles are synthetic calibration layers;
- arithmetic Betti/Moonshine work is diagnostic, not a Moonshine proof.

Unsupported statements:

- HAOS predicts external network recovery better than conventional baselines;
- HBP PB-02/PB-03/PB-04 are valid predictive bridges;
- HAOS derives EFT quantum gravity or a physical field theory;
- synthetic calibration establishes external geometry or physical mechanism;
- Betti recoverability proves Monstrous Moonshine or physical topology.

The root README remains release-state stale: it calls Release 66.5 the active
focus while also linking post-67.1 work and omitting newer active experiment
families. This is a navigation problem, not evidence against frozen results.

## 12. Frozen-Boundary Review

Read and preserved:

- `AGENTS.md`
- `API_CONTRACT.md`
- `ROOT_EXECUTION_GUIDE.md`
- `MIGRATION_MAP.md`
- `THEORY_NOTICE.md`
- `TERMINOLOGY.md`
- `haos_core/core.py`
- `haos_core/io.py`
- `telemetry/frozen_metrics.py`
- `run_phase.py`
- `phaseIII_phaseIV_authoritative_manifest.json`
- phase manifests and README/checker contracts through Phase XIX
- prediction ledger and experiment precommitments inspected for active bundles

No frozen API, metric definition, threshold, phase manifest, canonical result,
or historical experiment output was modified. Frozen telemetry edge-case bugs
were documented rather than patched.

## 13. Safe Fixes Implemented

1. Phase check fail-closed guard

   - Reproduced from dispatch logic: missing checkers fell back to builders.
   - Registered existing Phase 3-5 checkers.
   - Added `resolve_target()` and refused builder fallback in check mode.
   - Added three regression tests.

2. Clean dependency declaration

   - Reproduced in an isolated locked environment: PB-02 import failed because
     `h5py` was absent.
   - Added `h5py>=3.10` and regenerated `uv.lock`.
   - Clean environment now imports PB-02 with `h5py==3.16.0`.

3. Portable HBP data roots

   - Replaced hard-coded source paths in active HBP Python files.
   - Default is the workspace sibling `DATA/Powergraph`.
   - `HAOS_POWERGRAPH_ROOT` provides an explicit override.
   - Added default and override tests.

4. Claim and path documentation narrowing

   - Replaced machine-specific links in current API/milestone documents.
   - Marked PB-02 stored output as quarantined.
   - Added audit notes to PB-03/PB-04 pre-run readiness documents.
   - Corrected the HBP status snapshot so deterministic invalid artifacts are
     not described as unexecuted or operationally valid.

## 14. Files Changed

- `API_CONTRACT.md`
- `PROGRAM_STATE_MILESTONE_19.md`
- `ROOT_EXECUTION_GUIDE.md`
- `run_phase.py`
- `tests/test_run_phase.py`
- `pyproject.toml`
- `uv.lock`
- `experiments/hbp/data_paths.py`
- `experiments/hbp/hbp_status_snapshot.md`
- `experiments/hbp/pb02_external_network_recovery/README.md`
- `experiments/hbp/pb02_external_network_recovery/pb02_holdout.py`
- `experiments/hbp/pb03_temporal_recovery/execution_readiness.md`
- `experiments/hbp/pb03_temporal_recovery/loader.py`
- `experiments/hbp/pb04_cross_domain_transfer/execution_readiness.md`
- `experiments/hbp/pb04_cross_domain_transfer/loader.py`
- `experiments/hbp/pb04_cross_domain_transfer/run_pb04_cross_domain_transfer.py`
- `experiments/hbp/tests/test_data_paths.py`
- `experiments/hbp/tests/test_pb02_holdout.py`
- `experiments/hbp/tests/test_pb03_pb04_freeze.py`
- `docs/audits/HAOS_IIP_PROJECT_AUDIT_2026-07-11.md`

## 15. Remaining Prioritized Work

P1:

1. Freeze current PB-02/PB-03/PB-04 outputs as invalid historical artifacts;
   do not overwrite them.
2. Write a new HBP integrity precommitment before changing controls or model
   selection.
3. Centralize PB-02 feature families so scoring and prediction use the exact
   same arrays.
4. Select baselines on calibration only and open holdout once.
5. Implement control transformations before feature extraction and route them
   through the target model pipeline.
6. Make result checkers validate contracts, metrics, control behavior, and
   source hashes, not only result hashes.

P2:

1. Preserve missingness rather than imputing non-finite values as zeros.
2. Add read-only checkers for Phases 6-9.
3. Clarify current `haos_core` authority versus future `API_CONTRACT` module
   layout.
4. Consolidate release-state documentation into one generated status index.
5. Add schema validation for result JSON and precommitment contracts.
6. Add adapter-level tests for empty persistence grids and self-loop causal
   graphs without editing frozen telemetry.
7. Remove duplicate PB-02 test discovery while preserving a canonical command.

P3:

1. Establish an artifact retention policy for large trial tables, duplicate
   JSON blobs, PDFs, and video.
2. Configure a minimal lint/type/link-check layer for new code only.
3. Separate generated artifacts from source in a machine-readable manifest.

P4:

- Continue geometry, effective-operator, Bell, and physical bridge research
  only under new precommitments after the HBP instrument path is reliable.

## 16. Suggested Next Development Phase

Name: `HBP-IR-01 - Predictive Instrument Integrity Repair`

Objective:

Build one versioned, contract-enforced network-recovery benchmark that can
return `CONTROL_INVALID`, `PREDICTION_NOT_DISTINCT_FROM_BASELINES`, or
`PREDICTION_OUTPERFORMS_BASELINES` without target-dependent repair.

Subsystems:

- `experiments/hbp/benchmark_utils.py`
- `experiments/hbp/pb02_external_network_recovery/`
- `experiments/hbp/hbp_validation.py`
- `experiments/hbp/tests/`

Precommit before implementation:

- exact feature registry and dimensional schema;
- missingness policy;
- development/calibration/holdout selection rules;
- baseline selection on calibration only;
- one target/control execution interface;
- source-transform contracts for every control;
- actual independent repeat/seed policy;
- uncertainty on performance differences;
- result schema and checker gates;
- unchanged scientific thresholds.

Validation criteria:

- same feature registry drives fitting, scoring, and stored predictions;
- no holdout data influence model or feature selection;
- each control destroys its declared structure and preserves declared
  confounds;
- checker detects the historical PB-02 baseline mismatch and shuffled-label
  bug;
- missingness is reported separately from observed zero;
- two clean-environment runs produce the same hash;
- negative results remain valid outcomes.

Stopping condition:

Stop and return `CONTROL_INVALID` if any declared control fails, if holdout is
accessed before freeze, if target labels enter feature construction, or if the
same model is presented as both baseline and HAOS candidate. Do not tune around
the failure.

## 17. Explicit Non-Claims and Unresolved Questions

This audit does not establish:

- physical correspondence;
- a continuum limit;
- literal geometry, curvature, spacetime, or gravity;
- quantum mechanics or Bell correlations from HAOS-IIP;
- hydrogen or semiconductor derivation;
- EFT quantum gravity;
- a physical mechanism from HBP;
- Monstrous Moonshine from the Betti diagnostic;
- empirical validation of universal recoverability.

Unresolved questions:

- Can HBP add predictive value once controls and holdout discipline are valid?
- Are frozen telemetry edge cases relevant to any current external adapter
  input, or only latent risks?
- Can geometry transformation semantics survive independent topology and
  refinement holdouts?
- Can effective terms be extracted from HAOS artifacts rather than from a
  synthetic generator containing those terms?
- Which prediction-ledger rows have independent replication data rather than
  sidecar-seeded thresholds?

Self-audit:

- Invented repository assumptions: none used as findings. Remote state is based
  on local configured refs because no fetch was performed.
- Analogy-only mappings: geometry, field, EFT, Moonshine, and physics-facing
  language remain explicitly bounded and cannot support official prediction.
- Calibrated parameters masquerading as derivation: present as a general risk
  in synthetic geometry/effective-operator bundles; their documents currently
  maintain the calibration ceiling.
- Target-leakage paths: PB-02 holdout baseline selection is confirmed; field-
  facing prediction thresholds seeded from one sidecar remain a replication
  risk.
- Baselines HAOS may fail to beat: PowerGraph closeness and conventional
  spectral/topology models already dominate current HBP reads.
- Synthetic-calibration inflation risk: known-positive generators can validate
  instruments but cannot establish external mechanism ownership.
- Claims still requiring external evidence: every empirical physical bridge,
  continuum statement, and mechanism-level physics interpretation.

Final audit result: `PARTIALLY COHERENT / REPRODUCIBLE / HBP REPAIR REQUIRED`.
