# CST Toy Benchmark

Lifecycle status: `TERMINAL_NEGATIVE` for CST v0.2.2. The repaired instrument,
frozen results, and tests remain reproducible, but target-adjacent expansion is
not authorized. See the
[branch lifecycle registry](../../docs/branch_governance/branch_lifecycle_summary.md).

CST in this directory means **Closure Stability Test**. It is an experimental,
auditable numerical extension around existing HAOS-IIP frozen ledgers. It does
not assert new ontology and does not claim to replace quantum mechanics,
solve measurement, recover Born weights, or reproduce external physics
benchmarks.

## Scope

Implemented fact: the benchmark consumes frozen Phase 15-18 CSV artifacts and
writes deterministic CST telemetry under `experiments/cst/runs/`.

Design choice: a candidate branch is represented by `ClosureSignature`, with a
separate `CSTRunProvenance`. The branch identifier is a stable hash of a compact
`address_summary`; the run identifier hashes seed, config, source artifacts, and
control group.

Heuristic: non-front event times are rank-interpolated from Phase 18 recovery
time fields because the matching `n60/72/84` frozen slice does not expose all
Phase 16 event times directly.

Analogy: `generalized CST selectivity` is an engineering signal/noise ratio. It
is not a physical Q-factor or resonance claim.

Unverified hypothesis: repeated recovery of the same compact address under
matched perturbations may be a useful closure-stability diagnostic.

## Command

```bash
uv run python experiments/cst/run_cst_benchmark.py
```

Disable optional scalar compression:

```bash
uv run python experiments/cst/run_cst_benchmark.py --no-scalar
```

Check the generated CST bundle:

```bash
uv run python experiments/cst/diagnostics/check_cst_bundle.py
```

## Outputs

- `cst_runs.json`: observations, provenance, perturbations, and signatures.
- `closure_signatures.json`: branch signatures per observation.
- `closure_distance_matrix.csv`: scalar pairwise closure distances.
- `closure_distance_components.csv`: component-level distances.
- `recoverability_vectors.csv`: vector metrics; scalar is optional.
- `control_distributions.csv`: matched control distributions and missing controls.
- `seed_manifest.json`: deterministic seeds, config hash, artifact hashes.
- `benchmark_result.json`: PASS/OPEN/FAIL verdict with reasons.
- `benchmark_report.md`: concise technical report.

## Controls

Available controls use the same CST observation and distance path:

- shuffled-structure control
- randomized-edge control
- degraded-signal control
- label-permutation control
- generic graph/operator control
- seed-repeat control
- parameter-matched null control
- frozen periodic diagonal augmented control

Missing: perturbation-free control. It is reported as unavailable because the
selected frozen slice does not contain an unperturbed repeat with the same
signature contract.

## Verdict Semantics

`PASS` means all configured toy gates passed. `OPEN` means the benchmark produced
positive evidence but failed a scope or availability requirement. `FAIL` means
one or more configured gates failed. All three are valid outcomes.

The scalar distance uses explicit weights from
`configs/cst_toy_config.json`. Component-level values are authoritative for
review; scalar compression can be disabled.
