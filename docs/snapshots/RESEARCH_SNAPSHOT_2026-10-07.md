# Reviewed local research snapshot — 7 October 2026

Preservation commit: `2d0731e41458566b7668946fa7709897260f6384`.

Branch: `codex/research-snapshot-2026-10-07`.

Parent/base: `f4742c49f4a50b18f7b92b98b18f176cd1637987`.

The preservation commit captures the existing local research before status
reconciliation. It contains 63 formerly untracked files and seven existing
modified files, totaling 11,882,431 bytes across the 70 captured files. Every
captured file was compared with its committed blob using SHA-256. The exact
[inventory](research_snapshot_2026-10-07.json) records paths, sizes, hashes, and
whether each file was previously modified or untracked. Those hashes refer to
the preservation commit, not later edits to navigation documents.

## Included work

- Mesoscopic-field theory, implementation, tests, contract, 3,024 stored rows,
  summaries, plots, infographics, and pre-repair audit history.
- Adaptive-recoverability implementation, tests, contract, 3,888 stored rows,
  summaries, and hostile audit.
- PIX-7 finite-capacity exchange source, stored diagnostic result, trajectory
  arrays, comparison figure, construction, constraints, and candidate audit.
- The earlier blocked one-bit scale falsifier, including its non-frozen draft
  and M0 discovery record. It remains unexecuted historical evidence.
- Existing claim-governance schema, contract, validator, tests, foundation
  notes, and lifecycle updates needed to interpret these experiments.

All existing changes belonged to this research/governance set. Virtual
environments, Python caches, ignored reproduction outputs, and other repositories
were excluded. No frozen result, source implementation, contract, threshold,
or numerical row was edited during preservation. Original CSV line endings and
existing final blank lines were retained for byte-level fidelity.

## Review and verification

The review checked the change inventory and scope, result classifications,
source/artifact hashes, parseability of JSON and JSONL files, and common private
key/token patterns. This is not a comprehensive security or mathematical audit.

Executed successfully on 7 October 2026:

```bash
.venv/bin/python scripts/check_branch_lifecycle.py
.venv/bin/python scripts/check_claim_contracts.py
.venv/bin/python experiments/emergence_ladder/rung3_mesoscopic_field/check_bundle.py
.venv/bin/python experiments/antifragility/adaptive_recoverability_v1/check_bundle.py
.venv/bin/python -m unittest tests.test_branch_lifecycle tests.test_claim_contracts tests.test_rung3_mesoscopic_field experiments.antifragility.adaptive_recoverability_v1.tests.test_adaptive_recoverability experiments.hbp.integrity_repair_v2.tests.test_hbp_ir_01
.venv/bin/python examples/quick_reproduce.py
```

Results: 35 targeted tests passed; lifecycle and claim checks passed; both
experiment bundles passed their checks; the public artifact reproduction passed.
The mesoscopic classification remains `PARTIAL_RECOVERY_ONLY`, and the adaptive
classification remains `BOUNDED_TOY_ANTIFRAGILITY_SUPPORTED`.

For PIX-7, a read-only check matched `source_sha256` to `run_probe.py`, confirmed
the stored `PASS_SCOPED` flags, and loaded the trajectory with
`numpy.load(..., allow_pickle=False)`. All arrays were finite and agreed on
601 time samples: compact shape `(601, 2, 36, 3)`, tangent `(601, 4, 36)`.
The simulation was not rerun. This establishes stored-bundle consistency only.

No final seed partition was rerun, no independent replication was performed,
and no scientific claim was promoted. The later documentation reconciliation
does not change the lifecycle registry or its generated summary.

## Inspect and recover

Inspect the exact preservation commit:

```bash
git show --stat 2d0731e41458566b7668946fa7709897260f6384
git show 2d0731e41458566b7668946fa7709897260f6384:experiments/emergence_ladder/rung3_mesoscopic_field/final/summary.md
```

To inspect the original snapshot in another checkout without resetting current
work, choose an unused directory and run:

```bash
git worktree add --detach ../HAOS-IIP-snapshot-inspection 2d0731e41458566b7668946fa7709897260f6384
```

The snapshot is local. Nothing was pushed, tagged, or published, and `main`
remained at the base commit during this work. A local commit is recoverable Git
history on this drive; it is not an off-device backup.
