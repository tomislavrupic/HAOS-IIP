# EL-R3-MESOSCOPIC-FIELD-01

This experiment tests whether a uniformly quantized, fixed-budget `4 x 4`
mesoscopic field can restore the frozen low-frequency function across
microscopic refinement.

The payload is fixed at sixteen 8-bit values (`128` bits) at `8 x 8`,
`12 x 12`, and `16 x 16`. The decoder changes only each physical cell's mean
and leaves all within-cell residual structure untouched.

Authority and boundaries:

- theory: `docs/roadmaps/EL_R3_THEORY_REVISION_01.md`;
- contract: `precommitment_contract.json`;
- parent RT-01, RT-02, and RP-01 results remain frozen;
- no continuum, universality, physics, or general emergence claim;
- no commit, push, tag, or release is part of the experiment.

The contract was frozen before implementation or result inspection.

## Frozen execution

The representation stores sixteen signed, uniformly quantized 8-bit coarse-cell
means in canonical row-major order. The decoder applies one uniform correction
inside each fixed physical cell, restoring the stored quantized field exactly.
It never reads the independent functional score or its four phasors.

Run implementation checks and the smoke partition without consuming final seeds:

```bash
uv run python -m unittest tests.test_rung3_mesoscopic_field
uv run python experiments/emergence_ladder/rung3_mesoscopic_field/run_experiment.py \
  --partition smoke --output-root /tmp/el-r3-mesoscopic-smoke
```

The final partition is single-use at the selected output root:

```bash
uv run python experiments/emergence_ladder/rung3_mesoscopic_field/run_experiment.py \
  --partition final
uv run python experiments/emergence_ladder/rung3_mesoscopic_field/check_bundle.py
```

The checker recomputes aggregate decisions from raw rows, verifies every
predeclared schedule cell, payload and parent/source/artifact hash, and rejects
oracle admission or claim expansion.

## Frozen outcome

The final classification is `PARTIAL_RECOVERY_ONLY`. The field passed exact
quantized restoration, scale stability, uncertainty, variance, and destructive
control gates. It recovered the coarse-cell-offset and smooth-physical-twist
families at every tested level, but recovered none of the within-cell balanced
gradient family, leaving the per-level recovery rate at `2/3` below the frozen
`0.75` gate. Rung 3 therefore remains unsupported and this candidate is closed
without retuning. See `final/summary.md` and `final/result_manifest.json`.
