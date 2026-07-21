# EL-R3-RP-01 - Distributed Parity Recovery

Status: `VALIDATION_GATE_FAILED`  
Final evaluation: not authorized; zero final seeds consumed.

RP-01 tests whether 32 distributed parity bits over 48 local order relations
can restore identity, dynamics, and an independently measured low-frequency
response after perturbation of a 64-value continuous state.

The implementation is deliberately bounded:

- each decoder sees one three-symbol block and its two parity bits;
- the memory contains no continuous values, future states, function labels, or
  recovery scores;
- the correction bridge changes the continuous state through local inequality
  feedback rather than state reset;
- RT-02 is invoked unchanged as a frozen baseline;
- the full-state oracle and centralized decoder are non-admissible controls.

## Frozen Result

Within the declared one-error-per-block radius, the decoder recovered the
correct relational symbols in all 16 validation target cases and restored the
relation/operator region. It did not restore the independent function:
median functional restoration was `0.0`, and median trajectory-corridor gain
was about `0.0265`.

Above the radius, local decoding commonly eliminated the syndrome while
converging to the wrong codeword. Across all primary target validation runs,
the full recovery rate was `0.0`. The target did not beat passive relaxation,
the frozen RT-02 mechanism, random parity, or memory-budget-matched random
bits. Parity deletion therefore could not demonstrate degradation of a
functional recovery signal that was absent.

This is not a coding-theory failure. The local code exhibits its declared
distance-three behavior. It is a negative result for the stronger emergence
claim: this relational alphabet plus parity architecture does not contain
enough function-relevant information to support Rung 3.

## Run

```bash
UV_CACHE_DIR=/tmp/haos-uv-cache uv run python experiments/emergence_ladder/rung3_distributed_parity/run_experiment.py
UV_CACHE_DIR=/tmp/haos-uv-cache uv run python experiments/emergence_ladder/rung3_distributed_parity/check_bundle.py
```

Do not rerun under `EL-R3-RP-01` with changed symbols, checks, thresholds,
decoder, or bridge. A new mechanism requires a new identifier and contract.

