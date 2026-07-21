# EL-R3-RP-01 Classification

Classification: `VALIDATION_GATE_FAILED`

The validation gate failed 6 of 10 mandatory conditions. Final evaluation was
not authorized and no final seed was consumed.

## Gate Results

| Gate | Result |
| --- | --- |
| Nonzero functional recovery somewhere | PASS |
| Target exceeds passive with CI above zero | FAIL |
| Target exceeds frozen RT-02 with CI above zero | FAIL |
| Target exceeds random parity and random memory | FAIL |
| Identity threshold | PASS |
| Trivial-attractor control rejected by full recovery | FAIL |
| No leakage | PASS |
| Parity deletion degrades full recovery | FAIL |
| Above-radius recovery is lower than within-radius recovery | FAIL |
| Memory cannot reconstruct continuous state | PASS |

The last two recovery-rate gates fail because both within- and above-radius
full recovery rates are zero. This must not be compressed into an "almost"
result.

Rung 3 remains `NEGATIVE_RESULT`. Rung 4 stays dependency-blocked.
