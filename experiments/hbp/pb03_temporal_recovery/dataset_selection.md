# PB-03 Dataset Selection

PB-03 is frozen as a contract-first branch. This note freezes the current
dataset candidate direction so the next implementation step can be written
against a real source family instead of a placeholder.

## Selected Candidate Family

- PowerGraph cascades under `/DATA/Powergraph/dataset_cascades`
- infrastructure evolution style temporal recovery

## Candidate Networks

- `ieee24`
- `ieee39`
- `ieee118`
- `uk`

## Why this candidate

- It is temporal.
- It has explicit disturbance / cascade structure.
- It is closer to the PB-03 question than static topology-only graphs.
- It avoids reopening the already-frozen PB-02 holdout path.

## Limits

- This is still a candidate selection note, not a runner.
- The dataset still needs a frozen split manifest and benchmark implementation.
- The branch remains negative-until-proven.

## Freeze Statement

PB-03 will use PowerGraph cascades as the first implementation candidate unless a
later explicit precommitment supersedes this note.

