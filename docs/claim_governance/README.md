# Claim Governance

This directory is the machine-readable authority for public HAOS-IIP claim
boundaries. It supplements the narrative claim-gating and falsification
documents; it does not validate a scientific claim or alter a frozen result.

## Required contract

Every new public claim must have a `*.claim.json` file conforming to
`claim_contract.schema.json`. The contract makes the dependencies hidden by a
short claim visible:

- selection history and disclosure completeness;
- the declared perturbation domain;
- primary telemetry and reasonable alternatives;
- controls used, statistics matched, and dangerous controls still omitted;
- artifact, execution, and independent-replication status;
- allowed language, non-claims, downgrade conditions, and kill conditions.

Validate all committed contracts with:

```bash
uv run python scripts/check_claim_contracts.py
```

The validator checks required fields, controlled vocabulary, referenced local
artifacts, and cross-field claim ceilings. Passing it means only that the claim
boundary is complete and internally consistent.

## Evidence levels

| Level | Meaning | What it does not establish |
| --- | --- | --- |
| `ARTIFACT_VERIFIED` | Committed outputs match committed expected artifacts | Fresh numerical recomputation |
| `EXECUTION_REPRODUCED` | The declared computation was rerun successfully | Independent interpretation or external truth |
| `INDEPENDENTLY_REPLICATED` | An unaffiliated party reproduced outputs and classifications from the declared protocol | Physical correspondence unless separately tested |

The levels are cumulative only when the contract provides the evidence paths.
An absent higher level must be recorded as `NOT_DEMONSTRATED`, not inferred from
repository organization or document volume.

## Current contract

- `scale_bridge_66_5.claim.json` records the current scale-bridge claim as
  `OPEN`. It preserves the existing bounded numerical evidence while exposing
  incomplete selection-history disclosure, untested metric alternatives,
  omitted adversarial controls, absent independent replication, and absent
  physical correspondence.

