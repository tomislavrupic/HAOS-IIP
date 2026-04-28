# HAOS Minimal Dynamics

This is the stripped-back HAOS dynamics falsifier.

It tests one claim only:

> Does a relational address-preserving scalar update recover a perturbed field
> on the HAOS geometry graph better than controls?

No reaction-diffusion. No visual pattern search. No parameter garden.

## Rule

```text
x[t+1] = x[t]
       - diffusion * L x[t]
       + address_gain * address_restoration
       + invariant_gain * shell_variance_restoration
```

`address_restoration` tries to recover each node's weighted neighbor-difference
signature from the frozen reference field.

`shell_variance_restoration` tries to preserve each node's frozen branch-local
contrast energy: the weighted variance of its neighbor offsets. This adds one
explicit invariant-preserving pull without reintroducing reaction-diffusion
parameter complexity.

Controls:

- node label shuffle
- edge weight shuffle
- degree-preserving rewiring
- topology randomization

The strict specificity gate also includes a degree/shell-stat preserving null:
final field values are shuffled only inside matched weighted-degree and frozen
shell-variance buckets. This breaks branch identity while preserving the cheap
local statistics that made earlier controls too competitive.

## Run

```bash
uv run --with numpy --with matplotlib python experiments/dynamics/haos_minimal_dynamics/run_minimal_dynamics.py
```

## Outputs

Generated outputs are ignored:

- `bridge_status.json`
- `minimal_dynamics_report.md`
- `minimal_timeseries.csv`
- `recoverability.png`
- `address_specificity.png`
- `invariant_retention.png`

## Status Semantics

- `PASS`: observed recovers, strict branch-identity specificity passes, and controls do not.
- `MARGINAL`: observed recovers or is specific, but controls are too close.
- `FAIL`: recoverability or specificity does not survive.

This is a sidecar diagnostic, not a physics claim.
