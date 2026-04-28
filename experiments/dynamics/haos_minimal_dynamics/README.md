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

`address_restoration` now defaults to harmonic spectral projection: the field is
pulled toward the frozen low-frequency eigenmode coefficients of the branch
operator. Local neighbor-difference restoration remains available as
`--address-mode local`; `--address-mode hybrid` blends spectral and local pulls.
The default address gain is `0.45`, promoted from the fixed address-gain sweep.

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

`--null-level 2` adds a spectral-aware null on top of that gate:
candidate shuffles preserve degree/shell buckets and are selected to match the
final field's low-mode spectral energy before branch identity is scored.
By default, `--null-level 3` also matches one-hop local autocorrelation.
`--null-level 4` adds higher-order matching with 2-hop and triangle-supported
correlation proxies.

The run also writes a targeted ablation report. Ablations reuse the same initial
noise and perturbation noise, then remove one dynamics component at a time:

- no address restoration
- no shell-invariant restoration
- diffusion only
- randomized branch targets
- address only
- old local address restoration
- spectral address projection
- hybrid spectral/local address projection
- scalar phase address pull
- multi-scale shell-weighted address pull
- fixed address-gain sweep

## Run

```bash
uv run --with numpy --with matplotlib python experiments/dynamics/haos_minimal_dynamics/run_minimal_dynamics.py
```

Compare old local addressing:

```bash
uv run --with numpy --with matplotlib python experiments/dynamics/haos_minimal_dynamics/run_minimal_dynamics.py --address-mode local
```

To run only the address-focused ablation family:

```bash
uv run --with numpy --with matplotlib python experiments/dynamics/haos_minimal_dynamics/run_minimal_dynamics.py --focus-address
```

To compare null strictness:

```bash
uv run --with numpy --with matplotlib python experiments/dynamics/haos_minimal_dynamics/run_minimal_dynamics.py --null-level 1
uv run --with numpy --with matplotlib python experiments/dynamics/haos_minimal_dynamics/run_minimal_dynamics.py --null-level 4
```

## Outputs

Generated outputs are ignored:

- `bridge_status.json`
- `ablation_report.md`
- `minimal_dynamics_report.md`
- `minimal_timeseries.csv`
- `ablation_branch_identity.png`
- `null_erosion_spectral_vs_local.png`
- `recoverability.png`
- `address_specificity.png`
- `invariant_retention.png`

## Status Semantics

- `PASS`: observed recovers, strict branch-identity specificity passes, and controls do not.
- `MARGINAL`: observed recovers or is specific, but controls are too close.
- `FAIL`: recoverability or specificity does not survive.

This is a sidecar diagnostic, not a physics claim.
