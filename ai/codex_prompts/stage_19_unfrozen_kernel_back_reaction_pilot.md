# Stage 19 — Unfrozen Kernel Back-Reaction Pilot

## Objective
Test whether minimal local substrate back-reaction enables composite persistence or geometry-like capture on top of the frozen HAOS-IIP architecture.

This is the opening move of Phase II.
It is not a rescue patch for the frozen program.
It is a new architectural regime in which the interaction kernel itself becomes dynamical.

## Architectural change
The frozen interaction kernel `K0_ij` becomes a slow dynamical edge-weight field.

Canonical Stage 19A choice:

`K_ij(t) = K0_ij * (1 + eps * S_ij(t))`

with a slow stress-memory variable

`S_ij(t + dt) = (1 - dt / tau_relax) * S_ij(t) + rho_ij(t) - rho_bar_ij(t)`

where:
- `rho_ij(t)` is the canonical local occupation or envelope-energy source
- `rho_bar_ij(t)` is a zero-centered normalization term or capped local baseline used to prevent secular kernel inflation

## Canonical local source
Use amplitude or energy first, not phase.

Canonical source:

`rho_ij(t) = |A_ij(t)|^2`

where `A_ij(t)` is the local packet-edge amplitude or local envelope-energy proxy already available in the codebase.

Control source, only if Stage 19A shows real structure:

`rho_ij(t) = |phi_i(t) - phi_j(t)|^2`

Phase is not used in the first pilot.

## Hard constraints
- `support(K) = support(K0)`
- `K_ij = K_ji`
- bounded deviation: `|K_ij - K0_ij| <= kappa_max`
- no graph-topology changes
- no thresholding
- no additional nonlinear branches
- no sigma modulation
- no combined amplitude plus phase source

Higher local occupation is taken to mean stronger local effective coupling.
The opposite sign is deferred to a later sub-stage if needed.

## Timescale discipline
- `eps << 1`
- `tau_relax` must exceed packet crossing time
- packet propagation remains the fast sector
- substrate evolution remains the slow sector

Use only two relaxation choices in the first pilot:
- one near crossing time
- one clearly slower than crossing time

No broader sweep is allowed.

## Representatives
Keep the same three Phase I representatives:
1. clustered composite anchor
2. phase-ordered symmetric triad
3. counter-propagating corridor

## Conditions
For each representative run:
- null control: `eps = 0`
- very weak back-reaction
- weak back-reaction

Total base runs: `9`

Only if the base pilot earns promotion may the slower `tau_relax` variant or `12 -> 24` refinement be run.

## Diagnostics
Reuse the frozen ladder:
- composite lifetime
- binding persistence
- coarse basin persistence
- morphology classification label
- coarse classification label

Add Stage 19A-specific diagnostics:
- kernel deviation norm
- max local kernel deformation
- substrate-relaxation memory score
- local stress persistence

## Early-stop rule
Terminate Stage 19A if all base runs show:
- no composite lifetime gain
- no basin persistence gain
- no new stabilized ordering or morphology class

Correlation shifts or small kernel deformations alone do not count.

## Positive signal definition
Only count as structural emergence if at least one representative shows:
- statistically stable composite persistence increase
or
- clear basin stabilization
or
- a new morphology class that survives refinement

If positive appears:
- map only the minimal parameter window
- immediately run the scale-robustness check before interpretation

## Interpretation boundary
If Stage 19A remains null:
- conclude the binding deficit exceeds frozen-substrate memorylessness alone
- motivate deeper architectural hypotheses only after this bounded negative is frozen

If Stage 19A is positive:
- describe it only as substrate back-reaction enabling collective persistence
- do not claim spacetime or force law emergence

## Compression
Frozen pre-geometry can organize, but it does not bind.

Stage 19 tests whether bounded substrate back-reaction is the minimal missing ingredient.
