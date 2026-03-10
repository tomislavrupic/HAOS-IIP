# Stage 5 Prompt — Prototype: Dirac-Type Operator Lift (DH)

You are working inside the HAOS-IIP repository.

All existing operators and numerical infrastructure are frozen.

Do not modify:
- kernel definition
- graph construction
- node Laplacian `L_0`
- edge/Hodge operator `L_1`
- restricted transverse sector tests
- gauge phase implementation

Stage 5 introduces one new operator lift built on the same substrate.

The goal is not to claim fermions or gauge dynamics.

The goal is only to test whether a Dirac-type lift exists that is operator-consistent with the scalar branch.

## Core hypothesis

Construct an operator `D_H` such that

`D_H^\dagger D_H \approx L_0`

up to expected discretization error.

If this relation fails, the lift is rejected.

## Hilbert Space

Spinor field on nodes.

Two prototype choices:

Minimal prototype:

`\psi_i \in \mathbb{C}^2`

Full prototype:

`\psi_i \in \mathbb{C}^4`

Total operator size:

`spinor_dimension x N_nodes`

## Operator Construction

Reuse existing objects:
- Gaussian kernel weights `w_{ij}`
- node embedding `x_i`
- link phases `U_{ij}`

Define link direction:

`n_{ij} = (x_j - x_i) / |x_j - x_i|`

Define gamma matrices:

Standard Hermitian set satisfying

`{\gamma^\mu, \gamma^\nu} = 2\delta^{\mu\nu}`

Prototype operator:

`(D_H \psi)_i = \sum_j w_{ij} (\gamma \cdot n_{ij}) U_{ij} \psi_j - \sum_j w_{ij} (\gamma \cdot n_{ij})^\dagger U_{ji}^\dagger \psi_i`

Multiply by `i` if necessary to enforce Hermiticity.

## Experiment Block 1 — Hermiticity Test

Create:

`numerics/simulations/DH_hermiticity_test.py`

Verify:

`||D_H - D_H^\dagger||`

Report numerical deviation.

Expected result:

machine-precision Hermiticity.

## Experiment Block 2 — Square-Root Test

Create:

`numerics/simulations/DH_square_root_test.py`

Compute:

`D_H^\dagger D_H`

Compare spectrum against scalar operator:

`L0`

Diagnostics:
- eigenvalue correlation
- spectral deviation
- convergence with grid size

Plot:

`DH2_vs_L0_spectrum.png`

Interpretation:

Does `D_H` behave as a square root of the scalar operator?

## Experiment Block 3 — Spectral Structure

Create:

`numerics/simulations/DH_spectrum_structure.py`

Compute lowest 50 eigenvalues.

Check:
- `+/-` spectral pairing
- degeneracy structure
- localization (`IPR`)

Plots:

`DH_spectrum.png`
`DH_ipr_vs_mode.png`

Interpretation:

Does the spectrum resemble a Dirac-like operator?

## Experiment Block 4 — Flux Response (Optional)

Reuse existing torus + Landau gauge infrastructure.

Insert flux quanta through torus cycles.

Test:

how `D_H` spectrum responds to background flux.

Plots:

`DH_flux_spectral_flow.png`

Interpretation:

Does the spectrum shift smoothly with flux?

## Minimal Prototype Grid

Start with a small system:

`6 x 6` periodic torus (2D test)

Then:

`8^3` periodic lattice (3D prototype)

Do not scale larger until tests pass.

## Required Diagnostics

Every `D_H` experiment must output:
- eigenvalues
- `IPR`
- divergence-curl diagnostics if applicable
- spectral comparison with `L0`

## Output Format

Produce:

`data/YYYYMMDD_DH_test.json`
`plots/YYYYMMDD_plotname.png`
`experiments/spinor_sector/Experiment_Name_v1.md`

Update:

`EXPERIMENT_LOG.md`

Do not generate a paper.

## Interpretation Boundary

Stage 5 may support statements such as:
- a Dirac-type lift exists on the same kernel substrate
- `D_H^2` approximates the scalar operator
- the `D_H` spectrum shows expected symmetry structure

Stage 5 must not claim:
- fermions
- quantum field theory
- gauge dynamics
- particle interpretation

## Deliverables

1. `D_H` operator prototype
2. Hermiticity validation
3. square-root test vs `L0`
4. spectral structure analysis
5. optional flux response test

## Why Stage 5 matters

Stages 1–4 explore scalar and vector sectors.

Stage 5 asks a structural question:

Does the same kernel substrate admit a consistent spinor lift?

That completes the operator hierarchy:

`scalar -> vector -> spinor`
