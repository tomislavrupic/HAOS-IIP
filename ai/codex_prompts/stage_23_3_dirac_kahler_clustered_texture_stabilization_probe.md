# Stage 23.3 - Dirac-Kahler Clustered Texture Stabilization Probe

## Context

Previous Phase III runs established:
- Dirac-Kahler propagation introduces structured grade exchange during collisions.
- A single refinement-stable collision-texture family was identified:
  - `tight_clustered_pair :: in_phase`
- This family shows weak persistence gain, but does not yet produce bounded persistence or composite capture.

Stage 23.3 is therefore not a discovery scan.
It is a stabilization-mechanism probe targeted specifically at this collision texture.

The goal is to test whether persistence failure is due to:
- insufficient internal grade coupling structure
- absence of self-induced effective mass / curvature
- lack of phase-locked multi-component coherence
- or simple resolution-dependent dispersion

## Objective

Determine whether the refinement-stable clustered DK encounter can transition from:

`collision texture -> metastable composite -> bounded persistence regime`

under minimal structural extensions that remain within the Dirac-Kahler operator sector.

## Architecture (frozen baseline)
- Domain: 2D periodic complex (same as Stage 23.2)
- Operator:
  - `D_DK = d + delta`
- Propagator: identical numerical scheme as Stage 23.2
- Kernel: static Gaussian interaction kernel
- Packet initialization: reuse the exact tight clustered in-phase configuration

No substrate reinforcement or memory laws are introduced in this stage.

## New controlled probes

Each probe activates one additional structural degree of freedom only.

### Probe A - Internal grade-coupling strength
Introduce controlled mixing term:

`partial_t Psi = D_DK Psi + alpha * M(Psi)`

where `M` redistributes amplitude between 0- and 1-forms locally.

Test levels:
- `alpha = 0`
- `alpha = 0.01`
- `alpha = 0.02`

Goal:
Test whether stronger internal multicomponent coherence increases persistence.

### Probe B - Phase-locked composite constraint
Add weak phase-alignment term acting only when packet separation `< threshold`:

`partial_t Psi_i -> partial_t Psi_i + beta * (phi_j - phi_i)`

with smooth gating.

Levels:
- `beta = 0`
- `beta = 0.01`
- `beta = 0.02`

Goal:
Test whether persistence failure is due to relative phase drift inside clustered encounters.

### Probe C - Self-induced effective curvature proxy
Modify kernel locally by graded energy density:

`K(x,y) -> K(x,y) * (1 + gamma * rho_graded(x,t))`

Levels:
- `gamma = 0`
- `gamma = 0.01`
- `gamma = 0.02`

Goal:
Test whether clustered DK textures require self-induced trapping geometry to stabilize.

## Run design

Minimal targeted matrix:
- `3` probes
- `3` strength levels each
- single representative collision texture

Total base runs: `9`

Refinement rule:

Promote only if:
- composite lifetime increases `>=` predefined threshold
- and separation trace shows bounded oscillatory regime
- and effect survives `n = 12 -> 24`

## Required outputs

Per run:
- separation vs time
- composite lifetime metric
- grade-weight evolution
- phase-difference trace
- collision classification label

Summary:
- persistence-gain comparison panel
- stability classification table
- grade-coherence evolution map

## Stage interpretation criteria

### Negative freeze
- all runs show only weak persistence gain
- no bounded composite regime
- refinement destroys apparent stabilization

### Exploratory positive
- sustained clustered oscillatory regime
- stable grade-exchange pattern
- persistence advantage vs control

### Regime-level positive
- clear bounded composite class
- robustness under refinement
- selective appearance (not universal across probes)

## Conceptual purpose

Stage 23.3 tests a precise hypothesis:

Dirac-Kahler collision textures fail to become persistent not because the fermionic sector is irrelevant, but because the clustered encounter lacks an internal stabilizing constraint.

This stage therefore searches for the minimal structural trigger that converts grade-exchange texture into a true composite regime.
