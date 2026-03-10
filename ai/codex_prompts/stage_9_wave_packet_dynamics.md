STAGE 9 MASTER CODEX PROMPT

Wave-Packet Dynamics on the Validated Operator Architecture

Goal

Test whether the validated kernel-induced cochain operator architecture supports coherent propagating localized excitations.

This stage introduces time evolution only.
No new operators or architectural modifications are allowed.

The objective is to determine whether disturbances governed by the existing operators propagate as structured wave packets rather than diffusing incoherently.

Architectural Constraints

The following components are frozen and must not be modified:
- Gaussian kernel substrate
- node Laplacian L0
- edge Hodge operator L1
- restricted transverse projector
- Dirac–Kähler operator D_DK
- weighting conventions
- cochain complex structure

All dynamics must be constructed from these operators only.

Dynamics Model

Introduce second-order wave dynamics using existing operators.

Scalar sector:

\partial_t^2 \phi = -L_0 \phi

Transverse edge sector:

\partial_t^2 A = -L_1 A

Optionally, evolution may be restricted to the transverse subspace using the existing projector.

The evolution scheme must be explicit and numerically stable.

Recommended method:
- leapfrog
- symplectic Verlet
- or other energy-preserving integrator.

Experiments

Stage 9 contains three diagnostic blocks.

BLOCK 1 — Scalar Wave Packet Propagation

Purpose

Test propagation of a localized packet governed by the scalar Laplacian.

Initial condition

Construct a localized Gaussian packet:

\phi(x,0) = \exp(-|x-x_0|^2/\sigma^2)

Initial velocity:

\partial_t \phi(x,0) = 0

Optional second run:
- add momentum phase to produce directed motion.

Measurements

Track:
- packet center position vs time
- packet width vs time
- norm conservation
- energy conservation
- anisotropy of spreading

Outputs

Data:

*_scalar_packet_dynamics.json

Plots:

*_scalar_packet_position.png
*_scalar_packet_width.png
*_scalar_packet_energy.png

Note file:

Scalar_Wave_Packet_Test_v1.md

BLOCK 2 — Transverse Edge Packet Propagation

Purpose

Test propagation of localized disturbances in the transverse edge sector.

Initial condition

Create a localized edge packet satisfying the transverse constraint:

d_0^* A = 0

Suggested construction:
- localized edge Gaussian
- project into the restricted transverse subspace.

Evolution

Use:

\partial_t^2 A = -L_1 A

with transverse projection maintained.

Measurements

Track:
- packet center motion
- packet width evolution
- constraint preservation
- norm conservation
- anisotropy

Compute:

||d0* A(t)||

at every time step.

Outputs

Data:

*_transverse_packet_dynamics.json

Plots:

*_transverse_packet_position.png
*_transverse_packet_width.png
*_transverse_packet_constraint.png
*_transverse_packet_energy.png

Note file:

Transverse_Wave_Packet_Test_v1.md

BLOCK 3 — Constraint Stability During Evolution

Purpose

Verify that the operator constraints remain preserved during dynamics.

Check:
1. longitudinal redundancy consistency
2. divergence constraint preservation

Diagnostics:

||d0* A(t)||
||P_transverse A(t) - A(t)||

where P_transverse is the restricted projector.

Expected behavior:
- residuals remain near numerical precision.

Outputs:

*_constraint_stability.json
*_constraint_residuals.png

Note file:

Constraint_Stability_Test_v1.md

Numerical Guidelines

Use periodic domains identical to Stage 8 tests.

Recommended lattice sizes:

n = 16
n = 20
n = 24

Time step must satisfy stability condition:

dt < 2/\sqrt{\lambda_{max}}

Compute λ_max from the operator spectrum.

Interpretation Boundary

Stage 9 tests propagation behavior only.

Allowed statements:
- coherent wave-packet propagation
- dispersion behavior
- constraint preservation during evolution
- anisotropy measurements

Forbidden statements:
- particles
- electrons
- photons
- quantum mechanics
- gauge dynamics
- relativistic invariance

Those require additional theoretical layers not tested here.

Expected Outcomes

Positive outcome:
- packets propagate smoothly
- spreading remains structured
- operator constraints remain preserved
- propagation speed correlates with transverse stiffness.

Negative outcome:
- packet fragmentation
- immediate diffusion
- strong lattice anisotropy
- constraint breakdown.

Either outcome is scientifically valuable.

Repository Integration

Create scripts under:

numerics/simulations/
stage9_wave_packet_scalar.py
stage9_wave_packet_transverse.py
stage9_constraint_stability.py

Notes:

experiments/wave_dynamics/

Artifacts:

data/
plots/

Use timestamped artifact naming identical to previous stages.

End Goal

Determine whether the validated kernel-induced operator architecture supports coherent propagating excitations governed by its existing operators.

This stage transitions the project from static operator diagnostics to dynamic behavior analysis.
