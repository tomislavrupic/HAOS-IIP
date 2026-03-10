STAGE 9B MASTER CODEX PROMPT

Dirac–Kähler Wave-Packet Dynamics on the Validated Cochain Architecture

Goal

Test whether the validated Dirac–Kähler operator

D_{DK} = d + \delta

supports coherent first-order packet propagation on the frozen cochain architecture.

This stage extends Stage 9 from scalar/transverse dynamics to the graded Dirac–Kähler sector.

The aim is strictly operator-level:
- coherent first-order propagation
- stability of the graded packet
- relation between DK propagation and the already validated Laplacian/Hodge structure

Do not introduce new operators.

Do not claim:
- fermions
- particles
- electrons
- quantum field theory
- gauge dynamics

Architectural Constraints

All previously validated structures are frozen:
- Gaussian kernel substrate
- cochain complex
- exterior derivatives d_k
- codifferentials \delta_k
- Hodge Laplacian \Delta_H
- Dirac–Kähler operator D_{DK}

No new discretization is allowed.

Dynamics Model

Use the validated first-order operator directly.

Primary evolution equation:

\partial_t \Psi = i D_{DK} \Psi

where \Psi is a graded cochain state:

in 2D:
\Psi \in \Omega^0 \oplus \Omega^1 \oplus \Omega^2

in 3D:
\Psi \in \Omega^0 \oplus \Omega^1 \oplus \Omega^2 \oplus \Omega^3

Optionally compare with the induced second-order evolution

\partial_t^2 \Psi = - D_{DK}^2 \Psi = - \Delta_H \Psi

to confirm consistency.

Use a unitary or norm-preserving scheme if possible (e.g. exponential integrator, Crank–Nicolson, or sufficiently stable small-step explicit integrator).

Experiments

Stage 9B contains four blocks.

BLOCK 1 — 2D DK Packet Propagation

Purpose

Test coherent propagation of a localized graded packet in the validated 2D DK system.

Initial condition

Construct a localized packet centered at x_0, with support distributed over grades.

Two variants:
1. single-grade packet
- start in \Omega^0 only
- then in \Omega^1 only
2. graded mixed packet
- small Gaussian packet distributed across all grades

Optional momentum kick:
multiply by a smooth phase to induce directed motion.

Measurements

Track over time:
- packet center
- packet width
- norm conservation
- energy expectation
- grade weight distribution
- sign-paired spectral decomposition

Compute:

w_k(t) = \| P_{\Omega^k} \Psi(t) \|^2

for each grade k.

Outputs

Data:
- *_DK_2D_packet_dynamics.json

Plots:
- *_DK_2D_packet_position.png
- *_DK_2D_packet_width.png
- *_DK_2D_packet_norm.png
- *_DK_2D_packet_grade_weights.png

Note:
- DK_2D_Packet_Dynamics_v1.md

BLOCK 2 — 3D DK Packet Propagation

Purpose

Repeat the same test on the validated 3D DK architecture.

Initial condition

Use localized graded packets on

\Omega^0 \oplus \Omega^1 \oplus \Omega^2 \oplus \Omega^3

Start with:
- one packet concentrated in \Omega^1
- one packet mixed across all grades

Measurements

Same as Block 1, plus:
- anisotropy of propagation
- sensitivity to flat holonomy
- low-mode participation

Outputs

Data:
- *_DK_3D_packet_dynamics.json

Plots:
- *_DK_3D_packet_position.png
- *_DK_3D_packet_width.png
- *_DK_3D_packet_norm.png
- *_DK_3D_packet_grade_weights.png
- *_DK_3D_packet_anisotropy.png

Note:
- DK_3D_Packet_Dynamics_v1.md

BLOCK 3 — First-Order vs Second-Order Consistency

Purpose

Verify that first-order DK propagation and second-order Laplacian propagation remain consistent at the packet level.

Procedure

For the same initial packet, evolve:
1. with first-order DK dynamics

\partial_t \Psi = i D_{DK} \Psi

2. with second-order induced dynamics

\partial_t^2 \Psi = - \Delta_H \Psi

Compare:
- packet center
- packet width
- total norm / energy
- spectral support

Outputs

Data:
- *_DK_first_vs_second_order.json

Plots:
- *_DK_first_vs_second_order_comparison.png

Note:
- DK_First_vs_Second_Order_v1.md

Interpretation

A clean match supports that the first-order propagation really reflects the same validated graded operator structure.

BLOCK 4 — Flat-Holonomy DK Packet Response

Purpose

Test whether coherent DK packet propagation survives flat cycle holonomy, consistent with the Stage 6–7 spectral results.

Procedure

Run the 2D or 3D packet with several flat torus-cycle phase choices.

Measure:
- propagation speed
- packet spreading
- norm preservation
- grade redistribution

Outputs

Data:
- *_DK_holonomy_packet_response.json

Plots:
- *_DK_holonomy_packet_paths.png
- *_DK_holonomy_packet_widths.png
- *_DK_holonomy_grade_weights.png

Note:
- DK_Holonomy_Packet_Response_v1.md

Diagnostics to Preserve

For all blocks, record:
- norm conservation
- spectral decomposition symmetry
- packet center trajectory
- packet width evolution
- grade weight evolution
- anisotropy where relevant

Success Criteria

A successful Stage 9B result would show:
1. stable norm-preserving or approximately norm-preserving evolution
2. coherent packet propagation over multiple time steps
3. controlled spreading rather than immediate fragmentation
4. structured grade mixing rather than random grade leakage
5. consistency between DK first-order evolution and its second-order induced form

A failed result would show:
- immediate incoherent diffusion
- loss of norm / instability
- random grade scrambling
- strong unphysical anisotropy
- mismatch with the induced second-order propagation

Either result is scientifically useful.

Interpretation Boundary

Allowed statements:
- the Dirac–Kähler sector supports coherent first-order propagation
- packet evolution reflects the graded cochain structure
- the first-order and second-order forms are dynamically consistent
- flat holonomy modifies propagation smoothly

Not allowed:
- electrons
- fermions as particles
- quantum mechanics
- field quantization
- matter emergence

This stage is still strictly operator-level and kinematic.

Repository Integration

Create scripts under:

numerics/simulations/
    stage9b_DK_packet_2D.py
    stage9b_DK_packet_3D.py
    stage9b_DK_first_vs_second_order.py
    stage9b_DK_holonomy_packet.py

Create notes under:

experiments/wave_dynamics/

Store timestamped artifacts under:

data/
plots/

Use the same timestamped naming convention as all earlier stages.

End Goal

Determine whether the validated Dirac–Kähler operator on the frozen architecture supports coherent graded packet propagation.

This would extend the project from
- static DK factorization
to
- dynamic DK behavior.
