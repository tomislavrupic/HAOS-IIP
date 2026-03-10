# Codex Prompt — Stage 6: Dirac–Kähler Lift on the Existing Cochain Architecture

You are working inside the HAOS-IIP repository.

The existing operator architecture is frozen.

Do not modify:
- kernel definition
- scalar operator L0
- edge/Hodge operator L1
- restricted transverse sector
- Stage 1–5 experiments or interpretations

Stage 5 tested a node-spinor directional kernel lift and rejected it because the square-root criterion failed:

D_H^\dagger D_H \approx L0

Stage 6 switches to a Dirac–Kähler lift built directly from the validated discrete cochain operators.

The goal is to test whether a graded cochain Dirac–Kähler operator can be constructed such that its square reproduces the discrete Hodge Laplacian exactly or to numerical precision.

This is an operator-level test only.

Do not claim:
- fermions
- quantum field theory
- particle interpretation
- gauge dynamics

## Core Stage 6 hypothesis

On a graded cochain space, define

D_DK = d + \delta

where:
- d is the discrete exterior derivative
- \delta is the codifferential (adjoint of d under the weighted inner product)

Then test whether

D_DK^2 = \Delta_H

where \Delta_H is the full graded Hodge Laplacian.

## Stage 6 scope order

First prototype: 2D periodic cell complex only

Use a 2D periodic square cell complex (discrete 2-torus).

Graded cochain space:

\Omega^0 \oplus \Omega^1 \oplus \Omega^2

Required operators:
- d0 : \Omega0 -> \Omega1
- d1 : \Omega1 -> \Omega2
- \delta1 = d0*
- \delta2 = d1*

Dirac–Kähler operator:

D_DK = d0 + d1 + \delta1 + \delta2

acting on the full graded space.

Only after this passes should any 3D extension be attempted.

Do not implement 3D in this stage unless the 2D prototype is clearly validated.

## Experiment Block 1 — 2D Cochain Infrastructure Check

Create:

numerics/simulations/DK_stage6_common.py

Build the weighted periodic 2D cell complex using the same kernel-derived weight logic already used in the repo.

Required outputs:
- node set
- oriented edges
- oriented plaquettes
- incidence matrices
- weighted differentials d0, d1
- adjoints \delta1, \delta2

Verify cochain identities:
- d1 d0 = 0
- \delta1 \delta2 = 0

to numerical precision.

Output note:

experiments/spinor_sector/DK_2D_Cochain_Infrastructure_v1.md

## Experiment Block 2 — Dirac–Kähler Square Test

Create:

numerics/simulations/DK_square_test.py

Construct the graded block operator

D_DK = d + \delta

on

\Omega0 \oplus \Omega1 \oplus \Omega2

Construct the graded Hodge Laplacian

\Delta_H = diag(\delta1 d0, d0 \delta1 + \delta2 d1, d1 \delta2)

Test numerically:
1. operator residual

|| D_DK^2 - \Delta_H ||
2. relative Frobenius error
3. eigenvalue comparison between D_DK^2 and \Delta_H
4. blockwise agreement on 0-, 1-, and 2-form sectors

Expected result:

The square relation should hold exactly up to numerical precision if the implementation is correct.

Outputs:
- JSON artifact
- plot: DK_square_vs_Hodge.png
- note: DK_Square_Test_v1.md

## Experiment Block 3 — Spectrum of the Graded Dirac–Kähler Operator

Create:

numerics/simulations/DK_spectrum_structure.py

Compute the low spectrum of D_DK.

Measure:
- +/- spectral pairing
- low-mode degeneracy structure
- distribution of mode weight across \Omega0, \Omega1, \Omega2
- IPR of graded modes

Goal:

Determine whether the graded operator shows clean Dirac-like spectral structure on the validated cochain complex.

Outputs:
- DK_spectrum.png
- DK_mode_grading.png
- DK_ipr_vs_mode.png
- note: DK_Spectrum_Structure_v1.md

Interpretation boundary:

Only report spectral structure on the graded cochain space. Do not call these fermions.

## Experiment Block 4 — Flux Response in 2D

Create:

numerics/simulations/DK_flux_response_test.py

Reuse the existing periodic gauge / Landau-type flux infrastructure in 2D if compatible.

Inject flux through the 2D torus cycles while keeping the cochain architecture fixed.

Measure:
- low spectral flow of D_DK
- stability of the square relation under flux
- mode grading under flux

Outputs:
- DK_flux_spectral_flow.png
- DK_flux_square_residual.png
- note: DK_Flux_Response_Test_v1.md

Interpretation boundary:

Only report whether the graded operator responds smoothly to background flux while preserving the operator square relation.

## Experiment Block 5 — Optional 2D Continuum Comparison

Create:

numerics/simulations/DK_2D_continuum_comparison.py

Compare low graded spectra against the expected 2D Hodge/Laplacian structure on the torus.

Goal:

Check whether the low graded spectrum of D_DK is consistent with the known continuum cochain sector structure in 2D.

Outputs:
- DK_2D_continuum_comparison.png
- note: DK_2D_Continuum_Comparison_v1.md

This block is optional if Blocks 1–4 already give a decisive validation.

## Standard diagnostics

Every Stage 6 experiment must produce:
- stamped JSON output
- stamped plots
- experiment note in experiments/spinor_sector/
- log entry in EXPERIMENT_LOG.md

Use naming pattern:

data/YYYYMMDD_experiment_name.json
plots/YYYYMMDD_plotname.png

## Output requirements

Create only new files.

Do not overwrite earlier stage artifacts.

Do not generate a paper yet.

## Interpretation boundary

Stage 6 may support statements like:
- “a Dirac–Kähler lift exists on the validated cochain architecture”
- “the graded square relation holds numerically”
- “the graded low spectrum exhibits stable pairing and cochain-sector structure”
- “the operator responds smoothly to flux while preserving the square identity”

Stage 6 must not claim:
- physical fermions
- spin-1/2 particles
- quantum field theory
- emergent matter

## Success criterion

The stage is considered successful if:
1. D_DK^2 = \Delta_H holds to numerical precision
2. the low graded spectrum is stable and structured
3. the square relation survives mild flux insertion

If these fail, reject the prototype and document the failure exactly as in Stage 5.

## Deliverables summary

Produce:
1. 2D cochain infrastructure validation
2. Dirac–Kähler square test
3. graded spectrum analysis
4. flux response test
5. optional continuum comparison

No 3D extension in this stage unless 2D is clearly validated first.

This is the right Stage 6 because it stays fully aligned with what the repo already knows how to do:
- validated d
- validated \delta
- validated Hodge structure
- validated scalar and transverse sectors
