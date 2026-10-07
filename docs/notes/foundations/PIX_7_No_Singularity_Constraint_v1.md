# PIX-7 — No-Singularity Constraint v1

Date: 2026-09-08. Status: OPEN HYPOTHESIS / RESEARCH SPECIFICATION.
Origin: user-proposed NSP and C1–C10; formal qualifications below are research-design additions.
This document establishes no new physical result, gravitational correspondence, or change to existing HAOS-IIP evidence classifications. It is an internal specification, not a validated public claim.

## Hypothesis and scope

**NSP — No-Singularity Principle:** No physically realizable state contains divergent local observables. Singularities in an effective description indicate breakdown of that description rather than an infinite physical quantity.

NSP is a hypothesis. Universal upper bounds and complete evolution are stronger requirements imposed by this programme; neither follows merely from pointwise finiteness. Saturating interactions are one candidate mechanism, not part of NSP's definition. A successful nonsaturating mechanism would challenge that mechanism preference without refuting NSP.

The target is an interaction law L, with explicit degrees of freedom, matter coupling, units, locality or nonlocality, constraints, initial-data domain D, observable definitions, and continuation rule. A schematic correction tensor alone is not a candidate law.

## Mathematical contract

Use c = 1 unless restored explicitly. Before evaluating a candidate, specify D independently of whether its solutions remain regular. Declare topology, boundary conditions, regularity, matter models, parameter ranges, and the operational class of observers or probes.

For each declared observable O_a and evolution U_L(t;d), seek an a priori estimate

\[
\sup_{d\in D}\sup_{t\in J_d}\sup_{p\in M_d}
|O_a[U_L(t;d)](p)|\le B_a<\infty.
\]

The quantifier over D is the strong, uniform-bound version. A bound B_a(d) depending on initial data is a weaker result and must be labelled separately. Finiteness at each finite time does not imply a uniform bound as time tends to infinity. Bounds carry observable-specific dimensions; powers of a bounded observable do not share its numerical bound.

For Lorentzian curvature, require absolute values, for example

\[
|R|\le B_R,\quad |R_{\mu\nu}R^{\mu\nu}|\le B_2,\quad
|R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}|\le B_4.
\]

These contractions are not positive-definite tensor norms. Bounding this short list does not establish regularity of all tidal measurements. Include parallel-propagated tidal components along specified freely falling probes, matter observables, and whatever derivative estimates the continuation theorem needs. Density must specify a matter rest frame or observer: \(\rho_u=T_{\mu\nu}u^\mu u^\nu\). A universal cap over arbitrarily boosted observers is a different, generally incompatible demand from finite measurements along each physical observer trajectory.

For a classical metric domain, require all inextendible timelike and null geodesics to have unbounded proper and affine parameter respectively in the maximal physical extension. Coordinate patches do not determine completeness. Completeness, uniqueness, and continuous dependence of dynamical evolution are separate obligations: extension past a Cauchy horizon alone does not settle predictability. If geometry ceases to exist, give a fundamental continuation rule and a derived account of the effective endpoint; label this **replacement of the geometric completeness criterion**, not proof of geodesic completeness.

## Canonical C1–C10 gates

The final ten-item numbering supersedes the shorter preliminary numbering in the proposal.

| Gate | Required evidence | Rejection condition |
| --- | --- | --- |
| C1 Bounded observables | Declared observable/probe domain and a priori bounds | An admissible evolution violates a bound; a scalar-only test misses a physical divergence |
| C2 Complete evolution | Global continuation, causal-geodesic audit, and deterministic evolution where claimed | Finite physical endpoint, or unresolved loss of predictability |
| C3 GR recovery | Controlled expansion and quantitative agreement in tested regimes | Wrong leading equations or disagreement beyond declared observational tolerances |
| C4 Newton recovery | Weak-field, slow-motion derivation of Poisson dynamics with calibrated G | Wrong force law or uncontrolled limit |
| C5 Conservation closure | Covariant identity plus matter/extra-sector equations and constraint propagation | Unbalanced exchange or inconsistent constraints |
| C6 Pathology exclusion | Degree-of-freedom count, physical kinetic analysis, characteristic and well-posedness analysis | Physical ghost, uncontrolled instability, or unaccounted causal failure within the claimed domain |
| C7 Dynamical regularization | Explicit dynamics and estimates explaining arrest, bounce, reorganization, or another finite outcome | Clipping, numerical cutoff dependence, or merely assumed regular metric |
| C8 Independent admissibility | D fixed before testing and justified by constraints and preparation physics | Post hoc removal of failing data or regularity used to define admissibility |
| C9 Distinguishing prediction | Parameter-to-observable map, GR comparator, uncertainties and discriminating measurement | No calculated distinction, or complete nuisance degeneracy in the declared test |
| C10 Origin of apparent singularity | Controlled comparison to the truncated GR description for matched admissible data | No demonstration of where and why truncation fails |

Initial data need not include already divergent, ill-defined configurations. C8 requires explaining admissibility without assuming the desired conclusion. A symmetric solution is a scoped result; it cannot pass gates for generic inhomogeneous, anisotropic, or rotating evolution.

## Recovery and conservation

Write a possible architecture as

\[
G_{\mu\nu}+\Lambda g_{\mu\nu}+H_{\mu\nu}[g,\psi]
=8\pi G T_{\mu\nu}.
\]

Specify all dimensionless expansion parameters, including relevant curvature scales, gradients, frequencies and extra fields. Small R alone cannot define the GR regime: vacuum GR has R = 0 while its Weyl curvature can be large. Demonstrate suppression of corrections in solutions and observables, not merely H(0) = 0 in an equation. In restored units the Einstein source coefficient is 8πG/c⁴.

For separately conserved minimally coupled matter, the full equations must imply ∇^μH_{μν} = 0. With exchanging sectors, derive total conservation and the exchange law. A diffeomorphism-invariant action is a useful construction route, with the identity evaluated consistently with all field equations. It does not automatically establish stability. Higher derivatives require constraint and mode analysis rather than automatic rejection; an effective theory must also state its cutoff.

## Why saturation alone fails

Consider a schematic scalar response, not a covariant gravitational theory:

\[
F(K)=F_*\tanh(K/K_*),\qquad F(K)=\kappa\rho.
\]

Although F is bounded, its inverse is

\[
K=K_*\operatorname{artanh}(\kappa\rho/F_*).
\]

As the source approaches F*/κ from below, K diverges; above it there is no real solution. Thus bounded response can coexist with divergent curvature proxy or loss of solvability. Conversely, writing K = K*tanh(ρ/ρ*) bounds K but leaves ρ unbounded and supplies no conservation or evolution law. Neither passes C1–C7.

The required theorem concerns solutions of coupled dynamics, not the range of an isolated constitutive function. Any proposed evasion of a singularity theorem must identify the theorem's assumptions and which cease to hold.

## Observation and first research milestone

G1: weak-field dynamics and precision gravity. G2: inspiral, polarization and propagation. G3: exterior geometry and ringdown within actual measurement uncertainties; observational agreement does not require exact Kerr everywhere. G4: finite extreme-interior or early-universe predictions. An inaccessible interior is a theoretical consistency test until connected to data. G5: a discriminating observable with amplitude, parameter dependence, uncertainty, and feasibility.

The first milestone is a **candidate audit**, not a new equation: choose one fully specified action or evolution system; freeze D and its GR comparator; derive constraint propagation, recovery limits and one collapse-sector solution; test perturbations and identify the route to generic evolution. Numerical checks must include convergence and sensitivity to regulators. Finite runs provide bounded numerical evidence, never a proof of universal completeness. Missing proof means OPEN or NOT_DEMONSTRATED, not a failed empirical prediction. One admissible counterexample rejects the candidate's universal claim, not NSP for every possible theory.

C9 remains required by this programme, but logical distinguishability, present detectability, and uniqueness relative to other modified theories must be reported separately. No candidate has been established **in this document** to satisfy all ten gates; this is not an exhaustive literature verdict.

## Navier–Stokes boundary

The separate physical target is microscopic dynamics → kinetic description → hydrodynamic limit. Specify coarse-graining, Knudsen number, relaxation scales and convergence assumptions. Breakdown can occur when gradients or relaxation times invalidate local equilibrium, not only at a molecular length.

Subluminal physical three-velocity belongs to a relativistic description; appending |u| < c to classical incompressible Navier–Stokes changes its mathematical problem. Finite velocity alone supplies neither gradient estimates nor a global smoothness proof. A molecular or kinetic replacement does not settle the Clay problem for the original PDE. Density, pressure, temperature and gradients require operational definitions and separately justified bounds.

## HAOS-IIP boundary and sources

HAOS-IIP must independently supply a bridge from interaction variables to physical geometry, matter, units and observables, followed by GR recovery. Finite graph telemetry and bounded toy recoverability do not establish that bridge or NSP. Existing frozen APIs, results and evidence classifications remain authoritative.

Primary research anchors, consulted 2026-09-08; these are starting points, not a systematic review:

- [Chamseddine and Mukhanov, Resolving Cosmological Singularities](https://arxiv.org/abs/1612.05860): limiting-curvature construction with particular cosmological bounce solutions; not evidence for all C1–C10.
- [Chamseddine and Mukhanov, Nonsingular Black Hole](https://arxiv.org/abs/1612.05861): a related black-hole construction to audit rather than assume generally complete.
- [LIGO–Virgo–KAGRA, Tests of General Relativity with GWTC-3](https://arxiv.org/abs/2112.06861) and [associated data release](https://dcc-lho.ligo.org/LIGO-P2100456/public): a concrete observational benchmark, not asserted to be the latest catalogue.
- [Clay Mathematics Institute, Navier–Stokes Equation](https://www.claymath.org/millennium/Navier-Stokes-Equation/): the independent mathematical existence and smoothness problem.

Classification: NSP OPEN; bounded-interaction mechanism OPEN; GR/Newton recovery REQUIRED; candidate passing C1–C10 NOT_DEMONSTRATED; specific falsifiable prediction NOT_YET_SPECIFIED. The specification is test-oriented, while empirical testability requires an actual candidate and measurement map.
