# HAOS-IIP Celestial Holography Audit Roadmap

Status: queued after Phase 56 closure  
Scope: future physics-bridge sidecar, not HAOS core  
Version: v1  

## Purpose

This note saves the planned celestial-holography stress test for HAOS-IIP.

The purpose is not to claim that HAOS-IIP currently reproduces celestial
holography, BMS symmetry, Virasoro structure, gravitational memory, or an
S-matrix dictionary. The purpose is the opposite: use celestial holography as a
high-standard boundary audit for HAOS-IIP physics-facing language.

The current fix is a claim-boundary fix, not a core rewrite.

Correct boundary:

> HAOS-IIP currently supplies a reproducible pre-geometric interaction filter.
> Any contact with established physics must pass separate symmetry, boundary,
> scattering, and observable tests.

## Why This Matters

Celestial holography is top-down, symmetry-driven, and tied to flat-space
scattering structure: asymptotic symmetries, soft limits, memory effects,
celestial sphere operators, Mellin-basis amplitudes, and CFT-like correlators at
null infinity.

HAOS-IIP is bottom-up and deliberately stripped: frozen branch-local
cochain-Laplacian hierarchies, spectral feasibility, persistence,
recoverability, proto-geometry, cautious continuum proxies, and external
physics-bridge observables.

The two programs are not interchangeable. The celestial lens exposes exactly
where HAOS-IIP must avoid overclaiming:

- spectral smoothness is not conformal symmetry;
- power-law scaling is not a celestial amplitude;
- proto-geometry is not null infinity;
- recoverability is not unitarity;
- physics-facing proxies are not S-matrix observables;
- BMS/Virasoro-like structure must be tested, not inferred.

## Non-Goals

This roadmap does not attempt to:

- retrofit celestial holography into HAOS core;
- claim BMS, Virasoro, soft-theorem, or S-matrix recovery;
- derive gravity, spacetime, or a celestial CFT from HAOS;
- replace the existing physics bridge;
- weaken the frozen reproducibility discipline.

If a future test fails, the failure should remain visible.

## Entry Conditions

Start this roadmap only after Phase 56 is frozen enough that the biology/FMO
branch is no longer consuming active instrumentation attention.

Minimum entry conditions:

1. Phase 55.2 microtubule robustness remains reproducible.
2. Phase 56 FMO status is frozen honestly, whether PASS, MARGINAL, or FAIL.
3. No HAOS core files need modification.
4. The physics bridge remains explicitly external post-processing.

## Phase 57: Celestial Boundary Audit

Folder target:

```text
experiments/physics_bridge/celestial_boundary_audit/
```

Goal:

Formalize what HAOS-IIP does not yet prove when judged by celestial-holography
standards.

Deliverables:

- `README.md`
- `celestial_boundary_audit.md`
- `haos_vs_celestial_requirements.csv`
- `bridge_status.json`

Audit questions:

- Does HAOS recover an S2-like asymptotic boundary, or only local spectral
  envelopes?
- Does it recover conformal weights, or only graph-Laplacian scaling proxies?
- Does it preserve soft/collinear factorization, or only generic perturbation
  recoverability?
- Does it represent BMS-like charges, or only low-mode persistence?
- Does any current physics-bridge row overstate contact with real scattering
  observables?

Expected status:

`OPEN`, unless explicit evidence closes a narrower row.

Success condition:

The audit produces a clear boundary table and reduces overclaim risk.

## Phase 58: Spherical Harmonic Control Probe

Goal:

Build a known-target control benchmark before comparing HAOS structures to
celestial-style boundary data.

Experiment:

- construct a discrete sphere graph;
- compute graph Laplacian modes;
- compare against spherical-harmonic-like mode organization;
- measure refinement drift, mode mixing, degeneracy splitting, and spectral
  leakage;
- compare against HAOS frozen hierarchies only as a sidecar.

Required distinction:

Spherical-harmonic-like mode recovery is not celestial holography. It is only a
minimal boundary-geometry sanity check.

Success condition:

HAOS telemetry distinguishes real spherical boundary structure from generic
spectral smoothness under controls.

Failure condition:

Generic graph smoothness passes the same tests. In that case, the proxy is too
weak and must not be promoted.

## Phase 59: Soft-Theorem Proxy Test

Goal:

Create a toy amplitude-like benchmark with known soft/factorization behavior and
test whether HAOS-style telemetry can detect the relevant failure modes.

Experiment:

- generate toy amplitude tables with controlled soft poles/residues;
- add controls that preserve spectrum while breaking soft structure;
- measure residue stability, factorization drift, pole-location drift, and null
  specificity;
- compare standard telemetry against soft-structure-specific metrics.

Success condition:

HAOS-style telemetry, or a clearly named extension of it, detects soft-structure
breakage better than generic spectral controls.

Failure condition:

Telemetry preserves only broad smoothness or recoverability while missing
residue/factorization errors. That would keep the celestial bridge `OPEN`.

## Phase 60: Claim-Gated Physics Bridge Update

Goal:

Update the physics bridge language so future readers cannot confuse
pre-geometric proxies with established flat-space holography.

Allowed language:

- pre-geometric telemetry;
- recoverability and persistence;
- spectral emergence;
- continuum-sketch diagnostics;
- symmetry stress tests;
- external physics-facing proxy rows.

Disallowed until separately closed:

- S-matrix reconstruction;
- BMS recovery;
- Virasoro/celestial CFT recovery;
- exact soft-theorem recovery;
- gravitational memory prediction;
- flat-space gravity closure;
- source-code-of-physics language.

Deliverables:

- updated physics bridge note;
- explicit celestial boundary table;
- `PASS` / `OPEN` / `WATCH` rows for any celestial-facing proxy;
- no core changes unless a future phase explicitly justifies them.

## Decision Rule

Do not fix HAOS core for celestial holography.

Fix only the interface between HAOS and physics-facing claims:

1. If a celestial-facing proxy passes, mark the specific proxy `PASS`.
2. If it fails, mark it `OPEN`.
3. If generic controls also pass, mark it `MARGINAL` or `CONTROL_MATCH`.
4. If the test requires imposed asymptotic structure, state that the structure is
   borrowed, not emergent.

## Working Summary

Celestial holography should be used as a boundary validator for HAOS-IIP.

It does not invalidate HAOS-IIP as a reproducible interaction-invariance lab.
It does invalidate any premature claim that current HAOS physics-bridge proxies
already close on flat-space scattering, celestial operators, or asymptotic
gravity.

The future work is therefore disciplined:

> audit first, benchmark second, soft-structure proxy third, claim-gate last.
