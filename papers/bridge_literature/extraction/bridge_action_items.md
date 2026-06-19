# Bridge Action Items from the Literature Extraction

This file converts the extracted literature into concrete HAOS-IIP bridge work.

## A. Geometry bridge hardening

Use spectral clustering, diffusion maps, Hodgelets, graph wavelets, and
curvature diagnostics to harden geometry recovery.

### Add or verify

- normalized Laplacian variants: unnormalized, symmetric normalized, random-walk
  normalized
- multieigenvector spectral embedding
- Fiedler sign stability
- Cheeger sweep-cut estimate
- diffusion-map coordinates for holdout relation prediction
- spectral graph wavelet / Hodgelet-style multiscale localization
- Ollivier-Ricci curvature summary
- Forman curvature baseline
- curvature-flow perturbation probe
- MMD and sliced Wasserstein comparisons between target and controls

### Required controls

- topology destruction
- degree-preserving rewiring
- weight shuffling
- label permutation
- refinement-broken control
- kernel / bandwidth sensitivity
- normalized-vs-unnormalized Laplacian ablation

### Verdict discipline

Do not promote geometry recovery unless:

- holdout transfer passes
- control degradation is component-specific
- normalization sensitivity is bounded
- target improvement survives best conventional baselines

## B. Hodge / cochain bridge note

Use Hodge decomposition and persistent Laplacian literature to tighten the
existing HAOS-IIP cochain language.

### Extract into a technical note

- cochain space definition
- boundary / coboundary operator mapping
- Hodge Laplacian mapping
- exact / coexact / harmonic sectors
- harmonic residue as recoverable global signature
- persistent Laplacian as refinement / filtration diagnostic
- Hodgelets as localized sector probes

### Required claim boundary

The note must state:

- Hodge sector split is structural
- it is not a physical field equation
- it does not define a gauge group
- it does not derive fermions or continuum dynamics

## C. Curvature bridge note

Use Ollivier-Ricci, Forman-Ricci, and curvature-flow papers to define a
claim-gated curvature diagnostic layer.

### Extract into a technical note

- neighborhood-overlap curvature
- bridge / bottleneck interpretation
- Forman curvature as low-cost baseline
- curvature flow as diagnostic reweighting
- random geometric convergence as high-bar caution

### Required claim boundary

The note must state:

- graph curvature is not gravity
- curvature flow is not spacetime dynamics
- smooth-limit or Ricci-curvature claims require additional assumptions

## D. Transport / comparison metric layer

Use OT, Sinkhorn, sliced Wasserstein, barycenters, and MMD as metric tools.

### Add or verify

- exact or approximate Wasserstein on small supports
- Sinkhorn with frozen epsilon
- sliced Wasserstein as a cheap surrogate
- MMD with predeclared kernels and bandwidth policy
- transport/mmd ablations against existing closure-distance components
- raw and normalized metric reporting

### Required controls

- shuffled labels
- topology-destroyed features
- distribution-preserving null
- bandwidth / epsilon sensitivity
- best conventional baseline comparison

## E. Discrete field-theory analogy guardrail

Use Dirac-Kahler and geometric discretization papers as analogy anchors only.

### Extract into a note

- differential forms and cochains
- Dirac-Kahler operator as comparison object
- lattice/geometric discretization precedent
- distinction between analogy and derivation

### Required claim boundary

The note must state:

- HAOS-IIP cochain structure is not automatically a lattice field theory
- no fermion, gauge, or continuum claim follows from notation alone
- any future physical bridge must generate observables independently

## F. Immediate next bridge artifact

The clean next artifact is:

`papers/bridge_literature/extraction/bridge_methods_spine.md`

Recommended sections:

1. spectral identity
2. Hodge sectors
3. curvature
4. transport metrics
5. kernel distances
6. discrete field analogy
7. controls and failure gates

That spine can then feed future geometry-bridge or HBP work without creating a
new physical claim.
