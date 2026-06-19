# HAOS-IIP Bridge Methods Spine

This is the method-facing extraction from the bridge corpus.

It is written as a bridge toolkit, not a theory paper.

## 1. Spectral identity

### Methods to keep

- unnormalized graph Laplacian: `L = D - A`
- symmetric normalized Laplacian: `L_sym = I - D^(-1/2) A D^(-1/2)`
- random-walk normalized Laplacian: `L_rw = I - D^(-1) A`
- low-mode eigenvectors as stable identity coordinates
- Fiedler vector as the second nontrivial low mode
- spectral gap as a stability / separation diagnostic
- diffusion maps as slow relational coordinates
- graph spectral wavelets as multiscale localized probes

### HAOS-IIP use

- define or audit spectral address
- compare target/control low-mode stability
- test transformation recovery bottlenecks
- evaluate refinement persistence
- build multiscale holdout diagnostics

### Required output fields

- Laplacian type
- normalization
- eigenvalue window
- eigenvector sign/alignment policy
- spectral gap
- Fiedler sign stability
- control separation
- holdout transfer

## 2. Hodge sectors

### Methods to keep

- cochain spaces as data carriers
- boundary / coboundary maps as structure operators
- Hodge Laplacian on degree `k` cochains
- exact / coexact / harmonic sector split
- persistent Laplacians across filtrations or refinements
- localized Hodge spectra / Hodgelets

### HAOS-IIP use

- separate current closure into sector-level signals
- report harmonic residue as candidate stable global address
- distinguish local gradient-like response from circulation-like response
- avoid compressing sector behavior into one scalar

### Required output fields

- cochain degree
- operator family
- sector energy / mass by component
- harmonic dimension or near-null count
- exact/coexact/harmonic residuals
- refinement persistence
- missing-sector policy

## 3. Curvature and flow

### Methods to keep

- Ollivier-Ricci curvature as transport-overlap curvature
- Forman curvature as fast combinatorial curvature
- Cheeger / conductance sweep cuts as bottleneck diagnostics
- curvature-flow-style graph reweighting as diagnostic perturbation

### HAOS-IIP use

- identify local-to-global geometry stress points
- test whether bridge signals collapse under topology destruction
- compare target graph geometry against matched controls
- add curvature-specific controls to geometry bridge hardening

### Required output fields

- curvature type
- graph distance / cost definition
- lazy-random-walk or neighborhood measure policy
- transport approximation policy
- curvature mean / variance / quantiles
- negative-edge fraction
- conductance / sweep-cut estimate
- control degradation

## 4. Transport metrics

### Methods to keep

- exact Wasserstein on small supports
- Sinkhorn approximation for larger supports
- sliced Wasserstein as projection-based approximation
- Wasserstein barycenter when consensus geometry is needed

### HAOS-IIP use

- compare spectral embedding distributions
- compare curvature distributions
- compare closure signatures as feature distributions
- test whether controls preserve or destroy the relevant distribution

### Required output fields

- support definition
- cost matrix
- units or dimensionless policy
- regularization `epsilon` for Sinkhorn
- projection count for sliced Wasserstein
- raw value
- normalized value
- sensitivity to epsilon / projection count

## 5. Kernel distances

### Methods to keep

- MMD as a kernel two-sample statistic
- Gaussian/RBF kernel with frozen bandwidth policy
- Laplace or rational-quadratic kernel for robustness checks
- permutation or bootstrap uncertainty as uncertainty only, not synthetic
  evidence

### HAOS-IIP use

- target-control feature distribution testing
- holdout degradation testing
- detecting whether a proposed bridge signal separates at distribution level

### Required output fields

- kernel family
- bandwidth policy
- sample sizes
- MMD value
- uncertainty method
- shuffled-label null
- topology-destroyed null

## 6. Discrete field-theory analogy

### Methods to keep

- differential-form / cochain analogy
- Dirac-Kahler operator comparison
- lattice/geometric discretization as precedent

### HAOS-IIP use

- make DK-language comparisons precise
- constrain speculative field-theory notes
- prevent notation from doing unsupported physical work

### Required output fields

- analogy status
- mapped HAOS-IIP object
- mapped literature object
- missing physical ingredients
- disallowed claim language

## 7. Shared controls

Every bridge candidate using this spine should include:

- shuffled labels
- topology destruction
- degree-preserving rewiring
- weight shuffling
- parameter-matched null
- refinement-broken control
- perturbation-free baseline where available
- best conventional baseline
- leakage positive control where applicable

## 8. Shared failure gates

A bridge candidate should fail or remain open if:

- target improvement does not beat best conventional baseline
- controls preserve the measured signal
- normalization changes the verdict
- holdout transfer fails
- one metric silently dominates
- literature analogy is used as a derivation
- physical units / observation map are missing
- generated observables import the target answer

## 9. Minimal bridge package

The minimum defensible HAOS-IIP bridge package should report:

- source artifact
- mapping status: `DERIVED`, `CALIBRATED`, `HEURISTIC`, `ANALOGICAL`, or
  `UNAVAILABLE`
- method family
- output units or dimensionless policy
- baselines
- controls
- holdout split
- uncertainty
- verdict logic
- claim ceiling

## 10. Claim ceiling

This methods spine supports bridge instrumentation.

It does not establish:

- physical ontology
- continuum spacetime
- quantum correlations
- gravitational dynamics
- external empirical validation
