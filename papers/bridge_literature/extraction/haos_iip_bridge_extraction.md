# HAOS-IIP Bridge Extraction

This note extracts the bridge-relevant material from the PDF corpus.

It is organized around what HAOS-IIP can actually use:

1. spectral identity and embedding
2. Hodge / cochain sector structure
3. graph curvature and local-to-global geometry
4. optimal transport and kernel comparison metrics
5. discrete field-theory analogies
6. caution gates

## Extraction basis

- Source folder: `papers/bridge_literature/pdfs/`
- PDFs scanned: 47
- Text extraction status: 47/47 extractable
- Machine manifest: `pdf_extraction_manifest.csv`

The manifest is evidence for coverage and keyword emphasis only. It is not a
semantic proof.

## 1. Spectral identity and embedding

### What to keep

From spectral clustering, diffusion maps, graph wavelets, and normalized-cut
papers, HAOS-IIP needs:

- graph Laplacian and normalized Laplacian definitions
- low-mode eigenvectors as coherent structure coordinates
- Fiedler vector / second eigenmode as global partition diagnostic
- spectral gap as a stability / separation signal
- diffusion-map coordinates as slow relational coordinates
- spectral wavelets as multiscale localized structure probes
- normalized-cut / Cheeger language for bottlenecks

### HAOS-IIP mapping

| Literature object | HAOS-IIP bridge object | Status |
| --- | --- | --- |
| Laplacian eigenmode | spectral address / recoverable spectral identity | `DIRECT` |
| Low-mode embedding | proto-geometric coordinate system | `DIRECT` |
| Fiedler vector | transformation / partition diagnostic | `DIRECT` |
| Spectral gap | stability and separation margin | `DIRECT` |
| Diffusion coordinate | relational distance / slow-mode coordinate | `NEAR` |
| Spectral graph wavelet | multiscale localized address | `NEAR` |

### Legitimate bridge use

Use these papers to justify spectral diagnostics inside HAOS-IIP:

- spectral address as a low-mode identity coordinate
- geometry bridge features based on normalized Laplacian embeddings
- transformation bottleneck diagnostics based on Fiedler / multieigenvector
  clustering
- multiscale persistence checks using spectral wavelets or diffusion embeddings

### Claim boundary

Spectral graph structure does not imply physical spacetime.
It supports an operational geometry diagnostic, not an ontology.

## 2. Hodge / cochain sector structure

### What to keep

From Hodge decomposition, Hodgelets, weighted hypergraph Hodge theory,
persistent Laplacians, and Dirac-Kahler work, HAOS-IIP needs:

- cochain spaces as data-bearing state spaces
- boundary / coboundary operators as structural differentials
- Hodge Laplacian as the natural operator on cochains
- exact / coexact / harmonic decomposition
- harmonic representatives as topology-carrying residue
- localized Hodge spectral elements for flow diagnostics
- persistent Laplacians as filtration-sensitive spectral invariants
- Dirac-Kahler language as physics-adjacent, not physics-derived

### HAOS-IIP mapping

| Literature object | HAOS-IIP bridge object | Status |
| --- | --- | --- |
| Cochain complex | branch-local operator hierarchy | `DIRECT` |
| Hodge Laplacian | cochain-Laplacian operator family | `DIRECT` |
| Exact component | potential-like / gradient-like sector | `DIRECT` |
| Coexact component | circulation-like / transverse sector | `DIRECT` |
| Harmonic component | persistent global address / residue | `DIRECT` |
| Persistent Laplacian | refinement / filtration persistence diagnostic | `DIRECT` |
| Dirac-Kahler form language | discrete field-theory analogy | `NEAR` |

### Legitimate bridge use

Use these papers to sharpen:

- spectral address as a harmonic / low-mode identity candidate
- current closure as exact/coexact balance
- phase contracts that distinguish sector structure from scalar scores
- continuity / refinement checks using persistent Laplacian language
- speculative field-theory notes, with strict analogy labels

### Claim boundary

Hodge structure supports a mathematically serious sector split.
It does not by itself derive gauge fields, fermions, or continuum physics.

## 3. Graph curvature and local-to-global geometry

### What to keep

From Ollivier-Ricci, Forman-Ricci, graph curvature flow, and random geometric
graph convergence papers, HAOS-IIP needs:

- neighborhood-overlap curvature as local geometry signal
- negative curvature / low conductance as bridge or bottleneck marker
- Forman curvature as cheap combinatorial curvature baseline
- curvature flow as a graph reweighting / smoothing diagnostic
- convergence papers as high-standard caution for smooth-limit claims

### HAOS-IIP mapping

| Literature object | HAOS-IIP bridge object | Status |
| --- | --- | --- |
| Ollivier-Ricci curvature | transport-overlap curvature diagnostic | `DIRECT` |
| Forman curvature | cheap graph-local curvature baseline | `DIRECT` |
| Curvature flow | flow-like geometry hardening probe | `DIRECT` |
| Conductance / bottleneck | Cheeger / recovery obstruction diagnostic | `DIRECT` |
| Random geometric convergence | possible continuum limit comparison | `CAUTION` |

### Legitimate bridge use

Use these papers to justify:

- curvature diagnostics in geometry bridge hardening
- topology-destroyed and degree-preserving controls
- curvature-flow probes as diagnostic amplifiers
- explicit refusal to claim smooth Ricci curvature unless random-geometric
  or convergence assumptions are actually met

### Claim boundary

Graph curvature can diagnose relational geometry.
It does not equal gravity unless a full mapping to metric, dynamics, stress
source, and observational interface exists.

## 4. Optimal transport and kernel comparison metrics

### What to keep

From Sinkhorn, computational OT, sliced Wasserstein, barycenters, and MMD,
HAOS-IIP needs:

- Wasserstein distance as mass-transport geometry
- Sinkhorn regularization as tractable approximate OT
- sliced Wasserstein as a cheaper projection-based surrogate
- barycenters as consensus / average-shape geometry
- MMD as a kernel two-sample test for feature-distribution separation
- algorithm-performance papers as numerical caution

### HAOS-IIP mapping

| Literature object | HAOS-IIP bridge object | Status |
| --- | --- | --- |
| Wasserstein distance | closure / feature distribution distance | `DIRECT` |
| Sinkhorn distance | scalable transport diagnostic | `DIRECT` |
| Sliced Wasserstein | cheap surrogate transport score | `NEAR` |
| Wasserstein barycenter | consensus closure / average relational state | `NEAR` |
| MMD | nonparametric target-control feature test | `DIRECT` |
| Algorithm benchmarks | numerical robustness guard | `CAUTION` |

### Legitimate bridge use

Use these papers to build:

- target-vs-control distribution tests
- curvature distribution comparisons
- spectral embedding distribution comparisons
- holdout-transfer metrics that do not depend on a single scalar
- diagnostic alternatives when exact OT is too expensive

### Claim boundary

Transport distances can improve measurement.
They are not physical transport unless the transported mass, units, dynamics,
and observation map are independently defined.

## 5. Discrete field-theory analogies

### What to keep

From Dirac-Kahler and geometric discretization papers, HAOS-IIP needs:

- differential forms as lattice/discrete field carriers
- Dirac-Kahler factorization as a physics-adjacent comparison point
- geometric discretization as a serious precedent for cochain-like physics
  language
- exact lattice supersymmetry as evidence that discretized operators can carry
  nontrivial field-theory structure

### HAOS-IIP mapping

| Literature object | HAOS-IIP bridge object | Status |
| --- | --- | --- |
| Differential forms | cochain sectors | `NEAR` |
| Dirac-Kahler operator | HAOS-IIP DK analogies | `NEAR` |
| Lattice fermions | speculative discrete field bridge | `ADJACENT` |
| Gauge-like redundancy | transformation bottleneck analogy | `CAUTION` |

### Legitimate bridge use

Use these papers to:

- make the HAOS-IIP Dirac-Kahler language more precise
- compare cochain-sector splitting against known lattice constructions
- constrain speculative field-theory notes
- prevent false equivalence between HAOS cochains and physical fermions

### Claim boundary

The literature makes the analogy mathematically respectable.
It does not turn HAOS-IIP into a lattice field theory.

## 6. Caution gates

These are the main failure checks extracted from the literature.

### C1. Spectral geometry is not physical geometry

Spectral embeddings can recover useful coordinates, but the coordinates depend
on graph construction, weights, kernels, and sampling.

Bridge requirement:

- declare graph construction and normalization before scoring
- test degree-preserving and topology-destroyed controls
- report holdout transfer

### C2. Hodge sector split is not a field equation

Exact / coexact / harmonic decomposition is structural. It does not supply a
Lagrangian, gauge group, coupling constants, or physical observables.

Bridge requirement:

- keep sector diagnostics separate from physical interpretation
- require independent observation maps before physics claims

### C3. Graph curvature is not gravity

Ollivier/Forman curvature can behave geometrically, but gravity requires a
specific metric-dynamics-source-observation chain.

Bridge requirement:

- use curvature as a diagnostic unless a full generative model exists
- treat random-geometric convergence papers as high bar, not as permission

### C4. Transport metrics are not physical transport

Wasserstein/Sinkhorn/MMD compare distributions. They do not define physical
movement unless the underlying mass and cost have units and meaning.

Bridge requirement:

- label transport features as diagnostic or calibrated
- do not use them as derived physical quantities without a mapping contract

### C5. Discrete field analogies are not derivations

Dirac-Kahler and lattice papers show that cochain-like structures can appear
in serious physics, but HAOS-IIP still needs its own dynamics and observation
map.

Bridge requirement:

- use DK/lattice papers as analogy anchors
- require independent generation of physical spectra, correlations, or
  observables before any physical bridge claim

## 7. Minimal bridge toolkit extracted from the corpus

If HAOS-IIP uses this literature well, the next bridge work should preserve
these components:

1. normalized graph / Hodge Laplacian definitions
2. low-mode spectral address extraction
3. exact / coexact / harmonic sector reporting
4. persistent Laplacian or refinement persistence checks
5. Fiedler / Cheeger / spectral clustering diagnostics
6. Ollivier-Ricci and Forman curvature controls
7. Sinkhorn / sliced Wasserstein / MMD distribution comparisons
8. target-control-holdout separation
9. explicit mapping statuses: `DIRECT`, `NEAR`, `ADJACENT`, `CAUTION`
10. claim gates that separate diagnostic resemblance from physical derivation

## 8. What this enables

This extraction supports three next artifacts:

1. a bridge-methods technical note
2. a geometry bridge hardening plan
3. a claim-gated physics-bridge precommitment template

It does not authorize:

- Bell derivation claims
- hydrogen derivation claims
- gravity derivation claims
- consciousness claims
- continuum-limit claims
