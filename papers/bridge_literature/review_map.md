# HAOS-IIP Bridge Review Map

This is the reading guide for the bridge corpus.

Use it to answer three questions for each paper:

1. What does the paper actually contribute?
2. Where does it touch HAOS-IIP?
3. What kind of bridge use is legitimate?

The point is not to force a universal story. The point is to find the real
connection, if there is one, and keep the claim boundary honest.

## Legend

- `DIRECT` - structurally close to an existing HAOS-IIP operator or diagnostic.
- `NEAR` - useful mathematical neighbor, likely to sharpen a bridge note or test.
- `ADJACENT` - informative background, but not a direct bridge ingredient.
- `CAUTION` - valuable mostly as a warning that a plausible-looking mapping can fail.

## 1. Spectral clustering and graph Laplacians

| Paper | Why it matters | HAOS-IIP connection | Best use |
| --- | --- | --- | --- |
| von Luxburg, *A Tutorial on Spectral Clustering* | Canonical spectral clustering / normalized cuts / Fiedler vector reference | Directly supports HAOS spectral address, Fiedler diagnostics, Cheeger-style reasoning | `DIRECT` for spectral address and geometry-bridge notes |
| Nadler et al., *Diffusion Maps, Spectral Clustering and Eigenfunctions of Fokker-Planck operators* | Turns Laplacian eigenvectors into diffusion geometry / reaction coordinates | Good neighbor for the spectral-address and geometric-embedding story | `DIRECT` for diffusion-geometry interpretation |
| Nadler et al., *Diffusion maps, spectral clustering and reaction coordinates of dynamical systems* | Connects graph spectra to slow variables and reaction coordinates | Useful if HAOS-IIP is read as a dynamics-to-embedding pipeline | `NEAR` for transition / dynamics prose |
| Coifman and Hirn, *Diffusion maps for changing data* | Embeddings that move with parameters / data changes | Useful for holdout-transfer and change-tracking language | `NEAR` for refinement / parameter tracking |
| Gallier, *Notes on Elementary Spectral Graph Theory* | Clean normalized-cuts exposition with proofs | Good explanatory companion for internal docs | `ADJACENT` |
| Hammond, Vandergheynst, Gribonval, *Wavelets on Graphs via Spectral Graph Theory* | Spectral graph wavelets / multiscale localization | Helps bridge spectral address to multiscale diagnostics | `DIRECT` for multiscale spectral probes |
| Kuncheva and Montana, *Multi-scale Community Detection in Temporal Networks Using Spectral Graph Wavelets* | Multi-scale clustering on temporal graphs | Useful where HAOS-IIP cares about persistent structure across scales | `NEAR` for temporal / multiscale review |

## 2. Hodge decomposition, simplicial complexes, and higher-order spectra

| Paper | Why it matters | HAOS-IIP connection | Best use |
| --- | --- | --- | --- |
| Black and Nayyeri, *Hodge Decomposition and General Laplacian Solvers for Embedded Simplicial Complexes* | Higher-order Hodge decomposition and solvers | Strongly matches cochain-Laplacian language in HAOS-IIP | `DIRECT` for Hodge decomposition notes |
| Roddenberry et al., *Hodgelets* | Localized spectral representations on simplicial complexes | Natural neighbor for localized flows and sector-separated diagnostics | `DIRECT` for higher-order spectral diagnostics |
| Ren, Wu, Wu, *Hodge Decompositions for Weighted Hypergraphs* | Extends Hodge decomposition to weighted hypergraphs | Helps if the bridge is rewritten in hypergraph or weighted-complex terms | `NEAR` |
| Wei and Wei, *Persistent Topological Laplacians -- a Survey* | Persistence plus Laplacian spectra | Good for the “frozen structure over filtration” language | `DIRECT` for persistence / multiscale topology |
| Liu, Li, Wu, *The algebraic stability for persistent Laplacians* | Stability of persistent Laplacians | Useful when the bridge needs a stability guarantee rather than just a decomposition | `NEAR` |
| Farsi et al., *Wavelets and graph C*-algebras* | Graph wavelets / higher-rank graph structure | Bridge-adjacent, useful mostly as a conceptual cousin | `ADJACENT` |
| Catterall, *Dirac-Kähler fermions and exact lattice supersymmetry* | Lattice discretization with cochain-like fermions | One of the strongest physics-adjacent matches to cochain language | `DIRECT` for discrete field-theory analogy |
| Kruglov, *Dirac-KÄhler Equation* and *Review* | Dirac-Kähler formalism, symmetry, gauge language | Good for the “forms as spinors” analogy and lattice geometry discussion | `NEAR` |
| de Beauce, Sen, Sexton, *Chiral Dirac fermions on the lattice using Geometric Discretisation* | Discrete forms and lattice fermions | Useful if the bridge wants a more physics-specific lattice discretization note | `NEAR` |

## 3. Discrete curvature and curvature flow on graphs

| Paper | Why it matters | HAOS-IIP connection | Best use |
| --- | --- | --- | --- |
| Jost and Liu, *Ollivier's Ricci curvature, local clustering and curvature dimension inequalities on graphs* | Core graph curvature / neighborhood-overlap reference | Direct match to the curvature diagnostics and transport bottlenecks used in the geometry bridge | `DIRECT` |
| Sreejith et al., *Forman curvature for complex networks* | Fast combinatorial curvature | Good lightweight curvature baseline for HAOS-style graph diagnostics | `DIRECT` |
| Weber, Saucan, Jost, *Characterizing Complex Networks with Forman-Ricci Curvature and Associated Geometric Flows* | Curvature flow as network analysis | Useful for the flow-like update story in the bridge layer | `DIRECT` |
| van der Hoorn et al., *Ollivier curvature convergence in random geometric graphs* | Continuum convergence from random graphs | Strong cautionary bridge if someone starts talking about smooth limits too quickly | `CAUTION` |
| van der Hoorn et al., *Ollivier curvature of random geometric graphs converges to Ricci curvature of their Riemannian manifolds* | Stronger convergence link | Important if the bridge ever needs a “why curvature might matter” mathematical anchor | `CAUTION` / `NEAR` |

## 4. Optimal transport, Sinkhorn, and sliced Wasserstein

| Paper | Why it matters | HAOS-IIP connection | Best use |
| --- | --- | --- | --- |
| Cuturi, *Sinkhorn Distances* | Canonical entropic OT algorithm | Useful computationally for graph-neighborhood transport and curvature approximations | `DIRECT` |
| Altschuler et al., *Near-linear time approximation algorithms for optimal transport via Sinkhorn iteration* | Complexity and approximation guarantees | Good when the bridge needs tractable OT instead of exact transport | `NEAR` |
| Di Marino and Gerolin, *Optimal Transport losses and Sinkhorn algorithm with general convex regularization* | Generalized Sinkhorn / regularized OT | Helpful for robustness and alternative regularizations | `NEAR` |
| Dong et al., *A Study of Performance of Optimal Transport* | Practical algorithm comparison | Good reality check if a bridge note starts assuming Sinkhorn is always the right numerical answer | `CAUTION` |
| Wu et al., *Sliced Wasserstein Generative Models* | Sliced Wasserstein approximation | Useful if the bridge wants a cheaper transport surrogate | `NEAR` |
| Peyré and Cuturi, *Computational Optimal Transport* | Canonical OT reference | Best broad reference for the transport side of the bridge | `DIRECT` |
| Caicedo Torres et al., *A Survey on Optimal Transport for Machine Learning* | OT survey | Broad positioning reference | `ADJACENT` |
| Montesuma et al., *Recent Advances in Optimal Transport for Machine Learning* | Recent survey | Helpful as a modern overview, not a core bridge paper | `ADJACENT` |
| Kim and Pass, *Wasserstein Barycenters over Riemannian manifolds* | OT on manifolds | Strong geometric OT context if the bridge starts talking about barycenters or consensus geometry | `NEAR` |
| Claici et al., *Stochastic Wasserstein Barycenters* | Practical barycenters | Useful if one wants a stochastic approximation path | `ADJACENT` |
| Anderes et al., *Discrete Wasserstein Barycenters* | Discrete barycenter theory | Good if the data are finite and discrete, not continuous | `ADJACENT` |

## 5. Kernel distances and two-sample testing

| Paper | Why it matters | HAOS-IIP connection | Best use |
| --- | --- | --- | --- |
| Gretton et al., *A Kernel Two-Sample Test* | Canonical MMD reference | Good for comparing bridge-distributions on spectral / curvature features | `DIRECT` |

## 6. Physics-adjacent discrete geometry and lattice analogies

| Paper | Why it matters | HAOS-IIP connection | Best use |
| --- | --- | --- | --- |
| Catterall, *Dirac-Kähler fermions and exact lattice supersymmetry* | Lattice field theory from cochain-like structures | The cleanest physics-adjacent discrete-geometry bridge in the pack | `DIRECT` |
| Kruglov, *Dirac-KÄhler Equation* / *Review* | Symmetry and gauge language with forms | Useful if the bridge needs a careful statement of what the analogy is and is not | `NEAR` |
| de Beauce et al., *Chiral Dirac fermions on the lattice using Geometric Discretisation* | Geometric discretization of lattice fermions | A good bridge note if you want to compare HAOS discretization language to lattice methods | `NEAR` |
| Lederman and Toader, *On Manifold Learning in Plato's Cave* | Warning about geometry mismatch | Important cautionary reference for the bridge project itself | `CAUTION` |

## 7. Repository-local bridge papers

| Paper | Why it matters | HAOS-IIP connection | Best use |
| --- | --- | --- | --- |
| `14.1` Emergent Relational Geometry Probe | Early relational geometry diagnostics | Good starting point for the geometry bridge storyline | `DIRECT` |
| `18.1` Limits of Frozen Pre-Geometry | Clear boundary on what frozen operators do not give | The main “do not overclaim” anchor | `DIRECT` |
| `46.1` Master Synthesis | Program-level compression of the paper spine | Useful to orient the bridge corpus inside the archive | `NEAR` |
| `47.2` Cluster-Scale Dependence of Coupled Transport-Geometry Signatures | Transport / geometry coupling | Good bridge note for spectral transport and curvature adjacency | `DIRECT` |
| `48.1` Scalar Continuum Control... | Continuum-bridge caution | Useful for the “what we do not yet have” lane | `DIRECT` |
| `49.1` From Recoverable Coherence to Harmonic Operator Structure | Derivation roadmap with boundaries | Helpful as a scoped philosophy document for bridge language | `DIRECT` |
| `53.0` Physics Bridge Observables | Bridge observables release | Good external-facing gateway, but still claim-gated | `DIRECT` |
| `63.1` Celestial Boundary Toy Probes | Toy probes and boundary language | Useful as a cautionary example of bounded probes | `NEAR` |
| `64.1` / `65.1` Materials bridge sidecars | External data sidecars | Good as examples of claim-bounded external application | `NEAR` |
| `66.4` / `66.5` | Canonical entry and legitimacy audit | Core protocol references for the review discipline | `DIRECT` |
| `67.1` HBP recovery benchmarks | The most recent cross-domain bridge benchmark | Good for comparing how bridge claims are operationalized | `DIRECT` |

## 8. Practical review order

If you want the most useful reading sequence for bridge work, do this:

1. `von Luxburg` + `Nadler et al.` for spectral / diffusion geometry
2. `Black and Nayyeri` + `Hodgelets` for Hodge structure
3. `Jost and Liu` + `Forman curvature` + `Weber et al.` for curvature and flow
4. `Cuturi` + `Peyré and Cuturi` for transport
5. `Gretton et al.` for kernel-distance comparisons
6. `Catterall` + `Lederman/Toader` for the physics-adjacent cautionary layer
7. then the repo-local papers in `66.4`, `66.5`, `67.1`, and `18.1` order

## 9. How to use the pack

- If you are writing a new bridge note, start from the `DIRECT` papers.
- If you are trying to keep the bridge honest, read the `CAUTION` papers.
- If you want the bridge to stay connected to the repository, read the
  repo-local papers last.

## 10. Boundary

This map is still only a review guide.
It helps you see where the connection with HAOS-IIP is plausible, bounded, or
fragile.
It does not itself establish a physical bridge.
