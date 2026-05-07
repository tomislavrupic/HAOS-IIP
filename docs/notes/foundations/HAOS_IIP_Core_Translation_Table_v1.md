# HAOS-IIP Core Translation Table v1

Status:

- Phase 66 artifact
- terminology-reduction layer
- no new result
- no stronger physical claim

Purpose:

Reduce terminology debt and expose HAOS-IIP claims to standard-language
critique before internal vocabulary is used as shorthand.

## Rule

Standard-language meaning first.

HAOS-IIP term second.

Claim boundary always explicit.

Core enforcement line:

```text
A HAOS-IIP term is not allowed to function as an explanation until it has been
translated into operational language.
```

## How to Use This Table

For public-facing summaries:

1. Start with the standard-language equivalent.
2. State what is operationally measured or constructed.
3. State what is numerically tested.
4. State what is not claimed.
5. Only then introduce the HAOS-IIP term as shorthand.

For papers and releases:

- each claim box should use this table to translate internal terms
- each new term should either map to an existing row or add a new row
- no term should be used as proof of itself

## Core Table

| HAOS-IIP term | Standard-language equivalent | What it means operationally | What is numerically tested | What is not claimed | Nearby established concepts | Failure condition |
| --- | --- | --- | --- | --- | --- | --- |
| Harmonic address | Stable spectral / operator identity class | A mode, branch, peak, or diagnostic label that can be reidentified across perturbation, distance, delay, or refinement | Frequency / branch stability, amplitude retention, identity matching, control separation | A universal resonance, hidden essence, or proof of physical ontology | Eigenmodes, normal modes, spectral clusters, coherent modes | The address cannot be reidentified, drifts into controls, or is indistinguishable from arbitrary labeling |
| Spectral address | Specific spectral coordinate, band, or peak window | A frequency or eigenvalue-region used as an operational target | Peak detection, drift, linewidth proxy, peak-to-background, band retention | A guaranteed physical particle, field, or fundamental frequency | Resonance peaks, spectral bands, Fourier / STFT features | The target peak is absent, unstable, broadband, or no better than a null window |
| Recoverable coherence | Persistence of diagnostic structure under bounded perturbation | A baseline structure remains measurable, returns, or degrades continuously under declared stress | Recoverability score, delta persistence, safety margin, collapse threshold | Consciousness, mystical coherence, or final physical unity | Stability analysis, robustness, hysteresis, persistence metrics | The structure collapses under mild perturbation or fails matched controls |
| Recoverability metric | Bounded score comparing baseline and perturbed diagnostic structure | A scalar or vector score summarizes how much of the declared diagnostic survives a stressor | Baseline-normalized retention, drift, broadening, closure residual, delta persistence | A universal truth score or physical law by itself | Robustness metrics, resilience scores, stability margins | The score is threshold-tuned, non-reproducible, or fails null/model controls |
| k_star / collapse index | First declared perturbation level where recoverability shows sustained collapse | A condition index is marked when the metric remains below threshold for the required sustain window | Collapse threshold, sustain steps, safety margin, visible-failure lag | A physical singularity, ontic collapse, or universal critical point | Critical thresholds, failure onset, early-warning transition markers | k_star changes under harmless ordering, threshold, seed, or control variation |
| Interaction-native operator system | Discrete substrate with operators defined directly on interactions | A graph, complex, or relation set with weighted adjacency, incidence, Laplacian, or Hodge-like operators | Operator identities, spectra, branch behavior, perturbation response | Spacetime, quantum gravity, or a physical substrate by assertion | Graph Laplacians, finite complexes, DEC, kernel operators | Results depend on hidden geometric assumptions or disappear under equivalent operator formulations |
| Frozen regime | Fixed computational contract for a phase | Substrate class, operators, telemetry, thresholds, controls, and reproduction commands are held fixed | Whether declared tests reproduce under the frozen contract | General validity outside the declared regime | Benchmark protocols, preregistered tests, fixed baselines | Rules change after results are known or reproduction commands no longer match outputs |
| Branch-local structure | Structure restricted to a selected spectral or operator branch | A result is claimed only on a defined mode window, sector, seed family, or branch | Branch identity, leakage, overlap, local stability | Global structure across the whole operator system | Spectral subspaces, invariant subspaces, local charts | The branch cannot be separated, leaks into other sectors, or relabels under perturbation |
| Active branch closure | Stable tracking of the selected physically relevant branch | A low active sector remains identifiable across refinement or perturbation | Transported overlaps, principal angles, scaled-eigen drift, purity | Full continuum physics or universal branch identity | Coexact sectors, spectral-flow tracking, subspace continuation | The active branch reshuffles, mixes, or loses identity under refinement |
| Clean-baseline identity | Branch identity in the unperturbed reference case | The same branch can be tracked before adding disorder or stabilizing perturbations | Clean-background overlaps, baseline drift, degeneracy sensitivity | That disorder-stabilized results automatically hold in the clean case | Degenerate perturbation problems, symmetry baselines | The clean baseline is mixing-prone or cannot maintain stable identity |
| Holonomy-split family | Spectral family separated by cycle / twist response | A branch responds distinguishably to imposed periodic twists or cycle phases | Flux / twist sensitivity, family separation, branch splitting | Gauge theory recovered in the physical sense | Bloch phases, flat connections, holonomy parameters, boundary twists | Splitting vanishes, is numerical artifact, or matches scalar leakage controls |
| Scalar carrier | Scalar-valued field-like data on an interaction graph | Node-level scalar values propagated or tested through graph/operator response | Scalar Laplacian behavior, Green response, recoverability, perturbation response | A real physical scalar field without bridge validation | Graph signals, scalar fields on meshes, potential functions | Scalar behavior fails convergence, controls, or stability tests |
| Scalar-carrier geometry | Metric-like behavior inferred from scalar operator response | Scalar response produces distance-like or geometry-like diagnostics | Effective dimension, Green scaling, local metric surrogate, robustness | Actual spacetime geometry or Regge geometry | Spectral geometry, diffusion geometry, graph embeddings | Metric surrogate fails convergence, positivity, isotropy, or null controls |
| Local metric field | Local distance / metric surrogate derived from operator response | A per-location or pairwise metric-like diagnostic constructed from Green/eigen/operator response | Stability under perturbation, locality, isotropy, refinement behavior | Physical metric tensor by default | Diffusion distance, resistance distance, graph metric learning | Surrogate becomes nonlocal, unstable, anisotropic without claim support, or non-convergent |
| Recoverability gradient | Change in recoverability across location, scale, or condition | A directional or radial trend in recovery score under controlled variation | Slope, monotonicity, shell behavior, condition response | Physical force, field equation, or causal gradient | Sensitivity gradients, response surfaces, stability landscapes | Gradient is random, seed-specific, sign-unstable, or control-matched |
| Current closure | Flux-like consistency under selected operator constraints | A conservation-like balance or divergence/flux relation holds in the tested regime | Relative error, shell balance, closure residual, matched controls | Electromagnetic or physical current without bridge validation | Conservation laws, discrete divergence theorem, flux balance | Closure disappears under refinement, controls, or alternate discretization |
| Shell-native current | Flux-like balance measured on graph shells rather than imported coordinates | Shell-indexed response behaves consistently under the native graph distance/response structure | Shell kappa drift, radial balance, shell residuals | Physical radial current or inverse-square law by itself | Spherical shell averages, graph geodesic shells, radial flux tests | Shell signal is binning artifact or fails non-shell controls |
| Smooth-inhomogeneity current closure | Flux-like consistency under smooth coefficient variation | Closure test remains stable when the substrate has controlled smooth inhomogeneity | Closure residuals across smooth eta branches, p90 error, radial response | General inhomogeneous-field law | Variable-coefficient elliptic operators, effective media | Smooth branch closes only by tuning or fails matched bump/disorder controls |
| Collapse detection | Detection of loss of recoverability | A recoverability metric crosses a declared collapse threshold | Collapse threshold, sustained steps, delta persistence, visible failure lag | Ontic collapse, quantum measurement collapse, or irreversible doom | Early-warning signals, stability loss, bifurcation indicators | Collapse flag depends on arbitrary threshold, evaluation order, or noise |
| Collapse-ordering invariance | Order-independent detection of recoverability loss | The sequence or presence of collapse is not an artifact of how diagnostics are evaluated | Permuted-order checks, ranking stability, collapse-time consistency | Fundamental time ordering or causal law | Robust ranking, order-invariance tests, sensitivity analysis | Collapse changes under harmless evaluation reordering |
| Localized bump response | Operator response to a controlled localized perturbation | A bump perturbation produces measured local/global response under the frozen metric | Recovery loss, locality, response threshold, control separation | Real defect physics unless externally bridged | Impulse response, localized forcing, Green response | Bump response is indistinguishable from random perturbation or smoothing controls |
| Power-law scaling boundary | Boundary where response resembles a declared power-law regime | A fitted scaling relation holds only inside a claim-gated regime | Slope, fit quality, residuals, null comparison, range sensitivity | Universal physical law from fit alone | Scaling laws, critical exponents, Green-function asymptotics | Slope changes under range, seed, null model, or alternate estimator |
| Disorder flux | Propagation of controlled disorder through the operator system | A perturbation or disorder field produces measurable structured response | Seed universality, flux stability, residuals, disorder response | Physical energy flux or entropy law | Disorder response, random media, perturbation propagation | Signal is seed-specific, control-matched, or disappears under repeated trials |
| Physics bridge observable | HAOS-style metric mapped to a physics-adjacent measurable proxy | Internal diagnostic is connected to an external or toy observable with a claim gate | Toy/proxy score, control failure, external-data parse, reproducibility | Established physics recovery unless independently validated | Model validation, observable proxies, benchmark tasks | Proxy fails controls or cannot be distinguished from generic graph behavior |
| Biology telemetry | Recoverability-style metrics applied to biological or bio-adjacent data/models | Biological structure is represented through bounded diagnostic observables | Null ladders, robustness, external dataset bridge, control comparison | Life, consciousness, or biological mechanism proof | Systems biology metrics, network robustness, time-series biomarkers | Telemetry gives no advantage over standard biological null models |
| Materials bridge | Recoverability-style metrics applied to public materials data | External materials data are downloaded/loaded, parsed, and audited under bounded metrics | Provenance, hashes, parse status, peak stability, propagation, k_star logic | Material proves HAOS-IIP or embodies the framework | Spectroscopy, pump-probe analysis, materials informatics | Data unavailable, parse fails, metrics indistinguishable from standard analysis, or no predictive gain |
| Claim-gated bridge | Bridge test with explicit allowed and disallowed claims | Every PASS is scoped to toy/proxy/internal/external status before interpretation | Claim table, controls, status JSON, allowed-language checks | Unbounded physical inference | Benchmark governance, validation protocols, preregistered claim scope | Claim language exceeds evidence class or ignores failed controls |
| Internal PASS | A successful result inside a declared HAOS-IIP frozen regime | The code ran, thresholds were met, outputs reproduced, and controls behaved as specified | Reproduction, thresholds, summaries, validation flags | External physical validation or theory proof | Unit tests, benchmark PASS, numerical validation | Output cannot be reproduced or threshold/control logic is ambiguous |
| External bridge claim | Claim that HAOS-style telemetry has been applied to external data or predicts an external observable | A public dataset or independent observable is linked to a bounded HAOS-style metric | Data provenance, parsing, null comparison, recoverability, prediction status | Proof of HAOS-IIP or superiority over established models | External validation, cross-domain benchmark, predictive test | No real data, no provenance, no null model, or no advantage over standard explanations |

## Enforcement Checklist

Before a HAOS-IIP term is used in a public summary:

- [ ] standard-language equivalent is stated first
- [ ] operational construction is specified
- [ ] numerical test is named
- [ ] non-claim boundary is stated
- [ ] nearby known concepts are acknowledged
- [ ] failure condition is explicit

## Authority Boundary

This table translates vocabulary. It does not validate the translated claims.

Translation makes a claim easier to inspect. It does not make the claim true.

## Next Work

- turn this table into a compact public appendix for the canonical entry paper
- add term-specific links to the most relevant experiment, note, or paper
- use the Phase 66 falsification engine to attach failure gates to translated terms
- add a release-review rule that blocks new terms without translation rows
