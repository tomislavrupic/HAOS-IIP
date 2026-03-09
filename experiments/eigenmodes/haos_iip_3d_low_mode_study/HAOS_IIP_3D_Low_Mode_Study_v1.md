# HAOS / IIP 3D Low-Mode Study

March 9, 2026

Files produced in this study:

- `haos_iip_3d_low_modes.py`
- `haos_iip_results.json`
- `haos_iip_mode_plots/`

Guiding question:

> What do the first few eigenmodes actually look like? For a 3D network or continuum limit with tuned epsilon, do you get massless vectors, Dirac fermions, or curvature terms? Any explicit eigenvalue calc or plot showing matches to photon/electron masses or couplings?

Short answer:

- for the bare 3D node Laplacian, the low modes are mainly scalar diffusion / geometry-like modes
- no explicit massless vector sector appears
- no fermion-like or Dirac-like sector appears
- a phase-dressed branch does show connection sensitivity and circulating complex modes, but this is still a scalar field in a background connection, not a derived gauge-field sector

## 1. Purpose

The goal was to build the first explicit 3D low-spectrum study for the HAOS/IIP kernel

$$
A_{ij} = \exp\!\left(-\frac{|x_i-x_j|^2}{2\epsilon}\right),
\qquad
L = D - A,
$$

and to test whether the first low-lying eigenmodes of a concrete 3D substrate show any structure that plausibly bridges toward:

- geometry / curvature-like sector
- gauge / vector-like sector
- fermion-like or defect-bound excitation classes

The study was numerical and falsifiable. No sector was counted as present unless it appeared explicitly in the eigenvalues, eigenvectors, or current diagnostics.

## 2. Current Bottleneck

The current HAOS/IIP notes already support the scalar chain

$$
\text{kernel} \to \text{Laplacian} \to \text{continuum scalar operator} \to \text{geometry-like spectrum}.
$$

The bottleneck is that this does not yet produce:

- an autonomous gauge-field sector
- a spinorial or fermionic sector
- physical quantum numbers or mass/coupling matches

The practical question is therefore narrower:

Do the lowest modes of a concrete 3D substrate show anything beyond scalar geometry-like behavior?

## 3. Substrate Definition

Three numerical branches were used.

### 3.1 Bare 3D perturbed cubic substrate

- substrate: open cubic lattice in `[0,1]^3`
- grid: `6 x 6 x 6`
- node count: `N = 216`
- perturbation: interior nodes perturbed by `0.06 h`, with `h = 1 / 5 = 0.2`
- boundary condition: open
- kernel truncation: edges kept for `|x_i-x_j| <= 2.5 sqrt(epsilon)`
- reference epsilon:
  - mean nearest-neighbor squared distance `~ 0.037998`
  - epsilon sweep: `0.037998`, `0.056997`, `0.075996`
  - reference run: `epsilon = 0.056997`

### 3.2 Bare 3D random geometric substrate

- substrate: random geometric graph in `[0,1]^3`
- node count: `N = 256`
- sampling: uniform random point cloud
- boundary condition: open
- kernel truncation: edges kept for `|x_i-x_j| <= 2.5 sqrt(epsilon)`
- reference epsilon:
  - mean nearest-neighbor squared distance `~ 0.010290`
  - epsilon sweep: `0.015434`, `0.020579`, `0.025724`
  - reference run: `epsilon = 0.020579`

### 3.3 Optional phase-dressed branch

This branch was used only after the bare study worked.

- substrate: regular open cubic lattice in `[0,1]^3`
- grid: `6 x 6 x 6`
- node count: `N = 216`
- neighborhood: nearest-neighbor only via cutoff `1.05 h`
- epsilon: `epsilon = h^2 = 0.04`
- link dressing:

$$
A_{ij} \mapsto A_{ij} e^{i\theta_{ij}}
$$

with a Landau-gauge-like background chosen so that the mean `xy` plaquette phase is

$$
\Phi_{xy} \approx -\frac{\pi}{3}.
$$

This is a test of connection sensitivity, not a derivation of a gauge field.

## 4. Kernel and Laplacian Construction

For each substrate, the primary operator was the ordinary graph Laplacian:

$$
A_{ij} = \exp\!\left(-\frac{|x_i-x_j|^2}{2\epsilon}\right),
\qquad
D_{ii} = \sum_j A_{ij},
\qquad
L = D - A.
$$

Raw, unnormalized eigenvalues were used as the primary diagnostic. This matters because the overall scale of the ordinary Laplacian changes strongly with `epsilon`.

For the phase-dressed branch, the covariant node operator was

$$
(L_\theta \psi)_i
= \sum_j A_{ij} \big(\psi_i - e^{i\theta_{ij}}\psi_j\big).
$$

Edge currents were measured by

$$
J_{ij} = 2 A_{ij}\,\mathrm{Im}\!\left(\psi_i^* e^{i\theta_{ij}}\psi_j\right).
$$

## 5. Numerical Setup

Computed diagnostics:

- first nontrivial eigenvalues
- first low-lying eigenvectors
- inverse participation ratio (IPR)
- participation ratio
- near-degeneracy structure
- coordinate correlation for low modes
- epsilon sweep
- limited `N` sweep
- for the phase branch: plaquette holonomy and current ratio

Plots created:

- `low_spectrum_comparison.png`
- `cubic_open_perturbed_modes.png`
- `random_geometric_open_modes.png`
- `phase_dressed_lattice_modes.png`
- `localization_vs_mode_index.png`
- `epsilon_sweep.png`
- `n_sweep_spectral_gap.png`

## 6. Low-Spectrum Results

### 6.1 Bare cubic substrate

Reference parameters:

- `N = 216`
- `epsilon = 0.056997`
- cutoff radius `~ 0.597`
- connected components: `1`

Bottom of the spectrum:

| mode | eigenvalue | IPR | max coord corr | classification |
|---|---:|---:|---:|---|
| 1 | 2.410405 | 0.013429 | 0.7661 | geometry-like coordinate mode |
| 2 | 2.415782 | 0.014594 | 0.6308 | mixed scalar mode |
| 3 | 2.422536 | 0.012982 | 0.7193 | mixed scalar mode |
| 4 | 3.833708 | 0.030437 | 0.0072 | mixed scalar mode |
| 5 | 3.850625 | 0.033272 | 0.0067 | mixed scalar mode |
| 6 | 3.859648 | 0.033938 | 0.0054 | mixed scalar mode |

Observed structure:

- the first three nontrivial modes form a near-degenerate family
- the first family is smooth and strongly aligned with coordinate directions of the cube
- the visual plots show one-axis nodal variation, consistent with the first scalar box modes

This is the cleanest geometry-like result in the study.

### 6.2 Bare random geometric substrate

Reference parameters:

- `N = 256`
- `epsilon = 0.020579`
- cutoff radius `~ 0.359`
- connected components: `1`

Bottom of the spectrum:

| mode | eigenvalue | IPR | max coord corr | classification |
|---|---:|---:|---:|---|
| 1 | 0.461208 | 0.090158 | 0.5959 | localized candidate |
| 2 | 0.496521 | 0.010438 | 0.7230 | mixed scalar mode |
| 3 | 0.519401 | 0.062834 | 0.5309 | localized candidate |
| 4 | 0.745402 | 0.018514 | 0.6876 | mixed scalar mode |
| 5 | 0.800106 | 0.045599 | 0.5243 | localized candidate |
| 6 | 0.994653 | 0.037871 | 0.0869 | localized candidate |

Observed structure:

- some of the lowest modes are not clean long-wavelength box modes
- several low modes show strong boundary/corner concentration in the 3D plots
- the irregular point cloud plus ordinary unnormalized Laplacian produces low-mode localization before any credible non-scalar sector appears

Interpretation:

- this is not evidence for a fermion-like sector
- it is more plausibly a boundary / degree-irregularity effect of the open random geometric graph

### 6.3 Phase-dressed lattice branch

Reference parameters:

- regular `6 x 6 x 6` open lattice
- nearest-neighbor branch
- `epsilon = 0.04`
- mean `xy` plaquette phase `= -1.047198` with negligible spread

Bare nearest-neighbor cubic spectrum:

$$
0,\ 0.162519,\ 0.162519,\ 0.162519,\ 0.325039,\ 0.325039,\ 0.325039,\dots
$$

Phase-dressed spectrum:

$$
0.238821,\ 0.253441,\ 0.289311,\ 0.326528,\ 0.401341,\ 0.415960,\dots
$$

First low modes:

| mode | eigenvalue | IPR | current ratio |
|---|---:|---:|---:|
| 1 | 0.238821 | 0.006544 | 0.4300 |
| 2 | 0.253441 | 0.007862 | 0.2528 |
| 3 | 0.289311 | 0.005898 | 0.4582 |
| 4 | 0.326528 | 0.009685 | 0.3429 |
| 5 | 0.401341 | 0.009815 | 0.4806 |
| 6 | 0.415960 | 0.011793 | 0.2845 |

Observed structure:

- the exact bare degeneracies are split by the phase dressing
- the low modes become genuinely complex
- current ratios are substantial, around `0.25` to `0.48`
- the phase plots show organized phase winding / circulation patterns

Interpretation:

- this is explicit connection sensitivity
- it is not yet a vector boson sector
- these are still scalar modes in a background connection

## 7. Mode Classification

### Geometry / curvature-like sector

What appeared:

- on the perturbed cubic substrate, the first nontrivial family is smooth, near-degenerate, and coordinate-aligned
- this is exactly the structure expected from a scalar Laplacian on a 3D box

What did not appear:

- no explicit curvature term was extracted from these flat 3D substrates
- there was no mode family that required interpretation beyond scalar geometry-like behavior

Verdict:

- geometry-like scalar sector: yes
- explicit curvature sector: not yet

### Gauge / vector-like sector

What appeared:

- under phase dressing, low modes acquire nontrivial phase structure and circulating currents
- phase dressing lifts low degeneracies and changes the bottom of the spectrum in a controlled way

What did not appear:

- no autonomous transverse vector field was diagonalized
- no separate gauge-field mode family emerged from the bare node Laplacian
- no Maxwell-like operator or edge/Hodge sector was computed here

Verdict:

- connection-like sensitivity: yes
- credible massless vector sector: no

### Fermion-like / defect-like sector

What appeared:

- some low modes on the open random geometric graph are localized

What did not appear:

- no chiral doubling
- no spinorial pairing
- no Dirac-like linear branch
- no defect-bound family with robust topological quantum numbers

Verdict:

- localized excitations: yes, but presently boundary/irregularity dominated
- fermion-like sector: no clear analogue yet

## 8. Continuum Interpretation

The cleanest continuum bridge is the perturbed cubic substrate.

Its first low modes behave like the first scalar modes of a Laplace operator on a bounded 3D region:

- smooth
- long-wavelength
- near-degenerate by axis family
- strongly coordinate-aligned

This is consistent with the claim that the HAOS/IIP Gaussian kernel can reproduce a scalar geometry-like sector in a concrete 3D substrate.

The random geometric graph is less clean:

- low modes mix genuine smooth partitions with boundary/density artifacts
- the ordinary unnormalized Laplacian is sensitive to degree irregularity
- therefore this branch does not yet provide clean evidence for anything beyond scalar diffusion

The phase-dressed branch shows something real but limited:

- the dressed operator is sensitive to background holonomy
- low modes acquire circulation and complex phase structure
- however this is still a scalar connection-Laplacian test, not a derivation of a dynamical gauge field

Direct answer to the continuum question:

- the lowest bare modes do resemble ordinary scalar Laplacian modes on the clean cubic substrate
- there is no explicit sign here of an additional vector or fermion sector emerging from the bare 3D node Laplacian
- `epsilon` primarily changes the scale and mixing of the scalar family; it does not, by itself, create a new particle or force sector in this study

## 9. What Is Still Missing for Actual Particle / Force Claims

This study does **not** show:

- photon masses or couplings
- electron masses or couplings
- Dirac fermions
- spin structure
- chiral sectors
- a true gauge-boson spectrum

Reasons:

- the primary operator is a node Laplacian, which is a scalar operator
- no edge/Hodge operator was diagonalized
- no Dirac-type operator was constructed
- no controlled topological defect background was inserted
- no continuum renormalization was done to compare raw eigenvalues across `epsilon`

So the honest conclusion is:

- the bare HAOS/IIP 3D kernel currently shows a scalar geometry-like sector
- the phase-dressed extension shows background connection sensitivity
- actual particle / force claims still require new operator structure

## 10. Next-Step Program

The next steps should be narrower and harder, not broader.

1. Add a normalized-Laplacian comparison on the random geometric graph

- separate genuine low geometry from boundary/degree artifacts

2. Move from node operators to edge operators

- Hodge 1-Laplacian
- discrete curl / divergence decomposition
- explicit test for transverse low modes

3. Build a Dirac-type branch

- bipartite or chiral lattice
- discrete Dirac operator instead of only `D-A`
- test for linear low-lying branches and chiral doubling

4. Insert controlled defects

- vacancy defects
- phase vortices
- dislocation-like graph perturbations
- test whether localized modes persist under perturbation

5. Keep the phase-dressed branch, but do not oversell it

- it is useful for testing holonomy and current
- it is not yet a gauge-field derivation

6. Add a weakly curved 3D substrate

- if curvature-like terms are the goal, use a genuinely weakly curved point cloud rather than only flat boxes

## Useful Source Fragments

The local source search was broad, but only a small part was mathematically useful for this task:

- `HAOS_IIP kernel.docx`: Gaussian kernel and spectral-geometry route
- `HAOS_IIP Spectral Geometry Note.docx`: scalar geometry chain
- `Interaction-invariant Physics Full.docx`: phase structure, global phase gauge, and the observation that first-derivative terms need a background field
- `Rebuilding Duda’s picture in HAOS_IIP.docx`: topological charge and localized excitation placeholders
- `HAOS-IIP Emergent Gauge Sector Note.md`: covariant graph Laplacian and holonomy language

Most broader manifesto-style or metaphysical text was ignored.

## Current verdict

> **Current verdict**
>
> - what worked: the 3D cubic substrate produced explicit smooth low modes consistent with a scalar geometry-like sector, and the phase-dressed branch produced clear holonomy-sensitive circulating modes
> - what did not appear: no credible massless vector sector, no Dirac/fermion sector, and no observation-level mass/coupling match
> - what must be tried next: edge/Hodge operators, Dirac-type operators, and controlled defect backgrounds
