# TERMINOLOGY

This file standardizes symbols and names used in canonical HAOS-IIP documents.

## Canonical symbols

| Symbol | Meaning | Canonical use |
|---|---|---|
| `x_i` | node position or embedded coordinate | substrate definition |
| `epsilon_k` | interaction kernel width | Gaussian kernel bandwidth |
| `A_ij` | interaction kernel weight | weighted adjacency |
| `D` | degree matrix | `D_ii = sum_j A_ij` |
| `L0` | node Laplacian / scalar graph Laplacian | `L0 = D - A` |
| `L` | generic Laplacian symbol | only when sector is obvious |
| `mathcal{L}` | continuum operator | sector-specific continuum limit |
| `B` | node-edge incidence matrix | discrete exterior derivative on 0-forms |
| `C` | edge-face incidence matrix | discrete curl / 1-to-2 form lift |
| `W` | edge-weight diagonal matrix | weighted scalar/node sector |
| `W0` | edge-weight matrix in `L1` | 0->1 part of Hodge sector |
| `W2` | face-weight matrix in `L1` | 1->2 part of Hodge sector |
| `L1` | discrete edge / Hodge Laplacian | gauge / vector sector |
| `D_H` | discrete Dirac-type operator | fermion / spinor toy branch |
| `epsilon_c` | coherence threshold | recoverability inequality |
| `lambda_n` | eigenvalues of `L0`, `L1`, or `mathcal{L}` | spectral quantities |
| `lambda_1` | first nonzero eigenvalue | spectral gap / recoverability scale |
| `V(x)` | scalar potential term | continuum scalar operator |
| `theta_ij` | edge phase | phase-dressed kernel |
| `U_ij` | link variable `exp(i theta_ij)` | covariant transport |
| `L_theta` | covariant scalar graph Laplacian | phase-dressed node operator |
| `a_mu` | continuum connection 1-form | gauge-program notation |
| `F_mu_nu` | field-strength candidate | gauge-program notation |
| `psi` | complex field on nodes or continuum | mode / excitation field |
| `J_ij` | edge current | phase-dressed continuity diagnostics |
| `R(x)` | scalar curvature | heat-kernel geometry |
| `S_spec` | spectral action | gravity/gauge effective action |
| `Lambda` | spectral cutoff scale | spectral action |

## Canonical terms

### Recoverability

Operational meaning:

- stability under bounded perturbation over a finite horizon

Canonical inequality:

```text
exp(-lambda_1 T) delta < epsilon_c
```

### Geometry-like sector

Use when:

- low modes of `L0` are smooth, long-wavelength, and consistent with scalar Laplacian intuition

Do not use as a synonym for:

- full curvature derivation
- Einstein dynamics

### Gauge-like / connection-like

Use when:

- edge or connection data in `L1` or `L_theta` change the operator covariantly
- loop holonomy or edge-current structure is explicit

Do not use as a synonym for:

- photon
- derived Maxwell field

### Fermion-like / defect-bound

Use only when:

- localization survives controlled perturbation or defect insertion
- first-order or sign-structured operator data are present
- the structure is not explainable as an ordinary scalar boundary artifact

Do not use as a synonym for:

- electron
- spinor
- Dirac sector

## Symbol conflicts found in legacy material

### Conflict 1: `epsilon`

Legacy usage mixed:

- kernel width
- coherence threshold

Canonical resolution:

- `epsilon_k`: kernel width
- `epsilon_c`: coherence threshold

### Conflict 2: `L`

Legacy usage mixed:

- graph Laplacian
- Lagrangian-like object

Canonical resolution:

- `L0`: node Laplacian
- `L1`: edge/Hodge Laplacian
- `L`: generic Laplacian only when sector is obvious
- `mathcal{L}`: continuum operator
- `mathscr{L}`: Lagrangian density, only if needed

### Conflict 3: charge

Legacy usage mixed:

- winding number
- coupling asymmetry
- placeholder gauge quantity

Canonical resolution:

- use `charge candidate` or `topological charge candidate` unless a derivation exists

## Naming rules

- source drafts: `_v1`, `_v2`, `_v3`
- stabilized canonical files: kept in `docs/canon/`
- do not create new symbols in canonical files unless added here first

## Operator-sector naming

- `0-form`: node scalar field
- `1-form`: edge field / transport field
- `spinor field`: doubled or framed field acted on by `D_H`
