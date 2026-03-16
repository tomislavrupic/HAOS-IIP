# HAOS_IIP_CORE_THEORY

Canonical boundary:

- the frozen HAOS-IIP core is defined in [HAOS_IIP_CORE_AXIOMS.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/HAOS_IIP_CORE_AXIOMS.md)
- this file is an operator-theoretic realization layer beneath that core
- kernel choices, continuum bridges, and sector lifts in this file are model statements, not the axioms themselves

Status labels:

- `[E]` established in the current notes or explicit experiments
- `[P]` plausible extension consistent with the current notes
- `[O]` open

## 1. Recoverability axiom

`[E]` A structure is viable only if perturbations decay fast enough to restore coherence within a finite horizon.

Minimal operator form:

$$
x(t) = e^{-L_0 t}x(0), \qquad e^{-\lambda_1 T}\delta < \epsilon_c.
$$

`[E]` In the current operator realization, the low end of the interaction spectrum carries the recoverability burden.

For

$$
L_0 = D-A,
\qquad
0=\lambda_0 \le \lambda_1 \le \lambda_2 \le \dots,
$$

the primary recoverability diagnostic is the first nonzero eigenvalue `lambda_1`, with the low-mode band used as a refinement when one gap is not enough.

Current reading:

- if `lambda_1` is too small, perturbations linger, split clusters, or fail to re-enter the viable coherence band
- if `lambda_1` is large enough, bounded perturbations decay below the coherence threshold within the allowed horizon

This is the cleanest mathematical compression of the HAOS recoverability test on the scalar interaction branch:

$$
\text{recoverable coherence} \iff e^{-\lambda_1 T}\delta < \epsilon_c.
$$

`[P]` A useful state-dependent coherence functional on this branch is the normalized Rayleigh diagnostic

$$
C_L(x) = 1 - \frac{\langle x, L_0 x \rangle}{\lambda_{\max}(L_0)\,\|x\|^2},
\qquad
0 \le C_L(x) \le 1.
$$

Low-spectrum states satisfy `C_L(x) -> 1`, while high-frequency or noisy states drive `C_L(x) -> 0`.

This gives a viable band

$$
x \in \mathcal V_{\epsilon_c}
\iff
C_L(x) \ge \epsilon_c,
$$

with recoverability interpreted as dynamics that return perturbed states into `mathcal V_{epsilon_c}` inside the horizon `T`.

## 2. Interaction kernel

`[E]` The baseline kernel is Gaussian:

$$
A_{ij} = \exp\!\left(-\frac{|x_i-x_j|^2}{2\epsilon_k}\right), \qquad A_{ii}=0.
$$

`[E]` The kernel is shared across sectors. The mature reading is now:

```text
same substrate
same kernel
different operator lifts
```

## 3. Scalar / geometry operator

`[E]` The scalar sector acts on node `0-forms`.

Define

$$
D_{ii} = \sum_j A_{ij}, \qquad L_0 = D - A.
$$

`[E]` Equivalently, with node-edge incidence matrix `B` and edge-weight matrix `W`,

$$
L_0 = B^T W B.
$$

The spectral gap `lambda_1` of `L0` is the primary recoverability diagnostic in the current program.

`[P]` More cautiously: the gap is the minimal summary statistic, while the low-spectrum band carries the more faithful picture whenever several near-zero modes compete.

## 4. Continuum limit of `L0`

`[P]` Under dense sampling and standard normalization, the node Laplacian tends toward a scalar continuum operator

$$
L_0 \longrightarrow \mathcal{L}_0 = -\Delta_g + V(x).
$$

`[E]` This is the current mathematical bridge used by the geometry notes.

## 5. Spectral geometry

`[P]` The heat-kernel expansion of `mathcal{L}_0` contains curvature information:

$$
K(t,x,x) \sim (4\pi t)^{-d/2}\left(1 + \frac{t}{6}R(x) + O(t^2)\right).
$$

`[E]` The program already treats geometry as the spectral shadow of recoverable interaction.

## 6. Emergent gravity

`[P]` Built from the scalar continuum operator, the spectral action

$$
S_{\mathrm{spec}} = \mathrm{Tr}\,f(\mathcal{L}_0/\Lambda^2)
$$

is expected to reproduce Einstein-Hilbert-type terms in the low-energy expansion.

`[E]` The current notes have a coherent route to this statement.

`[O]` Full dynamical recovery of GR on a 3+1D substrate is not yet established.

## 7. Gauge / vector operator target

`[P]` Gauge structure should not be read off the node Laplacian alone. The next serious sector acts on edge `1-forms`.

With node-edge incidence matrix `B`, and eventually an edge-face incidence matrix `C`, the discrete Hodge candidate is

$$
L_1 = B W_0 B^T + C^T W_2 C.
$$

`[P]` Before a face structure is implemented, the roughest starter form is

$$
L_1 \approx B B^T.
$$

This is the correct location for:

- circulation
- flux
- transverse/longitudinal decomposition
- Maxwell-like candidate structure

## 8. Phase-dressed scalar bridge

`[P]` A useful bridge to the gauge program is obtained by dressing edges in the scalar sector:

$$
A_{ij} \mapsto A_{ij} e^{i\theta_{ij}}, \qquad U_{ij}=e^{i\theta_{ij}}.
$$

The covariant scalar Laplacian is

$$
(L_\theta \psi)_i = \sum_j A_{ij}(\psi_i - U_{ij}\psi_j).
$$

`[P]` This is useful, but it is still a scalar field in a background connection unless lifted into the edge/Hodge sector.

## 9. Fermion / spinor operator target

`[P]` A fermion sector requires a first-order operator, not only a Laplacian.

The minimal discrete Dirac-type toy architecture is

$$
D_H =
\begin{bmatrix}
0 & B^T \\
B & 0
\end{bmatrix},
$$

or in weighted / phase-dressed form

$$
D_H =
\begin{bmatrix}
0 & B^T U^T \\
U B & 0
\end{bmatrix}.
$$

This is not yet a physical Dirac operator, but it is the right algebraic direction for:

- first-order propagation
- sign-symmetric spectrum
- chiral or doubled structure
- defect-bound candidate fermion modes

## 10. Open problems

`[O]` The main unresolved steps are:

- prove the covariant continuum limit
- build and test `L1` on the current weighted substrates
- identify whether a genuine vector sector exists beyond scalar matter in a background connection
- build and test a minimal `D_H` branch
- derive or reject any particle-like sector with controlled operator structure

Current minimal reading:

```text
kernel -> L0  -> scalar / geometry sector is credible
kernel -> L1  -> gauge / vector sector is the next target
kernel -> D_H -> fermion / spinor sector is still open
```
