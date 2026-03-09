# GAUGE_PROGRAM

## Purpose

The gauge program asks whether the same weighted substrate that supports `L0` can be lifted to an edge/Hodge branch `L1` that supports genuine vector-sector structure.

## Minimal proposal

`[P]` Start from the same kernel-derived weighted geometry, then move from node fields to edge fields.

If `B` is the node-edge incidence matrix and `C` the edge-face incidence matrix, then the target operator is

$$
L_1 = B W_0 B^T + C^T W_2 C.
$$

Before explicit faces are built, the rough starter branch is

$$
L_1 \approx B B^T.
$$

## Phase-dressed scalar bridge

`[P]` The phase-dressed node operator remains useful as a bridge:

$$
A_{ij} \mapsto A_{ij}e^{i\theta_{ij}},
\qquad
U_{ij}=e^{i\theta_{ij}},
$$

with

$$
(L_\theta \psi)_i = \sum_j A_{ij}(\psi_i - U_{ij}\psi_j).
$$

This tests connection sensitivity, but it is still scalar.

## Local gauge law

`[P]`

$$
\psi_i \mapsto e^{i\alpha_i}\psi_i,
\qquad
U_{ij} \mapsto e^{i\alpha_i}U_{ij}e^{-i\alpha_j}.
$$

## Diagnostics

`[P]` The current practical diagnostics are:

- loop holonomy:

$$
W(C)=\prod_{(ij)\in C} U_{ij}
$$

- edge current:

$$
J_{ij}=2A_{ij}\,\mathrm{Im}(\psi_i^*U_{ij}\psi_j)
$$

- eventual decomposition of `L1` low modes into gradient-like vs transverse-like sectors

## Current evidence

`[E]` The current 3D phase-dressed lattice experiment shows:

- degeneracy splitting under nontrivial plaquette phase
- complex low modes with organized circulation
- substantial edge-current ratios in the lowest modes

`[E]` This is explicit connection sensitivity.

`[E]` The periodic/twisted `L1` experiment now shows:

- open-box low `L1` modes collapse into exact gradients
- periodic low `L1` modes include a mostly coexact or divergence-free family
- trivial-flux torus cycles lift under nontrivial plaquette holonomy while remaining mostly coexact
- the same qualitative behavior appears for the tested sizes `n = 4, 5`

`[E]` This clears the main falsification test for the current edge operator: the low periodic branch is not reducible to scalar-gradient leakage.

`[E]` The wider periodic flux scan now shows:

- across `m = 0, 1, 2, 3, 4`, the low periodic branch stays mostly coexact for the tested sizes `n = 4, 5`
- the branch remains clearly flux-responsive
- the `n = 5` flow is smoother than the `n = 4` flow
- the branch still does not organize into one clean smooth transverse band across the full scan

## What is still missing

`[O]` The current repository does not yet derive:

- a dynamical gauge field from the kernel alone
- an autonomous massless vector sector
- a Maxwell-like action term from the same numerics
- a clean propagating transverse band distinct from harmonic or topological edge modes

## Immediate next tests

`[O]`

1. widen the flux scan on the periodic/twisted `L1` branch
2. compare harmonic-cycle modes to puncture or defect-supported circulation modes
3. compare coexact-projected low modes across larger lattice sizes
4. continuum-limit checks for `L_theta` and `L1`
5. test whether any low family remains robustly transverse-like away from the torus-cycle sector
