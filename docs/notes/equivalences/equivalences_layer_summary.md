# Equivalences Layer Summary

This page summarizes the current equivalence notes inside HAOS-IIP.

## Scope

The equivalence layer records structural correspondences inside frozen discrete operator regimes.
It does not claim physical ontology, continuum derivation, or external physics validation.

## Current notes

1. **Spectral Address and Laplacian Eigenmodes**
   - spectral address behaves like a recoverable eigenmode identity coordinate
   - Hodge-language adjacency: spectral / harmonic identity in the Laplacian spectrum

2. **Current Closure and Discrete Conservation Balance**
   - current closure behaves like a bounded flux-like or divergence-like balance
   - Hodge-language adjacency: exact / coexact bookkeeping under operator constraints

3. **Hodge Decomposition and Discrete Sector Splitting**
   - the repo’s exact / coexact / harmonic sector split tracks discrete Hodge decomposition
   - Hodge-language adjacency: sector decomposition on cochain complexes

4. **Hodge Laplacian Extensions and Higher-Order Sector Recovery**
   - higher-order cochains and clique-complex Laplacians extend the sector split to k-simplices
   - Hodge-language adjacency: up / down Laplacians on higher-order discrete signals

5. **Discrete Curvature Models and Local-to-Global Diagnostics**
   - local curvature acts as a companion diagnostic for bottlenecks, bridges, and hubs
   - Hodge-language adjacency: local-to-global structure complementing spectral gap and sector split

6. **Wasserstein Distance as a Transport-Cost Diagnostic**
   - minimal transport cost acts as the primitive for comparing frozen neighborhood measures
   - Hodge-language adjacency: cost-based transport geometry underneath Ollivier-style probes

7. **Ollivier-Ricci Curvature as a Local Transport Diagnostic**
   - neighborhood transport overlap acts as a bridge-versus-community diagnostic
   - Hodge-language adjacency: local transport geometry complementing spectral and Cheeger probes

8. **Discrete Ricci Flow on Graphs as a Diagnostic Probe**
   - curvature-weighted reweighting acts as a short synthetic evolution pass
   - Hodge-language adjacency: local diagnostic evolution on discrete graph geometry

## How the pieces fit

In the current frozen language, the three notes form one layered picture:

- **spectral address** names recoverable identity in the harmonic / eigenmode neighborhood;
- **current closure** names flux-like balance in the exact / coexact neighborhood;
- **Hodge decomposition** names the decomposition rule that organizes those neighborhoods on discrete cochain complexes.

That is a structural relation, not a physics derivation. Curvature and Ricci-flow-like updates now sit beside that layered picture as local diagnostics for why the global structure may be hard to recover.

## Adding new notes

Add a new equivalence note only if all of the following are true:

1. the comparison class is already visible in frozen repo language
2. the note can state an explicit claim ceiling
3. the note has evidence gates or failure conditions
4. the note stays inside discrete/operator regimes
5. the note does not reopen any frozen physics bridge

If a candidate is only an analogy, say so.

## Index rule

The folder index should list every note and keep the scope statement visible.
