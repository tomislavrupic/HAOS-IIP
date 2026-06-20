# Betti Diagram and HAOS-IIP Connection

## Status

`GEO-MM-01` treats the Betti diagram as an arithmetic/topological diagnostic
sidecar, not as a proof of Monstrous Moonshine and not as a physical bridge.

The stored SVG references a Lean-certified Betti diagram. HAOS-IIP does not
upgrade that claim unless the Lean source is added to the repository and checked
locally.

## What the Betti Diagram Records

A Betti diagram summarizes homological features of a graph, simplicial complex,
or related algebraic object. In the simplest graph-level case:

- `Betti_0` records connected components.
- Higher Betti numbers record cycles, void-like structure, or higher-order
  topological features when the chosen complex supports them.

For this sidecar, the immediate HAOS-compatible target is `Betti_0`, because it
has a direct graph interpretation: component count.

## HAOS-IIP Interpretation

In HAOS-IIP language, a Betti diagram can serve as a closure signature:

```text
arithmetic constellation
-> declared relation / threshold graph
-> cochain or homology calculation
-> Betti vector / diagram
-> perturbation test
-> recoverable topological signature
```

This is useful because the signature is not a raw visual pattern. It is a
compressed invariant that can be tested under perturbation.

## Symmetry as Perturbation

For the Gaussian-prime and supersingular-prime framing, the natural first
perturbation family is the square-lattice symmetry group:

- rotations,
- reflections,
- coordinate sign changes,
- coordinate swaps.

Together these are commonly represented by the dihedral symmetry `D4` of the
square.

The operational claim to test is modest:

```text
valid D4 perturbation
-> same graph connectivity structure up to isomorphism
-> same Betti_0
```

This matches HAOS-IIP's recoverability principle:

```text
structure survives allowed interaction
-> invariant recovers
-> closure signature remains stable
```

## First Lean Target

The cleanest theorem target is not a full cochain-Laplacian result. It is the
graph-native component-count statement:

```lean
theorem d4_preserves_component_count
  (g : GaussianPrimeGraph R)
  (s : D4Symmetry) :
  componentCount (applyD4 s g) = componentCount g
```

This can later be connected to Laplacian language if the local Lean development
contains the required linear algebra:

```lean
theorem d4_preserves_laplacian_kernel_dim :
  kernelDim (laplacian (applyD4 s g)) =
  kernelDim (laplacian g)
```

The second statement is stronger because, for ordinary finite graphs, the
dimension of the kernel of the graph Laplacian corresponds to the number of
connected components. HAOS-IIP should not rely on that stronger form until the
needed definitions and proof obligations are present locally.

## Control Logic

A valid Betti diagnostic should separate symmetry-preserving perturbations from
structure-destroying controls.

Expected stable controls:

- D4 rotation/reflection of the lattice.
- Relabeling that preserves graph isomorphism.

Expected destructive controls:

- threshold changes that add/remove connectivity,
- topology-destroyed rewiring,
- relation shuffling,
- class swaps that alter the relation graph,
- prime-support replacement such as `71 -> 73` when the graph rule depends on
  supersingular support.

If destructive controls do not move the Betti diagram, the diagnostic is too
coarse for the tested structure.

If symmetry controls move the Betti diagram, the implementation or graph
construction is not symmetry-stable.

## Connection to Existing Outputs

The current executable `GEO-MM-01` diagnostic already checks bounded arithmetic
telemetry:

- supersingular prime support,
- Monster-order exponent address,
- low j-coefficient decomposition witnesses,
- Gaussian-prime residue classes.

The Betti SVG extends this as a visual/topological sidecar:

```text
arithmetic support telemetry
-> relation graph
-> Betti diagram
-> symmetry-recoverability reading
```

That is a natural HAOS-IIP bridge because it turns arithmetic structure into a
recoverable topological signature.

## Boundary

This note does not claim:

- a proof of Monstrous Moonshine,
- a construction of the Moonshine module,
- a complete Monster irrep analysis,
- a physical bridge,
- a continuum limit,
- quantum, gravity, or field-theory derivation.

It only states that Betti diagrams are a good HAOS-IIP diagnostic language for
testing whether an arithmetic constellation preserves topological structure
under declared perturbations.

