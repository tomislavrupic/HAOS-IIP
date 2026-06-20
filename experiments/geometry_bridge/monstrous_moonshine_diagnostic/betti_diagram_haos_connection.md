# Betti Diagram and HAOS-IIP Connection

## Status

`GEO-MM-01` treats the Betti diagram as an arithmetic/topological diagnostic
sidecar, not as a proof of Monstrous Moonshine and not as a physical bridge.

The executable local sidecar is:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_component_count.py
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_vector.py
```

The current declared relation graph has `Betti_0 = 4` and `33` relation edges.
Its result hash is `betti_component_1454dd195fe727984643d0dc`.

The Betti vector runner uses the same graph and reports:

```json
{
  "Betti_0": 4,
  "Betti_1": 22,
  "nodes": 15,
  "edges": 33
}
```

The `Betti_1` value is the finite-graph cycle count:

```text
Betti_1 = E - V + C = 33 - 15 + 4 = 22
```

Its result hash is `betti_vector_1d5373fa671fd513fe2d5ac0`.

The local robustness sweep is:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_threshold_sweep.py
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_null_ensemble.py
```

It reports:

- `Betti_0` stability band: `[8.0, 11.5]`
- exact edge-signature band: `[8.0, 8.5]`
- edge-neighborhood band: `[8.0, 10.0]`
- deterministic null false-positive rate: `0.110000`

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
has a direct graph interpretation: component count. That target is now
implemented locally by `run_betti_component_count.py`.

The first local extension is `Betti_1`, implemented by
`run_betti_vector.py`. This does not add a higher-dimensional Moonshine
topology claim. It only records independent cycles in the declared finite
undirected relation graph.

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

## Current Executable Target

The local executable target is:

```text
declared Gaussian-prime representative graph
-> D4 / relabel perturbation
-> same component count and edge signature
```

Current stable controls:

- `d4_rotate_90`
- `d4_reflect_real`
- `isomorphism_relabel`

Current destructive controls:

- `threshold_mutation_control`
- `topology_destroyed_control`
- `support_replacement_control`

Current output:

```text
known_positive: Betti_0 = 4
d4_rotate_90: Betti_0 = 4
d4_reflect_real: Betti_0 = 4
isomorphism_relabel: Betti_0 = 4
threshold_mutation_control: Betti_0 = 6
topology_destroyed_control: Betti_0 = 2
support_replacement_control: Betti_0 = 3
```

Current Betti-vector output:

```text
known_positive: Betti_0 = 4, Betti_1 = 22
d4_rotate_90: Betti_0 = 4, Betti_1 = 22
d4_reflect_real: Betti_0 = 4, Betti_1 = 22
isomorphism_relabel: Betti_0 = 4, Betti_1 = 22
threshold_mutation_control: Betti_0 = 6, Betti_1 = 9
topology_destroyed_control: Betti_0 = 2, Betti_1 = 20
support_replacement_control: Betti_0 = 3, Betti_1 = 28
```

## Future Lean Target

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

The Betti component-count sidecar extends this as a local topological
diagnostic:

```text
arithmetic support telemetry
-> relation graph
-> Betti diagram
-> symmetry-recoverability reading
```

The SVG remains a visual companion. The local executable authority for Betti_0
is `run_betti_component_count.py`. The local executable authority for the
finite-graph Betti vector is `run_betti_vector.py`.

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
