# Gauge Tests

This folder is reserved for explicit edge/Hodge and holonomy experiments.

## Current concrete target

Build the weighted edge/Hodge branch from the current 3D substrate:

1. construct incidence matrix `B`
2. construct weighted edge-face incidence `C`
3. construct weighted edge operator `L1`
4. compute its low spectrum
5. separate exact gradient-like from coexact transverse-like low modes
6. visualize circulation structure
7. compare pure-gauge and nontrivial-flux backgrounds

This is now an active numerical branch in the repository. The next hard question is whether any transverse-like `L1` family survives beyond one clean cubic-lattice setup.
