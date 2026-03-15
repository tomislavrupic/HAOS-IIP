# Stage C0.1 - Pure Combinatorial Kernel (No-Distance)

## Purpose

Design and execute the first austerity-grade kernel experiment in which:
- no coordinate metric
- no embedded distance oracle
- no Gaussian or continuous falloff
- no implicit geometry

are allowed.

The goal is to test whether structured packet transport / exchange topology can arise purely from graph-incidence combinatorics and internal state rules.

This marks the entry into Combinatorial Phase C0 of the Pre-Geometry Atlas.

## Architecture Constraints

The simulation must obey the following hard constraints:

1. Graph substrate
- finite node set
- edges defined only by adjacency list
- optional heterogeneous degree
- no coordinates stored anywhere

2. Kernel definition
Interaction weights must be constructed only from:
- node degree
- shared-neighbor count
- discrete graph distance (shortest path length as integer)
- local motif membership (triangle / square / cluster flags)

Continuous functions of spatial distance are forbidden.

3. Allowed kernel examples
The experiment should support selectable kernel modes such as:
- uniform adjacency kernel
- inverse-degree kernel
- shared-neighbor reinforcement kernel
- motif-gated kernel
- bounded hop-radius kernel

4. State variables
Each packet carries:
- amplitude
- internal phase
- optional vector / DK multi-component label

5. Dynamics
Update rule must be:
- local
- discrete-time
- linear or weakly nonlinear
- invariant under graph relabeling

## Experiment Objective

Test whether the following phenomena can appear without metric structure:
- persistent packet transport
- exchange-like crossing
- braid-like topology
- smeared transfer band
- localization or trapping

The key research question:

Can topology-class behavior emerge purely from combinatorial interaction structure?

## Measurement Layer

For each run record:
- packet centroid trajectory in node index space
- overlap integral vs time
- exchange topology classifier
- coherence order parameter
- transport efficiency metric

Topology classes:
- braid_like_exchange
- transfer_smeared
- trapped
- dispersive

## Variation Axes

Stage C0.1 should scan minimal axes:
- kernel mode
- graph irregularity level
- packet internal phase offset
- reinforcement strength parameter

Each run must vary only one dominant axis.

## Output Requirements

The runner must produce:
- JSON result file per run
- topology classification summary table
- packet flow visualization on graph
- coherence heatmap

## Interpretation Rule

If structured exchange or corridor-like topology persists:

-> evidence that metric distance is not required for early transport phenomenology.

If all runs reduce to diffusion:

-> metric-like weighting may be necessary for braid protection.

## Continuation Hint

Natural next steps after this branch include:
- C0.2 - degree-spectrum protection test
- C0.3 - dynamic edge rewiring memory probe
