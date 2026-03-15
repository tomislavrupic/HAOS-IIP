# Stage C0.6 Incidence-Seeded Mismatch Accumulation and Protected Arrow Forcing

Purpose:
- test whether the no-distance protected sequencing branch can generate positive incidence-mismatch ledger entries once external combinatorial forcing is introduced
- keep the C0.5 protected topology-stability constraint active throughout
- use only graph-internal perturbations already audited in C0.3 and C0.4

Frozen architecture:
- same graph size and projected-transverse seed family as Paper 24.1
- same adjacency and graph-shell kernel pair as C0.5
- same representative encounter families: clustered composite anchor, counter-propagating corridor, phase-ordered symmetric triad
- same no-distance edge-graph control branch

Forcing protocols:
- noise-seeded protected sequencing: step through the exact C0.3 incidence-noise ladder from ultra-weak rewiring to connectivity-edge stress
- motif-seeded protected sequencing: step through the exact C0.4 motif set of triangle, diamond, and hub micro-star
- hybrid seeded protected sequencing: compare the scheduled noise and motif candidates at each step and select the protected candidate with the stronger mismatch forcing signal

Primary observables:
- persistence-class stabilization depth
- incidence-mismatch ledger density
- kernel-update consistency failures
- protected arrow step
- positive mismatch density
- protected arrow step distribution

Honesty boundary:
- if positive mismatch density remains zero and protected arrow forcing appears only as an isolated event, the branch still behaves as a rigid persistence-class control
- if positive mismatch density becomes nonzero and protected arrow forcing recurs across multiple runs, the branch supports an incidence-driven sequencing mechanism worth lifting into kernel back-reaction
