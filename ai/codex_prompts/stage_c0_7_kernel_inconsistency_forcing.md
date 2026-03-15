# Stage C0.7 Kernel-Update Inconsistency Forcing and Explicit Irreversibility

Purpose:
- test whether explicit adjacency to graph-shell update inconsistency can mediate protected ordering on the no-distance branch
- keep the corrected C0.6 positive-mismatch forcing schedules and the protected topology-stability constraint fixed
- isolate the missing ingredient before any Phase II kernel back-reaction claim

Frozen architecture:
- same graph size, representative encounter families, projected-transverse sector, and shell-weight vocabulary as Paper 24.1 and C0.6
- same C0.3 incidence-noise and C0.4 motif forcing menus
- same no-distance edge-graph control branch

New forcing rule:
- once positive mismatch density has started, deliberately alternate adjacency and graph-shell survivor updates
- at each forced step compare the two second-order compositions: adjacency then graph-shell versus graph-shell then adjacency
- count explicit inconsistency as the survivor-partition disagreement between those two compositions on the same forced graph

Primary observables:
- positive mismatch density
- kernel update consistency failure count
- protected arrow step under inconsistency
- persistence depth under inconsistency

Honesty boundary:
- if explicit kernel-update inconsistency becomes positive and protected arrow forcing still appears across the matrix, the branch supports an irreversibility mechanism tied to the same incidence and kernel machinery as mismatch accumulation
- if explicit inconsistency stays absent or the protected arrow collapses when alternation is enforced, the graph-shell branch remains too rigid and Phase II is not yet justified
