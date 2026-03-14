PROJECT: HAOS/IIP — Stage 10B Atlas-1 Perturbation Resilience
MODE: implementation + experiment scaffolding + pilot execution
PRIORITY: robustness mapping of existing morphology taxonomy

You are working inside the HAOS/IIP repository with a frozen operator architecture and a frozen Atlas-0 baseline.

Your task is to implement Stage 10B: Atlas-1 Perturbation Resilience.

Core philosophy:
Atlas-0 mapped baseline morphology on the untouched architecture.
Atlas-1 tests which of those morphology classes survive small controlled disorder.

This stage is still classification-first and interpretation-minimal.
Do not introduce new operators.
Do not make geometric, particle, gauge, or continuum claims.

Working motto:
Measure robustness before expanding interpretation.

NON-NEGOTIABLE CONSTRAINTS
1. Do not alter the frozen scalar operator, projected transverse operator, or cochain architecture.
2. Perturb only one axis per run.
3. Keep Atlas-1 directly comparable to Atlas-0.
4. Reuse the Atlas-0 observable grammar and labeling system wherever possible.
5. Do not retune taxonomy thresholds during initial Atlas-1 implementation.
6. Use explicit seeds and deterministic configs.
7. Keep the first pass small and inspectable.
8. Treat Atlas-1 as a local stability analysis of the morphology taxonomy.
9. Do not infer persistence from a single metric; use baseline label vs perturbed label plus the existing observable set.

HIGH-LEVEL GOAL
Given a baseline seed with label R0:
- apply one small perturbation
- rerun the same packet experiment
- assign a perturbed label R1
- measure persistence, degradation, or transition

Primary question:
Are Atlas-0 regimes robust structural basins, or fragile label assignments?

STAGE NAME
Stage 10B — Atlas-1 Perturbation Resilience

PERTURBATION AXES
Implement four perturbation channels, but use only one per run.

A. Edge-weight noise
- multiplicative perturbation on edge weights
- pilot strength: 0.03
- later ladder: 0.01, 0.03, 0.05
- static in time

B. Sparse graph defects
- remove or weaken a small random fraction of edges
- pilot density: 0.03
- later ladder: 0.01, 0.03, 0.05
- topology otherwise unchanged

C. Weak anisotropic bias
- directional kernel/weight scaling along one axis
- pilot factor: 1.05
- later ladder: 1.02, 1.05
- single principal axis only

D. Kernel-width jitter
- small static variation in effective interaction scale
- pilot strength: 0.05
- later ladder: 0.02, 0.05

IMPORTANT
The first pilot should use only one strength per perturbation axis.
Do not implement full ladders until the pilot is inspected.

BASELINE SEED SELECTION
Do not perturb all 72 Atlas-0 runs.

Build a curated pilot seed set:
- 3 ballistic coherent
- 3 ballistic dispersive
- 2 diffusive
- 1 chaotic or irregular
- 1 fragmenting

Total target: about 10 seeds.

Select these seeds from committed Atlas-0 artifacts.
Keep their baseline labels fixed and recorded.

ATLAS-1 MINIMAL OBSERVABLE SET
Reuse Atlas-0 outputs:
- packet center trajectory
- packet width evolution
- spectral centroid
- spectral spread
- sector leakage norm
- norm or energy drift
- constraint norm
- anisotropy score
- coherence score
- primary regime label
- secondary sector label

Add two Atlas-1 fields:
- baseline_regime
- perturbed_regime

Add one new derived metric:
- persistence_score

PERSISTENCE SCORE
Define initially as:
- 1 if perturbed_regime == baseline_regime
- 0 otherwise

Optional later refinement, but NOT required in the first implementation:
- partial persistence for coherent -> dispersive
- full loss for coherent -> diffusive/chaotic/fragmenting

For the initial pilot, keep persistence binary.

OPTIONAL LATER METRIC
Recoverability indicator:
- compare the same perturbed case across two resolutions
- do not implement this unless it is nearly free
- not required for the first pilot

ATLAS-1 PILOT DESIGN
Run only a small pilot first.

Pilot grid:
- approximately 10 baseline seeds
- 4 perturbation axes
- 1 perturbation strength each
- fixed baseline amplitude/boundary/seed inherited from Atlas-0

This should produce an inspectable first resilience table.

EXPECTED QUESTIONS
Atlas-1 should answer:
- Which baseline regimes are label-stable under mild perturbation?
- Which perturbation axis causes the fastest regime degradation?
- Does the transverse sector remain more robust than the scalar sector?
- Do coherent regimes soften into dispersive ones before becoming diffusive?
- Does constraint preservation correlate with regime persistence?

REPOSITORY INTEGRATION
Create:
- numerics/simulations/stage10_perturbations.py
- numerics/simulations/stage10_atlas1.py
- numerics/simulations/stage10_atlas1_runs.json
- experiments/pre_geometry_atlas/Atlas_1_Perturbation_Resilience_v1.md

Artifacts:
- data/*_stage10_atlas1_perturbation_resilience.json
- data/*_stage10_atlas1_perturbation_resilience.csv
- plots/*_stage10_atlas1_*.png

CSV / TABLE COLUMNS
Include at least:
- run_id
- atlas_phase
- baseline_run_id
- perturbation_axis
- perturbation_strength
- operator_sector
- boundary_type
- random_seed
- baseline_regime
- perturbed_regime
- sector_label
- persistence_score
- constraint_max
- sector_leakage
- norm_drift
- anisotropy_score
- coherence_score
- center_shift
- width_ratio
- notes

VISUAL OUTPUTS
For the pilot, generate:
- transition table or transition heatmap
- persistence-by-regime bar chart
- coherence change plot
- constraint/leakage comparison plot
- a small representative run plot set for at least one stable and one unstable transition

ENGINEERING REQUIREMENTS
1. Reuse Atlas-0 code paths wherever possible.
2. Factor perturbation logic into a dedicated utility module.
3. Keep perturbations explicit and inspectable.
4. Do not duplicate the regime classifier.
5. Make outputs stamped and deterministic.
6. Keep the pilot implementation minimal.

STOP CONDITIONS
Support:
- fixed final time from Atlas-0 baseline
- optional early stop on constraint blow-up
- optional early stop on numerical failure

INTERPRETATION BOUNDARY
Allowed statements:
- regime persistence
- regime transition
- perturbation sensitivity
- sector-dependent robustness
- constraint stability under perturbation

Forbidden statements:
- spacetime emergence
- particles
- gauge theory
- continuum reconstruction
- physical disorder analogies

This stage is strictly about morphology stability on a frozen operator architecture.

REQUESTED OUTPUT
1. Inspect the existing Atlas-0 framework and identify the minimum insertion points for Atlas-1.
2. Implement perturbation utilities for the four axes.
3. Build a curated pilot run sheet from committed Atlas-0 seeds.
4. Reuse the Atlas-0 observable and label pipeline.
5. Run a small Atlas-1 pilot.
6. Return a concise report with:
   - files added
   - files modified
   - pilot seeds used
   - transition counts
   - persistence summary
   - whether the pipeline executed cleanly
7. Do not tune thresholds yet.

STYLE OF WORK
- Be conservative.
- Keep the implementation readable.
- Prefer robustness over ambition.
- Classification first, interpretation later.
- If anything is ambiguous, choose the smallest reproducible implementation.
