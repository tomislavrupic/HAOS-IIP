# Stage 24.4 / IV-A4 - Stable-Smeared Law-and-Ordering Extension Boundary Test

## Purpose
Run the first bounded extension test for the Phase IV law-and-ordering package on the frozen clustered Dirac-Kahler branch.

This stage carries forward exactly two frozen objects:
- `stable_closed_smeared = 1 iff flow_concentration_index <= 0.884308`
- `inverse_flow_survival_ordering`: inside the stable smeared subset, lower flow concentration orders with longer topology survival

## Scope discipline
- no operator change
- no threshold refit
- no braid rescue
- no family-wide generalization claim

The job is to identify the first honest extension boundary of the current Phase IV package under bounded structural perturbation.

## Extension classes
- `E1` clustered texture perturbation
- `E2` support anisotropy / support skew
- `E3` motif structure perturbation
- `E4` corridor geometry perturbation

Each class is tested at `baseline`, `mild`, and `moderate` levels on a small five-phase lattice.

## Required test order
1. Evaluate the fixed Stage 24.2 threshold law without retuning.
2. Evaluate the Stage 24.3 ordering relation only inside law-valid stable subsets.
3. Assign one of:
   - `law_survives_ordering_survives`
   - `law_survives_ordering_degrades`
   - `law_fails_ordering_not_applicable`
   - `branch_breakdown`

## Interpretation boundary
At most claim:
the Phase IV law-and-ordering package survives only within a bounded extension radius on the frozen clustered DK branch.

Do not broaden this into a family-wide transport law or a braid recovery story.
