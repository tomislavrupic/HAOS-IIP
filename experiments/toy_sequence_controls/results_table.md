# Toy Sequence Perturbation Family Comparison

- Timestamped JSON: `data/20260316_113225_toy_sequence_recoverability_controls.json`
- Timestamped CSV: `data/20260316_113225_toy_sequence_recoverability_controls.csv`
- Sequence length: `9`
- Protected positions: `6`
- Fixed candidate subset: `(0, 1, 2, 5, 6, 8)`
- Explicit HAOS protected-set selection rule available in workspace: `False`
- Default HAOS-story perturbation proxy: `anisotropic_local_disruptions`

## Family definitions
- `random_swaps`: three unconstrained swaps anywhere in the sequence
- `adjacent_swaps`: three adjacent swaps at arbitrary positions
- `corridor_limited_swaps`: three swaps restricted to one random contiguous corridor of width 4
- `burst_errors`: one random contiguous burst of width 4, internally permuted
- `anisotropic_local_disruptions` (story proxy): one directional local displacement built from three consecutive adjacent swaps

## Question

- Does the fixed candidate become consistently more special under `anisotropic_local_disruptions`, the chosen HAOS-story perturbation proxy?
- Answer: `yes`
- Rule: Answer yes only if the story family gives the fixed candidate the best or tied-best exhaustive rank among all 84 subsets in both baseline and bounded local repair, while also beating matched random and contiguous controls in both modes.

## Family comparison

| family | story proxy | baseline rank / 84 | local rank / 84 | baseline delta vs random | baseline delta vs contiguous | local delta vs random | local delta vs contiguous |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| anisotropic_local_disruptions | True | 7 | 3 | 0.027819 | 0.068755 | 0.022395 | 0.054574 |
| corridor_limited_swaps | False | 10 | 5 | 0.024842 | 0.073055 | 0.019554 | 0.060395 |
| burst_errors | False | 10 | 7 | 0.026900 | 0.078054 | 0.021792 | 0.067228 |
| random_swaps | False | 12 | 12 | 0.015827 | 0.051307 | 0.015361 | 0.049924 |
| adjacent_swaps | False | 17 | 19 | 0.006979 | 0.037220 | 0.004120 | 0.024918 |

## Note

- Full repair is omitted from the comparison table because it restores exact protected order by construction and therefore does not distinguish perturbation families.
- The exhaustive ranking remains the main nontrivial test: each family ranks the fixed candidate against all 84 possible 6-position subsets.
