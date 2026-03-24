# Toy Sequence Recoverability Controls

This directory holds the bounded multi-family control suite for the 9-event toy sequence claim.

What the script now tests:
- keeps the toy setup fixed at 9 events, 6 protected positions, and protected pairwise coherence as the score
- compares several perturbation families instead of one:
  - `random_swaps`
  - `adjacent_swaps`
  - `corridor_limited_swaps`
  - `burst_errors`
  - `anisotropic_local_disruptions`
- keeps the same exhaustive ranking over all `84` possible 6-position protected subsets for each perturbation family
- asks one decision question only:
  - does the fixed candidate become more special under the perturbation family used as the HAOS-story proxy?

Default command:

```bash
cd /Volumes/Samsung\ T5/2026/HAOS/HAOS\ DOCS/HAOS-IIP
python3 numerics/simulations/toy_sequence_recoverability_controls.py
```

Higher-stability run:

```bash
cd /Volumes/Samsung\ T5/2026/HAOS/HAOS\ DOCS/HAOS-IIP
python3 numerics/simulations/toy_sequence_recoverability_controls.py --trials 100000 --ranking-trials 30000
```

Optional story-family override:

```bash
cd /Volumes/Samsung\ T5/2026/HAOS/HAOS\ DOCS/HAOS-IIP
python3 numerics/simulations/toy_sequence_recoverability_controls.py --story-family corridor_limited_swaps
```

Outputs:
- timestamped JSON in `data/`
- timestamped CSV in `data/`
- stable markdown table at `experiments/toy_sequence_controls/results_table.md`
- memo at `experiments/toy_sequence_controls/memo.md`

Important constraints:
- no canonical protected-subset selection rule for the reported toy result was located in the workspace
- the default fixed candidate subset is therefore treated only as the subset under test, not as a derived HAOS rule
- no canonical HAOS perturbation-family definition was located in the workspace
- the default story-family proxy is therefore an explicit modeling choice: `anisotropic_local_disruptions`
- full repair remains trivial and is kept only in the machine-readable output; the family comparison focuses on baseline and bounded local repair
