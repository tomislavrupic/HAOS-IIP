# Phase XVII

Phase XVII consumed only frozen Phase XV–XVI artifacts and restricted the authority slice to `bias_onset`, `n_side = 60, 72, 84`, ensemble size `7`, branch seeds `1303` and `1302`, and the matched control `periodic_diagonal_augmented_control`.

Result: passed.

- Influence-edge reconstruction used only consecutive event-chain edges from the frozen Phase XVI ledgers.
- Branch mean edge reproducibility was `0.626263` with max adjacent edge-set distance `0.333333`; control was `0.585859` and `0.571429`.
- Branch acyclicity stayed in a compact band of width `0.101010`; control band width was `0.393939`.
- Branch mean causal depth was `1.714286` with drift `0.000000`; control drift was `0.714286`.
- Branch order-compatibility mismatch was `0.035906` with drift `0.001571` and remained below control `0.037018`.

Hard stop triggered: `False`.

Phase XVII establishes causal-closure feasibility for the frozen operator hierarchy.
