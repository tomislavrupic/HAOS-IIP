# Phase XVIII Summary

Objective: reuse the frozen Phase XV-XVII ledgers to test whether causal depth and propagation arrival order define a coherent distance surrogate without resimulating dynamics.

Frozen slice: `bias_onset`, `n_side = 60, 72, 84`, seeds `1303, 1302`, ensemble size `7`.

Gate results:
- Shell ordering coherence: `True` with branch max shell overlap `0.365779` and monotonic shell fraction `1.000000`.
- Refinement scaling stability: `True` with branch max adjacent slope drift `0.004762` and causal-depth drift `0.000000`.
- Branch/control separation: `True` with control max shell overlap `0.444716` and control max causal-depth drift `0.714286`.
- Triangle consistency: `True` with branch/control violation rates `0.000000` / `0.000000`.

Interpretation: the branch slice keeps shell-arrival ordering monotone, preserves a compact arrival-vs-depth slope band across refinement, and stays cleaner than the matched control under the same surrogate gates. No new dynamics were run; all values are derived from the frozen ledgers only.

Phase XVIII establishes proto-geometric distance-surrogate feasibility for the frozen operator hierarchy.
