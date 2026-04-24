# Biology Line B Numbered Paper v0.1

1. Title

Microtubule Lattice Demo: Recoverable Structure Under Lattice Perturbation.

2. Scope

This note records an experimental biology-layer artifact inside `experiments/biology/microtubule_lattice_demo/`. It is not a theory of consciousness, not Orch-OR, not quantum computation, and not a claim that HAOS explains microtubules.

3. Minimal Biological Background

Microtubules are cylindrical polymers made of tubulin dimers. A common eukaryotic microtubule has 13 protofilaments arranged laterally into a hollow tube. Local longitudinal and lateral interactions stabilize the lattice. This experiment is only a coarse-grained structural graph model.

4. Instrument Question

At what perturbation level does a microtubule-like lattice stop being recoverable under interaction?

5. Method

The demo builds a deterministic 13 x 24 cylindrical lattice with longitudinal, lateral, seam_or_diagonal, and weak_support edge classes. The primary sweep weakens lateral support edges over p in [0, 1]. Recoverability is a weighted proxy combining largest connected component fraction, weighted degree retention, and propagation retention.

6. Result

The baseline run reports `k_star_level=0.8333333333333334` and `first_visible_failure_level=0.8666666666666667` with `early_detection=True`. The robustness variation reports `early_detection=True` with `k_star_level=0.8333333333333334` and `first_visible_failure_level=0.8666666666666667`.

7. Limitation

This is a toy recoverability diagnostic on a coarse structural lattice. It is not molecular simulation, biological validation, consciousness theory, Orch-OR, quantum biology, or a modification of frozen HAOS-IIP core.
