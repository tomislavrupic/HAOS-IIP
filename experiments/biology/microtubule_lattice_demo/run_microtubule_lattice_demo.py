#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from microtubule_lattice_model import (
    MicrotubuleLattice,
    build_microtubule_lattice,
    run_perturbation_sweep,
)


EXPERIMENT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
RESULTS_CSV = OUTPUT_DIR / "results.csv"
ROBUSTNESS_CSV = OUTPUT_DIR / "robustness_check.csv"
VALIDATION_MD = EXPERIMENT_DIR / "microtubule_lattice_validation.md"
NUMBERED_PAPER_MD = EXPERIMENT_DIR / "Biology_Line_B_Numbered_Paper_v0.1.md"

PERTURBATION_LEVELS = np.linspace(0.0, 1.0, 31)
COLLAPSE_THRESHOLD = 0.70
SUSTAIN_STEPS = 2


def first_visible_failure_level(rows: list[dict[str, object]]) -> float | None:
    for row in rows:
        if row["visible_failure"]:
            return float(row["perturbation_level"])
    return None


def k_star_level(summary: dict[str, object]) -> float | None:
    k_star = summary["k_star"]
    return None if k_star is None else float(PERTURBATION_LEVELS[int(k_star)])


def baseline_stable(summary: dict[str, object]) -> bool:
    return (
        float(summary["baseline_recoverability"]) >= 0.98
        and float(summary["baseline_largest_component_fraction"]) >= 0.999
        and float(summary["baseline_weighted_degree_retention"]) >= 0.999
        and float(summary["baseline_propagation_retention"]) >= 0.999
    )


def recoverability_declines_gradually(rows: list[dict[str, object]]) -> bool:
    values = [float(row["recoverability"]) for row in rows]
    deltas = np.diff(values)
    if not deltas.size:
        return False
    total_decline = values[0] - values[-1]
    largest_drop = abs(float(np.min(deltas)))
    upward_bump = float(np.max(deltas))
    return total_decline > 0.30 and largest_drop < 0.20 and upward_bump <= 0.01


def early_detection(rows: list[dict[str, object]], summary: dict[str, object]) -> bool:
    collapse_level = k_star_level(summary)
    visible_level = first_visible_failure_level(rows)
    return collapse_level is not None and visible_level is not None and collapse_level < visible_level


def deterministic_repeat(first_rows: list[dict[str, object]], second_rows: list[dict[str, object]]) -> bool:
    for left, right in zip(first_rows, second_rows, strict=True):
        for field in (
            "perturbation_level",
            "recoverability",
            "delta_persistence",
            "safety_margin",
            "confidence",
            "largest_component_fraction",
            "weighted_degree_retention",
            "propagation_retention",
        ):
            if not np.isclose(float(left[field]), float(right[field]), atol=1.0e-12, rtol=0.0):
                return False
        if left["k_star"] != right["k_star"] or left["visible_failure"] != right["visible_failure"]:
            return False
    return True


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fieldnames})


def plot_lattice_overview(lattice: MicrotubuleLattice) -> None:
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for edge in lattice.edges:
        sp, sz, _ = lattice.coordinates(edge.source)
        tp, tz, _ = lattice.coordinates(edge.target)
        if abs(sp - tp) > 1 and {sp, tp} != {0, lattice.protofilaments - 1}:
            continue
        color = {
            "longitudinal": "#293241",
            "lateral": "#2a9d8f",
            "seam_or_diagonal": "#e76f51",
            "weak_support": "#8d99ae",
        }[edge.edge_class]
        alpha = 0.18 if edge.edge_class in ("seam_or_diagonal", "weak_support") else 0.35
        ax.plot([sp, tp], [sz, tz], color=color, linewidth=0.7, alpha=alpha)

    xs = []
    ys = []
    for node in range(lattice.node_count):
        p, z, _ = lattice.coordinates(node)
        xs.append(p)
        ys.append(z)
    ax.scatter(xs, ys, s=14, color="#111827", alpha=0.75)
    ax.set_title("Microtubule-inspired lattice overview")
    ax.set_xlabel("protofilament index")
    ax.set_ylabel("longitudinal dimer index")
    ax.set_xlim(-0.5, lattice.protofilaments - 0.5)
    ax.set_ylim(-0.5, lattice.dimers_per_protofilament - 0.5)
    ax.grid(True, alpha=0.18)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "lattice_overview.png", dpi=160)
    plt.close(fig)


def plot_metric(
    filename: str,
    rows: list[dict[str, object]],
    fields: tuple[str, ...],
    ylabel: str,
    title: str,
    *,
    threshold: float | None = None,
) -> None:
    x = [float(row["perturbation_level"]) for row in rows]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for field in fields:
        y = [float(row[field]) for row in rows]
        ax.plot(x, y, marker="o", linewidth=2.0, label=field)
    if threshold is not None:
        ax.axhline(threshold, color="black", linestyle="--", linewidth=1.0)
    if len(fields) > 1:
        ax.legend(fontsize=8)
    ax.set_title(title)
    ax.set_xlabel("perturbation level")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / filename, dpi=160)
    plt.close(fig)


def save_plots(lattice: MicrotubuleLattice, rows: list[dict[str, object]]) -> None:
    plot_lattice_overview(lattice)
    plot_metric(
        "recoverability_vs_perturbation.png",
        rows,
        ("recoverability",),
        "recoverability",
        "Recoverability vs perturbation",
        threshold=COLLAPSE_THRESHOLD,
    )
    plot_metric(
        "delta_persistence_vs_perturbation.png",
        rows,
        ("delta_persistence",),
        "delta persistence",
        "Delta persistence vs perturbation",
        threshold=0.0,
    )
    plot_metric(
        "safety_margin_vs_perturbation.png",
        rows,
        ("safety_margin",),
        "safety margin",
        "Safety margin vs perturbation",
        threshold=0.0,
    )
    plot_metric(
        "component_and_propagation_retention.png",
        rows,
        ("largest_component_fraction", "propagation_retention"),
        "retention",
        "Component and propagation retention",
        threshold=0.85,
    )


def run_primary() -> tuple[MicrotubuleLattice, list[dict[str, object]], dict[str, object]]:
    lattice = build_microtubule_lattice()
    rows, summary = run_perturbation_sweep(
        lattice,
        PERTURBATION_LEVELS,
        mode="lateral_weakening",
        collapse_threshold=COLLAPSE_THRESHOLD,
        sustain_steps=SUSTAIN_STEPS,
    )
    return lattice, rows, summary


def run_robustness_check() -> tuple[list[dict[str, object]], bool]:
    variation = build_microtubule_lattice(longitudinal_strength=1.02)
    rows, summary = run_perturbation_sweep(
        variation,
        PERTURBATION_LEVELS,
        mode="lateral_weakening",
        collapse_threshold=COLLAPSE_THRESHOLD,
        sustain_steps=SUSTAIN_STEPS,
    )
    result = {
        "variation": "longitudinal_strength_x_1_02",
        "baseline_stable": baseline_stable(summary),
        "k_star_level": k_star_level(summary),
        "first_visible_failure_level": first_visible_failure_level(rows),
        "early_detection": early_detection(rows, summary),
        "min_recoverability": summary["min_recoverability"],
    }
    robustness_rows = [result]
    write_csv(
        ROBUSTNESS_CSV,
        robustness_rows,
        [
            "variation",
            "baseline_stable",
            "k_star_level",
            "first_visible_failure_level",
            "early_detection",
            "min_recoverability",
        ],
    )
    return robustness_rows, bool(result["early_detection"])


def write_readme() -> None:
    (EXPERIMENT_DIR / "README.md").write_text(
        """# Microtubule Lattice Demo

This is Biology Line B.

It is an experimental application layer on top of frozen HAOS-IIP. It models a microtubule-inspired cylindrical lattice as a deterministic interaction graph and tests whether lightweight HAOS-style diagnostics detect early degradation under lattice perturbation.

The model uses the minimal grounded background that microtubules are cylindrical polymers made of tubulin dimers, and that a common eukaryotic microtubule has 13 protofilaments arranged laterally into a hollow tube. Local longitudinal and lateral interactions stabilize the coarse lattice.

This demo does not modify HAOS core. It makes no claims about consciousness, Orch-OR, quantum computation, or HAOS explaining microtubules.

Run from the repository root:

```bash
python experiments/biology/microtubule_lattice_demo/run_microtubule_lattice_demo.py
```

Outputs are written to `experiments/biology/microtubule_lattice_demo/outputs/`.
""",
        encoding="utf-8",
    )


def write_validation(
    lattice: MicrotubuleLattice,
    rows: list[dict[str, object]],
    summary: dict[str, object],
    checks: dict[str, bool],
    robustness_rows: list[dict[str, object]],
) -> None:
    collapse_level = k_star_level(summary)
    visible_level = first_visible_failure_level(rows)
    pass_status = all(checks.values())
    check_lines = "\n".join(f"- {name}: {'PASS' if passed else 'FAIL'}" for name, passed in checks.items())
    robustness = robustness_rows[0]
    edge_classes = sorted({edge.edge_class for edge in lattice.edges})

    VALIDATION_MD.write_text(
        f"""# Microtubule Lattice Demo Validation

## System
The system is a deterministic coarse-grained microtubule-inspired interaction lattice with {lattice.protofilaments} protofilaments and {lattice.dimers_per_protofilament} dimers per protofilament, for {lattice.node_count} total nodes.

Each node represents one tubulin dimer position. Edge classes are {", ".join(edge_classes)}. Longitudinal edges are the strongest local supports, lateral edges provide medium circumferential support, seam_or_diagonal edges provide weak deterministic geometric offsets, and weak_support edges provide very weak longer-range support.

The primary perturbation mode is lateral weakening: lateral edge strengths are multiplied by `1 - p` for p in [0, 1]. A deterministic defect patch mode is implemented in `microtubule_lattice_model.py` around protofilaments 4-6 and z indices 10-14, but it is not the primary validation sweep.

## Metrics
These are lightweight external proxy diagnostics aligned with HAOS-IIP stability logic. They are not frozen HAOS core metrics.

- recoverability: `0.05 * largest_component_fraction + 0.90 * weighted_degree_retention + 0.05 * propagation_retention`, clipped to [0, 1].
- delta_persistence: change in recoverability between consecutive perturbation levels.
- k_star: first perturbation index where recoverability remains below {COLLAPSE_THRESHOLD:.2f} for the current sample plus {SUSTAIN_STEPS} following samples.
- safety_margin: distance in perturbation parameter from current p to p at k_star. At and beyond k_star this value is zero or negative.
- confidence: magnitude of the strongest negative delta_persistence observed so far.
- visible_failure: True when `largest_component_fraction < 0.85` or `propagation_retention < 0.50`.

The edge score is `interaction_strength * locality_stability * perturbation_factor`, with `locality_stability = exp(-distance / lambda_locality)` and cylindrical distance used only as a coarse locality measure.

The filtration order sorts edge supports by this score, so strongest and most local supports enter before weaker or less local supports.

## Result
Pass status: {"PASS" if pass_status else "FAIL"}

{check_lines}

- baseline stable: {checks["baseline_stable"]}
- k_star perturbation level: {collapse_level}
- first visible failure level: {visible_level}
- early_detection: {checks["early_detection"]}
- robustness variation: {robustness["variation"]}
- robustness early_detection: {robustness["early_detection"]}
- minimum recoverability: {summary["min_recoverability"]:.6f}

## Limitations
This is a coarse-grained toy model.

It is not molecular dynamics.

It is not real biological validation.

It makes no consciousness claims.

It makes no Orch-OR claims.

It makes no quantum computation claims.

It does not modify or prove HAOS-IIP core.
""",
        encoding="utf-8",
    )


def write_numbered_paper(
    rows: list[dict[str, object]],
    summary: dict[str, object],
    robustness_rows: list[dict[str, object]],
) -> None:
    robustness = robustness_rows[0]
    NUMBERED_PAPER_MD.write_text(
        f"""# Biology Line B Numbered Paper v0.1

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

The baseline run reports `k_star_level={k_star_level(summary)}` and `first_visible_failure_level={first_visible_failure_level(rows)}` with `early_detection={early_detection(rows, summary)}`. The robustness variation reports `early_detection={robustness["early_detection"]}` with `k_star_level={robustness["k_star_level"]}` and `first_visible_failure_level={robustness["first_visible_failure_level"]}`.

7. Limitation

This is a toy recoverability diagnostic on a coarse structural lattice. It is not molecular simulation, biological validation, consciousness theory, Orch-OR, quantum biology, or a modification of frozen HAOS-IIP core.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_readme()

    lattice, rows, summary = run_primary()
    repeat_lattice, repeat_rows, repeat_summary = run_primary()
    del repeat_lattice, repeat_summary
    robustness_rows, robustness_early_detection = run_robustness_check()

    collapse_level = k_star_level(summary)
    visible_level = first_visible_failure_level(rows)
    checks = {
        "baseline_stable": baseline_stable(summary),
        "recoverability_declines_gradually": recoverability_declines_gradually(rows),
        "early_detection": collapse_level is not None and visible_level is not None and collapse_level < visible_level,
        "deterministic_repeated_runs": deterministic_repeat(rows, repeat_rows),
        "robustness_early_detection": robustness_early_detection,
    }

    write_csv(
        RESULTS_CSV,
        rows,
        [
            "perturbation_index",
            "perturbation_level",
            "recoverability",
            "delta_persistence",
            "k_star",
            "safety_margin",
            "confidence",
            "largest_component_fraction",
            "weighted_degree_retention",
            "propagation_retention",
            "visible_failure",
        ],
    )
    save_plots(lattice, rows)
    write_validation(lattice, rows, summary, checks, robustness_rows)
    write_numbered_paper(rows, summary, robustness_rows)

    print(f"baseline_stable: {checks['baseline_stable']}")
    print(f"k_star_level: {collapse_level}")
    print(f"first_visible_failure_level: {visible_level}")
    print(f"early_detection: {checks['early_detection']}")
    print(f"robustness_early_detection: {checks['robustness_early_detection']}")
    print("outputs_written: experiments/biology/microtubule_lattice_demo/outputs/")


if __name__ == "__main__":
    main()
