#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from morphogenesis_4d_model import (
    MorphogenesisModel,
    apply_developmental_drift,
    build_tissue_lattice,
    simulate_development,
    run_perturbation_sweep,
)


EXPERIMENT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
RESULTS_CSV = OUTPUT_DIR / "results.csv"
ROBUSTNESS_CSV = OUTPUT_DIR / "robustness_check.csv"
VALIDATION_MD = EXPERIMENT_DIR / "morphogenesis_4d_validation.md"

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
        float(summary["baseline_recoverability"]) >= 0.90
        and float(summary["baseline_final_morphology_match"]) >= 0.90
        and float(summary["baseline_trajectory_coherence"]) >= 0.90
        and float(summary["baseline_interaction_support_retention"]) >= 0.999
        and float(summary["baseline_spatial_continuity_retention"]) >= 0.999
    )


def recoverability_declines_gradually(rows: list[dict[str, object]]) -> bool:
    values = [float(row["recoverability"]) for row in rows]
    deltas = np.diff(values)
    if not deltas.size:
        return False
    total_decline = values[0] - values[-1]
    largest_drop = abs(float(np.min(deltas)))
    upward_bump = float(np.max(deltas))
    return total_decline > 0.25 and largest_drop < 0.16 and upward_bump <= 0.01


def deterministic_repeat(first_rows: list[dict[str, object]], second_rows: list[dict[str, object]]) -> bool:
    for left, right in zip(first_rows, second_rows, strict=True):
        for field in (
            "perturbation_level",
            "recoverability",
            "delta_persistence",
            "safety_margin",
            "confidence",
            "final_morphology_match",
            "trajectory_coherence",
            "interaction_support_retention",
            "spatial_continuity_retention",
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


def plot_snapshots(filename: str, model: MorphogenesisModel, trajectory: np.ndarray, title: str) -> None:
    x_size, y_size, z_size = model.grid_size
    mid_z = z_size // 2
    indices = (0, model.developmental_steps // 2, model.developmental_steps - 1)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    for ax, time_index in zip(axes, indices, strict=True):
        state = trajectory[time_index].reshape(model.grid_size)
        im = ax.imshow(state[:, :, mid_z].T, origin="lower", cmap="viridis", vmin=0.0, vmax=1.0)
        ax.set_title(f"t={time_index}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xticks(range(x_size))
        ax.set_yticks(range(y_size))
    fig.suptitle(title)
    fig.colorbar(im, ax=axes.ravel().tolist(), label="cell state", shrink=0.82)
    fig.savefig(OUTPUT_DIR / filename, dpi=160, bbox_inches="tight")
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


def save_plots(model: MorphogenesisModel, rows: list[dict[str, object]], summary: dict[str, object]) -> None:
    baseline_trajectory = simulate_development(model)
    collapse_index = summary["k_star"]
    selected_level = 1.0 if collapse_index is None else float(PERTURBATION_LEVELS[int(collapse_index)])
    perturbed_trajectory = simulate_development(apply_developmental_drift(model, selected_level))
    plot_snapshots(
        "developmental_snapshots_baseline.png",
        model,
        baseline_trajectory,
        "Baseline developmental snapshots",
    )
    plot_snapshots(
        "developmental_snapshots_perturbed.png",
        model,
        perturbed_trajectory,
        f"Perturbed developmental snapshots, p={selected_level:.3f}",
    )
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
        "morphology_and_trajectory_metrics.png",
        rows,
        ("final_morphology_match", "trajectory_coherence", "spatial_continuity_retention"),
        "retention",
        "Morphology and trajectory metrics",
        threshold=0.50,
    )


def run_primary() -> tuple[MorphogenesisModel, list[dict[str, object]], dict[str, object]]:
    model = build_tissue_lattice()
    rows, summary = run_perturbation_sweep(
        model,
        PERTURBATION_LEVELS,
        mode="developmental_drift",
        collapse_threshold=COLLAPSE_THRESHOLD,
        sustain_steps=SUSTAIN_STEPS,
    )
    return model, rows, summary


def run_robustness_check() -> tuple[list[dict[str, object]], bool]:
    model = build_tissue_lattice(target_weight=0.65 * 1.02)
    rows, summary = run_perturbation_sweep(
        model,
        PERTURBATION_LEVELS,
        mode="developmental_drift",
        collapse_threshold=COLLAPSE_THRESHOLD,
        sustain_steps=SUSTAIN_STEPS,
    )
    collapse_level = k_star_level(summary)
    visible_level = first_visible_failure_level(rows)
    result = {
        "variation": "target_weight_x_1_02",
        "baseline_stable": baseline_stable(summary),
        "k_star_level": collapse_level,
        "first_visible_failure_level": visible_level,
        "early_detection": collapse_level is not None and visible_level is not None and collapse_level < visible_level,
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
        """# 4D Morphogenesis Demo

This is Biology Line D.

It is an experimental application layer on top of frozen HAOS-IIP. It models 4D morphogenesis as a 3D tissue state evolving through developmental time and tests early recoverability loss under time-dependent developmental perturbation.

The model is a toy deterministic 8 x 8 x 8 cell lattice with 24 developmental steps. It does not modify HAOS core and makes no claims about real developmental biology.

Run from the repository root:

```bash
python experiments/biology/morphogenesis_4d_demo/run_morphogenesis_4d_demo.py
```

Outputs are written to `experiments/biology/morphogenesis_4d_demo/outputs/`.
""",
        encoding="utf-8",
    )


def write_validation(
    model: MorphogenesisModel,
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
    edge_classes = sorted({edge.edge_class for edge in model.edges})

    VALIDATION_MD.write_text(
        f"""# 4D Morphogenesis Demo Validation

## System
The system is a deterministic toy 4D developmental structure: a 3D grid with size {model.grid_size}, {model.node_count} cells, and {model.developmental_steps} developmental time steps.

Each cell has a scalar state in [0, 1]. The target morphology is deterministic: early development is a compact central form, mid development introduces an axis polarity gradient, and late development transitions toward a core-shell/anterior-posterior pattern.

Edge classes are {", ".join(edge_classes)}. The local_neighbor class gives 6-connected nearest-neighbor support, diagonal_neighbor gives weak local support, and developmental_signal gives deterministic longer coordination links.

The primary perturbation is developmental drift. It delays and weakens target pull in a deterministic region centered at {model.drift_center} during mid-to-late development, while weakening developmental_signal links in the same affected region. A local developmental lesion mode is implemented but is not the primary validation sweep.

## Metrics
These are lightweight external proxy diagnostics aligned with HAOS-IIP stability logic. They are not frozen HAOS core metrics.

- recoverability: `0.35 * final_morphology_match + 0.32 * trajectory_coherence + 0.25 * interaction_support_retention + 0.08 * spatial_continuity_retention`, clipped to [0, 1].
- delta_persistence: change in recoverability between consecutive perturbation levels.
- k_star: first perturbation index where recoverability remains below {COLLAPSE_THRESHOLD:.2f} for the current sample plus {SUSTAIN_STEPS} following samples.
- safety_margin: distance in perturbation parameter from current p to p at k_star. At and beyond k_star this value is zero or negative.
- confidence: magnitude of the strongest negative delta_persistence observed so far.
- visible_failure: True when `final_morphology_match < 0.50` or `spatial_continuity_retention < 0.50`.

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
This is a toy 4D developmental model.

It is not real morphogenesis.

It is not a biological validation claim.

It makes no consciousness claims.

It makes no quantum biology claims.

It does not modify or prove HAOS-IIP core.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_readme()

    model, rows, summary = run_primary()
    _, repeat_rows, _ = run_primary()
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
            "final_morphology_match",
            "trajectory_coherence",
            "interaction_support_retention",
            "spatial_continuity_retention",
            "visible_failure",
        ],
    )
    save_plots(model, rows, summary)
    write_validation(model, rows, summary, checks, robustness_rows)

    print(f"baseline_stable: {checks['baseline_stable']}")
    print(f"k_star_level: {collapse_level}")
    print(f"first_visible_failure_level: {visible_level}")
    print(f"early_detection: {checks['early_detection']}")
    print(f"robustness_early_detection: {checks['robustness_early_detection']}")
    print("outputs_written: experiments/biology/morphogenesis_4d_demo/outputs/")


if __name__ == "__main__":
    main()
