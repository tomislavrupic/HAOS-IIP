#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tissue_3d_model import TissueLattice, build_tissue_lattice, run_perturbation_sweep


EXPERIMENT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
RESULTS_CSV = OUTPUT_DIR / "results.csv"
ROBUSTNESS_CSV = OUTPUT_DIR / "robustness_check.csv"
VALIDATION_MD = EXPERIMENT_DIR / "tissue_3d_validation.md"

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
        and float(summary["baseline_spatial_coherence_retention"]) >= 0.999
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
    return total_decline > 0.25 and largest_drop < 0.18 and upward_bump <= 0.01


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
            "spatial_coherence_retention",
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


def plot_tissue_slice(lattice: TissueLattice) -> None:
    x_size, y_size, z_size = lattice.grid_size
    mid_z = z_size // 2
    lesion = np.zeros((x_size, y_size), dtype=float)
    center = np.array(lattice.lesion_center, dtype=float)
    for x in range(x_size):
        for y in range(y_size):
            coords = np.array((x, y, mid_z), dtype=float)
            distance = float(np.linalg.norm(coords - center))
            lesion[x, y] = max(0.0, 1.0 - distance / lattice.lesion_radius)

    fig, ax = plt.subplots(figsize=(5.5, 5.0))
    im = ax.imshow(lesion.T, origin="lower", cmap="magma", vmin=0.0, vmax=1.0)
    ax.set_title("3D tissue lesion slice overview")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xticks(range(x_size))
    ax.set_yticks(range(y_size))
    fig.colorbar(im, ax=ax, label="lesion exposure")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "tissue_slice_overview.png", dpi=160)
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


def save_plots(lattice: TissueLattice, rows: list[dict[str, object]]) -> None:
    plot_tissue_slice(lattice)
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
        "component_spatial_propagation_retention.png",
        rows,
        ("largest_component_fraction", "spatial_coherence_retention", "propagation_retention"),
        "retention",
        "Component, spatial, and propagation retention",
        threshold=0.85,
    )


def run_primary() -> tuple[TissueLattice, list[dict[str, object]], dict[str, object]]:
    lattice = build_tissue_lattice()
    rows, summary = run_perturbation_sweep(
        lattice,
        PERTURBATION_LEVELS,
        mode="local_lesion",
        collapse_threshold=COLLAPSE_THRESHOLD,
        sustain_steps=SUSTAIN_STEPS,
    )
    return lattice, rows, summary


def run_robustness_check() -> tuple[list[dict[str, object]], bool]:
    variation = build_tissue_lattice(local_neighbor_strength=1.02)
    rows, summary = run_perturbation_sweep(
        variation,
        PERTURBATION_LEVELS,
        mode="local_lesion",
        collapse_threshold=COLLAPSE_THRESHOLD,
        sustain_steps=SUSTAIN_STEPS,
    )
    collapse_level = k_star_level(summary)
    visible_level = first_visible_failure_level(rows)
    result = {
        "variation": "local_neighbor_strength_x_1_02",
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
        """# 3D Tissue Demo

This is Biology Line C.

It is an experimental application layer on top of frozen HAOS-IIP. It models a simple 3D tissue-like interaction structure as a deterministic cell graph and tests early recoverability loss under local lesion perturbation.

The model is a toy 8 x 8 x 8 multicellular lattice. It does not modify HAOS core and makes no claims about real tissue biology.

Run from the repository root:

```bash
python experiments/biology/tissue_3d_demo/run_tissue_3d_demo.py
```

Outputs are written to `experiments/biology/tissue_3d_demo/outputs/`.
""",
        encoding="utf-8",
    )


def write_validation(
    lattice: TissueLattice,
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
        f"""# 3D Tissue Demo Validation

## System
The system is a deterministic toy 3D tissue-like lattice with grid size {lattice.grid_size}, for {lattice.node_count} cells.

Each node is a cell with a scalar propagation state in [0, 1]. Edge classes are {", ".join(edge_classes)}. The local_neighbor class is 6-connected nearest-neighbor support, diagonal_neighbor provides weak local support, and long_range_signal provides sparse deterministic coordination links.

The primary perturbation is a local lesion centered at {lattice.lesion_center} with radius {lattice.lesion_radius}. Perturbation level p weakens edge support in that spherical region. A secondary signaling_weakening mode is implemented for long_range_signal edges but is not the primary validation sweep.

## Metrics
These are lightweight external proxy diagnostics aligned with HAOS-IIP stability logic. They are not frozen HAOS core metrics.

- recoverability: `0.15 * largest_component_fraction + 0.45 * weighted_degree_retention + 0.25 * spatial_coherence_retention + 0.15 * propagation_retention`, clipped to [0, 1].
- delta_persistence: change in recoverability between consecutive perturbation levels.
- k_star: first perturbation index where recoverability remains below {COLLAPSE_THRESHOLD:.2f} for the current sample plus {SUSTAIN_STEPS} following samples.
- safety_margin: distance in perturbation parameter from current p to p at k_star. At and beyond k_star this value is zero or negative.
- confidence: magnitude of the strongest negative delta_persistence observed so far.
- visible_failure: True when `largest_component_fraction < 0.85` or `spatial_coherence_retention < 0.50`.

## Result
Pass status: {"PASS" if pass_status else "FAIL"}

{check_lines}

- baseline_stable: {checks["baseline_stable"]}
- k_star perturbation level: {collapse_level}
- first visible failure level: {visible_level}
- early_detection: {checks["early_detection"]}
- robustness variation: {robustness["variation"]}
- robustness early_detection: {robustness["early_detection"]}
- minimum recoverability: {summary["min_recoverability"]:.6f}

## Limitations
This is a toy 3D model.

It is not real tissue simulation.

It is not morphogenesis.

It makes no biological claims.

It does not modify or prove HAOS-IIP core.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_readme()

    lattice, rows, summary = run_primary()
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
            "largest_component_fraction",
            "weighted_degree_retention",
            "spatial_coherence_retention",
            "propagation_retention",
            "visible_failure",
        ],
    )
    save_plots(lattice, rows)
    write_validation(lattice, rows, summary, checks, robustness_rows)

    print(f"baseline_stable: {checks['baseline_stable']}")
    print(f"k_star_level: {collapse_level}")
    print(f"first_visible_failure_level: {visible_level}")
    print(f"early_detection: {checks['early_detection']}")
    print(f"robustness_early_detection: {checks['robustness_early_detection']}")
    print("outputs_written: experiments/biology/tissue_3d_demo/outputs/")


if __name__ == "__main__":
    main()
