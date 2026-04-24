#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gene_network_model import (
    GeneNetwork,
    apply_edge_weakening,
    build_default_network,
    compute_effective_graph,
    compute_recoverability_metrics,
    final_window_mean,
    first_sustained_collapse,
    normalized_distance,
    simulate_network,
)


EXPERIMENT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
RESULTS_CSV = OUTPUT_DIR / "results.csv"
VALIDATION_MD = EXPERIMENT_DIR / "gene_network_validation.md"

COLLAPSE_THRESHOLD = 0.65
SUSTAIN_STEPS = 2
VISIBLE_FAILURE_THRESHOLD = 0.39
FINAL_WINDOW = 30
PERTURBATION_LEVELS = np.linspace(0.0, 1.0, 31)


def baseline_stable(trajectory: np.ndarray, *, window: int = 30, tolerance: float = 1.0e-6) -> bool:
    tail = trajectory[-window:]
    max_step_change = float(np.max(np.abs(np.diff(tail, axis=0))))
    return max_step_change < tolerance


def branch_recoverability(
    baseline_trajectory: np.ndarray,
    perturbed_trajectories: tuple[np.ndarray, ...],
    network: GeneNetwork,
) -> list[float]:
    idx = [network.index[gene] for gene in network.fragile_branch]
    baseline_state = final_window_mean(baseline_trajectory, FINAL_WINDOW)[idx]
    values: list[float] = []
    for trajectory in perturbed_trajectories:
        state = final_window_mean(trajectory, FINAL_WINDOW)[idx]
        values.append(float(np.clip(1.0 - normalized_distance(state, baseline_state), 0.0, 1.0)))
    return values


def run_primary_sweep() -> dict[str, object]:
    network = build_default_network()
    baseline_trajectory = simulate_network(network)
    perturbed_networks = tuple(apply_edge_weakening(network, p) for p in PERTURBATION_LEVELS)
    trajectories = tuple(simulate_network(item) for item in perturbed_networks)
    effective_graphs = tuple(
        compute_effective_graph(item, trajectory, final_window=FINAL_WINDOW)
        for item, trajectory in zip(perturbed_networks, trajectories, strict=True)
    )
    rows, summary = compute_recoverability_metrics(
        baseline_trajectory,
        trajectories,
        PERTURBATION_LEVELS,
        collapse_threshold=COLLAPSE_THRESHOLD,
        sustain_steps=SUSTAIN_STEPS,
        visible_failure_threshold=VISIBLE_FAILURE_THRESHOLD,
        final_window=FINAL_WINDOW,
    )
    branch_values = branch_recoverability(baseline_trajectory, trajectories, network)
    branch_k_star = first_sustained_collapse(
        branch_values,
        collapse_threshold=COLLAPSE_THRESHOLD,
        sustain_steps=SUSTAIN_STEPS,
    )
    return {
        "network": network,
        "baseline_trajectory": baseline_trajectory,
        "trajectories": trajectories,
        "effective_graphs": effective_graphs,
        "rows": rows,
        "summary": summary,
        "branch_recoverability": branch_values,
        "branch_k_star": branch_k_star,
    }


def first_visible_failure_level(rows: list[dict[str, object]]) -> float | None:
    for row in rows:
        if row["visible_failure"]:
            return float(row["perturbation_level"])
    return None


def recoverability_declines_gradually(rows: list[dict[str, object]]) -> bool:
    recoverability = [float(row["recoverability"]) for row in rows]
    deltas = np.diff(recoverability)
    total_decline = recoverability[0] - recoverability[-1]
    upward_bump = float(np.max(deltas)) if deltas.size else 0.0
    largest_drop = abs(float(np.min(deltas))) if deltas.size else 0.0
    return total_decline > 0.35 and upward_bump <= 0.01 and largest_drop < 0.18


def deterministic_repeat(first: dict[str, object], second: dict[str, object]) -> bool:
    first_rows = first["rows"]
    second_rows = second["rows"]
    assert isinstance(first_rows, list)
    assert isinstance(second_rows, list)
    for left, right in zip(first_rows, second_rows, strict=True):
        for key in ("perturbation_level", "recoverability", "delta_persistence", "safety_margin", "confidence"):
            if not np.isclose(float(left[key]), float(right[key]), atol=1.0e-12, rtol=0.0):
                return False
        if left["k_star"] != right["k_star"] or left["visible_failure"] != right["visible_failure"]:
            return False
    return True


def write_results_csv(rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "perturbation_index",
        "perturbation_level",
        "recoverability",
        "delta_persistence",
        "k_star",
        "safety_margin",
        "confidence",
        "visible_failure",
    ]
    with RESULTS_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_gene_activity(path: Path, trajectory: np.ndarray, network: GeneNetwork, title: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    for idx, gene in enumerate(network.genes):
        linewidth = 2.2 if gene in network.fragile_branch or gene == network.hub_gene else 1.0
        alpha = 0.95 if gene in network.fragile_branch or gene == network.hub_gene else 0.55
        ax.plot(trajectory[:, idx], label=gene, linewidth=linewidth, alpha=alpha)
    ax.set_title(title)
    ax.set_xlabel("time step")
    ax.set_ylabel("gene activity")
    ax.set_ylim(0.0, 1.0)
    ax.legend(ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_metric(
    path: Path,
    rows: list[dict[str, object]],
    field: str,
    ylabel: str,
    title: str,
    *,
    threshold: float | None = None,
) -> None:
    x = [float(row["perturbation_level"]) for row in rows]
    y = [float(row[field]) for row in rows]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x, y, marker="o", linewidth=2.0)
    if threshold is not None:
        ax.axhline(threshold, color="black", linestyle="--", linewidth=1.0)
    ax.set_title(title)
    ax.set_xlabel("perturbation level")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_plots(result: dict[str, object]) -> None:
    network = result["network"]
    baseline_trajectory = result["baseline_trajectory"]
    trajectories = result["trajectories"]
    rows = result["rows"]
    summary = result["summary"]
    assert isinstance(network, GeneNetwork)
    assert isinstance(baseline_trajectory, np.ndarray)
    assert isinstance(trajectories, tuple)
    assert isinstance(rows, list)
    assert isinstance(summary, dict)

    k_star = summary["k_star"]
    perturbed_index = int(k_star) if k_star is not None else len(trajectories) - 1
    perturbed_level = PERTURBATION_LEVELS[perturbed_index]
    plot_gene_activity(OUTPUT_DIR / "gene_activity_baseline.png", baseline_trajectory, network, "Baseline gene activity")
    plot_gene_activity(
        OUTPUT_DIR / "gene_activity_perturbed.png",
        trajectories[perturbed_index],
        network,
        f"Perturbed gene activity, p={perturbed_level:.3f}",
    )
    plot_metric(
        OUTPUT_DIR / "recoverability_vs_perturbation.png",
        rows,
        "recoverability",
        "recoverability",
        "Recoverability vs perturbation",
        threshold=COLLAPSE_THRESHOLD,
    )
    plot_metric(
        OUTPUT_DIR / "delta_persistence_vs_perturbation.png",
        rows,
        "delta_persistence",
        "delta persistence",
        "Delta persistence vs perturbation",
        threshold=0.0,
    )
    plot_metric(
        OUTPUT_DIR / "safety_margin_vs_perturbation.png",
        rows,
        "safety_margin",
        "safety margin",
        "Safety margin vs perturbation",
        threshold=0.0,
    )


def write_readme() -> None:
    readme = EXPERIMENT_DIR / "README.md"
    readme.write_text(
        """# Gene Network Demo

This is the first HAOS-IIP biology-line experiment.

It is an experimental application layer only. It does not modify frozen HAOS-IIP core code, theory, telemetry, phase artifacts, authority bundles, canonical documents, paper files, or already frozen experiment outputs.

The demo tests whether lightweight HAOS-style stability diagnostics can detect early degradation in a synthetic gene regulatory network under controlled perturbation. It does not claim to explain biology, validate a biological theory, or prove HAOS-IIP core behavior.

Run from the repository root:

```bash
python experiments/biology/gene_network_demo/run_gene_network_demo.py
```

Primary outputs are written to `experiments/biology/gene_network_demo/outputs/`.
""",
        encoding="utf-8",
    )


def write_validation(result: dict[str, object], checks: dict[str, bool]) -> None:
    network = result["network"]
    rows = result["rows"]
    summary = result["summary"]
    branch_k_star = result["branch_k_star"]
    branch_values = result["branch_recoverability"]
    assert isinstance(network, GeneNetwork)
    assert isinstance(rows, list)
    assert isinstance(summary, dict)
    assert isinstance(branch_values, list)

    k_star = summary["k_star"]
    k_star_level = None if k_star is None else float(PERTURBATION_LEVELS[int(k_star)])
    visible_level = first_visible_failure_level(rows)
    branch_k_star_level = None if branch_k_star is None else float(PERTURBATION_LEVELS[int(branch_k_star)])
    pass_status = all(checks.values())
    if checks["fragile_branch_degrades_earlier"]:
        origin_sentence = (
            "The collapse origin for this toy run is the fragile downstream branch, "
            "because its branch-level recoverability crosses the sustained-collapse "
            "threshold before the whole-network recoverability does."
        )
    else:
        origin_sentence = (
            "The collapse origin condition was not established by this run; the fragile "
            "branch did not show sustained threshold crossing before the whole network."
        )
    check_lines = "\n".join(
        f"- {name}: {'PASS' if passed else 'FAIL'}" for name, passed in checks.items()
    )

    VALIDATION_MD.write_text(
        f"""# Gene Network Demo Validation

## System
The system is a deterministic toy gene regulatory network with {len(network.genes)} genes: {", ".join(network.genes)}.

G0 is the hub regulator. The network includes a negative feedback loop where G1 activates G2 and G2 inhibits G1 and G0. It includes a positive feedback motif where G3 and G4 mutually activate. The fragile downstream branch is {", ".join(network.fragile_branch)}; its support depends on weakened inputs from G0 to G8, G4 to G8, and G5 to G9.

## Perturbation
The primary v0.1 experiment uses edge weakening. The perturbation parameter p is swept from 0.0 to 1.0 in {len(PERTURBATION_LEVELS)} deterministic samples. At p = 0.0 no fragile-branch support edges are weakened. At p = 1.0 those selected support edges are fully removed.

Hub knockdown is implemented in `gene_network_model.py` as an available perturbation mode, but it is not the primary validation path for this artifact.

## Metrics
These are lightweight external proxy diagnostics aligned with HAOS-IIP stability logic. They are not frozen HAOS core metrics.

- recoverability: `1 - normalized_distance(perturbed_final_state, baseline_final_state)`, clipped to [0, 1]. The final state is represented by the mean activity over the final {FINAL_WINDOW} time steps.
- delta_persistence: change in recoverability between consecutive perturbation levels.
- k_star: first perturbation index where recoverability stays below {COLLAPSE_THRESHOLD:.2f} for the current sample plus {SUSTAIN_STEPS} following samples.
- safety_margin: distance in perturbation parameter from current p to p at k_star. At and beyond k_star this value is zero or negative.
- confidence: magnitude of the strongest negative delta_persistence signal.

The effective graph uses `effective_weight(source, target) = abs(W[source, target]) * activity_factor`, where `activity_factor` is source-gene mean activity over the final trajectory window.

## Result
Pass status: {"PASS" if pass_status else "FAIL"}

{check_lines}

- baseline_stable: {checks["baseline_stable"]}
- k_star_level: {k_star_level}
- first_visible_failure_level: {visible_level}
- early_detection: {checks["early_detection"]}
- fragile_branch_k_star_level: {branch_k_star_level}
- minimum_fragile_branch_recoverability: {min(branch_values):.6f}
- minimum_whole_network_recoverability: {summary["min_recoverability"]:.6f}

{origin_sentence}

## Limitations
This is a toy model.

It is not real biological validation.

It is not a full systems biology model.

It does not modify or prove HAOS-IIP core.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_readme()

    result = run_primary_sweep()
    repeat = run_primary_sweep()
    rows = result["rows"]
    summary = result["summary"]
    baseline_trajectory = result["baseline_trajectory"]
    branch_k_star = result["branch_k_star"]
    assert isinstance(rows, list)
    assert isinstance(summary, dict)
    assert isinstance(baseline_trajectory, np.ndarray)

    k_star = summary["k_star"]
    visible_level = first_visible_failure_level(rows)
    k_star_level = None if k_star is None else float(PERTURBATION_LEVELS[int(k_star)])
    checks = {
        "baseline_stable": baseline_stable(baseline_trajectory),
        "recoverability_declines_gradually": recoverability_declines_gradually(rows),
        "early_detection": k_star_level is not None and visible_level is not None and k_star_level < visible_level,
        "deterministic_repeated_runs": deterministic_repeat(result, repeat),
        "fragile_branch_degrades_earlier": branch_k_star is not None and k_star is not None and int(branch_k_star) < int(k_star),
    }

    write_results_csv(rows)
    save_plots(result)
    write_validation(result, checks)

    print(f"baseline_stable: {checks['baseline_stable']}")
    print(f"k_star_level: {k_star_level}")
    print(f"first_visible_failure_level: {visible_level}")
    print(f"early_detection: {checks['early_detection']}")
    print("outputs_written: experiments/biology/gene_network_demo/outputs/")


if __name__ == "__main__":
    main()
