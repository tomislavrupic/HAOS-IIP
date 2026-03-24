#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = EXAMPLES_DIR / "output"
EXPECTED_TABLE = EXAMPLES_DIR / "expected_output.csv"
EXPECTED_PLOT = EXAMPLES_DIR / "expected_plot.svg"
GENERATED_TABLE = OUTPUT_DIR / "quick_reproduce_output.csv"
GENERATED_PLOT = OUTPUT_DIR / "quick_reproduce_plot.svg"
CONTINUUM_TABLE = ROOT / "continuum-sketch" / "continuum_sketch_table.csv"
CONTINUUM_PLOT = ROOT / "continuum-sketch" / "continuum_sketch_plot_family.svg"

PHASES = ("15", "16", "17", "18")
BRANCH = "frozen_branch"
CONTROL = "periodic_diagonal_augmented_control"


def run_python(args: list[str]) -> str:
    result = subprocess.run(
        [sys.executable, *args],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise SystemExit(result.stdout + result.stderr)
    return result.stdout.strip()


def run_phase_checks() -> dict[str, dict[str, object]]:
    outputs: dict[str, dict[str, object]] = {}
    for phase in PHASES:
        raw = run_python(["run_phase.py", phase, "--check"])
        outputs[phase] = json.loads(raw)
        if not outputs[phase].get("success"):
            raise SystemExit(f"Phase {phase} checker failed:\n{json.dumps(outputs[phase], indent=2)}")
    return outputs


def load_json(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def std(values: list[float]) -> float:
    if not values:
        return 0.0
    mu = mean(values)
    return (sum((value - mu) ** 2 for value in values) / len(values)) ** 0.5


def coeff_var(values: list[float]) -> float:
    mu = mean(values)
    return 0.0 if abs(mu) < 1.0e-12 else std(values) / abs(mu)


def max_adjacent_diff(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return max(abs(values[index] - values[index - 1]) for index in range(1, len(values)))


def select_series(rows: list[dict[str, str]], hierarchy_label: str, observable_name: str) -> list[float]:
    filtered = [
        row
        for row in rows
        if row["hierarchy_label"] == hierarchy_label and row["observable_name"] == observable_name
    ]
    filtered.sort(key=lambda row: int(row["n_side"]))
    return [float(row["value"]) for row in filtered]


def continuum_metrics() -> dict[str, float]:
    rows = read_csv(CONTINUUM_TABLE)
    branch_speed = select_series(rows, BRANCH, "effective_speed_band")
    control_speed = select_series(rows, CONTROL, "effective_speed_band")
    branch_shell = select_series(rows, BRANCH, "distance_surrogate_shell_slope")
    control_shell = select_series(rows, CONTROL, "distance_surrogate_shell_slope")
    return {
        "branch_speed_cv": coeff_var(branch_speed),
        "control_speed_cv": coeff_var(control_speed),
        "branch_shell_drift": max_adjacent_diff(branch_shell),
        "control_shell_drift": max_adjacent_diff(control_shell),
    }


def summary_rows() -> list[dict[str, str]]:
    phase15 = load_json(ROOT / "phase15-propagation" / "phase15_manifest.json")
    phase16 = load_json(ROOT / "phase16-temporal-ordering" / "phase16_manifest.json")
    phase17 = load_json(ROOT / "phase17-causal-closure" / "phase17_manifest.json")
    phase18 = load_json(ROOT / "phase18-distance-surrogate" / "phase18_manifest.json")
    continuum = continuum_metrics()

    rows = [
        {
            "narrative_step": "propagation_band",
            "source_artifact": "phase15-propagation/phase15_manifest.json",
            "branch_metric": "branch_effective_speed_drift",
            "control_metric": "control_effective_speed_drift",
            "branch_value": f"{phase15['aggregate_metrics']['branch_effective_speed_drift']:.6f}",
            "control_value": f"{phase15['aggregate_metrics']['control_effective_speed_drift']:.6f}",
            "expected_relation": "branch < control",
        },
        {
            "narrative_step": "temporal_ordering",
            "source_artifact": "phase16-temporal-ordering/phase16_manifest.json",
            "branch_metric": "branch_primary_robustness_distance",
            "control_metric": "control_primary_robustness_distance",
            "branch_value": f"{phase16['aggregate_metrics']['branch_primary_robustness_distance']:.6f}",
            "control_value": f"{phase16['aggregate_metrics']['control_primary_robustness_distance']:.6f}",
            "expected_relation": "branch < control",
        },
        {
            "narrative_step": "causal_closure",
            "source_artifact": "phase17-causal-closure/phase17_manifest.json",
            "branch_metric": "branch_causal_depth_drift",
            "control_metric": "control_causal_depth_drift",
            "branch_value": f"{phase17['aggregate_metrics']['branch_causal_depth_drift']:.6f}",
            "control_value": f"{phase17['aggregate_metrics']['control_causal_depth_drift']:.6f}",
            "expected_relation": "branch < control",
        },
        {
            "narrative_step": "distance_surrogate",
            "source_artifact": "phase18-distance-surrogate/phase18_manifest.json",
            "branch_metric": "branch_max_shell_overlap_fraction",
            "control_metric": "control_max_shell_overlap_fraction",
            "branch_value": f"{phase18['branch_metrics']['max_shell_overlap_fraction']:.6f}",
            "control_value": f"{phase18['control_metrics']['max_shell_overlap_fraction']:.6f}",
            "expected_relation": "branch < control",
        },
        {
            "narrative_step": "continuum_speed_band",
            "source_artifact": "continuum-sketch/continuum_sketch_table.csv",
            "branch_metric": "branch_speed_coefficient_of_variation",
            "control_metric": "control_speed_coefficient_of_variation",
            "branch_value": f"{continuum['branch_speed_cv']:.6f}",
            "control_value": f"{continuum['control_speed_cv']:.6f}",
            "expected_relation": "branch < control",
        },
        {
            "narrative_step": "continuum_shell_slope",
            "source_artifact": "continuum-sketch/continuum_sketch_table.csv",
            "branch_metric": "branch_shell_slope_max_adjacent_drift",
            "control_metric": "control_shell_slope_max_adjacent_drift",
            "branch_value": f"{continuum['branch_shell_drift']:.6f}",
            "control_value": f"{continuum['control_shell_drift']:.6f}",
            "expected_relation": "branch < control",
        },
    ]

    for row in rows:
        branch_value = float(row["branch_value"])
        control_value = float(row["control_value"])
        row["status"] = "pass" if branch_value < control_value else "fail"
    return rows


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "narrative_step",
                "source_artifact",
                "branch_metric",
                "control_metric",
                "branch_value",
                "control_value",
                "expected_relation",
                "status",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def compare_exact(path: Path, expected_path: Path) -> bool:
    return path.read_text(encoding="utf-8") == expected_path.read_text(encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    run_phase_checks()
    run_python(["continuum-sketch/build_continuum_sketch.py"])

    rows = summary_rows()
    if any(row["status"] != "pass" for row in rows):
        raise SystemExit("Generated summary row failed its branch/control relation check.")

    write_summary(GENERATED_TABLE, rows)
    shutil.copyfile(CONTINUUM_PLOT, GENERATED_PLOT)

    table_matches = compare_exact(GENERATED_TABLE, EXPECTED_TABLE)
    plot_matches = compare_exact(GENERATED_PLOT, EXPECTED_PLOT)
    if not table_matches or not plot_matches:
        raise SystemExit(
            "Frozen baseline mismatch.\n"
            f"table_matches={table_matches}\n"
            f"plot_matches={plot_matches}\n"
            f"generated_table={GENERATED_TABLE}\n"
            f"generated_plot={GENERATED_PLOT}"
        )

    print("Reproduction spine passed.")
    print(f"Table: {GENERATED_TABLE}")
    print(f"Plot: {GENERATED_PLOT}")
    print(
        "Statement: Frozen propagation -> ordering -> causality -> distance-surrogate "
        "reproduces from stored artifacts, and the artifact-only continuum sketch stays "
        "bounded and branch-distinct. No new dynamics or continuum ontology are introduced."
    )


if __name__ == "__main__":
    main()
