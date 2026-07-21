#!/usr/bin/env python3
"""Render descriptive tables from the frozen RP-01 validation rows."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "validation/validation_runs.csv"
PRIMARY = {
    "localized_block_shock", "sparse_node_corruption", "relational_twist",
    "clustered_symbol_corruption", "dispersed_equal_weight_corruption",
    "topology_preserving_identity_disruption",
}
FIELDS = (
    "recovered", "parity_decoding_success", "functional_restoration",
    "identity_preservation", "state_region_recovery_gain",
    "trajectory_corridor_recovery_gain", "operator_restoration",
    "corrective_cost", "wrong_codeword_convergence",
)


def numeric(value: str) -> float:
    if value == "True":
        return 1.0
    if value == "False":
        return 0.0
    return float(value)


def aggregate(rows: list[dict[str, str]], key: str) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row[key]].append(row)
    output = []
    for name, group in sorted(grouped.items()):
        item: dict[str, object] = {key: name, "run_count": len(group)}
        for field in FIELDS:
            values = [numeric(row[field]) for row in group if row[field] not in {"", "None"}]
            item[f"median_{field}"] = median(values) if values else None
        output.append(item)
    return output


def write(path: Path, rows: list[dict[str, object]]) -> None:
    keys = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    primary = [row for row in rows if row["perturbation_family"] in PRIMARY]
    target = [row for row in primary if row["condition"] == "target"]
    write(ROOT / "validation/control_results.csv", aggregate(primary, "condition"))
    write(ROOT / "validation/radius_results.csv", aggregate(target, "radius_class"))
    write(ROOT / "validation/perturbation_results.csv", aggregate(target, "perturbation_family"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

