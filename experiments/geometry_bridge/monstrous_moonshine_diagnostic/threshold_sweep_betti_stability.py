#!/usr/bin/env python3
"""Threshold-sweep and null-ensemble hardening for Betti_0 sidecar.

This script tests whether the selected Betti_0 result at threshold 8.0 is
isolated or sits inside a local robustness band. It also runs deterministic
support-replacement impostors to estimate how often the same Betti_0 and a
nearby edge count appear under matched null support.

It is still a topological diagnostic sidecar only: no Lean theorem, no
Moonshine proof, and no physical bridge.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run_betti_component_count import (  # noqa: E402
    BettiComponentConfig,
    SUPERSINGULAR_PRIMES,
    base_nodes,
    canonical_edge_signature,
    component_count,
    condition_nodes,
    jaccard_distance,
    relation_edges,
    stable_hash,
    stable_unit,
    write_csv,
    write_json,
)


SWEEP_CSV_PATH = ROOT / "betti_threshold_sweep.csv"
NULL_CSV_PATH = ROOT / "betti_null_ensemble.csv"
REPORT_PATH = ROOT / "betti_threshold_sweep_report.md"
RESULT_PATH = ROOT / "betti_threshold_sweep_result.json"

SWEEP_FIELDS = [
    "threshold",
    "edge_count",
    "betti0",
    "component_size_distribution",
    "largest_component_size",
    "edge_jaccard_distance_to_reference",
    "status_region",
]
NULL_FIELDS = [
    "seed",
    "replacement_count",
    "support",
    "component_count",
    "edge_count",
    "edge_count_delta",
    "edge_jaccard_distance_to_reference",
    "false_positive",
]

NONSUPERSINGULAR_PRIME_POOL = (37, 43, 53, 61, 67, 73, 79)


@dataclass(frozen=True)
class BettiSweepConfig:
    version: str = "betti-threshold-sweep-v0.1"
    threshold_min: float = 2.0
    threshold_max: float = 12.0
    threshold_step: float = 0.5
    reference_threshold: float = 8.0
    edge_count_tolerance: int = 2
    null_seed_count: int = 100
    null_min_replacements: int = 1
    null_max_replacements: int = 7


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Betti threshold-sweep hardening.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def threshold_values(config: BettiSweepConfig) -> list[float]:
    count = int(round((config.threshold_max - config.threshold_min) / config.threshold_step))
    return [round(config.threshold_min + index * config.threshold_step, 10) for index in range(count + 1)]


def connected_component_sizes(node_count: int, edges: list[tuple[int, int, float]]) -> list[int]:
    neighbors: list[list[int]] = [[] for _ in range(node_count)]
    for i, j, _ in edges:
        neighbors[i].append(j)
        neighbors[j].append(i)
    seen = [False] * node_count
    sizes: list[int] = []
    for start in range(node_count):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            for neighbor in neighbors[node]:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    stack.append(neighbor)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def status_region(
    threshold: float,
    edge_distance: float,
    edge_count: int,
    betti0: int,
    reference_edge_count: int,
    reference_betti0: int,
    config: BettiSweepConfig,
) -> str:
    if abs(threshold - config.reference_threshold) < 1.0e-12:
        return "REFERENCE_THRESHOLD"
    if betti0 != reference_betti0:
        return "BETTI_SHIFT"
    if edge_distance <= 1.0e-12:
        return "BETTI_AND_EDGE_STABLE"
    if abs(edge_count - reference_edge_count) <= config.edge_count_tolerance:
        return "BETTI_STABLE_EDGE_NEIGHBORHOOD"
    return "BETTI_STABLE_EDGE_DRIFT"


def sweep_rows(config: BettiSweepConfig) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    reference_nodes = condition_nodes("known_positive")
    reference_edges = relation_edges(reference_nodes, config.reference_threshold)
    reference_signature = canonical_edge_signature(reference_nodes, reference_edges)
    reference_betti0 = component_count(len(reference_nodes), reference_edges)
    reference_edge_count = len(reference_edges)
    rows: list[dict[str, Any]] = []
    stable_thresholds: list[float] = []
    exact_edge_thresholds: list[float] = []
    neighborhood_thresholds: list[float] = []

    for threshold in threshold_values(config):
        edges = relation_edges(reference_nodes, threshold)
        betti0 = component_count(len(reference_nodes), edges)
        sizes = connected_component_sizes(len(reference_nodes), edges)
        edge_signature = canonical_edge_signature(reference_nodes, edges)
        edge_distance = jaccard_distance(reference_signature, edge_signature)
        region = status_region(
            threshold,
            edge_distance,
            len(edges),
            betti0,
            reference_edge_count,
            reference_betti0,
            config,
        )
        if betti0 == reference_betti0:
            stable_thresholds.append(threshold)
        if edge_distance <= 1.0e-12:
            exact_edge_thresholds.append(threshold)
        if betti0 == reference_betti0 and abs(len(edges) - reference_edge_count) <= config.edge_count_tolerance:
            neighborhood_thresholds.append(threshold)
        rows.append(
            {
                "threshold": f"{threshold:.1f}",
                "edge_count": len(edges),
                "betti0": betti0,
                "component_size_distribution": json.dumps(sizes, separators=(",", ":")),
                "largest_component_size": sizes[0] if sizes else 0,
                "edge_jaccard_distance_to_reference": f"{edge_distance:.12g}",
                "status_region": region,
            }
        )

    summary = {
        "reference_betti0": reference_betti0,
        "reference_edge_count": reference_edge_count,
        "reference_threshold": config.reference_threshold,
        "betti0_stable_thresholds": stable_thresholds,
        "edge_exact_thresholds": exact_edge_thresholds,
        "edge_neighborhood_thresholds": neighborhood_thresholds,
        "betti0_stability_band": contiguous_band(stable_thresholds, config.reference_threshold),
        "edge_exact_band": contiguous_band(exact_edge_thresholds, config.reference_threshold),
        "edge_neighborhood_band": contiguous_band(neighborhood_thresholds, config.reference_threshold),
    }
    return rows, summary


def contiguous_band(values: list[float], anchor: float) -> list[float] | None:
    if anchor not in values:
        return None
    value_set = set(values)
    ordered = sorted(values)
    index = ordered.index(anchor)
    left = anchor
    right = anchor
    cursor = index - 1
    while cursor >= 0 and ordered[cursor] in value_set:
        if abs(ordered[cursor] - left) > 0.5000001:
            break
        left = ordered[cursor]
        cursor -= 1
    cursor = index + 1
    while cursor < len(ordered) and ordered[cursor] in value_set:
        if abs(ordered[cursor] - right) > 0.5000001:
            break
        right = ordered[cursor]
        cursor += 1
    return [left, right]


def deterministic_impostor_support(seed: int, config: BettiSweepConfig) -> list[int]:
    support = list(SUPERSINGULAR_PRIMES)
    replacement_count = config.null_min_replacements + (
        int(stable_unit("null_replacement_count", seed) * 10_000) % (config.null_max_replacements - config.null_min_replacements + 1)
    )
    positions = sorted(range(len(support)), key=lambda idx: stable_unit("null_position", seed, idx))[:replacement_count]
    replacement_order = sorted(NONSUPERSINGULAR_PRIME_POOL, key=lambda prime: stable_unit("null_prime", seed, prime))
    for index, position in enumerate(positions):
        support[position] = replacement_order[index % len(replacement_order)]
    return support


def null_rows(config: BettiSweepConfig, reference_summary: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    reference_nodes = condition_nodes("known_positive")
    reference_edges = relation_edges(reference_nodes, config.reference_threshold)
    reference_signature = canonical_edge_signature(reference_nodes, reference_edges)
    reference_betti0 = int(reference_summary["reference_betti0"])
    reference_edge_count = int(reference_summary["reference_edge_count"])
    rows: list[dict[str, Any]] = []
    false_positive_count = 0

    for seed in range(config.null_seed_count):
        support = deterministic_impostor_support(seed, config)
        nodes = base_nodes(support)
        edges = relation_edges(nodes, config.reference_threshold)
        betti0 = component_count(len(nodes), edges)
        edge_count = len(edges)
        edge_delta = abs(edge_count - reference_edge_count)
        edge_distance = jaccard_distance(reference_signature, canonical_edge_signature(nodes, edges))
        false_positive = betti0 == reference_betti0 and edge_delta <= config.edge_count_tolerance
        if false_positive:
            false_positive_count += 1
        rows.append(
            {
                "seed": seed,
                "replacement_count": len([prime for prime in support if prime not in SUPERSINGULAR_PRIMES]),
                "support": json.dumps(support, separators=(",", ":")),
                "component_count": betti0,
                "edge_count": edge_count,
                "edge_count_delta": edge_delta,
                "edge_jaccard_distance_to_reference": f"{edge_distance:.12g}",
                "false_positive": str(false_positive).lower(),
            }
        )

    summary = {
        "null_seed_count": config.null_seed_count,
        "false_positive_count": false_positive_count,
        "false_positive_rate": false_positive_count / config.null_seed_count if config.null_seed_count else 0.0,
        "false_positive_definition": "component_count equals reference Betti_0 and edge_count is within the configured edge-count tolerance",
        "nonsupersingular_prime_pool": list(NONSUPERSINGULAR_PRIME_POOL),
    }
    return rows, summary


def run_threshold_sweep(config: BettiSweepConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    sweep, sweep_summary = sweep_rows(config)
    nulls, null_summary = null_rows(config, sweep_summary)
    labels = [
        "BETTI_THRESHOLD_SWEEP_AVAILABLE",
        "BETTI_NULL_ENSEMBLE_AVAILABLE",
        "LEAN_THEOREM_NOT_INCLUDED",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        "CLAIM_GATED_TOPOLOGICAL_DIAGNOSTIC",
    ]
    if sweep_summary["betti0_stability_band"] is not None:
        labels.append("BETTI0_ROBUSTNESS_BAND_DETECTED")
    if null_summary["false_positive_rate"] > 0.0:
        labels.append("NULL_FALSE_POSITIVES_DETECTED")
    else:
        labels.append("NULL_FALSE_POSITIVES_NOT_DETECTED")

    result = {
        "bridge_id": "betti_threshold_sweep_v0_1",
        "status": "PASS",
        "classification": "TOPOLOGICAL_DIAGNOSTIC_SIDECAR",
        "labels": labels,
        "config": asdict(config),
        "sweep_summary": sweep_summary,
        "null_summary": null_summary,
        "claim_ceiling": "Threshold robustness hardening only; no Moonshine proof, local Lean theorem, or physical bridge claim.",
    }
    result["result_hash"] = stable_hash("betti_sweep", result)

    write_csv(output_dir / SWEEP_CSV_PATH.name, sweep, SWEEP_FIELDS)
    write_csv(output_dir / NULL_CSV_PATH.name, nulls, NULL_FIELDS)
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result)
    return result


def write_report(path: Path, result: dict[str, Any]) -> None:
    sweep = result["sweep_summary"]
    nulls = result["null_summary"]
    lines = [
        "# Betti Threshold-Sweep Report",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        f"- Reference threshold: `{sweep['reference_threshold']}`",
        f"- Reference `Betti_0`: `{sweep['reference_betti0']}`",
        f"- Reference edges: `{sweep['reference_edge_count']}`",
        f"- `Betti_0` stability band: `{sweep['betti0_stability_band']}`",
        f"- exact edge-signature band: `{sweep['edge_exact_band']}`",
        f"- edge-neighborhood band: `{sweep['edge_neighborhood_band']}`",
        f"- null false-positive rate: `{nulls['false_positive_rate']:.6f}`",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Boundary",
        "",
        "This sweep hardens the declared Betti_0 relation graph against threshold and null-support sensitivity.",
        "It does not prove Monstrous Moonshine, verify the external Lean SVG claim, or establish any physical bridge.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_and_run() -> None:
    args = parse_args()
    result = run_threshold_sweep(BettiSweepConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    parse_and_run()

