#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from fmo_spectral_telemetry import FMOTelemetryConfig, ROOT, run_fmo_spectral_telemetry


@dataclass(frozen=True)
class VariantSpec:
    name: str
    address_mode: str
    address_gain: float = 0.48
    local_address_gain: float = 0.18
    sink_gain: float = 0.0
    environment_assist_gain: float = 0.0


VARIANTS = [
    VariantSpec("spectral_baseline", "spectral"),
    VariantSpec("local_only", "local", address_gain=0.0, local_address_gain=0.36),
    VariantSpec("hybrid", "hybrid", address_gain=0.42, local_address_gain=0.24),
    VariantSpec("sink_spectral", "sink", address_gain=0.42, sink_gain=0.22),
    VariantSpec("hybrid_sink", "hybrid", address_gain=0.40, local_address_gain=0.22, sink_gain=0.22),
    VariantSpec("environment_assisted", "environment_assisted", address_gain=0.40, local_address_gain=0.20, sink_gain=0.18, environment_assist_gain=0.16),
]


def parse_float_list(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 56.2 FMO pathway-strengthening variant sweep.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "robustness_56_2_variants")
    parser.add_argument("--base-seed", type=int, default=20260429)
    parser.add_argument("--seed-count", type=int, default=10)
    parser.add_argument("--thermal-noise", default="0.00,0.04,0.08,0.12,0.16")
    parser.add_argument("--permutation-trials", type=int, default=32)
    parser.add_argument("--null-candidates", type=int, default=8)
    parser.add_argument("--null-level", type=int, choices=(1, 2, 4, 5), default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    noise_levels = parse_float_list(args.thermal_noise)
    for variant in VARIANTS:
        for noise in noise_levels:
            for seed_idx in range(args.seed_count):
                seed = args.base_seed + seed_idx
                run_dir = args.output_dir / "runs" / variant.name / f"noise_{noise:.3f}_seed_{seed}"
                config = FMOTelemetryConfig(
                    output_dir=run_dir,
                    seed=seed,
                    address_mode=variant.address_mode,
                    address_gain=variant.address_gain,
                    local_address_gain=variant.local_address_gain,
                    sink_gain=variant.sink_gain,
                    environment_assist_gain=variant.environment_assist_gain,
                    thermal_noise=noise,
                    null_level=args.null_level,
                    permutation_trials=args.permutation_trials,
                    null_candidates=args.null_candidates,
                )
                outcome = run_fmo_spectral_telemetry(config)
                observed = outcome["observed_summary"]
                rows.append(
                    {
                        "variant": variant.name,
                        "thermal_noise": noise,
                        "seed": seed,
                        "bridge_status": outcome["bridge_status"],
                        "recoverability_score": observed["recoverability_score"],
                        "site_identity_retention": observed["site_identity_retention"],
                        "pathway_identity_retention": observed["pathway_identity_retention"],
                        "delta_persistence": observed["delta_persistence"],
                        "active_null_z": observed["active_null_z"],
                        "control_strict_pass_count": outcome["control_strict_pass_count"],
                    }
                )
    summary = summarize(rows)
    write_csv(args.output_dir / "fmo_variant_runs.csv", rows)
    write_json(args.output_dir / "fmo_variant_summary.json", summary)
    write_report(args.output_dir / "fmo_variant_report.md", summary)
    write_plot(args.output_dir / "fmo_variant_summary.png", summary)
    best = max(summary, key=lambda item: (item["pass_rate"], item["pathway_identity_retention_mean"], item["active_null_z_mean"]))
    print(f"variants: {len(summary)}")
    print(f"runs: {len(rows)}")
    print(f"best_variant: {best['variant']}")
    print(f"best_pass_rate: {best['pass_rate']:.6f}")
    print(f"best_pathway_identity_mean: {best['pathway_identity_retention_mean']:.6f}")
    print(f"outputs: {args.output_dir}")


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for variant in sorted({str(row["variant"]) for row in rows}):
        subset = [row for row in rows if row["variant"] == variant]
        pass_count = sum(1 for row in subset if row["bridge_status"] == "PASS")
        control_passes = [int(row["control_strict_pass_count"]) for row in subset]
        item: dict[str, Any] = {
            "variant": variant,
            "runs": len(subset),
            "pass_rate": pass_count / max(len(subset), 1),
            "control_pass_max": max(control_passes) if control_passes else 0,
            "control_pass_mean": float(np.mean(control_passes)) if control_passes else 0.0,
        }
        for field in (
            "recoverability_score",
            "site_identity_retention",
            "pathway_identity_retention",
            "delta_persistence",
            "active_null_z",
        ):
            values = np.asarray([float(row[field]) for row in subset], dtype=float)
            item[f"{field}_mean"] = float(np.mean(values))
            item[f"{field}_std"] = float(np.std(values))
            item[f"{field}_min"] = float(np.min(values))
        out.append(item)
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "variant",
        "thermal_noise",
        "seed",
        "bridge_status",
        "recoverability_score",
        "site_identity_retention",
        "pathway_identity_retention",
        "delta_persistence",
        "active_null_z",
        "control_strict_pass_count",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, data: list[dict[str, Any]]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_report(path: Path, summary: list[dict[str, Any]]) -> None:
    lines = [
        "# Phase 56.2 FMO Pathway-Strengthening Variant Sweep",
        "",
        "This sweep tests bounded address-rule variants after the Phase 56.1 FMO falsifier.",
        "It should be read as a diagnostic sweep, not as biological validation.",
        "",
        "| variant | runs | pass_rate | recoverability | site_identity | pathway_identity | active_null_z | max_controls |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for item in sorted(summary, key=lambda row: row["pathway_identity_retention_mean"], reverse=True):
        lines.append(
            "| {variant} | {runs} | {pass_rate:.3f} | {rec:.6f} +/- {rec_std:.6f} | {site:.6f} +/- {site_std:.6f} | {path:.6f} +/- {path_std:.6f} | {z:.6f} +/- {z_std:.6f} | {controls} |".format(
                variant=item["variant"],
                runs=item["runs"],
                pass_rate=item["pass_rate"],
                rec=item["recoverability_score_mean"],
                rec_std=item["recoverability_score_std"],
                site=item["site_identity_retention_mean"],
                site_std=item["site_identity_retention_std"],
                path=item["pathway_identity_retention_mean"],
                path_std=item["pathway_identity_retention_std"],
                z=item["active_null_z_mean"],
                z_std=item["active_null_z_std"],
                controls=item["control_pass_max"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(path: Path, summary: list[dict[str, Any]]) -> None:
    import matplotlib.pyplot as plt

    ordered = sorted(summary, key=lambda row: row["pathway_identity_retention_mean"], reverse=True)
    labels = [row["variant"] for row in ordered]
    pathway = [row["pathway_identity_retention_mean"] for row in ordered]
    pass_rate = [row["pass_rate"] for row in ordered]
    z_values = [row["active_null_z_mean"] for row in ordered]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].bar(labels, pathway)
    axes[0].axhline(0.75, color="tab:red", linestyle="--", linewidth=1)
    axes[0].set_ylabel("pathway identity mean")
    axes[0].tick_params(axis="x", rotation=35)
    axes[1].bar(labels, pass_rate, color="tab:green")
    axes[1].axhline(0.60, color="tab:red", linestyle="--", linewidth=1)
    axes[1].set_ylabel("pass rate")
    axes[1].tick_params(axis="x", rotation=35)
    axes[2].bar(labels, z_values, color="tab:purple")
    axes[2].axhline(1.5, color="tab:red", linestyle="--", linewidth=1)
    axes[2].set_ylabel("active null z mean")
    axes[2].tick_params(axis="x", rotation=35)
    fig.suptitle("Phase 56.2 FMO variant sweep")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
