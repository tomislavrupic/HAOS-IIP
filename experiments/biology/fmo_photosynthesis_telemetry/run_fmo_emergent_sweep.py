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
class EmergentSpec:
    name: str
    address_gain: float
    sink_gain: float
    directed_bias_gain: float
    temporal_flux_gain: float
    environment_assist_gain: float = 0.0
    spectral_modes: int = 4


SPECS = [
    EmergentSpec("weak_directed", 0.44, 0.04, 0.12, 0.04),
    EmergentSpec("directed_sink", 0.42, 0.10, 0.14, 0.06),
    EmergentSpec("temporal_directed", 0.42, 0.06, 0.10, 0.14),
    EmergentSpec("assisted_directed", 0.40, 0.08, 0.12, 0.10, environment_assist_gain=0.06),
    EmergentSpec("modes3_directed", 0.42, 0.08, 0.12, 0.10, spectral_modes=3),
    EmergentSpec("modes5_directed", 0.42, 0.08, 0.12, 0.10, spectral_modes=5),
]


def parse_float_list(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 56.4 FMO intrinsic pathway dynamics sweep.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "robustness_56_4_emergent_pathway")
    parser.add_argument("--base-seed", type=int, default=20260429)
    parser.add_argument("--seed-count", type=int, default=10)
    parser.add_argument("--thermal-noise", default="0.00,0.04,0.08,0.12,0.16")
    parser.add_argument("--permutation-trials", type=int, default=32)
    parser.add_argument("--null-candidates", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    noise_levels = parse_float_list(args.thermal_noise)
    for spec in SPECS:
        for noise in noise_levels:
            for seed_idx in range(args.seed_count):
                seed = args.base_seed + seed_idx
                run_dir = args.output_dir / "runs" / spec.name / f"noise_{noise:.3f}_seed_{seed}"
                config = FMOTelemetryConfig(
                    output_dir=run_dir,
                    seed=seed,
                    address_mode="intrinsic_pathway",
                    address_gain=spec.address_gain,
                    sink_gain=spec.sink_gain,
                    directed_bias_gain=spec.directed_bias_gain,
                    temporal_flux_gain=spec.temporal_flux_gain,
                    environment_assist_gain=spec.environment_assist_gain,
                    spectral_modes=spec.spectral_modes,
                    thermal_noise=noise,
                    permutation_trials=args.permutation_trials,
                    null_candidates=args.null_candidates,
                    null_level=6,
                )
                outcome = run_fmo_spectral_telemetry(config)
                observed = outcome["observed_summary"]
                rows.append(
                    {
                        "spec": spec.name,
                        "thermal_noise": noise,
                        "seed": seed,
                        "bridge_status": outcome["bridge_status"],
                        "recoverability_score": observed["recoverability_score"],
                        "site_identity_retention": observed["site_identity_retention"],
                        "pathway_identity_retention": observed["pathway_identity_retention"],
                        "temporal_pathway_identity_retention": observed["temporal_pathway_identity_retention"],
                        "delta_persistence": observed["delta_persistence"],
                        "active_null_z": observed["active_null_z"],
                        "control_strict_pass_count": outcome["control_strict_pass_count"],
                    }
                )
    summary = summarize(rows)
    write_csv(args.output_dir / "fmo_emergent_runs.csv", rows)
    write_json(args.output_dir / "fmo_emergent_summary.json", summary)
    write_report(args.output_dir / "fmo_emergent_report.md", summary)
    write_plot(args.output_dir / "fmo_emergent_summary.png", summary)
    best = max(summary, key=lambda item: (item["pass_rate"], item["pathway_identity_retention_mean"], -item["control_pass_mean"]))
    print(f"specs: {len(summary)}")
    print(f"runs: {len(rows)}")
    print(f"best_spec: {best['spec']}")
    print(f"best_pass_rate: {best['pass_rate']:.6f}")
    print(f"best_pathway_identity_mean: {best['pathway_identity_retention_mean']:.6f}")
    print(f"best_control_pass_mean: {best['control_pass_mean']:.6f}")
    print(f"outputs: {args.output_dir}")


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for spec in sorted({str(row["spec"]) for row in rows}):
        subset = [row for row in rows if row["spec"] == spec]
        pass_count = sum(1 for row in subset if row["bridge_status"] == "PASS")
        control_passes = [int(row["control_strict_pass_count"]) for row in subset]
        item: dict[str, Any] = {
            "spec": spec,
            "runs": len(subset),
            "pass_rate": pass_count / max(len(subset), 1),
            "control_pass_max": max(control_passes) if control_passes else 0,
            "control_pass_mean": float(np.mean(control_passes)) if control_passes else 0.0,
        }
        for field in (
            "recoverability_score",
            "site_identity_retention",
            "pathway_identity_retention",
            "temporal_pathway_identity_retention",
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
        "spec",
        "thermal_noise",
        "seed",
        "bridge_status",
        "recoverability_score",
        "site_identity_retention",
        "pathway_identity_retention",
        "temporal_pathway_identity_retention",
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
        "# Phase 56.4 FMO Intrinsic Pathway Dynamics Sweep",
        "",
        "This sweep uses weak intrinsic directed/temporal pathway bias and the level-6 directed trajectory null.",
        "It is a diagnostic refinement, not a biological validation claim.",
        "",
        "| spec | runs | pass_rate | recoverability | pathway_identity | temporal_pathway | active_null_z | mean_controls | max_controls |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for item in sorted(summary, key=lambda row: (row["pass_rate"], row["pathway_identity_retention_mean"]), reverse=True):
        lines.append(
            "| {spec} | {runs} | {pass_rate:.3f} | {rec:.6f} +/- {rec_std:.6f} | {path:.6f} +/- {path_std:.6f} | {temporal:.6f} +/- {temporal_std:.6f} | {z:.6f} +/- {z_std:.6f} | {mean_controls:.3f} | {max_controls} |".format(
                spec=item["spec"],
                runs=item["runs"],
                pass_rate=item["pass_rate"],
                rec=item["recoverability_score_mean"],
                rec_std=item["recoverability_score_std"],
                path=item["pathway_identity_retention_mean"],
                path_std=item["pathway_identity_retention_std"],
                temporal=item["temporal_pathway_identity_retention_mean"],
                temporal_std=item["temporal_pathway_identity_retention_std"],
                z=item["active_null_z_mean"],
                z_std=item["active_null_z_std"],
                mean_controls=item["control_pass_mean"],
                max_controls=item["control_pass_max"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(path: Path, summary: list[dict[str, Any]]) -> None:
    import matplotlib.pyplot as plt

    ordered = sorted(summary, key=lambda row: (row["pass_rate"], row["pathway_identity_retention_mean"]), reverse=True)
    labels = [row["spec"] for row in ordered]
    pathway = [row["pathway_identity_retention_mean"] for row in ordered]
    pass_rate = [row["pass_rate"] for row in ordered]
    controls = [row["control_pass_mean"] for row in ordered]
    temporal = [row["temporal_pathway_identity_retention_mean"] for row in ordered]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    axes[0].bar(labels, pathway)
    axes[0].axhline(0.78, color="tab:red", linestyle="--", linewidth=1)
    axes[0].set_ylabel("pathway identity mean")
    axes[1].bar(labels, temporal, color="tab:cyan")
    axes[1].axhline(0.70, color="tab:red", linestyle="--", linewidth=1)
    axes[1].set_ylabel("temporal pathway mean")
    axes[2].bar(labels, pass_rate, color="tab:green")
    axes[2].axhline(0.50, color="tab:red", linestyle="--", linewidth=1)
    axes[2].set_ylabel("pass rate")
    axes[3].bar(labels, controls, color="tab:orange")
    axes[3].axhline(1.0, color="tab:red", linestyle="--", linewidth=1)
    axes[3].set_ylabel("mean controls passing")
    for ax in axes:
        ax.tick_params(axis="x", rotation=35)
    fig.suptitle("Phase 56.4 FMO intrinsic pathway dynamics")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
