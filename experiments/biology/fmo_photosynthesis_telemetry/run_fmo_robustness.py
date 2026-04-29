#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from fmo_spectral_telemetry import FMOTelemetryConfig, ROOT, run_fmo_spectral_telemetry


@dataclass(frozen=True)
class FMORobustnessRow:
    variant: str
    thermal_noise: float
    seed: int
    bridge_status: str
    recoverability_score: float
    site_identity_retention: float
    pathway_identity_retention: float
    delta_persistence: float
    safety_margin: float
    active_null_z: float
    control_strict_pass_count: int
    best_control_recoverability_score: float


def parse_float_list(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 56 FMO spectral robustness across seeds and noise.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "robustness_56")
    parser.add_argument("--base-seed", type=int, default=20260429)
    parser.add_argument("--seed-count", type=int, default=20)
    parser.add_argument("--thermal-noise", default="0.00,0.04,0.08,0.12,0.16")
    parser.add_argument("--variant", default="spectral", help="Named label for this robustness run.")
    parser.add_argument(
        "--address-mode",
        choices=("spectral", "local", "hybrid", "sink", "pathway_flux", "intrinsic_pathway", "environment_assisted"),
        default="spectral",
    )
    parser.add_argument("--address-gain", type=float, default=0.48)
    parser.add_argument("--local-address-gain", type=float, default=0.18)
    parser.add_argument("--sink-gain", type=float, default=0.0)
    parser.add_argument("--flux-gain", type=float, default=0.0)
    parser.add_argument("--directed-bias-gain", type=float, default=0.0)
    parser.add_argument("--temporal-flux-gain", type=float, default=0.0)
    parser.add_argument("--environment-assist-gain", type=float, default=0.0)
    parser.add_argument("--disorder-scale", type=float, default=0.18)
    parser.add_argument("--damage-scale", type=float, default=0.42)
    parser.add_argument("--null-level", type=int, choices=(1, 2, 4, 5, 6), default=5)
    parser.add_argument("--permutation-trials", type=int, default=32)
    parser.add_argument("--null-candidates", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[FMORobustnessRow] = []
    noise_levels = parse_float_list(args.thermal_noise)
    for noise in noise_levels:
        for seed_idx in range(args.seed_count):
            seed = args.base_seed + seed_idx
            run_dir = args.output_dir / "runs" / f"noise_{noise:.3f}_seed_{seed}"
            config = FMOTelemetryConfig(
                output_dir=run_dir,
                seed=seed,
                address_mode=args.address_mode,
                address_gain=args.address_gain,
                local_address_gain=args.local_address_gain,
                sink_gain=args.sink_gain,
                flux_gain=args.flux_gain,
                directed_bias_gain=args.directed_bias_gain,
                temporal_flux_gain=args.temporal_flux_gain,
                environment_assist_gain=args.environment_assist_gain,
                thermal_noise=noise,
                disorder_scale=args.disorder_scale,
                damage_scale=args.damage_scale,
                null_level=args.null_level,
                permutation_trials=args.permutation_trials,
                null_candidates=args.null_candidates,
            )
            rows.append(row_from_outcome(args.variant, noise, seed, run_fmo_spectral_telemetry(config)))
    summary = summarize(rows)
    write_csv(args.output_dir / "fmo_robustness_runs.csv", rows)
    write_json(args.output_dir / "fmo_robustness_summary.json", summary)
    write_report(args.output_dir / "fmo_robustness_report.md", rows, summary, args)
    write_plot(args.output_dir / "fmo_robustness_summary.png", rows, summary)
    print(f"runs: {len(rows)}")
    print(f"variant: {args.variant}")
    print(f"noise_levels: {','.join(str(level) for level in noise_levels)}")
    print(f"pass_rate: {summary['overall']['pass_rate']:.6f}")
    print(f"recoverability_mean: {summary['overall']['recoverability_score_mean']:.6f}")
    print(f"site_identity_mean: {summary['overall']['site_identity_retention_mean']:.6f}")
    print(f"pathway_identity_mean: {summary['overall']['pathway_identity_retention_mean']:.6f}")
    print(f"active_null_z_mean: {summary['overall']['active_null_z_mean']:.6f}")
    print(f"outputs: {args.output_dir}")


def row_from_outcome(variant: str, noise: float, seed: int, outcome: dict[str, Any]) -> FMORobustnessRow:
    observed = outcome["observed_summary"]
    return FMORobustnessRow(
        variant=variant,
        thermal_noise=noise,
        seed=seed,
        bridge_status=str(outcome["bridge_status"]),
        recoverability_score=float(observed["recoverability_score"]),
        site_identity_retention=float(observed["site_identity_retention"]),
        pathway_identity_retention=float(observed["pathway_identity_retention"]),
        delta_persistence=float(observed["delta_persistence"]),
        safety_margin=float(observed["safety_margin"]),
        active_null_z=float(observed["active_null_z"]),
        control_strict_pass_count=int(outcome["control_strict_pass_count"]),
        best_control_recoverability_score=float(outcome["best_control_recoverability_score"]),
    )


def summarize(rows: list[FMORobustnessRow]) -> dict[str, Any]:
    grouped: dict[float, list[FMORobustnessRow]] = defaultdict(list)
    for row in rows:
        grouped[row.thermal_noise].append(row)
    return {
        "overall": summarize_group(rows),
        "by_noise": {f"{noise:.3f}": summarize_group(grouped[noise]) for noise in sorted(grouped)},
    }


def summarize_group(rows: list[FMORobustnessRow]) -> dict[str, float | int]:
    pass_count = sum(1 for row in rows if row.bridge_status == "PASS")
    out: dict[str, float | int] = {
        "runs": len(rows),
        "pass_count": pass_count,
        "pass_rate": pass_count / max(len(rows), 1),
        "control_pass_mean": mean([row.control_strict_pass_count for row in rows]),
        "control_pass_max": max([row.control_strict_pass_count for row in rows], default=0),
    }
    for name in (
        "recoverability_score",
        "site_identity_retention",
        "pathway_identity_retention",
        "delta_persistence",
        "safety_margin",
        "active_null_z",
        "best_control_recoverability_score",
    ):
        values = [float(getattr(row, name)) for row in rows]
        out[f"{name}_mean"] = mean(values)
        out[f"{name}_std"] = std(values)
        out[f"{name}_min"] = min(values) if values else 0.0
    return out


def write_csv(path: Path, rows: list[FMORobustnessRow]) -> None:
    fields = list(FMORobustnessRow.__dataclass_fields__)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: getattr(row, field) for field in fields})


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_report(path: Path, rows: list[FMORobustnessRow], summary: dict[str, Any], args: argparse.Namespace) -> None:
    overall = summary["overall"]
    lines = [
        "# Phase 56 FMO Photosynthesis Telemetry Robustness",
        "",
        "Toy FMO-like 7-site weighted-network telemetry. This is not a molecular simulation or a quantum-biology claim.",
        "",
        "## Configuration",
        "",
        f"- seed_count: {args.seed_count}",
        f"- variant: {args.variant}",
        f"- address_mode: {args.address_mode}",
        f"- address_gain: {args.address_gain}",
        f"- local_address_gain: {args.local_address_gain}",
        f"- sink_gain: {args.sink_gain}",
        f"- flux_gain: {args.flux_gain}",
        f"- directed_bias_gain: {args.directed_bias_gain}",
        f"- temporal_flux_gain: {args.temporal_flux_gain}",
        f"- environment_assist_gain: {args.environment_assist_gain}",
        f"- thermal_noise: {args.thermal_noise}",
        f"- disorder_scale: {args.disorder_scale}",
        f"- damage_scale: {args.damage_scale}",
        f"- null_level: {args.null_level}",
        f"- permutation_trials: {args.permutation_trials}",
        f"- null_candidates: {args.null_candidates}",
        "",
        "## Overall Summary",
        "",
        f"- runs: {overall['runs']}",
        f"- pass_rate: {overall['pass_rate']:.6f}",
        f"- recoverability_score: {overall['recoverability_score_mean']:.6f} +/- {overall['recoverability_score_std']:.6f}",
        f"- site_identity_retention: {overall['site_identity_retention_mean']:.6f} +/- {overall['site_identity_retention_std']:.6f}",
        f"- pathway_identity_retention: {overall['pathway_identity_retention_mean']:.6f} +/- {overall['pathway_identity_retention_std']:.6f}",
        f"- active_null_z: {overall['active_null_z_mean']:.6f} +/- {overall['active_null_z_std']:.6f}",
        f"- control_pass_max: {overall['control_pass_max']}",
        "",
        "## By Noise Level",
        "",
        "| thermal_noise | runs | pass_rate | recoverability | site_identity | pathway_identity | active_null_z | max_control_passes |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for noise, group in summary["by_noise"].items():
        lines.append(
            "| {noise} | {runs} | {pass_rate:.3f} | {rec:.6f} +/- {rec_std:.6f} | {site:.6f} +/- {site_std:.6f} | {pathway:.6f} +/- {pathway_std:.6f} | {z:.6f} +/- {z_std:.6f} | {max_pass} |".format(
                noise=noise,
                runs=group["runs"],
                pass_rate=group["pass_rate"],
                rec=group["recoverability_score_mean"],
                rec_std=group["recoverability_score_std"],
                site=group["site_identity_retention_mean"],
                site_std=group["site_identity_retention_std"],
                pathway=group["pathway_identity_retention_mean"],
                pathway_std=group["pathway_identity_retention_std"],
                z=group["active_null_z_mean"],
                z_std=group["active_null_z_std"],
                max_pass=group["control_pass_max"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(path: Path, rows: list[FMORobustnessRow], summary: dict[str, Any]) -> None:
    import matplotlib.pyplot as plt

    noise_levels = sorted({row.thermal_noise for row in rows})
    rec = [summary["by_noise"][f"{noise:.3f}"]["recoverability_score_mean"] for noise in noise_levels]
    site = [summary["by_noise"][f"{noise:.3f}"]["site_identity_retention_mean"] for noise in noise_levels]
    pathway = [summary["by_noise"][f"{noise:.3f}"]["pathway_identity_retention_mean"] for noise in noise_levels]
    z_data = [[row.active_null_z for row in rows if row.thermal_noise == noise] for noise in noise_levels]
    pass_rates = [summary["by_noise"][f"{noise:.3f}"]["pass_rate"] for noise in noise_levels]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    axes[0, 0].plot(noise_levels, rec, marker="o", label="recoverability")
    axes[0, 0].plot(noise_levels, site, marker="s", label="site identity")
    axes[0, 0].plot(noise_levels, pathway, marker="^", label="pathway identity")
    axes[0, 0].set_ylim(0.0, 1.05)
    axes[0, 0].set_xlabel("thermal noise")
    axes[0, 0].set_ylabel("score")
    axes[0, 0].set_title("FMO telemetry score stability")
    axes[0, 0].grid(True, alpha=0.25)
    axes[0, 0].legend()

    axes[0, 1].boxplot(z_data, tick_labels=[f"{noise:.2f}" for noise in noise_levels], showmeans=True)
    axes[0, 1].axhline(1.5, color="tab:red", linestyle="--", linewidth=1, label="strict z gate")
    axes[0, 1].set_title("Level-5 active-null specificity")
    axes[0, 1].set_xlabel("thermal noise")
    axes[0, 1].set_ylabel("active null z")
    axes[0, 1].grid(True, alpha=0.25)
    axes[0, 1].legend()

    axes[1, 0].bar([f"{noise:.2f}" for noise in noise_levels], pass_rates, color="tab:green")
    axes[1, 0].set_title("PASS rate by noise level")
    axes[1, 0].set_ylim(0.0, 1.05)
    axes[1, 0].set_xlabel("thermal noise")
    axes[1, 0].set_ylabel("PASS rate")
    axes[1, 0].grid(True, axis="y", alpha=0.25)

    overall = summary["overall"]
    axes[1, 1].axis("off")
    table_rows = [
        ["runs", f"{overall['runs']}"],
        ["pass rate", f"{overall['pass_rate']:.3f}"],
        ["recoverability", f"{overall['recoverability_score_mean']:.3f} +/- {overall['recoverability_score_std']:.3f}"],
        ["site identity", f"{overall['site_identity_retention_mean']:.3f} +/- {overall['site_identity_retention_std']:.3f}"],
        ["pathway identity", f"{overall['pathway_identity_retention_mean']:.3f} +/- {overall['pathway_identity_retention_std']:.3f}"],
        ["active null z", f"{overall['active_null_z_mean']:.2f} +/- {overall['active_null_z_std']:.2f}"],
        ["max controls passing", f"{overall['control_pass_max']}"],
    ]
    table = axes[1, 1].table(cellText=table_rows, colLabels=["Metric", "Value"], loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.25)
    axes[1, 1].set_title("Phase 56 summary")
    fig.suptitle("HAOS-IIP Phase 56: FMO Photosynthesis Telemetry", fontsize=14)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def mean(values: list[float | int]) -> float:
    return float(np.mean(values)) if values else 0.0


def std(values: list[float | int]) -> float:
    return float(np.std(values)) if values else 0.0


if __name__ == "__main__":
    main()
