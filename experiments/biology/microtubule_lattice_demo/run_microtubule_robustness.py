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

from microtubule_spectral_telemetry import SpectralTelemetryConfig, run_microtubule_spectral_telemetry


ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class RobustnessRow:
    thermal_noise: float
    seed: int
    bridge_status: str
    recoverability_score: float
    protofilament_identity_retention: float
    delta_persistence: float
    safety_margin: float
    branch_null_z: float
    spectral_null_z: float
    higher_null_z: float
    trajectory_null_z: float
    active_null_z: float
    control_strict_pass_count: int
    best_control_recoverability_score: float


def parse_float_list(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Phase 55.2 microtubule spectral robustness across seeds and thermal-noise levels."
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "spectral_outputs" / "robustness_55_2")
    parser.add_argument("--base-seed", type=int, default=20260429, help="First deterministic seed.")
    parser.add_argument("--seed-count", type=int, default=10, help="Number of seeds per noise level.")
    parser.add_argument("--thermal-noise", default="0.00,0.09,0.18,0.27,0.36", help="Comma-separated noise levels.")
    parser.add_argument("--null-level", type=int, choices=(1, 2, 4, 5), default=5, help="Active null level.")
    parser.add_argument("--permutation-trials", type=int, default=32, help="Null trials per run.")
    parser.add_argument("--null-candidates", type=int, default=4, help="Candidate shuffles per null draw.")
    parser.add_argument("--steps", type=int, default=96, help="Dynamics steps per run.")
    parser.add_argument("--damage-scale", type=float, default=0.72, help="Cap-region damage amplitude.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[RobustnessRow] = []
    noise_levels = parse_float_list(args.thermal_noise)

    for noise in noise_levels:
        for seed_idx in range(args.seed_count):
            seed = args.base_seed + seed_idx
            run_dir = args.output_dir / "runs" / f"noise_{noise:.3f}_seed_{seed}"
            config = SpectralTelemetryConfig(
                output_dir=run_dir,
                seed=seed,
                steps=args.steps,
                thermal_noise=noise,
                damage_scale=args.damage_scale,
                null_level=args.null_level,
                permutation_trials=args.permutation_trials,
                null_candidates=args.null_candidates,
            )
            outcome = run_microtubule_spectral_telemetry(config)
            rows.append(row_from_outcome(noise, seed, outcome))

    summary = summarize(rows)
    write_csv(args.output_dir / "robustness_runs.csv", rows)
    write_json(args.output_dir / "robustness_summary.json", summary)
    write_report(args.output_dir / "robustness_report.md", rows, summary, args)
    write_distribution_plots(args.output_dir, rows)
    write_visual_summary(args.output_dir / "robustness_visual_summary.png", rows, summary, args)

    print(f"runs: {len(rows)}")
    print(f"noise_levels: {','.join(str(level) for level in noise_levels)}")
    print(f"pass_rate: {summary['overall']['pass_rate']:.6f}")
    print(f"recoverability_mean: {summary['overall']['recoverability_score_mean']:.6f}")
    print(f"protofilament_identity_mean: {summary['overall']['protofilament_identity_retention_mean']:.6f}")
    print(f"active_null_z_mean: {summary['overall']['active_null_z_mean']:.6f}")
    print(f"outputs: {args.output_dir}")


def row_from_outcome(noise: float, seed: int, outcome: dict[str, Any]) -> RobustnessRow:
    observed = outcome["observed_summary"]
    return RobustnessRow(
        thermal_noise=noise,
        seed=seed,
        bridge_status=str(outcome["bridge_status"]),
        recoverability_score=float(observed["recoverability_score"]),
        protofilament_identity_retention=float(observed["protofilament_identity_retention"]),
        delta_persistence=float(observed["delta_persistence"]),
        safety_margin=float(observed["safety_margin"]),
        branch_null_z=float(observed["branch_null_z"]),
        spectral_null_z=float(observed["spectral_null_z"]),
        higher_null_z=float(observed["higher_null_z"]),
        trajectory_null_z=float(observed["trajectory_null_z"]),
        active_null_z=float(observed["active_null_z"]),
        control_strict_pass_count=int(outcome["control_strict_pass_count"]),
        best_control_recoverability_score=float(outcome["best_control_recoverability_score"]),
    )


def summarize(rows: list[RobustnessRow]) -> dict[str, Any]:
    by_noise: dict[str, Any] = {}
    grouped: dict[float, list[RobustnessRow]] = defaultdict(list)
    for row in rows:
        grouped[row.thermal_noise].append(row)
    for noise in sorted(grouped):
        by_noise[f"{noise:.3f}"] = summarize_group(grouped[noise])
    return {
        "overall": summarize_group(rows),
        "by_noise": by_noise,
    }


def summarize_group(rows: list[RobustnessRow]) -> dict[str, float | int]:
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
        "protofilament_identity_retention",
        "delta_persistence",
        "safety_margin",
        "branch_null_z",
        "spectral_null_z",
        "higher_null_z",
        "trajectory_null_z",
        "active_null_z",
        "best_control_recoverability_score",
    ):
        values = [float(getattr(row, name)) for row in rows]
        out[f"{name}_mean"] = mean(values)
        out[f"{name}_std"] = std(values)
        out[f"{name}_min"] = min(values) if values else 0.0
    return out


def write_csv(path: Path, rows: list[RobustnessRow]) -> None:
    fields = list(RobustnessRow.__dataclass_fields__)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: getattr(row, field) for field in fields})


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_report(path: Path, rows: list[RobustnessRow], summary: dict[str, Any], args: argparse.Namespace) -> None:
    overall = summary["overall"]
    lines = [
        "# Phase 55.2 Microtubule Spectral Robustness",
        "",
        "This run expands the single-pass microtubule telemetry result into a seed/noise robustness check.",
        "It remains a toy cylindrical-lattice telemetry probe, not a biological or molecular simulation.",
        "",
        "## Configuration",
        "",
        f"- seed_count: {args.seed_count}",
        f"- thermal_noise: {args.thermal_noise}",
        f"- null_level: {args.null_level}",
        f"- permutation_trials: {args.permutation_trials}",
        f"- null_candidates: {args.null_candidates}",
        f"- steps: {args.steps}",
        f"- damage_scale: {args.damage_scale}",
        "",
        "## Overall Summary",
        "",
        f"- runs: {overall['runs']}",
        f"- pass_rate: {overall['pass_rate']:.6f}",
        f"- recoverability_score: {overall['recoverability_score_mean']:.6f} +/- {overall['recoverability_score_std']:.6f}",
        f"- protofilament_identity_retention: {overall['protofilament_identity_retention_mean']:.6f} +/- {overall['protofilament_identity_retention_std']:.6f}",
        f"- active_null_z: {overall['active_null_z_mean']:.6f} +/- {overall['active_null_z_std']:.6f}",
        f"- control_pass_mean: {overall['control_pass_mean']:.6f}",
        f"- control_pass_max: {overall['control_pass_max']}",
        "",
        "## By Noise Level",
        "",
        "| thermal_noise | runs | pass_rate | recoverability | identity | delta_persistence | active_null_z | max_control_passes |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for noise, group in summary["by_noise"].items():
        lines.append(
            "| {noise} | {runs} | {pass_rate:.3f} | {rec:.6f} +/- {rec_std:.6f} | {identity:.6f} +/- {identity_std:.6f} | {delta:.6f} +/- {delta_std:.6f} | {z:.6f} +/- {z_std:.6f} | {max_pass} |".format(
                noise=noise,
                runs=group["runs"],
                pass_rate=group["pass_rate"],
                rec=group["recoverability_score_mean"],
                rec_std=group["recoverability_score_std"],
                identity=group["protofilament_identity_retention_mean"],
                identity_std=group["protofilament_identity_retention_std"],
                delta=group["delta_persistence_mean"],
                delta_std=group["delta_persistence_std"],
                z=group["active_null_z_mean"],
                z_std=group["active_null_z_std"],
                max_pass=group["control_pass_max"],
            )
        )
    lines.extend(
        [
            "",
            "## Raw Runs",
            "",
            "| noise | seed | status | recoverability | identity | active_null_z | control_passes |",
            "| ---: | ---: | --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in rows:
        lines.append(
            "| {noise:.3f} | {seed} | {status} | {rec:.6f} | {identity:.6f} | {z:.6f} | {controls} |".format(
                noise=row.thermal_noise,
                seed=row.seed,
                status=row.bridge_status,
                rec=row.recoverability_score,
                identity=row.protofilament_identity_retention,
                z=row.active_null_z,
                controls=row.control_strict_pass_count,
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_distribution_plots(output_dir: Path, rows: list[RobustnessRow]) -> None:
    import matplotlib.pyplot as plt

    noise_levels = sorted({row.thermal_noise for row in rows})
    metrics = [
        ("recoverability_score", "Recoverability"),
        ("protofilament_identity_retention", "Protofilament Identity"),
        ("active_null_z", "Active Null z"),
        ("delta_persistence", "Delta Persistence"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, (field, title) in zip(axes.ravel(), metrics):
        data = [[float(getattr(row, field)) for row in rows if row.thermal_noise == noise] for noise in noise_levels]
        ax.boxplot(data, tick_labels=[f"{noise:.2f}" for noise in noise_levels], showmeans=True)
        ax.set_title(title)
        ax.set_xlabel("thermal noise")
        ax.grid(True, alpha=0.25)
    fig.suptitle("Phase 55.2 microtubule robustness distributions")
    fig.tight_layout()
    fig.savefig(output_dir / "robustness_distributions.png", dpi=170)
    plt.close(fig)


def write_visual_summary(path: Path, rows: list[RobustnessRow], summary: dict[str, Any], args: argparse.Namespace) -> None:
    import matplotlib.pyplot as plt

    noise_levels = sorted({row.thermal_noise for row in rows})
    rec_mean = [summary["by_noise"][f"{noise:.3f}"]["recoverability_score_mean"] for noise in noise_levels]
    rec_std = [summary["by_noise"][f"{noise:.3f}"]["recoverability_score_std"] for noise in noise_levels]
    ident_mean = [summary["by_noise"][f"{noise:.3f}"]["protofilament_identity_retention_mean"] for noise in noise_levels]
    ident_std = [summary["by_noise"][f"{noise:.3f}"]["protofilament_identity_retention_std"] for noise in noise_levels]
    pass_rates = [summary["by_noise"][f"{noise:.3f}"]["pass_rate"] for noise in noise_levels]
    z_values = [[row.active_null_z for row in rows if row.thermal_noise == noise] for noise in noise_levels]

    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.9])
    ax_scores = fig.add_subplot(gs[0, 0])
    ax_null = fig.add_subplot(gs[0, 1])
    ax_pass = fig.add_subplot(gs[1, 0])
    ax_table = fig.add_subplot(gs[1, 1])

    ax_scores.errorbar(noise_levels, rec_mean, yerr=rec_std, marker="o", label="recoverability")
    ax_scores.errorbar(noise_levels, ident_mean, yerr=ident_std, marker="s", label="protofilament identity")
    ax_scores.set_title("Robustness under thermal-noise proxy")
    ax_scores.set_xlabel("thermal noise")
    ax_scores.set_ylabel("score")
    ax_scores.set_ylim(0.0, 1.05)
    ax_scores.grid(True, alpha=0.25)
    ax_scores.legend()

    ax_null.boxplot(z_values, tick_labels=[f"{noise:.2f}" for noise in noise_levels], showmeans=True)
    ax_null.axhline(1.5, color="tab:red", linestyle="--", linewidth=1, label="strict z gate")
    ax_null.set_title("Level-5 active-null specificity")
    ax_null.set_xlabel("thermal noise")
    ax_null.set_ylabel("active null z")
    ax_null.grid(True, alpha=0.25)
    ax_null.legend()

    ax_pass.bar([f"{noise:.2f}" for noise in noise_levels], pass_rates, color="tab:green")
    ax_pass.set_title("PASS rate by noise level")
    ax_pass.set_xlabel("thermal noise")
    ax_pass.set_ylabel("PASS rate")
    ax_pass.set_ylim(0.0, 1.05)
    ax_pass.grid(True, axis="y", alpha=0.25)

    overall = summary["overall"]
    ax_table.axis("off")
    table_rows = [
        ["runs", f"{overall['runs']}"],
        ["pass rate", f"{overall['pass_rate']:.3f}"],
        ["recoverability", f"{overall['recoverability_score_mean']:.3f} +/- {overall['recoverability_score_std']:.3f}"],
        ["protofilament identity", f"{overall['protofilament_identity_retention_mean']:.3f} +/- {overall['protofilament_identity_retention_std']:.3f}"],
        ["active null z", f"{overall['active_null_z_mean']:.2f} +/- {overall['active_null_z_std']:.2f}"],
        ["max controls passing", f"{overall['control_pass_max']}"],
        ["null level", f"{args.null_level}"],
    ]
    table = ax_table.table(cellText=table_rows, colLabels=["Metric", "Value"], loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.35)
    ax_table.set_title("Phase 55.2 summary")

    fig.suptitle("HAOS-IIP Phase 55.2: Microtubule Spectral Robustness", fontsize=14)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def mean(values: list[float | int]) -> float:
    return float(np.mean(values)) if values else 0.0


def std(values: list[float | int]) -> float:
    return float(np.std(values)) if values else 0.0


if __name__ == "__main__":
    main()
