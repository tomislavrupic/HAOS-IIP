#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from microtubule_spectral_telemetry import SpectralTelemetryConfig, run_microtubule_spectral_telemetry


ROOT = Path(__file__).resolve().parent


def parse_float_list(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run microtubule spectral telemetry seed/noise sweep.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "spectral_outputs", help="Output directory.")
    parser.add_argument("--base-seed", type=int, default=20260429, help="First deterministic seed.")
    parser.add_argument("--seed-count", type=int, default=5, help="Number of seeds per noise level.")
    parser.add_argument("--thermal-noise", default="0.00,0.18,0.36", help="Comma-separated thermal noise levels.")
    parser.add_argument("--null-level", type=int, choices=(1, 2, 4, 5), default=4, help="Active null level.")
    parser.add_argument("--permutation-trials", type=int, default=16, help="Null trials per run.")
    parser.add_argument("--null-candidates", type=int, default=2, help="Candidate shuffles per null draw.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    noise_levels = parse_float_list(args.thermal_noise)
    for noise in noise_levels:
        for seed_idx in range(args.seed_count):
            seed = args.base_seed + seed_idx
            run_dir = args.output_dir / "sweep_runs" / f"noise_{noise:.3f}_seed_{seed}"
            config = SpectralTelemetryConfig(
                output_dir=run_dir,
                seed=seed,
                thermal_noise=noise,
                null_level=args.null_level,
                permutation_trials=args.permutation_trials,
                null_candidates=args.null_candidates,
            )
            outcome = run_microtubule_spectral_telemetry(config)
            observed = outcome["observed_summary"]
            rows.append(
                {
                    "thermal_noise": noise,
                    "seed": seed,
                    "bridge_status": outcome["bridge_status"],
                    "recoverability_score": observed["recoverability_score"],
                    "protofilament_identity_retention": observed["protofilament_identity_retention"],
                    "delta_persistence": observed["delta_persistence"],
                    "safety_margin": observed["safety_margin"],
                    "active_null_z": observed["active_null_z"],
                    "control_strict_pass_count": outcome["control_strict_pass_count"],
                }
            )
    summary = summarize(rows)
    write_csv(args.output_dir / "spectral_noise_sweep.csv", rows)
    write_summary(args.output_dir / "spectral_noise_sweep.md", rows, summary)
    write_plot(args.output_dir / "spectral_noise_sweep.png", rows)
    print(f"runs: {len(rows)}")
    print(f"noise_levels: {','.join(str(level) for level in noise_levels)}")
    print(f"pass_rate: {summary['pass_rate']:.6f}")
    print(f"outputs: {args.output_dir}")


def summarize(rows: list[dict[str, object]]) -> dict[str, float]:
    pass_count = sum(1 for row in rows if row["bridge_status"] == "PASS")
    return {
        "pass_rate": pass_count / max(len(rows), 1),
        "recoverability_mean": float(np.mean([float(row["recoverability_score"]) for row in rows])),
        "recoverability_std": float(np.std([float(row["recoverability_score"]) for row in rows])),
        "identity_mean": float(np.mean([float(row["protofilament_identity_retention"]) for row in rows])),
        "identity_std": float(np.std([float(row["protofilament_identity_retention"]) for row in rows])),
        "active_null_z_mean": float(np.mean([float(row["active_null_z"]) for row in rows])),
        "active_null_z_std": float(np.std([float(row["active_null_z"]) for row in rows])),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "thermal_noise",
        "seed",
        "bridge_status",
        "recoverability_score",
        "protofilament_identity_retention",
        "delta_persistence",
        "safety_margin",
        "active_null_z",
        "control_strict_pass_count",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, object]], summary: dict[str, float]) -> None:
    lines = [
        "# Microtubule Spectral Noise Sweep",
        "",
        f"- runs: {len(rows)}",
        f"- pass_rate: {summary['pass_rate']:.6f}",
        f"- recoverability_mean: {summary['recoverability_mean']:.6f}",
        f"- recoverability_std: {summary['recoverability_std']:.6f}",
        f"- protofilament_identity_mean: {summary['identity_mean']:.6f}",
        f"- protofilament_identity_std: {summary['identity_std']:.6f}",
        f"- active_null_z_mean: {summary['active_null_z_mean']:.6f}",
        f"- active_null_z_std: {summary['active_null_z_std']:.6f}",
        "",
        "| thermal_noise | seed | status | recoverability | protofilament_identity | active_null_z | controls_passed |",
        "| ---: | ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {noise:.3f} | {seed} | {status} | {rec:.6f} | {identity:.6f} | {z:.6f} | {controls} |".format(
                noise=float(row["thermal_noise"]),
                seed=int(row["seed"]),
                status=row["bridge_status"],
                rec=float(row["recoverability_score"]),
                identity=float(row["protofilament_identity_retention"]),
                z=float(row["active_null_z"]),
                controls=int(row["control_strict_pass_count"]),
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(path: Path, rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    noise_levels = sorted({float(row["thermal_noise"]) for row in rows})
    means = []
    stds = []
    identities = []
    for noise in noise_levels:
        subset = [row for row in rows if float(row["thermal_noise"]) == noise]
        values = np.array([float(row["recoverability_score"]) for row in subset], dtype=float)
        identity = np.array([float(row["protofilament_identity_retention"]) for row in subset], dtype=float)
        means.append(float(np.mean(values)))
        stds.append(float(np.std(values)))
        identities.append(float(np.mean(identity)))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.errorbar(noise_levels, means, yerr=stds, marker="o", label="recoverability")
    ax.plot(noise_levels, identities, marker="s", label="protofilament identity")
    ax.set_xlabel("thermal noise")
    ax.set_ylabel("score")
    ax.set_title("Microtubule spectral telemetry noise sweep")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
