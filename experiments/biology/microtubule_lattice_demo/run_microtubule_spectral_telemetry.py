#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from microtubule_spectral_telemetry import SpectralTelemetryConfig, run_microtubule_spectral_telemetry


ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run microtubule spectral telemetry prototype.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "spectral_outputs", help="Output directory.")
    parser.add_argument("--seed", type=int, default=20260429, help="Deterministic seed.")
    parser.add_argument("--protofilaments", type=int, default=13, help="Protofilament count.")
    parser.add_argument("--dimers", type=int, default=30, help="Dimers per protofilament.")
    parser.add_argument("--steps", type=int, default=96, help="Dynamics steps.")
    parser.add_argument("--thermal-noise", type=float, default=0.18, help="Small per-step thermal noise proxy.")
    parser.add_argument("--damage-scale", type=float, default=0.72, help="GTP-cap patch damage amplitude.")
    parser.add_argument("--null-level", type=int, choices=(1, 2, 4, 5), default=4, help="Active null level.")
    parser.add_argument("--permutation-trials", type=int, default=32, help="Null permutation trials.")
    parser.add_argument("--null-candidates", type=int, default=4, help="Candidate shuffles per null draw.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = SpectralTelemetryConfig(
        output_dir=args.output_dir,
        seed=args.seed,
        protofilaments=args.protofilaments,
        dimers_per_protofilament=args.dimers,
        steps=args.steps,
        thermal_noise=args.thermal_noise,
        damage_scale=args.damage_scale,
        null_level=args.null_level,
        permutation_trials=args.permutation_trials,
        null_candidates=args.null_candidates,
    )
    outcome = run_microtubule_spectral_telemetry(config)
    print(f"bridge_status: {outcome['bridge_status']}")
    print(f"failure_mode: {outcome['failure_mode']}")
    print(f"recoverability_score: {outcome['observed_summary']['recoverability_score']}")
    print(f"protofilament_identity_retention: {outcome['observed_summary']['protofilament_identity_retention']}")
    print(f"active_null_z: {outcome['observed_summary']['active_null_z']}")
    print(f"control_strict_pass_count: {outcome['control_strict_pass_count']}")
    print(f"outputs: {config.output_dir}")


if __name__ == "__main__":
    main()
