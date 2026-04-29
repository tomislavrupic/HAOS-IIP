#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from fmo_spectral_telemetry import FMOTelemetryConfig, ROOT, run_fmo_spectral_telemetry


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 56 FMO spectral telemetry prototype.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--seed", type=int, default=20260429)
    parser.add_argument("--steps", type=int, default=128)
    parser.add_argument("--address-mode", choices=("spectral", "local", "hybrid", "sink", "environment_assisted"), default="spectral")
    parser.add_argument("--address-gain", type=float, default=0.48)
    parser.add_argument("--local-address-gain", type=float, default=0.18)
    parser.add_argument("--sink-gain", type=float, default=0.0)
    parser.add_argument("--environment-assist-gain", type=float, default=0.0)
    parser.add_argument("--thermal-noise", type=float, default=0.08)
    parser.add_argument("--disorder-scale", type=float, default=0.18)
    parser.add_argument("--damage-scale", type=float, default=0.42)
    parser.add_argument("--null-level", type=int, choices=(1, 2, 4, 5), default=5)
    parser.add_argument("--permutation-trials", type=int, default=32)
    parser.add_argument("--null-candidates", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = FMOTelemetryConfig(
        output_dir=args.output_dir,
        seed=args.seed,
        steps=args.steps,
        address_mode=args.address_mode,
        address_gain=args.address_gain,
        local_address_gain=args.local_address_gain,
        sink_gain=args.sink_gain,
        environment_assist_gain=args.environment_assist_gain,
        thermal_noise=args.thermal_noise,
        disorder_scale=args.disorder_scale,
        damage_scale=args.damage_scale,
        null_level=args.null_level,
        permutation_trials=args.permutation_trials,
        null_candidates=args.null_candidates,
    )
    outcome = run_fmo_spectral_telemetry(config)
    print(f"bridge_status: {outcome['bridge_status']}")
    print(f"failure_mode: {outcome['failure_mode']}")
    print(f"recoverability_score: {outcome['observed_summary']['recoverability_score']}")
    print(f"site_identity_retention: {outcome['observed_summary']['site_identity_retention']}")
    print(f"pathway_identity_retention: {outcome['observed_summary']['pathway_identity_retention']}")
    print(f"active_null_z: {outcome['observed_summary']['active_null_z']}")
    print(f"control_strict_pass_count: {outcome['control_strict_pass_count']}")
    print(f"outputs: {config.output_dir}")


if __name__ == "__main__":
    main()
