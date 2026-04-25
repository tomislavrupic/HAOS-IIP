"""Run the stripped-back HAOS minimal dynamics probe."""

from __future__ import annotations

import argparse
from pathlib import Path

from haos_minimal_dynamics import MinimalConfig, run_minimal_probe


ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run minimal HAOS relational-address dynamics.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs", help="Output directory.")
    parser.add_argument("--seed", type=int, default=20260425, help="Deterministic random seed.")
    parser.add_argument("--steps", type=int, default=120, help="Dynamics steps.")
    parser.add_argument("--dt", type=float, default=0.08, help="Time step.")
    parser.add_argument("--diffusion", type=float, default=0.15, help="Diffusion strength.")
    parser.add_argument("--address-gain", type=float, default=0.05, help="Relational address restoration strength.")
    parser.add_argument("--perturbation-step", type=int, default=20, help="Perturbation step.")
    parser.add_argument("--perturbation-scale", type=float, default=0.85, help="Perturbation amplitude.")
    parser.add_argument("--permutation-trials", type=int, default=64, help="Address specificity permutation trials.")
    parser.add_argument("--max-nodes", type=int, default=96, help="Maximum graph nodes.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = MinimalConfig(
        output_dir=args.output_dir,
        seed=args.seed,
        steps=args.steps,
        dt=args.dt,
        diffusion=args.diffusion,
        address_gain=args.address_gain,
        perturbation_step=args.perturbation_step,
        perturbation_scale=args.perturbation_scale,
        permutation_trials=args.permutation_trials,
        max_nodes=args.max_nodes,
    )
    outcome = run_minimal_probe(config)
    print(f"bridge_status: {outcome['bridge_status']}")
    print(f"failure_mode: {outcome['failure_mode']}")
    print(f"recoverability_score: {outcome['observed_summary']['recoverability_score']}")
    print(f"address_specificity_z: {outcome['observed_summary']['address_specificity_z']}")
    print(f"control_specificity_pass_count: {outcome['control_specificity_pass_count']}")
    print(f"outputs: {config.output_dir}")


if __name__ == "__main__":
    main()
