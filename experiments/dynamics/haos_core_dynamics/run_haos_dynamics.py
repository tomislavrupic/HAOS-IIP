"""Launcher for the HAOS-IIP core dynamics sidecar probe."""

from __future__ import annotations

import argparse
from pathlib import Path

from haos_dynamics import DynamicsConfig, run_dynamics_probe


ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the HAOS-IIP core dynamics sidecar probe.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs", help="Output directory.")
    parser.add_argument("--seed", type=int, default=20260425, help="Deterministic random seed.")
    parser.add_argument("--steps", type=int, default=160, help="Discrete dynamics steps.")
    parser.add_argument("--mode", choices=("scalar", "rd"), default="scalar", help="Dynamics family.")
    parser.add_argument("--dt", type=float, default=0.08, help="Time-step size.")
    parser.add_argument("--diffusion", type=float, default=0.55, help="Local Laplacian-flow strength.")
    parser.add_argument("--recovery-gain", type=float, default=0.18, help="Reference-field recovery strength.")
    parser.add_argument("--rd-du", type=float, default=0.16, help="Gray-Scott U diffusion.")
    parser.add_argument("--rd-dv", type=float, default=0.08, help="Gray-Scott V diffusion.")
    parser.add_argument("--rd-feed", type=float, default=0.060, help="Gray-Scott feed rate.")
    parser.add_argument("--rd-kill", type=float, default=0.062, help="Gray-Scott kill rate.")
    parser.add_argument("--perturbation-step", type=int, default=24, help="Step where perturbation is applied.")
    parser.add_argument("--perturbation-scale", type=float, default=0.75, help="Perturbation amplitude.")
    parser.add_argument("--damage-fraction", type=float, default=0.0, help="Fraction of edges to remove before dynamics.")
    parser.add_argument("--max-nodes", type=int, default=96, help="Maximum graph nodes kept.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = DynamicsConfig(
        output_dir=args.output_dir,
        seed=args.seed,
        steps=args.steps,
        mode=args.mode,
        dt=args.dt,
        diffusion=args.diffusion,
        recovery_gain=args.recovery_gain,
        rd_du=args.rd_du,
        rd_dv=args.rd_dv,
        rd_feed=args.rd_feed,
        rd_kill=args.rd_kill,
        perturbation_step=args.perturbation_step,
        perturbation_scale=args.perturbation_scale,
        damage_fraction=args.damage_fraction,
        max_nodes=args.max_nodes,
    )
    outcome = run_dynamics_probe(config)
    print(f"bridge_status: {outcome.status}")
    print(f"graph_source: {outcome.graph_source}")
    print(f"nodes: {outcome.node_count}")
    print(f"edges: {outcome.edge_count}")
    print(f"recoverability_score: {outcome.summary.get('recoverability_score')}")
    print(f"control_contrast: {outcome.summary.get('control_contrast')}")
    print(f"outputs: {outcome.output_dir}")


if __name__ == "__main__":
    main()
