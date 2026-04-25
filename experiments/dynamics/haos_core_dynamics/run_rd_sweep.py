"""Run a small HAOS graph reaction-diffusion feed/kill sweep."""

from __future__ import annotations

import argparse
from pathlib import Path

from haos_dynamics import DynamicsConfig, run_rd_feed_kill_sweep


ROOT = Path(__file__).resolve().parent


def parse_float_list(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run an RD feed/kill sweep scored by pattern specificity.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "rd_sweep", help="Sweep output directory.")
    parser.add_argument("--seed", type=int, default=20260425, help="Deterministic random seed.")
    parser.add_argument("--steps", type=int, default=220, help="Discrete dynamics steps per run.")
    parser.add_argument("--dt", type=float, default=0.4, help="Time-step size.")
    parser.add_argument("--feeds", default="0.045,0.060,0.075", help="Comma-separated feed values.")
    parser.add_argument("--kills", default="0.055,0.062,0.069", help="Comma-separated kill values.")
    parser.add_argument("--pattern-permutation-trials", type=int, default=64, help="Pattern specificity permutation trials.")
    parser.add_argument("--max-nodes", type=int, default=96, help="Maximum graph nodes kept.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = DynamicsConfig(
        output_dir=args.output_dir,
        seed=args.seed,
        steps=args.steps,
        mode="rd",
        dt=args.dt,
        pattern_permutation_trials=args.pattern_permutation_trials,
        max_nodes=args.max_nodes,
    )
    rows = run_rd_feed_kill_sweep(config, parse_float_list(args.feeds), parse_float_list(args.kills))
    best = max(rows, key=lambda row: float(row.get("pattern_specificity_z") or -1.0e9))
    print(f"runs: {len(rows)}")
    print(f"best_feed: {best['feed']}")
    print(f"best_kill: {best['kill']}")
    print(f"best_status: {best['bridge_status']}")
    print(f"best_pattern_specificity_z: {best['pattern_specificity_z']}")
    print(f"best_pattern_specificity_p: {best['pattern_specificity_p']}")
    print(f"outputs: {args.output_dir}")


if __name__ == "__main__":
    main()
