"""Launcher for the HAOS-IIP Kuramoto oscillator bridge."""

from __future__ import annotations

import argparse
from pathlib import Path

from kuramoto_bridge import KuramotoConfig, run_bridge


ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the HAOS-IIP Kuramoto oscillator bridge.")
    parser.add_argument("--graph", type=Path, help="Optional graph JSON/CSV/NPZ path.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs", help="Output directory.")
    parser.add_argument("--seed", type=int, default=20260425, help="Deterministic random seed.")
    parser.add_argument("--max-nodes", type=int, default=96, help="Maximum nodes kept from large graphs.")
    parser.add_argument("--k-min", type=float, default=0.0, help="Minimum coupling K.")
    parser.add_argument("--k-max", type=float, default=4.0, help="Maximum coupling K.")
    parser.add_argument("--k-count", type=int, default=17, help="Number of K values to scan.")
    parser.add_argument("--steps", type=int, default=900, help="Integration steps per K value.")
    parser.add_argument("--dt", type=float, default=0.03, help="Euler integration step size.")
    parser.add_argument("--omega-mean", type=float, default=0.0, help="Gaussian natural-frequency mean.")
    parser.add_argument("--omega-std", type=float, default=1.0, help="Gaussian natural-frequency standard deviation.")
    parser.add_argument(
        "--frequency-mode",
        choices=("gaussian", "custom"),
        default="gaussian",
        help="Natural-frequency source.",
    )
    parser.add_argument("--frequency-file", type=Path, help="Numeric frequency file for custom mode.")
    parser.add_argument(
        "--proxy-mode",
        choices=("generic", "line_e_like"),
        default="generic",
        help="HAOS proxy family used for recoverability diagnostics.",
    )
    parser.add_argument("--higher-order", action="store_true", help="Enable deterministic triangle coupling.")
    parser.add_argument("--higher-order-strength", type=float, default=0.15, help="Triangle coupling strength.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = KuramotoConfig(
        explicit_graph_path=args.graph,
        output_dir=args.output_dir,
        seed=args.seed,
        max_nodes=args.max_nodes,
        k_min=args.k_min,
        k_max=args.k_max,
        k_count=args.k_count,
        steps=args.steps,
        dt=args.dt,
        omega_mean=args.omega_mean,
        omega_std=args.omega_std,
        frequency_mode=args.frequency_mode,
        frequency_file=args.frequency_file,
        proxy_mode=args.proxy_mode,
        higher_order=args.higher_order,
        higher_order_strength=args.higher_order_strength,
    )
    outcome = run_bridge(config)
    print(f"bridge_status: {outcome.status}")
    print(f"graph_source: {outcome.graph.source_label}")
    print(f"nodes: {outcome.graph.node_count}")
    print(f"edges: {outcome.graph.edge_count}")
    print(f"proxy_mode: {outcome.primary_summary.get('proxy_mode')}")
    print(f"k_star_time: {outcome.primary_summary.get('k_star_time')}")
    print(f"early_detection: {outcome.primary_summary.get('early_detection')}")
    print(f"outputs: {outcome.output_dir}")


if __name__ == "__main__":
    main()
