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
    parser.add_argument("--address-gain", type=float, default=0.45, help="Relational address restoration strength.")
    parser.add_argument("--address-mode", choices=("spectral", "local", "hybrid", "phase", "multi_scale"), default="spectral", help="Address restoration mode.")
    parser.add_argument("--spectral-modes", type=int, default=8, help="Low-frequency Laplacian modes for spectral address.")
    parser.add_argument("--hybrid-spectral-weight", type=float, default=0.7, help="Spectral weight for hybrid address mode.")
    parser.add_argument("--invariant-gain", type=float, default=0.025, help="Branch-local shell invariant restoration strength.")
    parser.add_argument("--perturbation-step", type=int, default=20, help="Perturbation step.")
    parser.add_argument("--perturbation-scale", type=float, default=0.85, help="Perturbation amplitude.")
    parser.add_argument("--permutation-trials", type=int, default=64, help="Address specificity permutation trials.")
    parser.add_argument("--identity-bins", type=int, default=4, help="Degree/shell bins for branch-identity specificity.")
    parser.add_argument("--null-level", type=int, choices=(1, 2, 3), default=2, help="Specificity null strictness: 1 degree/shell, 2 spectral-aware, 3 spectral+autocorrelation.")
    parser.add_argument("--spectral-null-candidates", type=int, default=8, help="Candidate shuffles sampled per spectral-aware null draw.")
    parser.add_argument("--focus-address", action="store_true", help="Only run address-focused ablations.")
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
        address_mode=args.address_mode,
        spectral_modes=args.spectral_modes,
        hybrid_spectral_weight=args.hybrid_spectral_weight,
        invariant_gain=args.invariant_gain,
        perturbation_step=args.perturbation_step,
        perturbation_scale=args.perturbation_scale,
        permutation_trials=args.permutation_trials,
        identity_bins=args.identity_bins,
        null_level=args.null_level,
        spectral_null_candidates=args.spectral_null_candidates,
        focus_address=args.focus_address,
        max_nodes=args.max_nodes,
    )
    outcome = run_minimal_probe(config)
    print(f"bridge_status: {outcome['bridge_status']}")
    print(f"failure_mode: {outcome['failure_mode']}")
    print(f"recoverability_score: {outcome['observed_summary']['recoverability_score']}")
    print(f"address_specificity_z: {outcome['observed_summary']['address_specificity_z']}")
    print(f"invariant_specificity_z: {outcome['observed_summary']['invariant_specificity_z']}")
    print(f"branch_identity_z: {outcome['observed_summary']['branch_identity_z']}")
    print(f"spectral_aware_z: {outcome['observed_summary']['spectral_aware_z']}")
    print(f"autocorr_aware_z: {outcome['observed_summary']['autocorr_aware_z']}")
    print(f"control_combined_specificity_pass_count: {outcome['control_combined_specificity_pass_count']}")
    print(f"control_strict_specificity_pass_count: {outcome['control_strict_specificity_pass_count']}")
    print(f"control_spectral_aware_pass_count: {outcome['control_spectral_aware_pass_count']}")
    print(f"control_autocorr_aware_pass_count: {outcome['control_autocorr_aware_pass_count']}")
    print(f"ablation_report: {config.output_dir / 'ablation_report.md'}")
    print(f"outputs: {config.output_dir}")


if __name__ == "__main__":
    main()
