#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent

PHASE_COMMANDS = {
    "3": [ROOT / "phase3-stability" / "runs" / "run_phase3.py"],
    "4": [ROOT / "phase4-sector-freeze" / "runs" / "run_phase4.py"],
    "5": [ROOT / "phase5-readout" / "runs" / "run_phase5.py"],
    "6": [ROOT / "phase6-operator" / "build_operator_family.py"],
    "7": [ROOT / "phase7-spectral" / "build_spectral_feasibility.py"],
    "8": [ROOT / "phase8-trace" / "build_trace_asymptotic_feasibility.py"],
    "9": [ROOT / "phase9-invariants" / "build_coefficient_stability.py"],
    "10": [ROOT / "phase10-bridge" / "build_continuum_bridge_feasibility.py"],
    "x": [ROOT / "phaseX-proto-particle" / "build_integrated_proto_particle_feasibility.py"],
    "11": [ROOT / "phase11-protection" / "build_localized_mode_protection.py"],
    "12": [ROOT / "phase12-interactions" / "build_two_mode_interaction.py"],
    "13": [ROOT / "phase13-sector-formation" / "build_multi_mode_sector.py"],
    "14": [ROOT / "phase14-collective-dynamics" / "build_collective_dynamics.py"],
    "15": [ROOT / "phase15-propagation" / "build_propagation_structure.py"],
    "16": [ROOT / "phase16-temporal-ordering" / "build_temporal_ordering.py"],
    "17": [ROOT / "phase17-causal-closure" / "build_causal_closure.py"],
    "18": [ROOT / "phase18-distance-surrogate" / "build_distance_surrogate.py"],
}

PHASE_CHECKERS = {
    "10": ROOT / "phase10-bridge" / "diagnostics" / "check_phase10_bundle.py",
    "11": ROOT / "phase11-protection" / "diagnostics" / "check_phase11_bundle.py",
    "12": ROOT / "phase12-interactions" / "diagnostics" / "check_phase12_bundle.py",
    "13": ROOT / "phase13-sector-formation" / "diagnostics" / "check_phase13_bundle.py",
    "14": ROOT / "phase14-collective-dynamics" / "diagnostics" / "check_phase14_bundle.py",
    "15": ROOT / "phase15-propagation" / "diagnostics" / "check_phase15_bundle.py",
    "16": ROOT / "phase16-temporal-ordering" / "diagnostics" / "check_phase16_bundle.py",
    "17": ROOT / "phase17-causal-closure" / "diagnostics" / "check_phase17_bundle.py",
    "18": ROOT / "phase18-distance-surrogate" / "diagnostics" / "check_phase18_bundle.py",
    "x": ROOT / "phaseX-proto-particle" / "diagnostics" / "check_phaseX_bundle.py",
}


def normalize_phase_id(value: str) -> str:
    phase = value.strip().lower()
    if phase == "proto":
        return "x"
    return phase


def main() -> None:
    parser = argparse.ArgumentParser(description="Run or check one frozen HAOS-IIP phase builder.")
    parser.add_argument("phase", help="Phase id, for example: 10, 16, 18, or X")
    parser.add_argument("--check", action="store_true", help="Run the phase checker instead of the builder")
    args = parser.parse_args()

    phase_id = normalize_phase_id(args.phase)
    target = PHASE_CHECKERS.get(phase_id) if args.check else None
    if target is None:
        command = PHASE_COMMANDS.get(phase_id)
        if command is None:
            valid = ", ".join(sorted(PHASE_COMMANDS))
            raise SystemExit(f"Unknown phase '{args.phase}'. Valid phase ids: {valid}")
        target = command[0]

    if not Path(target).exists():
        raise SystemExit(f"Missing target for phase '{args.phase}': {target}")

    subprocess.run([sys.executable, str(target)], check=True, cwd=str(ROOT))


if __name__ == "__main__":
    main()
