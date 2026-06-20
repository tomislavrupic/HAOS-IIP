#!/usr/bin/env python3
"""Compatibility entry point for Betti null-ensemble hardening."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from threshold_sweep_betti_stability import BettiSweepConfig, parse_args, run_threshold_sweep  # noqa: E402


def main() -> None:
    args = parse_args()
    result = run_threshold_sweep(BettiSweepConfig(), args.output_dir)
    null_summary = result["null_summary"]
    print(
        json.dumps(
            {
                "false_positive_rate": null_summary["false_positive_rate"],
                "result_hash": result["result_hash"],
                "status": result["status"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

