#!/usr/bin/env python3
"""Compatibility entry point for Betti threshold sweep hardening."""

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
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()

