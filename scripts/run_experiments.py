#!/usr/bin/env python3

from __future__ import annotations

import json
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from numerics.simulations.laplacian_modes import load_config, run_laplacian_test
from numerics.simulations.gauge_modes import run_gauge_test

RESULTS = ROOT / "data"
PLOTS = ROOT / "plots"
LOG = ROOT / "experiments" / "EXPERIMENT_LOG.md"
CONFIG = ROOT / "config.json"


def save_results(name: str, result: dict[str, Any], timestamp: str) -> Path:
    RESULTS.mkdir(exist_ok=True)
    stamped = RESULTS / f"{timestamp}_{name}.json"
    latest = RESULTS / f"{name}_latest.json"
    payload = json.dumps(result, indent=2)
    stamped.write_text(payload, encoding="utf-8")
    latest.write_text(payload, encoding="utf-8")
    return stamped


def copy_plots(name: str, result: dict[str, Any], timestamp: str) -> list[str]:
    saved: list[str] = []
    for rel_path in result.get("plots", []):
        src = ROOT / rel_path
        if not src.exists():
            continue
        dst = PLOTS / f"{timestamp}_{name}_{src.name}"
        shutil.copy2(src, dst)
        saved.append(str(dst.relative_to(ROOT)))
    return saved


def log_experiment(name: str, result_path: Path, result: dict[str, Any], saved_plots: list[str]) -> None:
    if not LOG.exists():
        LOG.write_text("# EXPERIMENT_LOG\n", encoding="utf-8")
    cfg = result.get("config", {})
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(f"\n## {name}\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(f"- Config: substrate={cfg.get('substrate')}, nodes={cfg.get('nodes')}, epsilon={cfg.get('epsilon')}, seed={cfg.get('random_seed')}\n")
        handle.write(f"- Results: `{result_path.relative_to(ROOT)}`\n")
        if saved_plots:
            handle.write(f"- Plots: {', '.join(f'`{p}`' for p in saved_plots)}\n")
        handle.write(f"- Observation: {result.get('observation', 'n/a')}\n")
        handle.write(f"- Conclusion: {result.get('conclusion', 'n/a')}\n")


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    PLOTS.mkdir(exist_ok=True)
    LOG.parent.mkdir(exist_ok=True)

    config = load_config(config_path=CONFIG)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    lap = run_laplacian_test(config)
    p1 = save_results("laplacian_modes", lap, timestamp)
    plots1 = copy_plots("laplacian_modes", lap, timestamp)
    log_experiment("Laplacian geometry test", p1, lap, plots1)

    gauge = run_gauge_test(config)
    p2 = save_results("gauge_modes", gauge, timestamp)
    plots2 = copy_plots("gauge_modes", gauge, timestamp)
    log_experiment("Gauge sector test", p2, gauge, plots2)

    print(f"Config: {CONFIG}")
    print(f"Saved: {p1.relative_to(ROOT)}")
    print(f"Saved: {p2.relative_to(ROOT)}")
    print(f"Log: {LOG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
