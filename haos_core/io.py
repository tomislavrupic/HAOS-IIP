from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]


def timestamp_slug() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_timestamped_json(output_dir: Path, stem: str, payload: dict[str, Any]) -> tuple[Path, Path]:
    ensure_dir(output_dir)
    stamped = output_dir / f"{timestamp_slug()}_{stem}.json"
    latest = output_dir / f"{stem}_latest.json"
    write_json(stamped, payload)
    write_json(latest, payload)
    return stamped, latest


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def relpath(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))
