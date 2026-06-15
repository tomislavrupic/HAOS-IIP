from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]


def canonicalize(value: Any, *, float_precision: int = 12) -> Any:
    if isinstance(value, float):
        return round(value, float_precision)
    if isinstance(value, dict):
        return {
            str(key): canonicalize(value[key], float_precision=float_precision)
            for key in sorted(value)
        }
    if isinstance(value, (list, tuple)):
        return [canonicalize(item, float_precision=float_precision) for item in value]
    return value


def canonical_json(value: Any, *, float_precision: int = 12) -> str:
    return json.dumps(
        canonicalize(value, float_precision=float_precision),
        sort_keys=True,
        separators=(",", ":"),
    )


def stable_hash(value: Any, *, prefix: str = "", float_precision: int = 12) -> str:
    data = canonical_json(value, float_precision=float_precision).encode("utf-8")
    digest = hashlib.sha256(data).hexdigest()
    return f"{prefix}{digest[:24]}" if prefix else digest


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repo_rel(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(str(path))
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def code_version() -> str:
    try:
        result = subprocess.run(
            ["git", "-c", "core.fsmonitor=false", "rev-parse", "--short", "HEAD"],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return "git-unavailable"
    if result.returncode != 0:
        return "git-unavailable"
    return result.stdout.strip() or "git-unavailable"


def artifact_hashes(paths: dict[str, Path]) -> dict[str, str]:
    return {name: file_sha256(path) for name, path in sorted(paths.items())}
