from __future__ import annotations

import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
POWERGRAPH_ROOT_ENV = "HAOS_POWERGRAPH_ROOT"


def powergraph_data_root() -> Path:
    configured = os.environ.get(POWERGRAPH_ROOT_ENV)
    if configured:
        return Path(configured).expanduser()
    return REPO_ROOT.parent / "DATA" / "Powergraph"


DEFAULT_POWERGRAPH_ROOT = powergraph_data_root()
