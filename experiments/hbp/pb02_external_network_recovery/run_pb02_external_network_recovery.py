#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from experiments.hbp.pb02_external_network_recovery.pb02_holdout import main


if __name__ == "__main__":
    raise SystemExit(main())
