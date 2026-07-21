from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

from experiments.hbp.data_paths import POWERGRAPH_ROOT_ENV, REPO_ROOT, powergraph_data_root


class HBPDataPathTests(unittest.TestCase):
    def test_default_powergraph_root_is_workspace_relative(self) -> None:
        with patch.dict("os.environ", {}, clear=True):
            self.assertEqual(powergraph_data_root(), REPO_ROOT.parent / "DATA" / "Powergraph")

    def test_environment_override_is_supported(self) -> None:
        configured = "/tmp/haos-powergraph-fixture"
        with patch.dict("os.environ", {POWERGRAPH_ROOT_ENV: configured}, clear=True):
            self.assertEqual(powergraph_data_root(), Path(configured))


if __name__ == "__main__":
    unittest.main()
