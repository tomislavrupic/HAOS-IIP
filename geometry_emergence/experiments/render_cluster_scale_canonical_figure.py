#!/usr/bin/env python3

from __future__ import annotations

import json
import os
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODULE_ROOT = Path(__file__).resolve().parents[1]
FULL_DATA_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_sweep_full.json"
SUMMARY_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_sweep_summary.json"
SVG_OUTPUT_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_transition_canonical.svg"
PDF_OUTPUT_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_transition_canonical.pdf"
PNG_OUTPUT_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_transition_canonical.png"
PDF_PAPER_OUTPUT_PATH = ROOT / "papers" / "figures" / "cluster_scale_transition_canonical.pdf"
SVG_PAPER_OUTPUT_PATH = ROOT / "papers" / "figures" / "cluster_scale_transition_canonical.svg"
PNG_PAPER_OUTPUT_PATH = ROOT / "papers" / "figures" / "cluster_scale_transition_canonical.png"

os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from geometry_emergence.figure_protocol import build_cluster_scale_phase_diagram, save_figure_bundle


def read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def build_plot() -> None:
    full_payload = read_json(FULL_DATA_PATH)
    summary_payload = read_json(SUMMARY_PATH)
    figure = build_cluster_scale_phase_diagram(
        full_payload=full_payload,
        summary_payload=summary_payload,
    )
    save_figure_bundle(
        figure,
        svg_path=SVG_OUTPUT_PATH,
        pdf_path=PDF_OUTPUT_PATH,
        png_path=PNG_OUTPUT_PATH,
    )
    save_figure_bundle(
        figure,
        svg_path=SVG_PAPER_OUTPUT_PATH,
        pdf_path=PDF_PAPER_OUTPUT_PATH,
        png_path=PNG_PAPER_OUTPUT_PATH,
    )
    figure.clf()


if __name__ == "__main__":
    build_plot()
