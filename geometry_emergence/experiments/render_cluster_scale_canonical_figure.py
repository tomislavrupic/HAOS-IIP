#!/usr/bin/env python3

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE_ROOT = Path(__file__).resolve().parents[1]
FULL_DATA_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_sweep_full.json"
SUMMARY_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_sweep_summary.json"
SVG_OUTPUT_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_transition_canonical.svg"
PDF_OUTPUT_PATH = ROOT / "papers" / "figures" / "cluster_scale_transition_canonical.pdf"
SVG_PAPER_OUTPUT_PATH = ROOT / "papers" / "figures" / "cluster_scale_transition_canonical.svg"

os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def round_float(value: float, digits: int = 6) -> float:
    return round(float(value), digits)


def build_plot() -> None:
    full_payload = read_json(FULL_DATA_PATH)
    summary_payload = read_json(SUMMARY_PATH)

    kernel_widths = np.array(full_payload["kernel_widths"], dtype=float)
    scale_runs = full_payload["scale_runs"]
    per_scale_summary = summary_payload["per_scale"]
    threshold = float(summary_payload["config"]["coupled_threshold"])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(9.5, 8.2),
        gridspec_kw={"height_ratios": [3.2, 1.7]},
        facecolor="white",
    )

    colors = plt.get_cmap("viridis")(np.linspace(0.15, 0.9, len(scale_runs)))
    scale_to_color = {}

    for color, scale_entry in zip(colors, scale_runs):
        scale = float(scale_entry["scale"])
        scale_to_color[scale] = color
        label = f"scale {scale:.2f}"
        curves = []
        for seed_run in scale_entry["seed_runs"]:
            coupled = np.array(
                seed_run["onset_detection"]["coupled_progress"]["coupled_diagnostic_score"],
                dtype=float,
            )
            curves.append(coupled)
            ax_top.plot(
                kernel_widths,
                coupled,
                color=color,
                linewidth=1.0,
                alpha=0.15,
            )
        mean_curve = np.mean(np.vstack(curves), axis=0)
        ax_top.plot(
            kernel_widths,
            mean_curve,
            color=color,
            linewidth=2.6,
            label=label,
        )

    ax_top.axhline(threshold, color="#666666", linestyle="--", linewidth=1.0)
    ax_top.set_title("Canonical Cluster-Scale Transition Signal")
    ax_top.set_xlabel("Kernel width")
    ax_top.set_ylabel("Coupled diagnostic score")
    ax_top.set_xlim(float(kernel_widths[0]), float(kernel_widths[-1]))
    ax_top.set_ylim(-0.02, 1.02)
    ax_top.grid(alpha=0.2, color="#777777")
    ax_top.legend(frameon=False, ncol=2, loc="lower right")

    scales = np.array([row["scale"] for row in per_scale_summary], dtype=float)
    means = np.array([row["mean_onset"] for row in per_scale_summary], dtype=float)
    mins = np.array([row["min_onset"] for row in per_scale_summary], dtype=float)
    maxs = np.array([row["max_onset"] for row in per_scale_summary], dtype=float)
    detected = np.array([row["detected_fraction"] for row in per_scale_summary], dtype=float)
    lower = means - mins
    upper = maxs - means

    for scale, mean, lo, hi, fraction in zip(scales, means, lower, upper, detected):
        color = scale_to_color[float(scale)]
        ax_bottom.errorbar(
            scale,
            mean,
            yerr=np.array([[lo], [hi]], dtype=float),
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=2.0,
            capsize=4,
            markersize=6,
        )
        ax_bottom.text(
            scale,
            mean + hi + 0.003,
            f"{int(round(fraction * 5))}/5",
            color=color,
            ha="center",
            va="bottom",
            fontsize=8.5,
        )

    ax_bottom.plot(scales, means, color="#303030", linewidth=1.2, alpha=0.7)
    ax_bottom.set_xlabel("Cluster scale")
    ax_bottom.set_ylabel("Onset kernel width")
    ax_bottom.set_xlim(float(scales[0]) - 0.05, float(scales[-1]) + 0.05)
    ax_bottom.set_ylim(float(min(mins)) - 0.006, float(max(maxs)) + 0.012)
    ax_bottom.grid(alpha=0.2, color="#777777")
    ax_bottom.set_title("Onset Location by Cluster Scale")

    fig.text(
        0.015,
        0.985,
        "A",
        ha="left",
        va="top",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.015,
        0.36,
        "B",
        ha="left",
        va="top",
        fontsize=14,
        fontweight="bold",
    )

    fig.tight_layout()
    SVG_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    PDF_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(SVG_OUTPUT_PATH, format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(PDF_OUTPUT_PATH, format="pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(SVG_PAPER_OUTPUT_PATH, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    build_plot()

