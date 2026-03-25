"""Frozen visual protocol for geometry-emergence figures."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator


CANONICAL_FIGSIZE = (6.5, 8.25)
CANONICAL_FOOTER = "HAOS-IIP geometry_emergence canonical sweep v1"
PNG_DPI = 450


def configure_geometry_figure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.8,
            "grid.linewidth": 0.4,
            "grid.color": "#d7d7d7",
            "grid.linestyle": ":",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure_bundle(
    fig: plt.Figure,
    *,
    svg_path: Path | None = None,
    pdf_path: Path | None = None,
    png_path: Path | None = None,
) -> None:
    for path in (svg_path, pdf_path, png_path):
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)

    if svg_path is not None:
        fig.savefig(svg_path, format="svg", bbox_inches="tight", facecolor="white")
    if pdf_path is not None:
        fig.savefig(pdf_path, format="pdf", bbox_inches="tight", facecolor="white")
    if png_path is not None:
        fig.savefig(
            png_path,
            format="png",
            dpi=PNG_DPI,
            bbox_inches="tight",
            facecolor="white",
        )


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        0.0,
        1.02,
        f"({label})",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
    )


def _apply_major_grid(ax: plt.Axes) -> None:
    ax.grid(True, which="major")


def _representative_seed(
    scale_runs: list[dict[str, Any]],
    per_scale_summary: list[dict[str, Any]],
) -> int:
    scale_means = {
        float(entry["scale"]): None
        if entry["mean_onset"] is None
        else float(entry["mean_onset"])
        for entry in per_scale_summary
    }
    candidate_seeds = sorted(
        set.intersection(
            *[
                {int(seed_run["seed"]) for seed_run in scale_entry["seed_runs"]}
                for scale_entry in scale_runs
            ]
        )
    )

    best_seed = candidate_seeds[0]
    best_score = float("inf")
    for seed in candidate_seeds:
        deviations: list[float] = []
        missing = 0
        for scale_entry in scale_runs:
            scale = float(scale_entry["scale"])
            mean_onset = scale_means.get(scale)
            seed_run = next(
                seed_run
                for seed_run in scale_entry["seed_runs"]
                if int(seed_run["seed"]) == seed
            )
            onset = seed_run["onset_kernel_width"]
            if mean_onset is None or onset is None:
                missing += 1
                continue
            deviations.append((float(onset) - mean_onset) ** 2)
        if missing or not deviations:
            score = float("inf")
        else:
            score = float(np.sqrt(np.mean(deviations)))
        if score < best_score or (np.isclose(score, best_score) and seed < best_seed):
            best_seed = seed
            best_score = score
    return int(best_seed)


def _detection_fraction_heatmap(
    scale_runs: list[dict[str, Any]],
    kernel_widths: np.ndarray,
    threshold: float,
) -> np.ndarray:
    heatmap = np.zeros((len(scale_runs), kernel_widths.size), dtype=float)
    for row_index, scale_entry in enumerate(scale_runs):
        curves = [
            np.array(
                seed_run["onset_detection"]["coupled_progress"]["coupled_diagnostic_score"],
                dtype=float,
            )
            for seed_run in scale_entry["seed_runs"]
        ]
        curve_stack = np.vstack(curves)
        heatmap[row_index, :] = np.mean(curve_stack >= float(threshold), axis=0)
    return heatmap


def build_cluster_scale_phase_diagram(
    full_payload: dict[str, Any],
    summary_payload: dict[str, Any],
) -> plt.Figure:
    configure_geometry_figure_style()

    kernel_widths = np.array(full_payload["kernel_widths"], dtype=float)
    scale_runs = sorted(full_payload["scale_runs"], key=lambda entry: float(entry["scale"]))
    per_scale_summary = sorted(
        summary_payload["per_scale"], key=lambda entry: float(entry["scale"])
    )
    threshold = float(summary_payload["config"]["coupled_threshold"])
    seed_count = len(summary_payload["config"]["seeds"])

    scales = np.array([float(entry["scale"]) for entry in per_scale_summary], dtype=float)
    mean_onsets = np.array([float(entry["mean_onset"]) for entry in per_scale_summary], dtype=float)
    std_onsets = np.array([float(entry["std_onset"]) for entry in per_scale_summary], dtype=float)

    figure = plt.figure(figsize=CANONICAL_FIGSIZE, facecolor="white")
    grid = figure.add_gridspec(3, 1, height_ratios=[1.0, 1.15, 1.0])
    ax_a = figure.add_subplot(grid[0, 0])
    ax_b = figure.add_subplot(grid[1, 0])
    ax_c = figure.add_subplot(grid[2, 0])

    marker_color = "#222222"
    ax_a.errorbar(
        scales,
        mean_onsets,
        yerr=std_onsets,
        fmt="o-",
        color=marker_color,
        linewidth=1.0,
        elinewidth=1.0,
        capsize=2.5,
        markersize=4,
        markerfacecolor=marker_color,
        markeredgecolor=marker_color,
    )
    ax_a.set_xlabel("Cluster scale parameter")
    ax_a.set_ylabel("Transition onset kernel width")
    ax_a.set_xlim(float(scales[0]) - 0.05, float(scales[-1]) + 0.05)
    y_min = float(np.min(mean_onsets - std_onsets))
    y_max = float(np.max(mean_onsets + std_onsets))
    ax_a.set_ylim(np.floor((y_min - 0.002) * 100.0) / 100.0, np.ceil((y_max + 0.002) * 100.0) / 100.0)
    ax_a.set_xticks(scales)
    ax_a.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax_a.yaxis.set_major_locator(MultipleLocator(0.02))
    ax_a.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    _apply_major_grid(ax_a)
    ax_a.text(
        0.03,
        0.95,
        f"Detection threshold: coupled_score >= {threshold:.2f}\nSeeds: N = {seed_count}",
        transform=ax_a.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    _panel_label(ax_a, "A")

    heatmap = _detection_fraction_heatmap(scale_runs, kernel_widths, threshold)
    kernel_step = float(summary_payload["config"]["kernel_step"])
    scale_step = float(np.min(np.diff(scales))) if scales.size > 1 else 0.1
    extent = (
        float(kernel_widths[0] - kernel_step / 2.0),
        float(kernel_widths[-1] + kernel_step / 2.0),
        float(scales[0] - scale_step / 2.0),
        float(scales[-1] + scale_step / 2.0),
    )
    image = ax_b.imshow(
        heatmap,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=extent,
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
    )
    ax_b.set_xlabel("Kernel interaction width")
    ax_b.set_ylabel("Cluster scale")
    ax_b.set_xlim(float(kernel_widths[0]), float(kernel_widths[-1]))
    ax_b.set_ylim(float(scales[0]) - scale_step / 2.0, float(scales[-1]) + scale_step / 2.0)
    ax_b.xaxis.set_major_locator(MultipleLocator(0.02))
    ax_b.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax_b.set_yticks(scales)
    ax_b.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    colorbar = figure.colorbar(image, ax=ax_b, fraction=0.045, pad=0.02)
    colorbar.set_label("Detection fraction", fontsize=8)
    colorbar.set_ticks([0.0, 0.5, 1.0])
    colorbar.ax.tick_params(labelsize=8)
    _panel_label(ax_b, "B")

    representative_seed = _representative_seed(scale_runs, per_scale_summary)
    grayscale = plt.get_cmap("Greys")(np.linspace(0.85, 0.50, len(scale_runs)))
    for color, scale_entry in zip(grayscale, scale_runs):
        seed_run = next(
            seed_run
            for seed_run in scale_entry["seed_runs"]
            if int(seed_run["seed"]) == representative_seed
        )
        coupled = np.array(
            seed_run["onset_detection"]["coupled_progress"]["coupled_diagnostic_score"],
            dtype=float,
        )
        ax_c.plot(
            kernel_widths,
            coupled,
            color=color,
            linewidth=1.2,
            label=f"scale {float(scale_entry['scale']):.2f}",
        )
    ax_c.axhline(threshold, color="#5f5f5f", linestyle="--", linewidth=1.0)
    ax_c.text(
        float(kernel_widths[-1]) - 0.001,
        threshold + 0.02,
        "Transition threshold",
        ha="right",
        va="bottom",
        fontsize=8,
        color="#4f4f4f",
    )
    ax_c.set_xlabel("Kernel width")
    ax_c.set_ylabel("Coupled transport stability score")
    ax_c.set_xlim(float(kernel_widths[0]), float(kernel_widths[-1]))
    ax_c.set_ylim(-0.02, 1.02)
    ax_c.xaxis.set_major_locator(MultipleLocator(0.02))
    ax_c.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax_c.yaxis.set_major_locator(MultipleLocator(0.2))
    _apply_major_grid(ax_c)
    ax_c.legend(frameon=False, loc="upper left", ncol=1, handlelength=2.2)
    _panel_label(ax_c, "C")

    figure.text(
        0.985,
        0.015,
        CANONICAL_FOOTER,
        ha="right",
        va="bottom",
        fontsize=6.5,
        color="#4d4d4d",
    )
    figure.subplots_adjust(left=0.12, right=0.88, top=0.98, bottom=0.07, hspace=0.42)
    return figure
