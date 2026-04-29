#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class ContrastRow:
    label: str
    runs: int
    pass_rate: float
    recoverability: float
    identity: float
    pathway_identity: float | None
    active_null_z: float
    control_pass_max: int
    status: str
    note: str


ROWS = [
    ContrastRow(
        label="Microtubule 55.2",
        runs=50,
        pass_rate=1.0,
        recoverability=0.916789,
        identity=0.871810,
        pathway_identity=None,
        active_null_z=28.705018,
        control_pass_max=0,
        status="ROBUST PASS",
        note="Structured cylindrical lattice preserves protofilament identity under level-5 controls.",
    ),
    ContrastRow(
        label="FMO 56.1 baseline",
        runs=50,
        pass_rate=0.04,
        recoverability=0.803724,
        identity=0.984934,
        pathway_identity=0.546597,
        active_null_z=2.472024,
        control_pass_max=2,
        status="DIAGNOSTIC FAIL",
        note="Site identity is strong, but pathway identity and control discrimination are weak.",
    ),
    ContrastRow(
        label="FMO 56.2 sink",
        runs=50,
        pass_rate=0.08,
        recoverability=0.827528,
        identity=0.988291,
        pathway_identity=0.622057,
        active_null_z=2.659739,
        control_pass_max=2,
        status="DIAGNOSTIC FAIL",
        note="Reaction-center sink improves pathway retention but does not break controls.",
    ),
    ContrastRow(
        label="FMO 56.3 flux",
        runs=50,
        pass_rate=0.0,
        recoverability=0.930908,
        identity=0.998190,
        pathway_identity=0.906155,
        active_null_z=3.227554,
        control_pass_max=2,
        status="CONTROL MATCH",
        note="Explicit flux restoration solves retention but is too easy for controls to exploit.",
    ),
    ContrastRow(
        label="FMO 56.4 intrinsic",
        runs=50,
        pass_rate=0.16,
        recoverability=0.848259,
        identity=0.991098,
        pathway_identity=0.671540,
        active_null_z=2.612894,
        control_pass_max=2,
        status="PRODUCTIVE FALSIFIER",
        note="Intrinsic directed dynamics improves pass rate slightly but loses pathway fidelity.",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write Phase 56.5 topology contrast summary.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir = args.output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    write_markdown(args.output_dir / "fmo_topology_contrast_56_5.md")
    write_plot(figure_dir / "topology_contrast_56_5.png")
    print(f"rows: {len(ROWS)}")
    print(f"report: {args.output_dir / 'fmo_topology_contrast_56_5.md'}")
    print(f"figure: {figure_dir / 'topology_contrast_56_5.png'}")


def write_markdown(path: Path) -> None:
    lines = [
        "# Phase 56.5 FMO Topology Contrast and Lessons",
        "",
        "This closes the initial FMO telemetry arc as a diagnostic module rather than forcing a PASS.",
        "The comparison is bounded: microtubule and FMO are toy telemetry substrates, not biological validation.",
        "",
        "## Summary Table",
        "",
        "| system | runs | pass_rate | recoverability | identity | pathway_identity | active_null_z | max_controls | status |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in ROWS:
        pathway = "n/a" if row.pathway_identity is None else f"{row.pathway_identity:.6f}"
        lines.append(
            f"| {row.label} | {row.runs} | {row.pass_rate:.6f} | {row.recoverability:.6f} | "
            f"{row.identity:.6f} | {pathway} | {row.active_null_z:.6f} | {row.control_pass_max} | {row.status} |"
        )
    lines.extend(
        [
            "",
            "## Topology Sensitivity",
            "",
            "The microtubule lattice and the FMO-like network expose different strengths and limits of the same telemetry instrument.",
            "",
            "- Microtubule 55.2 is a robust external toy PASS: structured cylindrical topology, repeated protofilament identity, and multi-scale local coupling remain specific under strict controls.",
            "- FMO preserves site identity strongly, but compact delocalized transfer topology makes pathway identity and control discrimination much harder.",
            "- Explicit flux restoration in 56.3 proves pathway retention can be forced, but it also shows that engineered pathway terms can be copied by controls.",
            "- Intrinsic pathway dynamics in 56.4 reduces the obviousness of the term and slightly improves pass rate, but sacrifices pathway fidelity.",
            "",
            "## FMO Lesson",
            "",
            "FMO should be kept as a productive falsifier. It does not invalidate the spectral dynamics instrument; it shows that the instrument is topology-sensitive.",
            "",
            "The current bottleneck is not simple recoverability or site identity. The bottleneck is branch-specific pathway discrimination on compact transfer graphs.",
            "",
            "## Closure Decision",
            "",
            "Do not keep tuning FMO until it turns green. Freeze the initial FMO arc as an honest diagnostic result:",
            "",
            "- 56.1: baseline falsifier;",
            "- 56.2: sink improves pathway retention but not specificity;",
            "- 56.3: explicit flux solves retention but creates control leakage;",
            "- 56.4: intrinsic directed dynamics improves pass rate modestly but loses retention;",
            "- 56.5: topology sensitivity closure.",
            "",
            "Future FMO work should require a new transfer-state model or a branch-specific control-discrimination metric, not just stronger scalar address terms.",
            "",
            "## Non-Claims",
            "",
            "- No claim is made about molecular FMO dynamics.",
            "- No claim is made about quantum biology or photosynthesis.",
            "- No claim is made that HAOS spectral dynamics universally solve biological telemetry.",
            "- The positive result remains the microtubule toy lattice; the FMO result remains diagnostic.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(path: Path) -> None:
    import matplotlib.pyplot as plt

    labels = [row.label.replace(" ", "\n") for row in ROWS]
    x = np.arange(len(ROWS))
    width = 0.35

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))

    axes[0, 0].bar(labels, [row.pass_rate for row in ROWS], color=["tab:green", "tab:orange", "tab:orange", "tab:orange", "tab:orange"])
    axes[0, 0].set_ylim(0.0, 1.05)
    axes[0, 0].set_ylabel("pass rate")
    axes[0, 0].set_title("Robust pass rate")
    axes[0, 0].tick_params(axis="x", labelsize=8)

    axes[0, 1].bar(x - width / 2, [row.recoverability for row in ROWS], width, label="recoverability")
    axes[0, 1].bar(x + width / 2, [row.identity for row in ROWS], width, label="identity")
    axes[0, 1].set_xticks(x, labels)
    axes[0, 1].set_ylim(0.0, 1.05)
    axes[0, 1].set_title("Recoverability and primary identity")
    axes[0, 1].legend()
    axes[0, 1].tick_params(axis="x", labelsize=8)

    fmo_rows = [row for row in ROWS if row.pathway_identity is not None]
    fmo_labels = [row.label.replace("FMO ", "").replace(" ", "\n") for row in fmo_rows]
    axes[1, 0].bar(fmo_labels, [float(row.pathway_identity) for row in fmo_rows], color="tab:blue")
    axes[1, 0].axhline(0.78, color="tab:red", linestyle="--", linewidth=1, label="target")
    axes[1, 0].set_ylim(0.0, 1.05)
    axes[1, 0].set_ylabel("pathway identity")
    axes[1, 0].set_title("FMO pathway-retention tradeoff")
    axes[1, 0].legend()
    axes[1, 0].tick_params(axis="x", labelsize=8)

    axes[1, 1].barh(labels, [row.active_null_z for row in ROWS], color="tab:purple")
    axes[1, 1].axvline(1.5, color="tab:red", linestyle="--", linewidth=1, label="strict z gate")
    axes[1, 1].set_xlabel("active null z")
    axes[1, 1].set_title("Null-specificity contrast")
    axes[1, 1].legend()

    fig.suptitle("Phase 56.5: Topology Sensitivity in HAOS-IIP Biology Telemetry", fontsize=14)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
