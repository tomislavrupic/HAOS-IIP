"""Run Biology Line E external dataset bridge."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "outputs"
VALIDATION_PATH = ROOT / "external_dataset_bridge_validation.md"
REPO_ROOT = ROOT.parent.parent.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run HAOS-IIP Biology Line E external dataset bridge.")
    parser.add_argument("--csv", type=Path, help="Path to expression time-series CSV.")
    parser.add_argument("--soft", type=Path, help="Path to GEO SOFT or SOFT.gz dataset.")
    parser.add_argument("--known-gene-set", type=Path, help="Optional newline-delimited known response gene set.")
    parser.add_argument("--max-features", type=int, default=2000, help="Keep top variable features for diagnostics.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Directory for generated outputs.")
    parser.add_argument("--validation-path", type=Path, default=VALIDATION_PATH, help="Validation markdown path.")
    return parser.parse_args()


def display_path(path: Path) -> str:
    try:
        return f"{path.resolve().relative_to(REPO_ROOT.resolve())}/"
    except ValueError:
        return str(path.resolve())


def write_status_json(path: Path, status: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(status, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_open_no_data_report() -> None:
    status = {
        "bridge_status": "OPEN_NO_DATA",
        "primary_target": "GDS112 / GSE18 yeast heat-shock time course",
        "early_detection": None,
        "scope_note": "No external data file was supplied; this is a bridge scaffold, not a validation result.",
    }
    write_status_json(OUTPUT_DIR / "bridge_status.json", status)
    VALIDATION_PATH.write_text(
        "\n".join(
            [
                "# External Dataset Bridge Validation",
                "",
                "## System",
                "Biology Line E is an external dataset bridge for real perturbation datasets.",
                "The v0.1 primary target is NCBI GEO GDS112, a yeast heat-shock expression time course from GSE18.",
                "",
                "## Metrics",
                "When data are supplied, the bridge computes expression_match, step_coherence, trajectory_identity_retention, support_retention, recoverability, delta_persistence, k_star, safety_margin, confidence, and visible_failure proxy.",
                "The supplied-data path also runs deterministic negative controls and standard expression-change summaries.",
                "Known-gene-set ranking is optional and remains open until a curated external gene set is supplied.",
                "",
                "## Result",
                "- bridge_status: OPEN_NO_DATA",
                "- baseline_stable: not evaluated",
                "- k_star_time: not evaluated",
                "- first_visible_failure_time: not evaluated",
                "- early_detection: not evaluated",
                "",
                "## Limitations",
                "- This is not a toy pass.",
                "- No external validation claim is made without an input dataset.",
                "- Visible failure is an expression-pattern proxy unless an external phenotype label is supplied.",
                "- A future probe pass must beat deterministic negative controls.",
                "- A biology-specific ranking claim must beat simple fold-change against a curated external gene set.",
                "- This does not modify or prove HAOS-IIP core.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def plot_results(rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    times = [float(row["time_value"]) for row in rows]

    plt.figure(figsize=(7, 4))
    plt.plot(times, [float(row["recoverability"]) for row in rows], marker="o")
    plt.axhline(0.70, color="firebrick", linestyle="--", linewidth=1, label="collapse threshold")
    plt.xlabel("time / ordered sample")
    plt.ylabel("recoverability")
    plt.title("Recoverability vs external perturbation time")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "recoverability_vs_time.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4))
    plt.plot(times, [float(row["delta_persistence"]) for row in rows], marker="o", color="darkorange")
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.xlabel("time / ordered sample")
    plt.ylabel("delta persistence")
    plt.title("Delta persistence vs external perturbation time")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "delta_persistence_vs_time.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4))
    plt.plot(times, [float(row["expression_match"]) for row in rows], marker="o", label="expression match")
    plt.plot(times, [float(row["step_coherence"]) for row in rows], marker="s", label="step coherence")
    plt.axhline(0.50, color="firebrick", linestyle="--", linewidth=1, label="visible proxy threshold")
    plt.xlabel("time / ordered sample")
    plt.ylabel("retention")
    plt.title("Expression match and step coherence")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "expression_match_and_step_coherence.png", dpi=160)
    plt.close()


def plot_control_summary(runs: list[object]) -> None:
    import matplotlib.pyplot as plt

    labels = [run.label for run in runs]
    k_values = [
        float(run.summary["k_star_time"]) if run.summary["k_star_time"] is not None else float("nan")
        for run in runs
    ]
    visible_values = [
        float(run.summary["first_visible_failure_time"])
        if run.summary["first_visible_failure_time"] is not None
        else float("nan")
        for run in runs
    ]

    plt.figure(figsize=(8, 4))
    x_values = range(len(labels))
    plt.scatter(x_values, k_values, label="k_star time", s=60)
    plt.scatter(x_values, visible_values, label="visible proxy time", s=60, marker="x")
    plt.xticks(list(x_values), labels, rotation=20, ha="right")
    plt.ylabel("time / ordered sample")
    plt.title("Observed bridge vs deterministic controls")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "observed_vs_controls.png", dpi=160)
    plt.close()


def plot_failure_decomposition(rows: list[dict[str, object]], standard_metrics: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    times = [float(row["time_value"]) for row in rows]
    plt.figure(figsize=(8, 4))
    plt.plot(times, [float(row["recoverability"]) for row in rows], marker="o", label="recoverability")
    plt.plot(times, [float(row["expression_match"]) for row in rows], marker="s", label="expression match")
    plt.plot(times, [float(row["step_coherence"]) for row in rows], marker="^", label="step coherence")
    plt.plot(
        times,
        [float(row["fraction_abs_change_ge_1_0"]) for row in standard_metrics],
        marker="x",
        label="fraction abs change >= 1.0",
    )
    plt.axhline(0.70, color="firebrick", linestyle="--", linewidth=1, label="collapse threshold")
    plt.axhline(0.50, color="gray", linestyle=":", linewidth=1, label="visible/comparator threshold")
    plt.xlabel("time / ordered sample")
    plt.ylabel("value")
    plt.title("Line E failure decomposition")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "failure_decomposition.png", dpi=160)
    plt.close()


def plot_threshold_sensitivity(rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    collapse_values = sorted({float(row["collapse_threshold"]) for row in rows})
    visible_values = sorted({float(row["visible_match_threshold"]) for row in rows})
    matrix = np.zeros((len(collapse_values), len(visible_values)))

    for row in rows:
        c_idx = collapse_values.index(float(row["collapse_threshold"]))
        v_idx = visible_values.index(float(row["visible_match_threshold"]))
        matrix[c_idx, v_idx] = 1.0 if bool(row["early_detection"]) else 0.0

    plt.figure(figsize=(6, 4))
    plt.imshow(matrix, cmap="Greys", vmin=0, vmax=1, aspect="auto")
    plt.xticks(range(len(visible_values)), [str(value) for value in visible_values])
    plt.yticks(range(len(collapse_values)), [str(value) for value in collapse_values])
    plt.xlabel("visible match threshold")
    plt.ylabel("collapse threshold")
    plt.title("Early-detection threshold sensitivity")
    plt.colorbar(label="early detection")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "threshold_sensitivity.png", dpi=160)
    plt.close()


def write_failure_analysis_markdown(path: Path, analysis: dict[str, object]) -> None:
    path.write_text(
        "\n".join(
            [
                "# External Bridge Failure Analysis",
                "",
                f"- failure_mode: {analysis['failure_mode']}",
                f"- interpretation: {analysis['interpretation']}",
                f"- k_star_index: {analysis['k_star_index']}",
                f"- k_star_time: {analysis['k_star_time']}",
                f"- first_visible_index: {analysis['first_visible_index']}",
                f"- first_visible_time: {analysis['first_visible_time']}",
                f"- first_standard_large_change_index: {analysis['first_standard_large_change_index']}",
                f"- first_standard_large_change_time: {analysis['first_standard_large_change_time']}",
                f"- recoverability_expression_match_correlation: {analysis['recoverability_expression_match_correlation']}",
                f"- recoverability_step_coherence_correlation: {analysis['recoverability_step_coherence_correlation']}",
                "",
                "## Reading",
                str(analysis["dominant_failure_note"]),
                "",
                "This report is diagnostic only. It does not change the failed probe result.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def write_validation(summary: dict[str, object], runs: list[object]) -> None:
    observed = runs[0]
    gene_set = summary.get("gene_set_evaluation", {"gene_set_status": "OPEN_GENE_SET_NOT_SUPPLIED"})
    failure_analysis = summary.get("failure_analysis", {})
    control_lines = [
        (
            f"- {run.label}: early_detection={run.summary['early_detection']}, "
            f"k_star_time={run.summary['k_star_time']}, "
            f"first_visible_failure_time={run.summary['first_visible_failure_time']}"
        )
        for run in runs[1:]
    ]

    VALIDATION_PATH.write_text(
        "\n".join(
            [
                "# External Dataset Bridge Validation",
                "",
                "## System",
                f"Source: {observed.summary['source_name']}",
                f"Features analyzed: {observed.summary['feature_count']}",
                f"Samples analyzed: {observed.summary['sample_count']}",
                "Primary v0.1 target: public yeast heat-shock expression time course if GDS112 is supplied.",
                "",
                "## Metrics",
                "recoverability is a weighted proxy combining expression_match, step_coherence, and support_retention.",
                "delta_persistence is the change in recoverability between ordered samples.",
                "k_star is the first sustained recoverability crossing below 0.70.",
                "safety_margin is distance in time or ordered-sample units to k_star.",
                "visible_failure is a strict expression-pattern or trajectory-identity proxy, not a biological failure label.",
                "standard expression metrics report ordinary absolute expression movement at 0.5 and 1.0 log-ratio thresholds.",
                "trajectory_identity_retention is a standard time-series comparator used to catch identity-breaking controls.",
                "feature_rankings.csv ranks genes by HAOS-weighted loss score and simple fold-change comparators.",
                "",
                "## Negative Controls",
                "The bridge runs deterministic time-shuffle, feature-shuffle, and row-permutation controls.",
                "A probe pass requires observed early detection and zero matching early-detection hits in these controls.",
                *control_lines,
                "",
                "## Known Gene-Set Comparator",
                f"- gene_set_status: {gene_set['gene_set_status']}",
                f"- matched_gene_count: {gene_set.get('matched_gene_count')}",
                f"- haos_average_precision: {gene_set.get('haos_average_precision')}",
                f"- fold_change_average_precision: {gene_set.get('fold_change_average_precision')}",
                f"- haos_beats_fold_change: {gene_set.get('haos_beats_fold_change')}",
                "",
                "## Result",
                f"- bridge_status: {summary['bridge_status']}",
                f"- baseline_stable: {observed.summary['baseline_stable']}",
                f"- k_star_time: {observed.summary['k_star_time']}",
                f"- first_visible_failure_time: {observed.summary['first_visible_failure_time']}",
                f"- early_detection: {observed.summary['early_detection']}",
                f"- contrast_pass: {summary['contrast_pass']}",
                f"- control_early_detection_count: {summary['control_early_detection_count']}",
                "",
                "## Failure Analysis",
                f"- failure_mode: {failure_analysis.get('failure_mode')}",
                f"- interpretation: {failure_analysis.get('interpretation')}",
                f"- recoverability_expression_match_correlation: {failure_analysis.get('recoverability_expression_match_correlation')}",
                f"- first_standard_large_change_time: {failure_analysis.get('first_standard_large_change_time')}",
                "",
                "## Explicit Failure Conditions",
                "- The bridge fails if observed data do not show early detection.",
                "- The bridge fails if deterministic controls show the same early-detection signature.",
                "- A biology-specific claim fails if HAOS feature ranking does not beat simple fold-change ranking against a supplied curated gene set.",
                "- The bridge remains open if visible_failure is only a proxy and no external phenotype label is supplied.",
                "",
                "## Limitations",
                "- External expression data do not by themselves provide visible biological failure labels.",
                "- This bridge tests diagnostic ordering, not biological mechanism.",
                "- This is not proof that HAOS-IIP explains biology.",
                "- This does not modify or prove HAOS-IIP core.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    global OUTPUT_DIR, VALIDATION_PATH

    args = parse_args()
    OUTPUT_DIR = args.output_dir
    VALIDATION_PATH = args.validation_path
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not args.csv and not args.soft:
        write_open_no_data_report()
        print("bridge_status: OPEN_NO_DATA")
        print("early_detection: not_evaluated")
        print(f"outputs_written: {display_path(OUTPUT_DIR)}")
        return 0

    if args.csv and args.soft:
        print("Use only one of --csv or --soft.", file=sys.stderr)
        return 2

    from external_dataset_bridge import (
        load_csv_expression,
        load_geo_soft_expression,
        compute_feature_rankings,
        analyze_failure_modes,
        evaluate_gene_set,
        load_gene_set,
        run_bridge_with_controls,
        select_top_variable_features,
        threshold_sensitivity_rows,
        write_control_summary_csv,
        write_feature_rankings_csv,
        write_metric_decomposition_csv,
        write_results_csv,
        write_standard_metrics_csv,
        write_threshold_sensitivity_csv,
    )

    dataset = load_csv_expression(args.csv) if args.csv else load_geo_soft_expression(args.soft)
    dataset = select_top_variable_features(dataset, max_features=args.max_features)
    bundle = run_bridge_with_controls(dataset)
    runs = bundle["runs"]
    observed = runs[0]
    status = {
        key: value
        for key, value in bundle.items()
        if key != "runs"
    }
    status["observed"] = observed.summary
    status["controls"] = [run.summary for run in runs[1:]]

    feature_rankings = compute_feature_rankings(observed.dataset, observed.rows)
    failure_analysis = analyze_failure_modes(observed.rows, observed.standard_metrics)
    threshold_rows = threshold_sensitivity_rows(observed.rows)
    status["failure_analysis"] = failure_analysis
    if args.known_gene_set:
        gene_set = load_gene_set(args.known_gene_set)
        status["gene_set_evaluation"] = evaluate_gene_set(feature_rankings, gene_set)
    else:
        status["gene_set_evaluation"] = {
            "gene_set_status": "OPEN_GENE_SET_NOT_SUPPLIED",
            "failure_rule": "A biology-specific ranking claim is not evaluated until a curated known-gene-set file is supplied.",
        }

    write_results_csv(OUTPUT_DIR / "results.csv", observed.rows)
    write_standard_metrics_csv(OUTPUT_DIR / "standard_expression_metrics.csv", observed.standard_metrics)
    write_metric_decomposition_csv(OUTPUT_DIR / "metric_decomposition.csv", observed.rows, observed.standard_metrics)
    write_threshold_sensitivity_csv(OUTPUT_DIR / "threshold_sensitivity.csv", threshold_rows)
    write_feature_rankings_csv(OUTPUT_DIR / "feature_rankings.csv", feature_rankings)
    write_control_summary_csv(OUTPUT_DIR / "control_summary.csv", runs)
    write_status_json(OUTPUT_DIR / "bridge_status.json", status)
    write_status_json(OUTPUT_DIR / "failure_analysis.json", failure_analysis)
    write_failure_analysis_markdown(OUTPUT_DIR / "failure_analysis.md", failure_analysis)
    plot_results(observed.rows)
    plot_control_summary(runs)
    plot_failure_decomposition(observed.rows, observed.standard_metrics)
    plot_threshold_sensitivity(threshold_rows)
    write_validation(status, runs)

    print(f"bridge_status: {status['bridge_status']}")
    print(f"baseline_stable: {observed.summary['baseline_stable']}")
    print(f"k_star_time: {observed.summary['k_star_time']}")
    print(f"first_visible_failure_time: {observed.summary['first_visible_failure_time']}")
    print(f"early_detection: {observed.summary['early_detection']}")
    print(f"contrast_pass: {status['contrast_pass']}")
    print(f"control_early_detection_count: {status['control_early_detection_count']}")
    print(f"gene_set_status: {status['gene_set_evaluation']['gene_set_status']}")
    print(f"failure_mode: {failure_analysis['failure_mode']}")
    print(f"outputs_written: {display_path(OUTPUT_DIR)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
