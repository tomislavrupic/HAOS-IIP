"""One-command reproduction and audit for Biology Line E.

This wrapper keeps the external bridge lab-facing:
- downloads only public, no-credential inputs when a full rerun is requested
- runs the bounded GDS112/GDS20 probes through the existing bridge script
- audits existing outputs with --check without requiring numpy, matplotlib, or network
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent.parent.parent
MANIFEST_PATH = ROOT / "input_manifest.json"
OUTPUT_ROOT = ROOT / "outputs"
RUNNER = ROOT / "run_external_dataset_bridge.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Reproduce or audit Biology Line E external bridge.")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Audit existing outputs only; do not download data or rerun probes.",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Require raw inputs to already exist under data/raw/ during a full rerun.",
    )
    parser.add_argument(
        "--only",
        choices=("all", "gds112", "gds20", "gds20_go_0006970"),
        default="all",
        help="Limit a full rerun or check to one bounded run.",
    )
    return parser.parse_args()


def load_manifest() -> dict[str, Any]:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def rel(path_text: str) -> Path:
    return ROOT / path_text


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "HAOS-IIP-Line-E-Reproducer/0.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        with destination.open("wb") as handle:
            shutil.copyfileobj(response, handle)


def input_maps(manifest: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    raw_inputs = {item["id"]: item for item in manifest["raw_inputs"]}
    references = {item["id"]: item for item in manifest["reference_inputs"]}
    return raw_inputs, references


def selected_runs(manifest: dict[str, Any], only: str) -> list[dict[str, Any]]:
    runs = list(manifest["bounded_runs"])
    if only == "all":
        return runs
    return [run for run in runs if run["id"] == only]


def ensure_raw_inputs(
    manifest: dict[str, Any],
    runs: list[dict[str, Any]],
    no_download: bool,
) -> list[dict[str, str]]:
    raw_inputs, references = input_maps(manifest)
    required_raw_ids = sorted({run["input_id"] for run in runs})
    required_reference_ids = sorted(
        {run["known_gene_set_id"] for run in runs if "known_gene_set_id" in run}
    )
    digests: list[dict[str, str]] = []

    for input_id in required_raw_ids:
        spec = raw_inputs[input_id]
        path = rel(spec["local_path"])
        if not path.exists():
            if no_download:
                raise FileNotFoundError(f"missing raw input and --no-download was set: {path}")
            print(f"downloading {input_id}: {spec['url']}")
            download_file(spec["url"], path)
        digest = sha256_file(path)
        expected = spec.get("sha256")
        if expected and digest != expected:
            raise ValueError(f"sha256 mismatch for {input_id}: expected {expected}, got {digest}")
        digests.append({"id": input_id, "path": display_path(path), "sha256": digest})

    for reference_id in required_reference_ids:
        spec = references[reference_id]
        path = rel(spec["local_path"])
        if not path.exists():
            raise FileNotFoundError(f"missing committed reference input: {path}")
        digest = sha256_file(path)
        expected = spec.get("sha256")
        if expected and digest != expected:
            raise ValueError(f"sha256 mismatch for {reference_id}: expected {expected}, got {digest}")
        digests.append({"id": reference_id, "path": display_path(path), "sha256": digest})

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "input_digests.json").write_text(
        json.dumps(digests, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return digests


def runner_command() -> list[str]:
    uv = shutil.which("uv")
    if uv:
        return [uv, "run", "--with", "numpy", "--with", "matplotlib", "python", str(RUNNER)]
    return [sys.executable, str(RUNNER)]


def run_probe(run: dict[str, Any], manifest: dict[str, Any]) -> None:
    raw_inputs, references = input_maps(manifest)
    input_spec = raw_inputs[run["input_id"]]
    command = [
        *runner_command(),
        "--soft",
        str(rel(input_spec["local_path"])),
        "--output-dir",
        str(rel(run["output_dir"])),
        "--validation-path",
        str(rel(run["validation_path"])),
    ]
    if "known_gene_set_id" in run:
        reference_spec = references[run["known_gene_set_id"]]
        command.extend(["--known-gene-set", str(rel(reference_spec["local_path"]))])

    print(f"running {run['id']}")
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def load_status(run: dict[str, Any]) -> dict[str, Any]:
    path = rel(run["output_dir"]) / "bridge_status.json"
    if not path.exists():
        raise FileNotFoundError(f"missing bridge status for {run['id']}: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def check(condition: bool, run_id: str, metric: str, detail: str) -> dict[str, Any]:
    return {
        "run_id": run_id,
        "metric": metric,
        "passed": bool(condition),
        "detail": detail,
    }


def audit_run(run: dict[str, Any], status: dict[str, Any]) -> list[dict[str, Any]]:
    expected = run["expected"]
    observed = status.get("observed", {})
    failure_analysis = status.get("failure_analysis", {})
    gene_set = status.get("gene_set_evaluation", {})
    checks: list[dict[str, Any]] = []
    run_id = run["id"]

    checks.append(
        check(
            status.get("bridge_status") == expected.get("bridge_status"),
            run_id,
            "bridge_status",
            f"expected {expected.get('bridge_status')}, got {status.get('bridge_status')}",
        )
    )
    checks.append(
        check(
            observed.get("baseline_stable") is expected.get("baseline_stable"),
            run_id,
            "baseline_stable",
            f"expected {expected.get('baseline_stable')}, got {observed.get('baseline_stable')}",
        )
    )
    checks.append(
        check(
            observed.get("early_detection") is expected.get("early_detection"),
            run_id,
            "early_detection",
            f"expected {expected.get('early_detection')}, got {observed.get('early_detection')}",
        )
    )
    checks.append(
        check(
            status.get("control_early_detection_count") == expected.get("control_early_detection_count"),
            run_id,
            "control_contrast",
            (
                f"expected control count {expected.get('control_early_detection_count')}, "
                f"got {status.get('control_early_detection_count')}"
            ),
        )
    )
    checks.append(
        check(
            failure_analysis.get("failure_mode") == expected.get("failure_mode"),
            run_id,
            "failure_mode",
            f"expected {expected.get('failure_mode')}, got {failure_analysis.get('failure_mode')}",
        )
    )

    if expected.get("k_star_before_visible"):
        k_idx = observed.get("k_star_index")
        visible_idx = observed.get("first_visible_failure_index")
        checks.append(
            check(
                k_idx is not None and visible_idx is not None and int(k_idx) < int(visible_idx),
                run_id,
                "k_star_before_visible",
                f"k_star_index={k_idx}, first_visible_failure_index={visible_idx}",
            )
        )

    if expected.get("gene_set_status"):
        checks.append(
            check(
                gene_set.get("gene_set_status") == expected["gene_set_status"],
                run_id,
                "gene_set_status",
                f"expected {expected['gene_set_status']}, got {gene_set.get('gene_set_status')}",
            )
        )
        checks.append(
            check(
                int(gene_set.get("matched_gene_count") or 0)
                >= int(expected.get("minimum_matched_gene_count", 0)),
                run_id,
                "matched_gene_count",
                (
                    f"expected at least {expected.get('minimum_matched_gene_count', 0)}, "
                    f"got {gene_set.get('matched_gene_count')}"
                ),
            )
        )
        margin = gene_set.get("average_precision_margin")
        minimum_margin = float(expected.get("minimum_average_precision_margin", 0.0))
        checks.append(
            check(
                bool(gene_set.get("haos_beats_fold_change"))
                and margin is not None
                and float(margin) >= minimum_margin,
                run_id,
                "gene_set_comparator",
                (
                    f"haos_beats_fold_change={gene_set.get('haos_beats_fold_change')}, "
                    f"margin={margin}, minimum={minimum_margin}"
                ),
            )
        )

    required_files = [
        "bridge_status.json",
        "results.csv",
        "standard_expression_metrics.csv",
        "metric_decomposition.csv",
        "control_summary.csv",
        "failure_analysis.md",
        "external_dataset_bridge_validation.md",
    ]
    output_dir = rel(run["output_dir"])
    for file_name in required_files:
        checks.append(
            check(
                (output_dir / file_name).exists(),
                run_id,
                f"artifact:{file_name}",
                display_path(output_dir / file_name),
            )
        )
    return checks


def write_probe_comparison(statuses: list[tuple[str, dict[str, Any]]]) -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_ROOT / "probe_comparison.csv"
    fieldnames = [
        "dataset",
        "bridge_status",
        "k_star_time",
        "first_visible_failure_time",
        "early_detection",
        "contrast_pass",
        "control_early_detection_count",
        "gene_set_status",
        "haos_average_precision",
        "fold_change_average_precision",
        "average_precision_margin",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for run_id, status in statuses:
            observed = status.get("observed", {})
            gene_set = status.get("gene_set_evaluation", {})
            writer.writerow(
                {
                    "dataset": run_id,
                    "bridge_status": status.get("bridge_status"),
                    "k_star_time": observed.get("k_star_time"),
                    "first_visible_failure_time": observed.get("first_visible_failure_time"),
                    "early_detection": observed.get("early_detection"),
                    "contrast_pass": status.get("contrast_pass"),
                    "control_early_detection_count": status.get("control_early_detection_count"),
                    "gene_set_status": gene_set.get("gene_set_status"),
                    "haos_average_precision": gene_set.get("haos_average_precision"),
                    "fold_change_average_precision": gene_set.get("fold_change_average_precision"),
                    "average_precision_margin": gene_set.get("average_precision_margin"),
                }
            )

    lines = [
        "# External Dataset Bridge Probe Comparison",
        "",
        "| dataset | bridge_status | k_star_time | first_visible_failure_time | early_detection | controls early | gene_set_status | AP margin |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for run_id, status in statuses:
        observed = status.get("observed", {})
        gene_set = status.get("gene_set_evaluation", {})
        lines.append(
            "| "
            + " | ".join(
                [
                    run_id,
                    str(status.get("bridge_status")),
                    str(observed.get("k_star_time")),
                    str(observed.get("first_visible_failure_time")),
                    str(observed.get("early_detection")),
                    str(status.get("control_early_detection_count")),
                    str(gene_set.get("gene_set_status")),
                    str(gene_set.get("average_precision_margin")),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "GDS112 is retained as an honest negative result. GDS20 is the current external pass with deterministic controls.",
            "The GO:0006970 row is a narrow sourced gene-set comparator, not a full osmotic-stress biology validation.",
        ]
    )
    (OUTPUT_ROOT / "probe_comparison.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_check_report(checks: list[dict[str, Any]], statuses: list[tuple[str, dict[str, Any]]]) -> bool:
    passed = all(bool(item["passed"]) for item in checks)
    report = {
        "check_status": "PASS" if passed else "FAIL",
        "checks_passed": sum(1 for item in checks if item["passed"]),
        "checks_total": len(checks),
        "checks": checks,
    }
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "quick_reproduce_check.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# Biology Line E Quick Reproduce Check",
        "",
        f"- check_status: {report['check_status']}",
        f"- checks_passed: {report['checks_passed']}",
        f"- checks_total: {report['checks_total']}",
        "",
        "## Dataset Status",
    ]
    for run_id, status in statuses:
        observed = status.get("observed", {})
        lines.extend(
            [
                f"- {run_id}: {status.get('bridge_status')}",
                f"  k_star_time: {observed.get('k_star_time')}",
                f"  first_visible_failure_time: {observed.get('first_visible_failure_time')}",
                f"  early_detection: {observed.get('early_detection')}",
            ]
        )
    lines.extend(["", "## Failed Checks"])
    failed = [item for item in checks if not item["passed"]]
    if failed:
        for item in failed:
            lines.append(f"- {item['run_id']} / {item['metric']}: {item['detail']}")
    else:
        lines.append("- none")
    (OUTPUT_ROOT / "quick_reproduce_check.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return passed


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def main() -> int:
    args = parse_args()
    manifest = load_manifest()
    runs = selected_runs(manifest, args.only)

    if not args.check:
        ensure_raw_inputs(manifest, runs, no_download=args.no_download)
        for run in runs:
            run_probe(run, manifest)

    statuses = [(run["id"], load_status(run)) for run in runs]
    if args.only == "all":
        write_probe_comparison(statuses)
    checks = [item for run, (_, status) in zip(runs, statuses) for item in audit_run(run, status)]
    passed = write_check_report(checks, statuses)

    print(f"quick_reproduce_check: {'PASS' if passed else 'FAIL'}")
    print(f"runs_checked: {len(runs)}")
    print(f"checks_passed: {sum(1 for item in checks if item['passed'])}/{len(checks)}")
    print(f"outputs_written: {display_path(OUTPUT_ROOT)}/")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
