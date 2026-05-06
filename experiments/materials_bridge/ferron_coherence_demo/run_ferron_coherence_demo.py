#!/usr/bin/env python3
"""Run Materials Bridge Line A - Ferron Coherence Recoverability Probe."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ferron_data_loader import (
    MANIFEST_PATH,
    OUTPUTS_DIR,
    PRIMARY_DATA_DOI,
    PRIMARY_DATA_URL,
    PRIMARY_PAPER_DOI,
    build_smoke_test_dataset_if_needed,
    download_or_locate_dataset,
    inspect_dataset_files,
    load_frequency_spectra,
    load_stft_amplitudes,
    load_time_traces,
)
from ferron_recoverability_model import run_ferron_sweep
from spectral_feature_audit import run_spectral_feature_audit
from stft_feature_audit import run_stft_feature_audit


ROOT = Path(__file__).resolve().parent
RESULTS_PATH = OUTPUTS_DIR / "results.csv"
SUMMARY_PATH = OUTPUTS_DIR / "summary.json"
INSPECTION_PATH = OUTPUTS_DIR / "data_inspection.json"
VALIDATION_PATH = ROOT / "ferron_coherence_validation.md"

FIELDNAMES = [
    "source_mode",
    "sample_id",
    "condition_index",
    "condition_type",
    "condition_value",
    "distance_um",
    "time_ps",
    "frequency_THz",
    "peak_amplitude",
    "peak_frequency_THz",
    "linewidth_proxy",
    "ferron_peak_recoverability",
    "propagation_recoverability",
    "temporal_recoverability",
    "combined_recoverability",
    "delta_persistence",
    "k_star",
    "safety_margin",
    "confidence",
    "visible_failure",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Skip public data download and use existing raw files or smoke-test fallback.",
    )
    parser.add_argument(
        "--force-smoke-test",
        action="store_true",
        help="Bypass external data and run the synthetic software smoke test.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, Any] | None = None
    inspection: dict[str, Any]
    time_traces: list[dict[str, Any]]
    frequency_spectra: list[dict[str, Any]]
    stft_amplitudes: list[dict[str, Any]]
    smoke_files: list[str] = []

    if args.force_smoke_test:
        smoke = build_smoke_test_dataset_if_needed(ROOT)
        inspection = smoke["inspection"]
        time_traces = smoke["time_traces"]
        frequency_spectra = smoke["frequency_spectra"]
        stft_amplitudes = smoke["stft_amplitudes"]
        smoke_files = smoke["files"]
        data_status = "DATA_UNAVAILABLE_SMOKE_TEST_ONLY"
        source_mode = "SMOKE_TEST"
    else:
        manifest = download_or_locate_dataset(ROOT, allow_download=not args.no_download)
        inspection = inspect_dataset_files(ROOT)
        time_traces = load_time_traces(ROOT, inspection)
        frequency_spectra = load_frequency_spectra(ROOT, inspection)
        stft_amplitudes = load_stft_amplitudes(ROOT, inspection)
        parsed_count = len(time_traces) + len(frequency_spectra) + len(stft_amplitudes)
        raw_file_count = len(manifest.get("files", [])) if manifest else 0

        if parsed_count > 0:
            data_status = "REAL_DATA_LOADED"
            source_mode = "REAL_DATA_LOADED"
        elif raw_file_count > 0:
            data_status = "PARSE_FAILED"
            source_mode = "REAL_DATA_UNPARSED"
        else:
            smoke = build_smoke_test_dataset_if_needed(ROOT)
            inspection = smoke["inspection"]
            time_traces = smoke["time_traces"]
            frequency_spectra = smoke["frequency_spectra"]
            stft_amplitudes = smoke["stft_amplitudes"]
            smoke_files = smoke["files"]
            data_status = "DATA_UNAVAILABLE_SMOKE_TEST_ONLY"
            source_mode = "SMOKE_TEST"

    INSPECTION_PATH.write_text(json.dumps(inspection, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    rows: list[dict[str, Any]] = []
    if data_status != "PARSE_FAILED":
        rows = run_ferron_sweep(
            time_traces=time_traces,
            frequency_spectra=frequency_spectra,
            stft_amplitudes=stft_amplitudes,
            source_mode=source_mode,
        )
        if not rows:
            data_status = "PARSE_FAILED" if source_mode == "REAL_DATA_LOADED" else "DATA_UNAVAILABLE_SMOKE_TEST_ONLY"

    validation_status = _validation_status(data_status, rows)
    summary = _build_summary(
        data_status=data_status,
        validation_status=validation_status,
        rows=rows,
        manifest=manifest,
        inspection=inspection,
        smoke_files=smoke_files,
        time_traces=time_traces,
        frequency_spectra=frequency_spectra,
        stft_amplitudes=stft_amplitudes,
    )

    _write_csv(RESULTS_PATH, rows)
    if data_status == "REAL_DATA_LOADED":
        _write_csv(OUTPUTS_DIR / "real_data_metrics.csv", rows)
    elif data_status == "DATA_UNAVAILABLE_SMOKE_TEST_ONLY":
        _write_csv(OUTPUTS_DIR / "smoke_test_metrics.csv", rows)

    _write_plots(
        rows=rows,
        time_traces=time_traces,
        frequency_spectra=frequency_spectra,
        smoke_test=data_status == "DATA_UNAVAILABLE_SMOKE_TEST_ONLY",
    )
    spectral_summary = run_spectral_feature_audit(
        frequency_spectra if data_status == "REAL_DATA_LOADED" else [],
        time_traces=time_traces if data_status == "REAL_DATA_LOADED" else [],
        output_dir=OUTPUTS_DIR,
        stft_count=len(stft_amplitudes),
    )
    summary["spectral_feature_audit"] = spectral_summary
    stft_summary = run_stft_feature_audit(root=ROOT, output_dir=OUTPUTS_DIR)
    summary["stft_time_frequency_audit"] = stft_summary
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_validation_markdown(summary, rows, manifest, inspection, spectral_summary, stft_summary)
    _print_console_summary(summary)
    return 0


def _validation_status(data_status: str, rows: list[dict[str, Any]]) -> str:
    if data_status == "DATA_UNAVAILABLE_SMOKE_TEST_ONLY":
        return "DATA_UNAVAILABLE_SMOKE_TEST_ONLY"
    if data_status == "REAL_DATA_LOADED" and rows:
        return "PASS"
    return "FAIL"


def _build_summary(
    *,
    data_status: str,
    validation_status: str,
    rows: list[dict[str, Any]],
    manifest: dict[str, Any] | None,
    inspection: dict[str, Any],
    smoke_files: list[str],
    time_traces: list[dict[str, Any]],
    frequency_spectra: list[dict[str, Any]],
    stft_amplitudes: list[dict[str, Any]],
) -> dict[str, Any]:
    baseline = _first_metric(rows, "combined_recoverability")
    k_star = _first_non_empty(rows, "k_star")
    first_visible_failure = _first_visible_failure(rows)
    early_detection = _early_detection(rows, k_star, first_visible_failure)
    files_loaded = _files_loaded(data_status, manifest, inspection, smoke_files)

    return {
        "experiment": "Materials Bridge Line A - Ferron Coherence Recoverability Probe",
        "data_status": data_status,
        "validation_status": validation_status,
        "source_doi": PRIMARY_DATA_DOI if data_status == "REAL_DATA_LOADED" else None,
        "source_url": PRIMARY_DATA_URL if data_status == "REAL_DATA_LOADED" else None,
        "primary_paper_doi": PRIMARY_PAPER_DOI,
        "files_loaded": files_loaded,
        "parsed_signals": {
            "time_traces": len(time_traces),
            "frequency_spectra_or_peak_records": len(frequency_spectra),
            "stft_amplitudes": len(stft_amplitudes),
        },
        "baseline_recoverability": baseline,
        "k_star_level": k_star,
        "first_visible_failure_level": first_visible_failure,
        "early_detection": early_detection,
        "confidence_summary": _confidence_summary(rows),
        "outputs_written": str(OUTPUTS_DIR),
        "claim_boundary": [
            "external-data probe only",
            "no proof of HAOS-IIP",
            "no claim that ferrons embody HAOS-IIP",
            "no new materials theory",
            "no consciousness claim",
        ],
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _csv_value(row.get(field)) for field in FIELDNAMES})


def _write_plots(
    *,
    rows: list[dict[str, Any]],
    time_traces: list[dict[str, Any]],
    frequency_spectra: list[dict[str, Any]],
    smoke_test: bool,
) -> None:
    if not rows:
        return
    prefix = "smoke_test_" if smoke_test else ""
    metric_plots = {
        "peak_amplitude_vs_condition.png": "peak_amplitude",
        "peak_frequency_vs_condition.png": "peak_frequency_THz",
        "linewidth_vs_condition.png": "linewidth_proxy",
        "recoverability_vs_condition.png": "combined_recoverability",
        "delta_persistence_vs_condition.png": "delta_persistence",
        "safety_margin_vs_condition.png": "safety_margin",
    }
    for filename, metric in metric_plots.items():
        _plot_metric(rows, metric, OUTPUTS_DIR / f"{prefix}{filename}")
    _plot_spectrum_examples(frequency_spectra, OUTPUTS_DIR / f"{prefix}spectrum_examples.png")
    _plot_time_trace_examples(time_traces, OUTPUTS_DIR / f"{prefix}time_trace_examples.png")


def _plot_metric(rows: list[dict[str, Any]], metric: str, path: Path) -> None:
    series = []
    for row in rows:
        y_value = _float_or_none(row.get(metric))
        if y_value is None:
            continue
        x_value = _float_or_none(row.get("condition_value"))
        if x_value is None:
            x_value = _float_or_none(row.get("condition_index"))
        if x_value is None:
            continue
        series.append((x_value, y_value, str(row.get("sample_id", "sample"))))
    if not series:
        return
    plt.figure(figsize=(7.2, 4.4))
    labels = sorted(set(label for _, _, label in series))
    for label in labels:
        subset = sorted((x, y) for x, y, item_label in series if item_label == label)
        plt.plot([x for x, _ in subset], [y for _, y in subset], marker="o", linewidth=1.6, label=label)
    plt.xlabel("condition value")
    plt.ylabel(metric)
    if len(labels) <= 6:
        plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def _plot_spectrum_examples(records: list[dict[str, Any]], path: Path) -> None:
    candidates = [record for record in records if record.get("kind") == "frequency_spectrum"]
    if not candidates:
        return
    plt.figure(figsize=(7.2, 4.4))
    for record in candidates[:5]:
        frequency = record.get("frequency_THz")
        amplitude = record.get("amplitude")
        if frequency is None or amplitude is None:
            continue
        label = f"{record.get('condition_type')}={record.get('condition_value')}"
        plt.plot(frequency, amplitude, linewidth=1.2, label=label)
    plt.xlabel("frequency THz")
    plt.ylabel("amplitude")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def _plot_time_trace_examples(records: list[dict[str, Any]], path: Path) -> None:
    if not records:
        return
    plt.figure(figsize=(7.2, 4.4))
    for record in records[:5]:
        label = f"{record.get('condition_type')}={record.get('condition_value')}"
        plt.plot(record.get("time_ps"), record.get("signal"), linewidth=1.0, label=label)
    plt.xlabel("time ps")
    plt.ylabel("signal")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def _write_validation_markdown(
    summary: dict[str, Any],
    rows: list[dict[str, Any]],
    manifest: dict[str, Any] | None,
    inspection: dict[str, Any],
    spectral_summary: dict[str, Any],
    stft_summary: dict[str, Any],
) -> None:
    provenance_lines = []
    if summary["data_status"] == "REAL_DATA_LOADED" and manifest:
        provenance_lines.append(f"- Source DOI: {PRIMARY_DATA_DOI}")
        provenance_lines.append(f"- Source URL: {PRIMARY_DATA_URL}")
        for file_record in manifest.get("files", []):
            provenance_lines.append(
                "- "
                f"{file_record.get('file_name')} | "
                f"sha256={file_record.get('sha256')} | "
                f"source={file_record.get('source_url')}"
            )
    elif summary["data_status"] == "PARSE_FAILED" and manifest:
        provenance_lines.append(f"- Source DOI: {PRIMARY_DATA_DOI}")
        provenance_lines.append(f"- Source URL: {PRIMARY_DATA_URL}")
        provenance_lines.append("- Raw files were present or downloaded, but no supported ferron signal parsed.")
        for file_record in manifest.get("files", []):
            provenance_lines.append(
                "- "
                f"{file_record.get('file_name')} | "
                f"sha256={file_record.get('sha256')}"
            )
    else:
        provenance_lines.append("- External dataset was not downloaded automatically.")
        provenance_lines.append("- Synthetic smoke-test files are software verification inputs only.")

    parsed = summary["parsed_signals"]
    parsed_lines = [
        f"- Time traces: {parsed['time_traces']}",
        f"- Frequency spectra or peak records: {parsed['frequency_spectra_or_peak_records']}",
        f"- STFT amplitudes: {parsed['stft_amplitudes']}",
    ]
    if not rows and inspection.get("files"):
        parsed_lines.append("- Parsed ferron-related signals: 0")
        parsed_lines.append("- Unsupported or unrecognized files:")
        for file_info in inspection.get("files", []):
            if file_info.get("parse_status") != "candidate":
                parsed_lines.append(f"  - {file_info.get('relative_path')} ({file_info.get('extension')})")
    spectral_lines = _spectral_markdown_lines(spectral_summary)
    stft_lines = _stft_markdown_lines(stft_summary)

    content = f"""# Ferron Coherence Validation

## Status
{summary["validation_status"]}

## Data Provenance
{chr(10).join(provenance_lines)}

## Parsed Signals
{chr(10).join(parsed_lines)}

## HAOS Mapping
- carrier: ferroelectric polarization order
- perturbation: laser pulse, fluence, distance, temperature, sample thickness, boundary/disorder effects
- recoverable structure: coherent ferron / polarization-wave mode
- observable readout: narrow-band THz emission, 3.13 THz amplitude, transient reflectance oscillation, Fourier/STFT peak
- visible failure: loss of narrow-band signal, spectral broadening, amplitude collapse, incoherent background dominance, or propagation failure

## Results
- baseline_recoverability: {_format_optional(summary["baseline_recoverability"])}
- k_star_level: {_format_optional(summary["k_star_level"])}
- first_visible_failure_level: {_format_optional(summary["first_visible_failure_level"])}
- early_detection: {summary["early_detection"]}
- confidence summary: {json.dumps(summary["confidence_summary"], sort_keys=True)}

## Spectral Feature Audit
{chr(10).join(spectral_lines)}

## STFT / Time-Frequency Audit
{chr(10).join(stft_lines)}

## Limitations
- external-data probe only
- no proof of HAOS-IIP
- no claim that ferrons embody HAOS-IIP
- no new materials theory
- no consciousness claim
- metric proxies are provisional
- results depend on available data quality
"""
    VALIDATION_PATH.write_text(content, encoding="utf-8")


def _print_console_summary(summary: dict[str, Any]) -> None:
    spectral = summary.get("spectral_feature_audit") or {}
    stft = summary.get("stft_time_frequency_audit") or {}
    print(f"data_status: {summary['data_status']}")
    print(f"source_doi: {_format_optional(summary['source_doi'])}")
    print(f"files_loaded: {len(summary['files_loaded'])}")
    print(f"baseline_recoverability: {_format_optional(summary['baseline_recoverability'])}")
    print(f"k_star_level: {_format_optional(summary['k_star_level'])}")
    print(f"first_visible_failure_level: {_format_optional(summary['first_visible_failure_level'])}")
    print(f"early_detection: {summary['early_detection']}")
    print(f"validation_status: {summary['validation_status']}")
    print(f"spectral_audit_status: {_format_optional(spectral.get('status'))}")
    print(f"spectral_records_audited: {_format_optional(spectral.get('records_audited'))}")
    print(f"target_peak_records_found: {_format_optional(spectral.get('peaks_found'))}")
    print(f"mean_peak_frequency_THz: {_format_optional(spectral.get('mean_peak_frequency_THz'))}")
    print(f"mean_spectral_recoverability: {_format_optional(spectral.get('mean_spectral_recoverability'))}")
    print(f"k_star_detected: {bool(spectral.get('k_star_detected'))}")
    print(f"bounded_interpretation: {_format_optional(spectral.get('bounded_interpretation'))}")
    print(f"stft_audit_status: {_format_optional(stft.get('status'))}")
    print(f"stft_maps_found: {_format_optional(stft.get('stft_maps_found', 0))}")
    print(f"stft_time_points_audited: {_format_optional(stft.get('time_points_audited', 0))}")
    print(f"stft_target_peak_records_found: {_format_optional(stft.get('target_peak_records_found', 0))}")
    print(f"mean_stft_recoverability: {_format_optional(stft.get('mean_stft_recoverability'))}")
    print(f"stft_k_star_detected: {bool(stft.get('k_star_detected'))}")
    print(f"stft_bounded_interpretation: {_format_optional(stft.get('bounded_interpretation'))}")
    print("outputs_written: experiments/materials_bridge/ferron_coherence_demo/outputs/")


def _files_loaded(
    data_status: str,
    manifest: dict[str, Any] | None,
    inspection: dict[str, Any],
    smoke_files: list[str],
) -> list[dict[str, Any]]:
    if data_status == "DATA_UNAVAILABLE_SMOKE_TEST_ONLY":
        return [{"file_name": Path(path).name, "path": path, "mode": "synthetic_smoke_test"} for path in smoke_files]
    if manifest and manifest.get("files"):
        return list(manifest["files"])
    return [
        {
            "file_name": file_info.get("file_name"),
            "relative_path": file_info.get("relative_path"),
            "parse_status": file_info.get("parse_status"),
        }
        for file_info in inspection.get("files", [])
    ]


def _first_metric(rows: list[dict[str, Any]], key: str) -> float | None:
    for row in rows:
        value = _float_or_none(row.get(key))
        if value is not None:
            return value
    return None


def _first_non_empty(rows: list[dict[str, Any]], key: str) -> Any | None:
    for row in rows:
        value = row.get(key)
        if value not in (None, ""):
            return value
    return None


def _first_visible_failure(rows: list[dict[str, Any]]) -> Any | None:
    for row in rows:
        if bool(row.get("visible_failure")):
            return row.get("condition_value")
    return None


def _early_detection(rows: list[dict[str, Any]], k_star: Any | None, first_visible_failure: Any | None) -> bool:
    if k_star is None:
        return False
    if first_visible_failure is None:
        return False
    k_index = _first_condition_index(rows, k_star)
    failure_index = _first_condition_index(rows, first_visible_failure)
    if k_index is None or failure_index is None:
        return False
    return k_index < failure_index


def _first_condition_index(rows: list[dict[str, Any]], condition_value: Any) -> int | None:
    for row in rows:
        if row.get("condition_value") == condition_value:
            return int(row.get("condition_index", 0))
    return None


def _confidence_summary(rows: list[dict[str, Any]]) -> dict[str, float | None]:
    values = [_float_or_none(row.get("confidence")) for row in rows]
    values = [value for value in values if value is not None]
    if not values:
        return {"min": None, "mean": None, "max": None}
    return {
        "min": min(values),
        "mean": sum(values) / len(values),
        "max": max(values),
    }


def _spectral_markdown_lines(spectral_summary: dict[str, Any]) -> list[str]:
    if not spectral_summary:
        return ["- Spectral audit was not run."]
    window = _format_optional(spectral_summary.get("search_window_THz"))
    lines = [
        f"- status: {_format_optional(spectral_summary.get('status'))}",
        f"- spectral records audited: {_format_optional(spectral_summary.get('records_audited'))}",
        f"- target peak search window: 3.13 THz +/- {window} THz",
        f"- target peak records found: {_format_optional(spectral_summary.get('peaks_found'))}",
        f"- strongest / most stable peak frequency summary: mean {_format_optional(spectral_summary.get('mean_peak_frequency_THz'))} THz",
        f"- mean spectral recoverability: {_format_optional(spectral_summary.get('mean_spectral_recoverability'))}",
        f"- k_star remains absent: {not bool(spectral_summary.get('k_star_detected'))}",
        f"- visible failure detected: {bool(spectral_summary.get('visible_failure_detected'))}",
    ]
    features = spectral_summary.get("standout_features") or []
    for feature in features:
        lines.append(f"- {feature}")
    missing = spectral_summary.get("missing_data_notes") or []
    if missing:
        lines.append("- missing data notes:")
        for note in missing[:8]:
            lines.append(f"  - {note}")
    interpretation = spectral_summary.get("bounded_interpretation")
    if interpretation:
        lines.append(f"- bounded interpretation: {interpretation}")
    return lines


def _stft_markdown_lines(stft_summary: dict[str, Any]) -> list[str]:
    if not stft_summary:
        return ["- STFT audit was not run."]
    lines = [
        f"- status: {_format_optional(stft_summary.get('status'))}",
        f"- raw files inspected: {_format_optional(stft_summary.get('raw_files_inspected'))}",
        f"- candidate sheets inspected: {_format_optional(stft_summary.get('candidate_sheets_inspected'))}",
        f"- usable STFT maps found: {_format_optional(stft_summary.get('stft_maps_found', 0))}",
        f"- target frequency search window: 3.13 THz +/- {_format_optional(stft_summary.get('search_window_THz', 0.2))} THz",
    ]
    accepted = stft_summary.get("accepted_maps") or []
    if accepted:
        lines.append("- accepted STFT maps:")
        for item in accepted:
            lines.append(
                "  - "
                f"{item.get('source_file')}#{item.get('sheet_name')} | "
                f"time_points={item.get('time_points')} | "
                f"frequency_points={item.get('frequency_points')}"
            )
        lines.extend(
            [
                f"- STFT time points audited: {_format_optional(stft_summary.get('time_points_audited'))}",
                f"- target peak records found: {_format_optional(stft_summary.get('target_peak_records_found'))}",
                f"- mean peak frequency: {_format_optional(stft_summary.get('mean_peak_frequency_THz'))} THz",
                f"- mean STFT recoverability: {_format_optional(stft_summary.get('mean_stft_recoverability'))}",
                f"- STFT k_star detected: {bool(stft_summary.get('k_star_detected'))}",
                f"- visible failure detected: {bool(stft_summary.get('visible_failure_detected'))}",
            ]
        )
    else:
        lines.append("- usable STFT data found: False")
        lines.append(
            "- No usable raw STFT/time-frequency map was found in the downloaded XLSX files. "
            "No STFT recoverability claim is made."
        )
    rejected = stft_summary.get("rejected_candidate_sheets") or []
    if rejected:
        lines.append("- rejected candidate sheets:")
        for item in rejected[:12]:
            lines.append(f"  - {item.get('source_file')}#{item.get('sheet_name')}: {item.get('reason')}")
    missing = stft_summary.get("missing_data_notes") or []
    if missing:
        lines.append("- missing data notes:")
        for note in missing[:12]:
            lines.append(f"  - {note}")
    interpretation = stft_summary.get("bounded_interpretation")
    if interpretation:
        lines.append(f"- bounded interpretation: {interpretation}")
    return lines


def _csv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.10g}"
    return value


def _float_or_none(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _format_optional(value: Any) -> str:
    if value is None:
        return "None"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "None"
        return f"{value:.6g}"
    return str(value)


if __name__ == "__main__":
    raise SystemExit(main())
