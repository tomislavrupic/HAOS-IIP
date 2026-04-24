#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RAW_GRADIENT = REPO_ROOT / "data" / "scalar_kernel_graph_recoverability_gradient_latest.json"

CSV_PATH = RESULTS_DIR / "raw_gradient_shape_audit.csv"
JSON_PATH = RESULTS_DIR / "raw_gradient_shape_audit.json"
SUMMARY_PATH = RESULTS_DIR / "raw_gradient_shape_audit_summary.md"

FIELDNAMES = [
    "readout_id",
    "scope",
    "source_artifact",
    "transform",
    "value",
    "target",
    "status",
    "scope_note",
]

DEFAULT_CONFIG: dict[str, Any] = {
    "raw_artifact": str(RAW_GRADIENT.relative_to(REPO_ROOT)),
    "thresholds": {
        "min_radial_alignment": 0.94,
        "max_power_deviation": 0.18,
        "min_power_fit_r2": 0.95,
        "max_flux_constancy_cv": 0.12,
        "max_refinement_profile_drift": 0.08,
    },
    "windowed_shape_readout": {
        "min_radius": 0.28,
        "smoothing_window": 5,
        "normalize_scaled_flux": False,
    },
    "shape_normalized_readout": {
        "min_radius": 0.28,
        "smoothing_window": 5,
        "normalize_scaled_flux": True,
    },
    "drift_samples": 16,
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing frozen scalar artifact: {path.relative_to(REPO_ROOT)}")
    return json.loads(path.read_text(encoding="utf-8"))


def artifact(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def fmt(value: float) -> str:
    return f"{value:.6f}"


def all_cases(data: dict[str, Any]) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for key in ("disorder_cases", "kernel_cases", "cases"):
        cases.extend(data.get(key, []))
    return cases


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def std(values: list[float]) -> float:
    if not values:
        return 0.0
    mu = mean(values)
    return math.sqrt(sum((value - mu) ** 2 for value in values) / len(values))


def fit_power(rows: list[dict[str, Any]]) -> dict[str, float]:
    pairs = [
        (float(row["radius"]), abs(float(row["radial_flux"])))
        for row in rows
        if float(row["radius"]) > 0.0 and abs(float(row["radial_flux"])) > 0.0
    ]
    log_x = [math.log(x) for x, _ in pairs]
    log_y = [math.log(y) for _, y in pairs]
    mean_x = mean(log_x)
    mean_y = mean(log_y)
    sxx = sum((x - mean_x) ** 2 for x in log_x)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_x, log_y))
    slope = sxy / max(sxx, 1.0e-18)
    intercept = mean_y - slope * mean_x
    predictions = [slope * x + intercept for x in log_x]
    ss_res = sum((y - pred) ** 2 for y, pred in zip(log_y, predictions))
    ss_tot = sum((y - mean_y) ** 2 for y in log_y)
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(1.0 - ss_res / max(ss_tot, 1.0e-18)),
    }


def smooth_rows(rows: list[dict[str, Any]], window: int) -> list[dict[str, Any]]:
    if window <= 1 or len(rows) < window:
        return rows
    half_window = window // 2
    smoothed: list[dict[str, Any]] = []
    for index in range(len(rows)):
        lo = max(0, index - half_window)
        hi = min(len(rows), index + half_window + 1)
        chunk = rows[lo:hi]
        smoothed.append(
            {
                "radius": mean([float(row["radius"]) for row in chunk]),
                "radial_flux": mean([float(row["radial_flux"]) for row in chunk]),
                "scaled_flux": mean([float(row["scaled_flux"]) for row in chunk]),
                "count": sum(int(row["count"]) for row in chunk),
            }
        )
    return smoothed


def normalize_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    scale = max(mean([abs(float(row["scaled_flux"])) for row in rows]), 1.0e-18)
    normalized: list[dict[str, Any]] = []
    for row in rows:
        normalized.append(
            {
                "radius": float(row["radius"]),
                "radial_flux": float(row["radial_flux"]) / scale,
                "scaled_flux": float(row["scaled_flux"]) / scale,
                "count": int(row["count"]),
            }
        )
    return normalized


def transform_rows(
    rows: list[dict[str, Any]],
    min_radius: float,
    smoothing_window: int,
    normalize_scaled_flux: bool,
) -> list[dict[str, Any]]:
    filtered = [row for row in rows if float(row["radius"]) >= float(min_radius)]
    filtered = smooth_rows(filtered, int(smoothing_window))
    if normalize_scaled_flux:
        filtered = normalize_rows(filtered)
    return filtered


def profile_metrics(rows: list[dict[str, Any]]) -> dict[str, float]:
    power = fit_power(rows)
    scaled_flux = [float(row["scaled_flux"]) for row in rows]
    return {
        "power_slope": float(power["slope"]),
        "power_fit_r2": float(power["r2"]),
        "flux_constancy_cv": float(std(scaled_flux) / max(mean([abs(value) for value in scaled_flux]), 1.0e-18)),
    }


def interpolate_scaled_flux(rows: list[dict[str, Any]], radius: float) -> float:
    if radius <= float(rows[0]["radius"]):
        return float(rows[0]["scaled_flux"])
    if radius >= float(rows[-1]["radius"]):
        return float(rows[-1]["scaled_flux"])
    for left, right in zip(rows, rows[1:]):
        left_radius = float(left["radius"])
        right_radius = float(right["radius"])
        if left_radius <= radius <= right_radius:
            weight = (radius - left_radius) / max(right_radius - left_radius, 1.0e-18)
            return float(left["scaled_flux"]) + weight * (
                float(right["scaled_flux"]) - float(left["scaled_flux"])
            )
    return float(rows[-1]["scaled_flux"])


def interpolation_drift(
    rows_a: list[dict[str, Any]],
    rows_b: list[dict[str, Any]],
    samples: int,
) -> float:
    lo = max(float(rows_a[0]["radius"]), float(rows_b[0]["radius"]))
    hi = min(float(rows_a[-1]["radius"]), float(rows_b[-1]["radius"]))
    if hi <= lo:
        return float("inf")
    grid = [lo + (hi - lo) * index / max(samples - 1, 1) for index in range(samples)]
    flux_a = [interpolate_scaled_flux(rows_a, radius) for radius in grid]
    flux_b = [interpolate_scaled_flux(rows_b, radius) for radius in grid]
    norm = max(mean([abs(value) for value in flux_a]), 1.0e-18)
    rms = math.sqrt(sum((left - right) ** 2 for left, right in zip(flux_a, flux_b)) / len(grid))
    return float(rms / norm)


def pair_labels(raw: dict[str, Any]) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    disorder_cases = raw.get("disorder_cases", [])
    kernel_cases = raw.get("kernel_cases", [])
    for jitter_fraction in sorted({float(case["jitter_fraction"]) for case in disorder_cases}):
        left = next(
            case["label"]
            for case in disorder_cases
            if case["label"].startswith("disorder|n=11") and float(case["jitter_fraction"]) == jitter_fraction
        )
        right = next(
            case["label"]
            for case in disorder_cases
            if case["label"].startswith("disorder|n=13") and float(case["jitter_fraction"]) == jitter_fraction
        )
        pairs.append((left, right))
    for family in sorted({str(case["kernel_family"]) for case in kernel_cases}):
        left = next(
            case["label"]
            for case in kernel_cases
            if case["label"].startswith("kernel|n=11") and str(case["kernel_family"]) == family
        )
        right = next(
            case["label"]
            for case in kernel_cases
            if case["label"].startswith("kernel|n=13") and str(case["kernel_family"]) == family
        )
        pairs.append((left, right))
    return pairs


def summarize_canonical(raw: dict[str, Any], thresholds: dict[str, float]) -> dict[str, Any]:
    cases = all_cases(raw)
    weakest = min(cases, key=lambda case: float(case["profile"]["power_fit_r2"]))
    min_power_r2 = min(float(case["profile"]["power_fit_r2"]) for case in cases)
    max_power_deviation = max(abs(float(case["profile"]["power_slope"]) + 2.0) for case in cases)
    max_flux_cv = max(float(case["profile"]["flux_constancy_cv"]) for case in cases)
    max_drift = float(raw["max_profile_drift"])
    passes = bool(
        min_power_r2 >= float(thresholds["min_power_fit_r2"])
        and max_power_deviation <= float(thresholds["max_power_deviation"])
        and max_flux_cv <= float(thresholds["max_flux_constancy_cv"])
        and max_drift <= float(thresholds["max_refinement_profile_drift"])
    )
    return {
        "readout": "canonical_raw_local_gradient",
        "scope": "full observable",
        "transform": "stored raw local-gradient profile from 53.3",
        "pass": passes,
        "case_count": len(cases),
        "min_power_fit_r2": min_power_r2,
        "max_power_deviation": max_power_deviation,
        "max_flux_constancy_cv": max_flux_cv,
        "max_refinement_profile_drift": max_drift,
        "weakest_case": {
            "label": str(weakest["label"]),
            "power_slope": float(weakest["profile"]["power_slope"]),
            "power_fit_r2": float(weakest["profile"]["power_fit_r2"]),
            "flux_constancy_cv": float(weakest["profile"]["flux_constancy_cv"]),
        },
        "scope_note": "Keeps the canonical 53.3 raw boundary intact: full raw local-gradient power scaling remains open.",
    }


def summarize_transformed(
    raw: dict[str, Any],
    thresholds: dict[str, float],
    readout_id: str,
    scope: str,
    transform: str,
    min_radius: float,
    smoothing_window: int,
    normalize_scaled_flux: bool,
    drift_samples: int,
    scope_note: str,
) -> dict[str, Any]:
    cases = all_cases(raw)
    rows_by_label: dict[str, list[dict[str, Any]]] = {}
    case_metrics: list[tuple[str, dict[str, float]]] = []
    for case in cases:
        rows = transform_rows(
            rows=case["profile"]["rows"],
            min_radius=min_radius,
            smoothing_window=smoothing_window,
            normalize_scaled_flux=normalize_scaled_flux,
        )
        if len(rows) < 5:
            raise ValueError(f"transformed profile too short for {case['label']}")
        rows_by_label[str(case["label"])] = rows
        case_metrics.append((str(case["label"]), profile_metrics(rows)))

    min_power_r2 = min(metric["power_fit_r2"] for _, metric in case_metrics)
    max_power_deviation = max(abs(metric["power_slope"] + 2.0) for _, metric in case_metrics)
    max_flux_cv = max(metric["flux_constancy_cv"] for _, metric in case_metrics)
    max_drift = max(
        interpolation_drift(rows_by_label[left], rows_by_label[right], samples=drift_samples)
        for left, right in pair_labels(raw)
    )
    weakest_label, weakest_metric = min(case_metrics, key=lambda item: item[1]["power_fit_r2"])
    passes = bool(
        min_power_r2 >= float(thresholds["min_power_fit_r2"])
        and max_power_deviation <= float(thresholds["max_power_deviation"])
        and max_flux_cv <= float(thresholds["max_flux_constancy_cv"])
        and max_drift <= float(thresholds["max_refinement_profile_drift"])
    )
    return {
        "readout": readout_id,
        "scope": scope,
        "transform": transform,
        "pass": passes,
        "case_count": len(cases),
        "min_power_fit_r2": min_power_r2,
        "max_power_deviation": max_power_deviation,
        "max_flux_constancy_cv": max_flux_cv,
        "max_refinement_profile_drift": max_drift,
        "weakest_case": {
            "label": weakest_label,
            "power_slope": float(weakest_metric["power_slope"]),
            "power_fit_r2": float(weakest_metric["power_fit_r2"]),
            "flux_constancy_cv": float(weakest_metric["flux_constancy_cv"]),
        },
        "scope_note": scope_note,
        "transform_config": {
            "min_radius": float(min_radius),
            "smoothing_window": int(smoothing_window),
            "normalize_scaled_flux": bool(normalize_scaled_flux),
            "drift_samples": int(drift_samples),
        },
    }


def row(
    readout_id: str,
    scope: str,
    source_artifact: str,
    transform: str,
    value: str,
    target: str,
    status: str,
    scope_note: str,
) -> dict[str, str]:
    return {
        "readout_id": readout_id,
        "scope": scope,
        "source_artifact": source_artifact,
        "transform": transform,
        "value": value,
        "target": target,
        "status": status,
        "scope_note": scope_note,
    }


def make_result(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))
    if config is not None:
        for key, value in config.items():
            if key not in ("thresholds", "windowed_shape_readout", "shape_normalized_readout"):
                cfg[key] = value
        for key in ("thresholds", "windowed_shape_readout", "shape_normalized_readout"):
            if isinstance(config.get(key), dict):
                cfg[key].update(config[key])

    raw = load_json(REPO_ROOT / str(cfg["raw_artifact"]))
    thresholds = cfg["thresholds"]
    canonical = summarize_canonical(raw, thresholds)
    windowed = summarize_transformed(
        raw=raw,
        thresholds=thresholds,
        readout_id="windowed_raw_local_gradient",
        scope="full observable",
        transform="source-core exclusion + local smoothing",
        min_radius=float(cfg["windowed_shape_readout"]["min_radius"]),
        smoothing_window=int(cfg["windowed_shape_readout"]["smoothing_window"]),
        normalize_scaled_flux=bool(cfg["windowed_shape_readout"]["normalize_scaled_flux"]),
        drift_samples=int(cfg["drift_samples"]),
        scope_note=(
            "Improves raw power-law shape metrics after source-core exclusion, "
            "but full-observable refinement drift remains open."
        ),
    )
    shape_normalized = summarize_transformed(
        raw=raw,
        thresholds=thresholds,
        readout_id="shape_normalized_windowed_raw_local_gradient",
        scope="shape-only observable",
        transform="source-core exclusion + local smoothing + per-case scaled-flux normalization",
        min_radius=float(cfg["shape_normalized_readout"]["min_radius"]),
        smoothing_window=int(cfg["shape_normalized_readout"]["smoothing_window"]),
        normalize_scaled_flux=bool(cfg["shape_normalized_readout"]["normalize_scaled_flux"]),
        drift_samples=int(cfg["drift_samples"]),
        scope_note=(
            "Closes a shape-only raw local-gradient proxy after source-core exclusion "
            "and per-case amplitude normalization; this does not replace the canonical raw bridge row."
        ),
    )

    rows = [
        row(
            readout_id=canonical["readout"],
            scope=canonical["scope"],
            source_artifact=str(cfg["raw_artifact"]),
            transform=canonical["transform"],
            value=(
                f"{fmt(canonical['min_power_fit_r2'])} / "
                f"{fmt(canonical['max_power_deviation'])} / "
                f"{fmt(canonical['max_flux_constancy_cv'])} / "
                f"{fmt(canonical['max_refinement_profile_drift'])}"
            ),
            target=(
                f">= {fmt(float(thresholds['min_power_fit_r2']))} / "
                f"<= {fmt(float(thresholds['max_power_deviation']))} / "
                f"<= {fmt(float(thresholds['max_flux_constancy_cv']))} / "
                f"<= {fmt(float(thresholds['max_refinement_profile_drift']))}"
            ),
            status="PASS" if canonical["pass"] else "OPEN",
            scope_note=canonical["scope_note"],
        ),
        row(
            readout_id=windowed["readout"],
            scope=windowed["scope"],
            source_artifact=str(cfg["raw_artifact"]),
            transform=windowed["transform"],
            value=(
                f"{fmt(windowed['min_power_fit_r2'])} / "
                f"{fmt(windowed['max_power_deviation'])} / "
                f"{fmt(windowed['max_flux_constancy_cv'])} / "
                f"{fmt(windowed['max_refinement_profile_drift'])}"
            ),
            target=(
                f">= {fmt(float(thresholds['min_power_fit_r2']))} / "
                f"<= {fmt(float(thresholds['max_power_deviation']))} / "
                f"<= {fmt(float(thresholds['max_flux_constancy_cv']))} / "
                f"<= {fmt(float(thresholds['max_refinement_profile_drift']))}"
            ),
            status="PASS" if windowed["pass"] else "OPEN",
            scope_note=windowed["scope_note"],
        ),
        row(
            readout_id=shape_normalized["readout"],
            scope=shape_normalized["scope"],
            source_artifact=str(cfg["raw_artifact"]),
            transform=shape_normalized["transform"],
            value=(
                f"{fmt(shape_normalized['min_power_fit_r2'])} / "
                f"{fmt(shape_normalized['max_power_deviation'])} / "
                f"{fmt(shape_normalized['max_flux_constancy_cv'])} / "
                f"{fmt(shape_normalized['max_refinement_profile_drift'])}"
            ),
            target=(
                f">= {fmt(float(thresholds['min_power_fit_r2']))} / "
                f"<= {fmt(float(thresholds['max_power_deviation']))} / "
                f"<= {fmt(float(thresholds['max_flux_constancy_cv']))} / "
                f"<= {fmt(float(thresholds['max_refinement_profile_drift']))}"
            ),
            status="PASS" if shape_normalized["pass"] else "OPEN",
            scope_note=shape_normalized["scope_note"],
        ),
    ]

    if shape_normalized["pass"]:
        observation = (
            "the raw local-gradient boundary splits into three readable layers: "
            "canonical raw remains open, source-core exclusion alone improves shape but not refinement drift, "
            "and a shape-normalized windowed raw read passes the existing raw thresholds as a shape-only proxy"
        )
        conclusion = (
            f"canonical raw local-gradient stays open with minimum R^2 {canonical['min_power_fit_r2']:.4f} "
            f"and drift {canonical['max_refinement_profile_drift']:.4f}; "
            f"source-core exclusion plus smoothing improves the raw fit to minimum R^2 {windowed['min_power_fit_r2']:.4f} "
            f"but leaves drift open at {windowed['max_refinement_profile_drift']:.4f}; "
            f"after the same windowing with per-case scaled-flux normalization, the shape-only read passes with minimum R^2 "
            f"{shape_normalized['min_power_fit_r2']:.4f}, maximum |slope+2| {shape_normalized['max_power_deviation']:.4f}, "
            f"maximum CV {shape_normalized['max_flux_constancy_cv']:.4f}, and drift {shape_normalized['max_refinement_profile_drift']:.4f}; "
            "this supports a bounded shape-only raw-gradient proxy, not full raw inverse-square closure"
        )
    else:
        observation = (
            "source-core exclusion and shape normalization improve the raw local-gradient read, "
            "but the refined raw audit remains outside at least one existing threshold"
        )
        conclusion = (
            f"the strongest audited raw-gradient variant remains open at minimum R^2 {shape_normalized['min_power_fit_r2']:.4f}, "
            f"maximum |slope+2| {shape_normalized['max_power_deviation']:.4f}, CV {shape_normalized['max_flux_constancy_cv']:.4f}, "
            f"and drift {shape_normalized['max_refinement_profile_drift']:.4f}"
        )

    return {
        "experiment": "raw_gradient_shape_audit",
        "config": cfg,
        "source_artifacts": {
            "raw_local_gradient": str(cfg["raw_artifact"]),
        },
        "summaries": {
            "canonical_raw_local_gradient": canonical,
            "windowed_raw_local_gradient": windowed,
            "shape_normalized_windowed_raw_local_gradient": shape_normalized,
        },
        "rows": rows,
        "observation": observation,
        "conclusion": conclusion,
    }


def write_csv(rows: list[dict[str, str]]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(result: dict[str, Any]) -> None:
    JSON_PATH.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


def write_summary(result: dict[str, Any]) -> None:
    rows = result["rows"]
    table_rows = "\n".join(
        "| {readout_id} | {scope} | {transform} | {value} | {target} | {status} |".format(**row)
        for row in rows
    )
    summary = f"""# HAOS-IIP Raw Gradient Shape Audit

## Scope

This file is generated by `experiments/physics_bridge/raw_gradient_shape_audit.py`.

It is an external post-processing audit over the frozen raw recoverability-gradient artifact. It does not modify HAOS core, change the canonical 53.0 bridge rows, or loosen any thresholds.

## Audit Table

| Readout | Scope | Transform | Value | Target | Status |
| --- | --- | --- | --- | --- | --- |
{table_rows}

## Observation

{result['observation']}

## Conclusion

{result['conclusion']}

## Boundary

- `canonical_raw_local_gradient` remains the authoritative raw bridge row and stays `OPEN`.
- `shape_normalized_windowed_raw_local_gradient` is a sidecar shape-only proxy and must not be promoted into full raw inverse-square closure.
"""
    SUMMARY_PATH.write_text(summary, encoding="utf-8")


def main() -> None:
    result = make_result()
    write_csv(result["rows"])
    write_json(result)
    write_summary(result)
    print("Raw gradient shape audit generated.")
    print(f"CSV: {CSV_PATH.relative_to(REPO_ROOT)}")
    print(f"JSON: {JSON_PATH.relative_to(REPO_ROOT)}")
    print(f"Summary: {SUMMARY_PATH.relative_to(REPO_ROOT)}")
    print(
        "Statement: This is a sidecar audit over frozen raw recoverability-gradient data; "
        "it introduces no HAOS core changes and no change to the canonical bridge rows."
    )


if __name__ == "__main__":
    main()
