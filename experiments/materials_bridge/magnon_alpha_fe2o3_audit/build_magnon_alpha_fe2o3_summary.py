from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"

RELEASE_DATE = "2026-05-07"

SUMMARY = {
    "bridge_line": "Materials Bridge Line B",
    "title": "Magnon Address Stability in alpha-Fe2O3",
    "status": "PASS_LITERATURE_USER_AUDIT_SEEDED",
    "data_status": "NO_RAW_PUBLIC_DATA_LOADED",
    "evidence_class": "literature_report_plus_user_provided_audit_anchor",
    "raw_stft_time_frequency_grid_status": "NO_PUBLIC_STFT_GRID_FOUND",
    "local_reproduction_status": "PENDING_RAW_DATA_ARTIFACT",
    "release_date": RELEASE_DATE,
    "primary_material": "alpha-Fe2O3",
    "carrier": "antiferromagnetic / altermagnetic spin order",
    "recoverable_structure": "coherent sub-THz exchange magnon packet",
    "observable_readout": "time-resolved magneto-optical magnon traces, spectra, and propagation delay",
    "perturbation_axis": "crystal orientation, wavevector direction, thickness/boundary condition, magnetic field, pump delay",
    "visible_failure": "broadband thermalized response, amplitude collapse, loss of directional propagation, or no recoverable narrowband magnon feature",
    "persistence_proxy": 0.902,
    "persistence_proxy_status": "USER_PROVIDED_PENDING_LOCAL_REPRODUCTION",
    "group_velocity_min_km_s": 17.31,
    "group_velocity_max_km_s": 24.43,
    "group_velocity_span_km_s": 7.12,
    "propagation_status": "micrometer_scale_directional_literature_reported",
    "spectral_status": "narrowband_sub_thz_exchange_magnons_literature_reported",
    "k_star_status": "NOT_CLAIMED_WITHOUT_RAW_SWEEP",
    "first_visible_failure_status": "NOT_CLAIMED_WITHOUT_RAW_SWEEP",
    "prediction_records": [
        "HAOS-PRED-010",
        "HAOS-PRED-011",
        "HAOS-PRED-012",
        "HAOS-PRED-013",
    ],
    "public_sources": [
        {
            "label": "Li et al., Anisotropic Coherent Propagation of Sub-Terahertz Magnons in Altermagnetic alpha-Fe2O3",
            "doi": "10.1002/adom.202503604",
            "url": "https://doi.org/10.1002/adom.202503604",
            "role": "primary literature anchor for group velocity and anisotropic coherent propagation",
        },
        {
            "label": "Yang et al., Spin-orbit torque manipulation of sub-terahertz magnons in antiferromagnetic alpha-Fe2O3",
            "doi": "10.1038/s41467-024-48431-w",
            "url": "https://doi.org/10.1038/s41467-024-48431-w",
            "role": "secondary alpha-Fe2O3 magnon control context",
        },
    ],
    "bounded_interpretation": (
        "Line B records a first-class magnon bridge artifact for alpha-Fe2O3, "
        "but its current metrics are literature/user-audit seeded. It supports a "
        "bounded spin-sector persistence note and prediction package, while leaving "
        "raw-data recoverability validation pending until open time-domain or "
        "time-frequency tables are committed and parsed locally."
    ),
    "non_claims": [
        "not a proof of HAOS-IIP",
        "not a new magnon theory",
        "not a raw-data external bridge pass",
        "not a claim that alpha-Fe2O3 embodies HAOS-IIP",
        "no k_star claim without raw sweep data",
        "no raw STFT/time-frequency claim without a public validated grid",
    ],
}

METRICS = [
    {
        "metric": "persistence_proxy",
        "value": 0.902,
        "unit": "unitless",
        "evidence_kind": "user_provided_audit_anchor",
        "source": "local conversation record, pending committed raw-data reproduction",
        "included_in_line_b_summary": True,
        "boundary_note": "Treat as seeded anchor until raw alpha-Fe2O3 audit tables are committed.",
    },
    {
        "metric": "group_velocity_min",
        "value": 17.31,
        "unit": "km/s",
        "evidence_kind": "literature_report",
        "source": "Li et al. 2026, DOI 10.1002/adom.202503604",
        "included_in_line_b_summary": True,
        "boundary_note": "Reported velocity range; no local trace refit is claimed here.",
    },
    {
        "metric": "group_velocity_max",
        "value": 24.43,
        "unit": "km/s",
        "evidence_kind": "literature_report",
        "source": "Li et al. 2026, DOI 10.1002/adom.202503604",
        "included_in_line_b_summary": True,
        "boundary_note": "Reported velocity range; no local trace refit is claimed here.",
    },
    {
        "metric": "raw_stft_time_frequency_grid_found",
        "value": 0.0,
        "unit": "boolean_0_or_1",
        "evidence_kind": "boundary_status",
        "source": "current repo audit state",
        "included_in_line_b_summary": True,
        "boundary_note": "NO_PUBLIC_STFT_GRID_FOUND; no raw STFT recoverability claim is made.",
    },
    {
        "metric": "local_raw_data_loaded",
        "value": 0.0,
        "unit": "boolean_0_or_1",
        "evidence_kind": "boundary_status",
        "source": "current repo audit state",
        "included_in_line_b_summary": True,
        "boundary_note": "No alpha-Fe2O3 raw data files are committed in this sidecar yet.",
    },
]


def write_json() -> Path:
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    path = OUTPUTS / "magnon_persistence_summary.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(SUMMARY, handle, indent=2)
        handle.write("\n")
    return path


def write_csv() -> Path:
    path = OUTPUTS / "magnon_audit_metrics.csv"
    fieldnames = [
        "metric",
        "value",
        "unit",
        "evidence_kind",
        "source",
        "included_in_line_b_summary",
        "boundary_note",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in METRICS:
            writer.writerow(row)
    return path


def write_markdown() -> Path:
    path = OUTPUTS / "magnon_address_stability_summary.md"
    text = f"""# Materials Bridge Line B - Magnon Address Stability in alpha-Fe2O3

## Status

`{SUMMARY["status"]}`

## Evidence Class

`{SUMMARY["evidence_class"]}`

This sidecar closes the repo bookkeeping loop for the alpha-Fe2O3 magnon
prediction layer. It is a first-class artifact, but it is not yet a locally
reproduced raw-data bridge. The persistence proxy is recorded as a user-provided
audit anchor until public raw tables or a dedicated local audit dataset are
committed and parsed.

## Metrics

- persistence_proxy: `{SUMMARY["persistence_proxy"]}`
- group_velocity_range_km_s: `{SUMMARY["group_velocity_min_km_s"]}-{SUMMARY["group_velocity_max_km_s"]}`
- raw_stft_time_frequency_grid_status: `{SUMMARY["raw_stft_time_frequency_grid_status"]}`
- data_status: `{SUMMARY["data_status"]}`
- k_star_status: `{SUMMARY["k_star_status"]}`

## HAOS Mapping

- carrier: {SUMMARY["carrier"]}
- perturbation: {SUMMARY["perturbation_axis"]}
- recoverable structure: {SUMMARY["recoverable_structure"]}
- observable readout: {SUMMARY["observable_readout"]}
- visible failure: {SUMMARY["visible_failure"]}

## Bounded Interpretation

{SUMMARY["bounded_interpretation"]}

## Non-Claims

""" + "\n".join(f"- {item}" for item in SUMMARY["non_claims"]) + "\n"
    path.write_text(text, encoding="utf-8")
    return path


def write_plot() -> Path:
    import matplotlib.pyplot as plt

    path = OUTPUTS / "magnon_persistence_components.png"

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.suptitle("Materials Bridge Line B: alpha-Fe2O3 Magnon Address Stability", fontsize=12)

    axes[0].bar(["persistence"], [SUMMARY["persistence_proxy"]], color="#3566a8")
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("proxy")
    axes[0].set_title("Seeded Persistence Proxy")
    axes[0].text(0, SUMMARY["persistence_proxy"] + 0.03, "0.902", ha="center", fontsize=10)

    axes[1].barh(["reported range"], [SUMMARY["group_velocity_max_km_s"] - SUMMARY["group_velocity_min_km_s"]],
                 left=[SUMMARY["group_velocity_min_km_s"]], color="#3a8f71")
    axes[1].set_xlim(0, 30)
    axes[1].set_xlabel("km/s")
    axes[1].set_title("Reported Group Velocity")
    axes[1].text(20.9, 0, "17.31-24.43 km/s", ha="center", va="center", color="white", fontsize=9)

    boundary_labels = ["raw data", "raw STFT grid", "k_star sweep"]
    boundary_values = [0, 0, 0]
    axes[2].bar(boundary_labels, boundary_values, color="#8a8f98")
    axes[2].set_ylim(0, 1)
    axes[2].set_title("Current Boundaries")
    axes[2].tick_params(axis="x", rotation=25)
    for idx, label in enumerate(["pending", "not found", "not claimed"]):
        axes[2].text(idx, 0.08, label, ha="center", fontsize=9)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main() -> None:
    written = [write_json(), write_csv(), write_markdown(), write_plot()]
    print("magnon_line_b_status:", SUMMARY["status"])
    print("data_status:", SUMMARY["data_status"])
    print("persistence_proxy:", SUMMARY["persistence_proxy"])
    print("group_velocity_range_km_s:", f"{SUMMARY['group_velocity_min_km_s']}-{SUMMARY['group_velocity_max_km_s']}")
    print("raw_stft_time_frequency_grid_status:", SUMMARY["raw_stft_time_frequency_grid_status"])
    print("outputs_written:")
    for path in written:
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
