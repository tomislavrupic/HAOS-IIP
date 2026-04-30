#!/usr/bin/env python3
"""Phase 60 claim-gated physics bridge update.

This generator consolidates the celestial-facing sidecars without promoting toy
or known-target PASS results into established physics claims.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(__file__).resolve().parent

SPHERICAL_STATUS = REPO_ROOT / "experiments/physics_bridge/spherical_harmonic_control_probe/bridge_status.json"
SOFT_STATUS = REPO_ROOT / "experiments/physics_bridge/soft_structure_proxy_test/bridge_status.json"
CELESTIAL_STATUS = REPO_ROOT / "experiments/physics_bridge/celestial_boundary_audit/bridge_status.json"
MEMORY_TOY_STATUS = REPO_ROOT / "experiments/physics_bridge/gravitational_memory_toy_probe/bridge_status.json"
MULTIPOLE_STATUS = REPO_ROOT / "experiments/physics_bridge/multipole_supertranslation_probe/bridge_status.json"
COMPOSITION_STATUS = REPO_ROOT / "experiments/physics_bridge/supertranslation_memory_composition_probe/bridge_status.json"
GW_ENTRY_STATUS = REPO_ROOT / "experiments/physics_bridge/gw_memory_entry_gate/bridge_status.json"
GW_HARDENING_STATUS = REPO_ROOT / "experiments/physics_bridge/gw_memory_entry_gate/event_window_hardening_status.json"

CSV_PATH = ROOT / "claim_gate_table.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "claim_gated_physics_bridge.md"

FIELDNAMES = [
    "claim_id",
    "claim_or_probe",
    "status",
    "scope",
    "supporting_artifact",
    "promotion_gate",
    "allowed_language",
    "disallowed_language",
    "notes",
]


@dataclass(frozen=True)
class ClaimGateRow:
    claim_id: str
    claim_or_probe: str
    status: str
    scope: str
    supporting_artifact: str
    promotion_gate: str
    allowed_language: str
    disallowed_language: str
    notes: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 60 claim-gated bridge update.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated CSV, status JSON, and markdown report.",
    )
    return parser.parse_args()


def load_status(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def artifact(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def build_rows() -> list[ClaimGateRow]:
    spherical = load_status(SPHERICAL_STATUS)
    soft = load_status(SOFT_STATUS)
    celestial = load_status(CELESTIAL_STATUS)
    memory_toy = load_status(MEMORY_TOY_STATUS)
    multipole = load_status(MULTIPOLE_STATUS)
    composition = load_status(COMPOSITION_STATUS)
    gw_entry = load_status(GW_ENTRY_STATUS)
    gw_hardening = load_status(GW_HARDENING_STATUS)
    return [
        ClaimGateRow(
            claim_id="spherical_harmonic_boundary_geometry",
            claim_or_probe="Spherical harmonic boundary geometry",
            status=str(spherical["bridge_status"]),
            scope=(
                "Known-target S2 mode recovery on a constructed discrete sphere; hardened "
                "against weight-shuffle, geodesic-bin shuffle, degree-rewire, and smooth-ring controls."
            ),
            supporting_artifact=artifact(SPHERICAL_STATUS),
            promotion_gate="Already closed only as a boundary-geometry sanity check.",
            allowed_language="known-target S2 mode sanity check; spherical-harmonic control probe",
            disallowed_language="null infinity recovered; celestial sphere derived; asymptotic boundary proven",
            notes=(
                f"latest_real_score={spherical['latest_real_score']:.6f}; "
                f"control_max={spherical['control_max_score']:.6f}; "
                f"best_control={spherical['control_max_sample']}"
            ),
        ),
        ClaimGateRow(
            claim_id="toy_soft_structure_proxy",
            claim_or_probe="Toy soft-structure proxy",
            status=str(soft["bridge_status"]),
            scope=(
                "Synthetic amplitude-like pole/residue/factorization benchmark with controls "
                "that preserve smoothness or spectra while breaking soft-specific structure."
            ),
            supporting_artifact=artifact(SOFT_STATUS),
            promotion_gate="Already closed only for toy data with known residues and factorization.",
            allowed_language="toy soft-pole/residue/factorization proxy; synthetic soft-structure sanity check",
            disallowed_language="real soft theorem recovery; celestial amplitude recovery; S-matrix evidence",
            notes=(
                f"target_soft_score={soft['target_soft_specificity_score']:.6f}; "
                f"control_max={soft['control_max_score']:.6f}; "
                f"best_control={soft['control_max_sample']}"
            ),
        ),
        ClaimGateRow(
            claim_id="toy_gravitational_memory_proxy",
            claim_or_probe="Toy gravitational-memory proxy",
            status=str(memory_toy["bridge_status"]),
            scope=(
                "Synthetic permanent-displacement benchmark on sampled S2 with known target "
                "mode, permanent step, temporal profile, and controls."
            ),
            supporting_artifact=artifact(MEMORY_TOY_STATUS),
            promotion_gate="Already closed only as a synthetic displacement-memory toy benchmark.",
            allowed_language="toy memory-like displacement proxy; synthetic memory benchmark",
            disallowed_language="real gravitational memory detected; memory observable recovered; BMS soft hair recovered",
            notes=(
                f"target_memory_score={memory_toy['target_memory_specificity_score']:.6f}; "
                f"control_max={memory_toy['control_max_score']:.6f}; "
                f"best_control={memory_toy['control_max_sample']}"
            ),
        ),
        ClaimGateRow(
            claim_id="toy_multipole_supertranslation_proxy",
            claim_or_probe="Toy multi-pole supertranslation proxy",
            status=str(multipole["bridge_status"]),
            scope=(
                "Synthetic ordered l=2,3,4 permanent mode-shift sequence on sampled S2; "
                "tests ordered address-binding, field identity, event order, and permanence."
            ),
            supporting_artifact=artifact(MULTIPOLE_STATUS),
            promotion_gate="Already closed only as an ordered synthetic S2 mode-shift benchmark.",
            allowed_language="toy multi-pole supertranslation-like proxy; ordered S2 mode-shift benchmark",
            disallowed_language="real BMS supertranslation recovered; BMS charge recovered; soft hair proven",
            notes=(
                f"target_supertranslation_score={multipole['target_supertranslation_specificity_score']:.6f}; "
                f"control_max={multipole['control_max_score']:.6f}; "
                f"best_control={multipole['control_max_sample']}"
            ),
        ),
        ClaimGateRow(
            claim_id="toy_supertranslation_memory_composition_proxy",
            claim_or_probe="Toy supertranslation + memory composition proxy",
            status=str(composition["bridge_status"]),
            scope=(
                "Synthetic two-stage benchmark that composes a supertranslation-like shift "
                "with a later induced memory-like response through a fixed toy response map."
            ),
            supporting_artifact=artifact(COMPOSITION_STATUS),
            promotion_gate="Already closed only for the synthetic response-map composition benchmark.",
            allowed_language="toy supertranslation-memory composition proxy; synthetic response relation",
            disallowed_language="real BMS-memory composition recovered; gravitational memory relation proven",
            notes=(
                f"target_composition_score={composition['target_composition_specificity_score']:.6f}; "
                f"control_max={composition['control_max_score']:.6f}; "
                f"best_control={composition['control_max_sample']}"
            ),
        ),
        ClaimGateRow(
            claim_id="gw_memory_entry_gate_proxy",
            claim_or_probe="GW memory entry-gate proxy",
            status=str(gw_entry["bridge_status"]),
            scope=(
                "Strain-derived proxy telemetry on a deterministic GW150914-like surrogate "
                "or optional local strain file with time-shuffle, phase-scramble, amplitude-preserving, "
                "event-window, off-event, and noise controls."
            ),
            supporting_artifact=artifact(GW_ENTRY_STATUS),
            promotion_gate=(
                "Closed only as a strain-derived entry gate. Promotion requires external "
                "gravitational-memory observable benchmarks and real-data replication."
            ),
            allowed_language="strain-derived proxy entry gate; structured strain-like event separated from controls",
            disallowed_language="real gravitational memory detected; BMS charge recovered; GW memory observed",
            notes=(
                f"target_strain_proxy_score={gw_entry['target_strain_proxy_score']:.6f}; "
                f"control_max={gw_entry['control_max_score']:.6f}; "
                f"best_control={gw_entry['control_max_sample']}"
            ),
        ),
        ClaimGateRow(
            claim_id="gw_event_window_hardening_proxy",
            claim_or_probe="GW event-window hardening proxy",
            status=str(gw_hardening["bridge_status"]),
            scope=(
                "Multi-event strain-derived proxy hardening against sliding-window, partial-overlap, "
                "chirp-reversal, envelope-locked phase, timing-randomization, and micro-chunk controls."
            ),
            supporting_artifact=artifact(GW_HARDENING_STATUS),
            promotion_gate=(
                "Requires mean margin >= 0.25 and zero event-window control leakage before being "
                "promoted from MARGINAL to bounded PASS."
            ),
            allowed_language="event-window hardening diagnostic; remaining strain-proxy leakage characterized",
            disallowed_language="real gravitational memory detected; event-window leak closed; GW memory observed",
            notes=(
                f"mean_margin={gw_hardening['mean_margin_over_best_control']:.6f}; "
                f"min_margin={gw_hardening['min_margin_over_best_control']:.6f}; "
                f"highest_control={gw_hardening['highest_control_sample']}"
            ),
        ),
        ClaimGateRow(
            claim_id="celestial_holography_recovery",
            claim_or_probe="Celestial holography recovery",
            status="OPEN",
            scope="Full flat-space holographic dictionary is not addressed by current HAOS sidecars.",
            supporting_artifact=artifact(CELESTIAL_STATUS),
            promotion_gate="Requires a real map from 4D scattering data to celestial correlator-like observables.",
            allowed_language="celestial boundary remains open; sidecar stress tests only",
            disallowed_language="HAOS recovers celestial holography; HAOS derives CCFT",
            notes=f"Phase 57 status remains {celestial['bridge_status']}.",
        ),
        ClaimGateRow(
            claim_id="bms_charge_recovery",
            claim_or_probe="BMS charge recovery",
            status="OPEN",
            scope="No asymptotic BMS charge algebra, Ward identity, or charge-conservation test has been closed.",
            supporting_artifact=artifact(CELESTIAL_STATUS),
            promotion_gate="Requires explicit supertranslation/superrotation charge construction and Ward tests.",
            allowed_language="BMS-facing requirement remains open",
            disallowed_language="BMS charges recovered; asymptotic symmetry proven",
            notes="Do not infer BMS structure from low-mode persistence or spherical harmonics.",
        ),
        ClaimGateRow(
            claim_id="virasoro_ccft_recovery",
            claim_or_probe="Virasoro / celestial CFT recovery",
            status="OPEN",
            scope="No Virasoro generators, stress-tensor action, CCFT OPE, or conformal primary dictionary is closed.",
            supporting_artifact=artifact(CELESTIAL_STATUS),
            promotion_gate="Requires controlled boundary operator algebra and conformal/OPE tests.",
            allowed_language="Virasoro/CCFT remains a comparison target",
            disallowed_language="Virasoro recovered; celestial CFT derived",
            notes="Spherical-harmonic mode recovery is not conformal field theory recovery.",
        ),
        ClaimGateRow(
            claim_id="s_matrix_reconstruction",
            claim_or_probe="S-matrix reconstruction",
            status="OPEN",
            scope="No real unitary scattering amplitudes or S-matrix dictionary are reconstructed.",
            supporting_artifact=artifact(CELESTIAL_STATUS),
            promotion_gate="Requires toy-unitary and real-amplitude tests before any S-matrix language is allowed.",
            allowed_language="no S-matrix reconstruction claim",
            disallowed_language="HAOS reconstructs amplitudes; HAOS contains the S-matrix",
            notes="Phase 59 uses synthetic amplitude-like tables only.",
        ),
        ClaimGateRow(
            claim_id="real_soft_theorem_recovery",
            claim_or_probe="Real soft theorem recovery",
            status="OPEN",
            scope="No gravitational, gauge-theory, or celestial soft theorem has been recovered from real amplitude data.",
            supporting_artifact=artifact(SOFT_STATUS),
            promotion_gate="Requires known physical soft theorem datasets with residues, Ward identities, and controls.",
            allowed_language="toy soft-structure proxy passed",
            disallowed_language="soft theorem recovered; soft graviton theorem closed",
            notes="The Phase 59 PASS is synthetic and row-local.",
        ),
        ClaimGateRow(
            claim_id="gravitational_memory_observable",
            claim_or_probe="Gravitational memory observable",
            status="OPEN",
            scope="No displacement memory, spin memory, detector observable, or memory/soft-charge equivalence test is closed.",
            supporting_artifact=artifact(GW_ENTRY_STATUS),
            promotion_gate="Requires an explicit memory observable benchmark and comparison to controls.",
            allowed_language="toy memory and strain-derived proxy rows exist; real observable remains open",
            disallowed_language="gravitational memory predicted; memory effect recovered",
            notes="Phases 61, 64, and 65 are proxy gates only; they do not close this row.",
        ),
        ClaimGateRow(
            claim_id="claim_boundary_enforcement",
            claim_or_probe="Claim boundary enforcement",
            status="PASS",
            scope="This Phase 60 table separates bounded proxy PASS rows from OPEN established-physics claims.",
            supporting_artifact="experiments/physics_bridge/claim_gated_bridge_update/claim_gate_table.csv",
            promotion_gate="Use this table when editing physics-facing README, reports, papers, or public summaries.",
            allowed_language="claim-gated physics bridge; bounded proxy rows",
            disallowed_language="source code of physics; theory of everything; real gravity closure",
            notes="A PASS here means claim hygiene, not physics closure.",
        ),
    ]


def build_status(rows: list[ClaimGateRow]) -> dict[str, Any]:
    counts = {
        "PASS": sum(1 for row in rows if row.status == "PASS"),
        "OPEN": sum(1 for row in rows if row.status == "OPEN"),
        "WATCH": sum(1 for row in rows if row.status == "WATCH"),
        "MARGINAL": sum(1 for row in rows if row.status == "MARGINAL"),
    }
    open_established_claims = [
        row.claim_id
        for row in rows
        if row.status == "OPEN" and row.claim_id not in {"claim_boundary_enforcement"}
    ]
    return {
        "phase": "60",
        "bridge_name": "claim_gated_physics_bridge_update",
        "bridge_status": "PASS",
        "failure_mode": "NONE",
        "success_definition": (
            "PASS means bounded sidecar claims are separated from established physics claims. "
            "It does not mean celestial holography, BMS, Virasoro, S-matrix, real soft theorem, "
            "or gravitational memory recovery has been closed."
        ),
        "status_counts": counts,
        "open_established_physics_claims": open_established_claims,
        "bounded_proxy_passes": [
            row.claim_id
            for row in rows
            if row.status == "PASS" and row.claim_id != "claim_boundary_enforcement"
        ],
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "66_event_window_specificity_or_real_gwosc_replicates",
        "non_claims": [
            "no celestial holography recovery claim",
            "no BMS charge recovery claim",
            "no Virasoro or celestial CFT recovery claim",
            "no S-matrix reconstruction claim",
            "no real soft theorem recovery claim",
            "no gravitational memory observable claim",
        ],
    }


def write_csv(path: Path, rows: list[ClaimGateRow]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def write_status(path: Path, status: dict[str, Any]) -> None:
    path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")


def write_report(path: Path, rows: list[ClaimGateRow], status: dict[str, Any]) -> None:
    table_rows = "\n".join(
        "| {claim_or_probe} | {status} | {scope} | {allowed_language} | {disallowed_language} |".format(
            **asdict(row)
        )
        for row in rows
    )
    open_rows = "\n".join(
        f"- `{row.claim_id}`: {row.promotion_gate}" for row in rows if row.status == "OPEN"
    )
    pass_rows = "\n".join(
        f"- `{row.claim_id}`: {row.notes}" for row in rows if row.status == "PASS"
    )
    report = f"""# Phase 60 Claim-Gated Physics Bridge Update

## Scope

This file is generated by `experiments/physics_bridge/claim_gated_bridge_update/run_claim_gated_bridge_update.py`.

Phase 60 consolidates the celestial-facing bridge after Phase 57 (`OPEN` boundary audit), Phase 58.1 (`PASS` known-target S2 mode probe), Phase 59 (`PASS` toy soft-structure proxy), Phase 61 (`PASS` toy memory proxy), Phase 62 (`PASS` toy multi-pole supertranslation proxy), Phase 63 (`PASS` toy composition proxy), Phase 64 (`PASS` strain-derived entry gate on the default surrogate), and Phase 65 (`MARGINAL` event-window hardening).

The purpose is claim control. It does not modify HAOS core and does not promote toy or known-target sidecar results into established physics claims.

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- `PASS`: {status['status_counts']['PASS']}
- `OPEN`: {status['status_counts']['OPEN']}
- `WATCH`: {status['status_counts']['WATCH']}
- `MARGINAL`: {status['status_counts']['MARGINAL']}
- core modified: `{str(status['core_modified']).lower()}`

## Claim Gate Table

| Claim / probe | Status | Scope | Allowed language | Disallowed language |
| --- | --- | --- | --- | --- |
{table_rows}

## Bounded PASS Rows

{pass_rows}

## Established Physics Claims Still OPEN

{open_rows}

## Use Rule

Use the PASS rows as bounded instrumentation results only:

- `spherical_harmonic_boundary_geometry`: known-target S2 boundary-geometry sanity check.
- `toy_soft_structure_proxy`: synthetic pole/residue/factorization sanity check.
- `toy_gravitational_memory_proxy`: synthetic permanent-displacement sanity check.
- `toy_multipole_supertranslation_proxy`: synthetic ordered S2 mode-shift benchmark.
- `toy_supertranslation_memory_composition_proxy`: synthetic response-map composition benchmark.
- `gw_memory_entry_gate_proxy`: strain-derived entry gate on a deterministic surrogate or optional local strain file.

The Phase 65 event-window hardening row is deliberately `MARGINAL`: the target remains strong, but envelope-locked and sliding/micro-window controls still compete under the stricter ladder.

Do not use any of these rows as evidence that HAOS-IIP has recovered celestial holography, BMS charges, Virasoro/CCFT structure, a real S-matrix, real soft theorems, or gravitational memory.

## Next

Phase 66 should either tighten event-window specificity directly or replicate the Phase 64/65 entry gate on curated local GWOSC strain files. In both cases, real gravitational-memory observable and BMS/soft-theorem rows remain `OPEN` until directly tested.

## Authority

- CSV: `experiments/physics_bridge/claim_gated_bridge_update/claim_gate_table.csv`
- JSON: `experiments/physics_bridge/claim_gated_bridge_update/bridge_status.json`
- Generator: `experiments/physics_bridge/claim_gated_bridge_update/run_claim_gated_bridge_update.py`
"""
    path.write_text(report, encoding="utf-8")


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    status = build_status(rows)
    write_csv(output_dir / CSV_PATH.name, rows)
    write_status(output_dir / STATUS_PATH.name, status)
    write_report(output_dir / REPORT_PATH.name, rows, status)

    print("Phase 60 claim-gated bridge update generated.")
    print(f"CSV: {relative(output_dir / CSV_PATH.name)}")
    print(f"JSON: {relative(output_dir / STATUS_PATH.name)}")
    print(f"Report: {relative(output_dir / REPORT_PATH.name)}")
    print(
        f"bridge_status: {status['bridge_status']} with "
        f"{len(status['open_established_physics_claims'])} established physics claims still OPEN"
    )


if __name__ == "__main__":
    main()
