#!/usr/bin/env python3
"""Phase 57 celestial-holography boundary audit for HAOS-IIP.

This sidecar is deliberately conservative. It does not try to simulate
celestial holography. It writes a reproducible requirements table stating which
celestial-facing structures HAOS-IIP currently has not closed.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
AUDIT_DIR = Path(__file__).resolve().parent

CSV_NAME = "haos_vs_celestial_requirements.csv"
STATUS_NAME = "bridge_status.json"
REPORT_NAME = "celestial_boundary_audit.md"

FIELDNAMES = [
    "audit_id",
    "celestial_requirement",
    "haos_current_evidence",
    "current_coverage",
    "status",
    "boundary_risk",
    "next_test",
    "claim_gate",
]


@dataclass(frozen=True)
class AuditRow:
    audit_id: str
    celestial_requirement: str
    haos_current_evidence: str
    current_coverage: str
    status: str
    boundary_risk: str
    next_test: str
    claim_gate: str


ROWS: tuple[AuditRow, ...] = (
    AuditRow(
        audit_id="asymptotic_s2_boundary",
        celestial_requirement="Recover an S2-like asymptotic boundary at null infinity.",
        haos_current_evidence=(
            "Phases 17-18 build proto-geometry and distance surrogates from frozen "
            "branch-local operator data."
        ),
        current_coverage="local proto-geometric proxies only",
        status="OPEN",
        boundary_risk="high",
        next_test="Phase 58 spherical-harmonic control probe on a known discrete sphere.",
        claim_gate="Do not describe HAOS proto-geometry as null infinity or a celestial sphere.",
    ),
    AuditRow(
        audit_id="conformal_mellin_basis",
        celestial_requirement="Represent conformal weights via Mellin-basis boost eigenstates.",
        haos_current_evidence=(
            "Phase VII and later bridges track graph-Laplacian spectra, low-mode "
            "organization, and power-law-like scaling."
        ),
        current_coverage="spectral scaling proxies only",
        status="OPEN",
        boundary_risk="high",
        next_test="Construct a known Mellin/conformal-weight toy table before any HAOS comparison.",
        claim_gate="Power-law spectral fits are not celestial conformal dimensions.",
    ),
    AuditRow(
        audit_id="bms_supertranslation_charges",
        celestial_requirement="Recover BMS supertranslation charges and their Ward identities.",
        haos_current_evidence=(
            "Current HAOS telemetry measures persistence, recoverability, causal closure, "
            "and scalar-carrier current proxies."
        ),
        current_coverage="no asymptotic charge algebra",
        status="OPEN",
        boundary_risk="high",
        next_test="Define an explicit charge-conservation benchmark with known Ward failures.",
        claim_gate="Do not infer BMS charges from low-mode persistence or recoverability.",
    ),
    AuditRow(
        audit_id="virasoro_superrotation_structure",
        celestial_requirement="Recover superrotation/Virasoro-like generators and stress-tensor action.",
        haos_current_evidence=(
            "HAOS has harmonic addresses and branch-local spectral envelopes, but no "
            "stress tensor or conformal generator construction."
        ),
        current_coverage="no Virasoro representation",
        status="OPEN",
        boundary_risk="high",
        next_test="Only test after a controlled boundary operator algebra exists.",
        claim_gate="Do not use Virasoro/CFT language except as an external comparison target.",
    ),
    AuditRow(
        audit_id="soft_theorem_factorization",
        celestial_requirement="Preserve known soft limits, residues, and factorization behavior.",
        haos_current_evidence=(
            "Physics bridge rows check scalar response laws, current closure, disorder flux, "
            "and localized-bump behavior."
        ),
        current_coverage="generic perturbation telemetry only",
        status="OPEN",
        boundary_risk="high",
        next_test="Phase 59 soft-theorem proxy test with known pole and residue structure.",
        claim_gate="Recoverability under perturbation is not a soft theorem.",
    ),
    AuditRow(
        audit_id="s_matrix_unitarity",
        celestial_requirement="Connect to unitary flat-space scattering data or an S-matrix dictionary.",
        haos_current_evidence=(
            "HAOS propagation and temporal-ordering phases are abstract operator diagnostics, "
            "not scattering amplitudes."
        ),
        current_coverage="no amplitude or unitarity map",
        status="OPEN",
        boundary_risk="high",
        next_test="Start with a toy unitary scattering table before real amplitude claims.",
        claim_gate="Do not call current HAOS trajectories S-matrix data.",
    ),
    AuditRow(
        audit_id="gravitational_memory_observable",
        celestial_requirement="Map to gravitational memory observables or memory/soft-charge equivalence.",
        haos_current_evidence=(
            "HAOS tracks persistence and recoverability after perturbations, including "
            "biology and physics-bridge stress tests."
        ),
        current_coverage="memory-like metaphor only",
        status="OPEN",
        boundary_risk="high",
        next_test="Require an explicit memory observable benchmark before promotion.",
        claim_gate="Persistence is not gravitational memory.",
    ),
    AuditRow(
        audit_id="celestial_ope_collinear_limits",
        celestial_requirement="Recover celestial OPE coefficients or collinear-limit constraints.",
        haos_current_evidence=(
            "No current HAOS artifact encodes collinear amplitude limits or celestial OPE data."
        ),
        current_coverage="absent",
        status="OPEN",
        boundary_risk="high",
        next_test="Use a toy amplitude dataset with known collinear limits.",
        claim_gate="Do not claim OPE recovery from branch-local hierarchy closure.",
    ),
    AuditRow(
        audit_id="ir_pole_residue_structure",
        celestial_requirement="Track infrared poles, analytic continuation, and residue stability.",
        haos_current_evidence=(
            "Current scalar artifacts are finite-grid numerical observables with explicit "
            "frozen manifests and telemetry."
        ),
        current_coverage="finite-grid telemetry only",
        status="OPEN",
        boundary_risk="high",
        next_test="Create a residue-recoverability benchmark with controls that preserve spectra.",
        claim_gate="Finite-grid nullspaces and spectral drift are not IR pole data.",
    ),
    AuditRow(
        audit_id="physics_bridge_proxy_boundary",
        celestial_requirement="Prevent current physics-facing proxies from being promoted into established gravity claims.",
        haos_current_evidence=(
            "Physics Bridge 53.x already uses PASS/OPEN/WATCH rows and explicitly separates "
            "shell-native proxies from raw local-gradient failures."
        ),
        current_coverage="bounded proxy dictionary",
        status="WATCH",
        boundary_risk="medium",
        next_test="Phase 60 claim-gated physics bridge language update.",
        claim_gate="Keep scalar-carrier observables framed as proxy tests, not ontology.",
    ),
    AuditRow(
        audit_id="reproducibility_and_falsification",
        celestial_requirement="Make failures auditable and rerunnable.",
        haos_current_evidence=(
            "HAOS maintains frozen baselines, generated tables, public reproduction commands, "
            "and explicit FAIL/MARGINAL/PASS semantics."
        ),
        current_coverage="strong operational coverage",
        status="PASS",
        boundary_risk="low",
        next_test="Keep this sidecar generated and versioned.",
        claim_gate="This passes reproducibility discipline, not celestial physics recovery.",
    ),
    AuditRow(
        audit_id="claim_boundary_enforcement",
        celestial_requirement="State clearly that HAOS-IIP has not closed celestial holography.",
        haos_current_evidence=(
            "The saved Phase 57 roadmap and this audit require OPEN status for untested "
            "BMS, Virasoro, S-matrix, soft-theorem, and memory claims."
        ),
        current_coverage="explicit boundary enforcement",
        status="PASS",
        boundary_risk="low",
        next_test="Use the generated table as the gate for future physics-bridge text.",
        claim_gate="No celestial-facing claim may be promoted without a row-specific test.",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Phase 57 celestial boundary audit.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=AUDIT_DIR,
        help="Directory where the generated CSV, JSON, and markdown audit are written.",
    )
    return parser.parse_args()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def status_counts(rows: tuple[AuditRow, ...]) -> dict[str, int]:
    statuses = ("PASS", "WATCH", "OPEN")
    return {status: sum(1 for row in rows if row.status == status) for status in statuses}


def build_status(rows: tuple[AuditRow, ...]) -> dict[str, Any]:
    counts = status_counts(rows)
    celestial_open = [
        row.audit_id
        for row in rows
        if row.status == "OPEN" and row.audit_id != "physics_bridge_proxy_boundary"
    ]
    return {
        "phase": "57",
        "bridge_name": "celestial_boundary_audit",
        "bridge_status": "OPEN" if celestial_open else "WATCH",
        "failure_mode": "CELESTIAL_REQUIREMENTS_NOT_CLOSED" if celestial_open else "PROXY_BOUNDARY_WATCH",
        "success_definition": (
            "Phase 57 succeeds by making the celestial boundary explicit. It does not "
            "attempt to close celestial holography."
        ),
        "status_counts": counts,
        "open_celestial_requirements": celestial_open,
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "58_spherical_harmonic_control_probe",
        "non_claims": [
            "no BMS recovery claim",
            "no Virasoro or celestial CFT recovery claim",
            "no S-matrix reconstruction claim",
            "no soft-theorem or collinear-factorization recovery claim",
            "no gravitational-memory observable claim",
            "no continuum spacetime ontology claim",
        ],
    }


def write_csv(path: Path, rows: tuple[AuditRow, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def write_status(path: Path, rows: tuple[AuditRow, ...]) -> dict[str, Any]:
    status = build_status(rows)
    path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def write_report(path: Path, rows: tuple[AuditRow, ...], status: dict[str, Any]) -> None:
    table_rows = "\n".join(
        "| {audit_id} | {status} | {current_coverage} | {boundary_risk} | {claim_gate} |".format(
            **asdict(row)
        )
        for row in rows
    )
    open_rows = "\n".join(
        f"- `{row.audit_id}`: {row.next_test}" for row in rows if row.status == "OPEN"
    )
    if not open_rows:
        open_rows = "- none"

    report = f"""# Phase 57 Celestial Boundary Audit

## Scope

This file is generated by `experiments/physics_bridge/celestial_boundary_audit/run_celestial_boundary_audit.py`.

This is an external physics-bridge sidecar. It does not modify HAOS core, re-run scalar simulations, or claim that HAOS-IIP recovers celestial holography.

The audit uses celestial holography as a high-standard boundary check: BMS symmetry, Virasoro/superrotation structure, Mellin-basis conformal weights, soft/collinear factorization, S-matrix unitarity, and gravitational memory are treated as row-specific requirements, not metaphors.

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- `PASS`: {status['status_counts']['PASS']}
- `WATCH`: {status['status_counts']['WATCH']}
- `OPEN`: {status['status_counts']['OPEN']}
- core modified: `{str(status['core_modified']).lower()}`

Phase 57 is therefore a successful boundary audit and an `OPEN` celestial bridge.

## Audit Table

| Audit row | Status | Current HAOS coverage | Risk | Claim gate |
| --- | --- | --- | --- | --- |
{table_rows}

## Open Requirements

{open_rows}

## Interpretation

HAOS-IIP remains a reproducible pre-geometric interaction filter with strong audit discipline. Under a celestial-holography lens, that is not enough to claim flat-space holography, BMS/Virasoro recovery, soft-theorem recovery, S-matrix reconstruction, or gravitational-memory observables.

The correct bridge language is:

> HAOS-IIP supplies reproducible interaction-invariance telemetry. Any contact with established physics must pass separate symmetry, boundary, scattering, and observable tests.

## Next

Phase 58 should build a known-target spherical-harmonic control probe before making any comparison between HAOS spectral structures and celestial-style boundary data. Passing that probe would still be only a boundary-geometry sanity check, not celestial holography.

## Authority

- CSV: `experiments/physics_bridge/celestial_boundary_audit/{CSV_NAME}`
- JSON: `experiments/physics_bridge/celestial_boundary_audit/{STATUS_NAME}`
- Generator: `experiments/physics_bridge/celestial_boundary_audit/run_celestial_boundary_audit.py`
"""
    path.write_text(report, encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / CSV_NAME
    status_path = output_dir / STATUS_NAME
    report_path = output_dir / REPORT_NAME

    write_csv(csv_path, ROWS)
    status = write_status(status_path, ROWS)
    write_report(report_path, ROWS, status)

    print("Phase 57 celestial boundary audit generated.")
    print(f"CSV: {relative(csv_path)}")
    print(f"JSON: {relative(status_path)}")
    print(f"Report: {relative(report_path)}")
    print(
        "Statement: bridge_status is OPEN by design; this audit protects the claim "
        "boundary rather than closing celestial holography."
    )


if __name__ == "__main__":
    main()
