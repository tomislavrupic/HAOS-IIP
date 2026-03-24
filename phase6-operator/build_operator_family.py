#!/usr/bin/env python3

from __future__ import annotations

import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_core import build_graph, read_json, relpath, write_csv_rows, write_json

PHASE_ROOT = ROOT / "phase6-operator"
PHASE4_BUNDLE_PATH = ROOT / "phase4-sector-freeze" / "runs" / "phase4_sector_freeze_bundle_latest.json"
PHASE4_CONFIG_PATH = ROOT / "phase4-sector-freeze" / "configs" / "stage24_4_stable_smeared_law_and_ordering_extension_boundary_runs.json"
PHASE5_LEDGER_PATH = ROOT / "phase5-readout" / "runs" / "phase5_runs_ledger_latest.json"

RUN_LEDGER_PATH = PHASE_ROOT / "phase6_runs.json"
REFINEMENT_LEDGER_PATH = PHASE_ROOT / "phase6_refinement_ledger.csv"
MANIFEST_PATH = PHASE_ROOT / "phase6_operator_manifest.json"
FREEZE_NOTE_PATH = PHASE_ROOT / "phase6_operator_freeze.md"

OPERATOR_CLASS = "cochain_laplacian"
OPERATOR_LABEL = "branch_local_periodic_dk2d_cochain_laplacian"
OPERATOR_RULE = (
    "For each refinement level, build the periodic DK2D complex with the frozen Phase IV "
    "base epsilon and set O_h to the full block cochain Laplacian delta_h."
)
REFINEMENT_MULTIPLIERS = (1, 2, 3, 4)
TRANSLATION_SHIFT = (1, 2)
SMALL_EIGEN_SAMPLE_COUNT = 12
LARGE_EIGEN_SAMPLE_COUNT = 6
RAYLEIGH_TRIAL_COUNT = 6
BASE_RNG_SEED = 20260321
HERMITICITY_TOL = 1.0e-12
SEMIBOUNDED_TOL = 1.0e-8
NULLSPACE_TOL = 1.0e-8
REORDER_TOL = 1.0e-8
CLAIM_BOUNDARY = (
    "Phase VI authority is limited to the frozen cochain-Laplacian family built on the "
    "periodic DK2D branch-local complex across the declared refinement hierarchy. It does not "
    "assert continuum Laplacians, heat-kernel asymptotics, geometric coefficients, spectral "
    "invariants, or continuum correspondence."
)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def timestamp_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def closure_statement(success: bool) -> str:
    if success:
        return "Phase VI now has a frozen operator family suitable for spectral feasibility testing."
    return "Phase VI does not yet have a sufficiently frozen operator family for spectral feasibility testing."


def load_frozen_inputs() -> dict[str, Any]:
    phase4_bundle = read_json(PHASE4_BUNDLE_PATH)
    if phase4_bundle.get("phase_name") != "phase4-sector-freeze":
        raise ValueError("Phase VI requires the frozen Phase IV bundle.")
    phase4_summary = dict(phase4_bundle.get("summary", {}))
    if not bool(phase4_summary.get("phaseIV_freeze_valid")):
        raise ValueError("Phase VI requires a valid frozen Phase IV bundle.")

    phase4_config = read_json(PHASE4_CONFIG_PATH)
    common_fields = dict(phase4_config["common_fields"])
    if str(common_fields["boundary_type"]) != "periodic":
        raise ValueError("Phase VI freezes only the periodic branch-local operator family.")

    phase5_ledger = read_json(PHASE5_LEDGER_PATH)
    if not bool(phase5_ledger.get("success")):
        raise ValueError("Phase VI requires the authoritative Phase V ledger to have succeeded.")

    sector_summary = dict(phase5_ledger["baseline_frozen_run"]["summary"])
    return {
        "phase4_bundle": phase4_bundle,
        "phase4_summary": phase4_summary,
        "phase4_config": phase4_config,
        "common_fields": common_fields,
        "phase5_ledger": phase5_ledger,
        "phase5_summary": sector_summary,
    }


def cochain_laplacian_for_level(n_side: int, epsilon: float) -> tuple[sp.csr_matrix, tuple[int, int, int]]:
    graph = build_graph(
        {
            "kind": "dk2d_periodic",
            "n_side": int(n_side),
            "epsilon": float(epsilon),
            "cycle_phase_x": 0.0,
            "cycle_phase_y": 0.0,
        }
    )
    operator = graph.delta_h.tocsr()
    imag_max = float(np.max(np.abs(operator.data.imag))) if operator.nnz else 0.0
    if imag_max > 1.0e-12:
        raise ValueError(f"unexpected imaginary component in frozen operator family: {imag_max}")
    return operator.real.astype(float).tocsr(), tuple(int(value) for value in graph.block_sizes)


def sample_eigenvalues(operator: sp.csr_matrix, k: int, which: str) -> np.ndarray:
    dim = int(operator.shape[0])
    if dim <= max(32, k + 2):
        dense_values = np.linalg.eigvalsh(operator.toarray())
        if which == "SA":
            return np.sort(dense_values)[:k]
        return np.sort(dense_values)[-k:]
    ncv = min(dim - 1, max(2 * k + 1, 24))
    values = spla.eigsh(operator, k=k, which=which, return_eigenvectors=False, tol=1.0e-10, ncv=ncv)
    return np.sort(np.real_if_close(values))


def translation_permutation(n_side: int, shift_i: int, shift_j: int) -> np.ndarray:
    n_nodes = n_side * n_side
    n_x = n_nodes
    n_y = n_nodes
    n_faces = n_nodes
    total = n_nodes + n_x + n_y + n_faces
    perm = np.empty(total, dtype=int)

    def node_index(i: int, j: int) -> int:
        return i * n_side + j

    def x_edge_index(i: int, j: int) -> int:
        return n_nodes + i * n_side + j

    def y_edge_index(i: int, j: int) -> int:
        return n_nodes + n_x + i * n_side + j

    def face_index(i: int, j: int) -> int:
        return n_nodes + n_x + n_y + i * n_side + j

    for i in range(n_side):
        for j in range(n_side):
            ni = (i + shift_i) % n_side
            nj = (j + shift_j) % n_side
            perm[node_index(i, j)] = node_index(ni, nj)
            perm[x_edge_index(i, j)] = x_edge_index(ni, nj)
            perm[y_edge_index(i, j)] = y_edge_index(ni, nj)
            perm[face_index(i, j)] = face_index(ni, nj)
    return perm


def random_rayleigh_samples(operator: sp.csr_matrix, trial_count: int, seed: int) -> list[float]:
    rng = np.random.default_rng(seed)
    dim = int(operator.shape[0])
    samples: list[float] = []
    for _ in range(trial_count):
        vec = rng.normal(size=dim)
        vec /= max(float(np.linalg.norm(vec)), 1.0e-12)
        value = float(vec @ (operator @ vec))
        samples.append(value)
    return samples


def round_float(value: float, digits: int = 12) -> float:
    return round(float(value), digits)


def compute_level_summary(
    level_index: int,
    multiplier: int,
    base_resolution: int,
    epsilon: float,
) -> dict[str, Any]:
    n_side = int(base_resolution * multiplier)
    h_value = 1.0 / float(n_side)
    operator, block_sizes = cochain_laplacian_for_level(n_side, epsilon)
    dim = int(operator.shape[0])
    nnz = int(operator.nnz)
    density = float(nnz / (dim * dim))
    sparsity_fraction = 1.0 - density

    hermitian_diff = operator - operator.transpose()
    hermitian_diff_norm = float(spla.norm(hermitian_diff, ord="fro"))
    operator_norm = max(float(spla.norm(operator, ord="fro")), 1.0e-12)
    hermiticity_defect = hermitian_diff_norm / operator_norm

    small_values = sample_eigenvalues(
        operator,
        min(SMALL_EIGEN_SAMPLE_COUNT, max(2, dim - 2)),
        "SA",
    )
    large_values = sample_eigenvalues(
        operator,
        min(LARGE_EIGEN_SAMPLE_COUNT, max(2, dim - 2)),
        "LA",
    )
    rayleigh_samples = random_rayleigh_samples(operator, RAYLEIGH_TRIAL_COUNT, BASE_RNG_SEED + n_side)
    minimum_rayleigh = min(rayleigh_samples)

    smallest_eigenvalue = float(small_values[0])
    nullspace_estimate = int(np.sum(np.abs(small_values) <= NULLSPACE_TOL))
    first_positive_eigenvalue = next(
        (float(value) for value in small_values if value > NULLSPACE_TOL),
        float("nan"),
    )
    spectral_radius = float(max(abs(float(value)) for value in large_values))

    permutation = translation_permutation(n_side, *TRANSLATION_SHIFT)
    permuted_operator = operator[permutation, :][:, permutation]
    inverse_permutation = np.argsort(permutation)
    restored_operator = permuted_operator[inverse_permutation, :][:, inverse_permutation]
    reorder_similarity_defect = float(spla.norm(restored_operator - operator, ord="fro")) / operator_norm

    admissible = (
        hermiticity_defect <= HERMITICITY_TOL
        and smallest_eigenvalue >= -SEMIBOUNDED_TOL
        and minimum_rayleigh >= -SEMIBOUNDED_TOL
        and reorder_similarity_defect <= REORDER_TOL
    )

    return {
        "level_id": f"R{level_index}",
        "refinement_multiplier": int(multiplier),
        "n_side": int(n_side),
        "h": round_float(h_value),
        "dimension": dim,
        "zero_form_size": int(block_sizes[0]),
        "one_form_size": int(block_sizes[1]),
        "two_form_size": int(block_sizes[2]),
        "nnz": nnz,
        "density": round_float(density),
        "sparsity_fraction": round_float(sparsity_fraction),
        "hermiticity_defect": round_float(hermiticity_defect),
        "minimum_random_rayleigh": round_float(minimum_rayleigh),
        "sampled_smallest_eigenvalue": round_float(smallest_eigenvalue),
        "sampled_first_positive_eigenvalue": round_float(first_positive_eigenvalue),
        "sampled_spectral_radius": round_float(spectral_radius),
        "nullspace_estimate": int(nullspace_estimate),
        "node_reordering_similarity_defect": round_float(reorder_similarity_defect),
        "admissible": bool(admissible),
        "sampled_smallest_eigenvalues": [round_float(value) for value in small_values.tolist()],
        "sampled_largest_eigenvalues": [round_float(value) for value in large_values.tolist()],
        "random_rayleigh_samples": [round_float(value) for value in rayleigh_samples],
    }


def aggregate_level_summaries(levels: list[dict[str, Any]]) -> dict[str, Any]:
    nullspace_values = [int(level["nullspace_estimate"]) for level in levels]
    spectral_radii = [float(level["sampled_spectral_radius"]) for level in levels]
    first_positive_values = [float(level["sampled_first_positive_eigenvalue"]) for level in levels]
    all_admissible = all(bool(level["admissible"]) for level in levels)
    spectral_radius_non_decreasing = all(
        spectral_radii[index + 1] >= spectral_radii[index] - 1.0e-10
        for index in range(len(spectral_radii) - 1)
    )
    nullspace_consistent = len(set(nullspace_values)) == 1
    reordering_consistent = all(
        float(level["node_reordering_similarity_defect"]) <= REORDER_TOL
        for level in levels
    )
    first_positive_strictly_positive = all(
        value > NULLSPACE_TOL for value in first_positive_values
    )
    return {
        "all_levels_admissible": bool(all_admissible),
        "spectral_radius_non_decreasing": bool(spectral_radius_non_decreasing),
        "nullspace_estimate_consistent": bool(nullspace_consistent),
        "node_reordering_consistent": bool(reordering_consistent),
        "first_positive_eigenvalue_positive": bool(first_positive_strictly_positive),
        "max_hermiticity_defect": round_float(max(float(level["hermiticity_defect"]) for level in levels)),
        "minimum_sampled_eigenvalue": round_float(min(float(level["sampled_smallest_eigenvalue"]) for level in levels)),
        "minimum_random_rayleigh": round_float(min(float(level["minimum_random_rayleigh"]) for level in levels)),
        "nullspace_estimate_values": nullspace_values,
        "spectral_radius_values": [round_float(value) for value in spectral_radii],
        "first_positive_eigenvalue_values": [round_float(value) for value in first_positive_values],
    }


def make_csv_rows(levels: list[dict[str, Any]]) -> tuple[list[str], list[dict[str, Any]]]:
    fieldnames = [
        "level_id",
        "refinement_multiplier",
        "n_side",
        "h",
        "dimension",
        "zero_form_size",
        "one_form_size",
        "two_form_size",
        "nnz",
        "density",
        "sparsity_fraction",
        "hermiticity_defect",
        "minimum_random_rayleigh",
        "sampled_smallest_eigenvalue",
        "sampled_first_positive_eigenvalue",
        "sampled_spectral_radius",
        "nullspace_estimate",
        "node_reordering_similarity_defect",
        "admissible",
    ]
    rows = [{name: level[name] for name in fieldnames} for level in levels]
    return fieldnames, rows


def build_manifest(
    timestamp: str,
    frozen_inputs: dict[str, Any],
    levels: list[dict[str, Any]],
    aggregate: dict[str, Any],
    overall_success: bool,
) -> dict[str, Any]:
    common_fields = dict(frozen_inputs["common_fields"])
    phase5_summary = dict(frozen_inputs["phase5_summary"])
    return {
        "timestamp": timestamp,
        "phase": 6,
        "phase_name": "phase6-operator",
        "status": "passed" if overall_success else "failed",
        "selected_operator_class": OPERATOR_CLASS,
        "operator_label": OPERATOR_LABEL,
        "operator_rule": OPERATOR_RULE,
        "frozen_inputs": {
            "phase4_bundle": relpath(PHASE4_BUNDLE_PATH),
            "phase4_config": relpath(PHASE4_CONFIG_PATH),
            "phase5_runs_ledger": relpath(PHASE5_LEDGER_PATH),
            "phase4_sector_identifier": str(phase5_summary["sector_identifier"]),
            "operator_sector": str(common_fields["operator_sector"]),
            "boundary_type": str(common_fields["boundary_type"]),
            "graph_type": str(common_fields["graph_type"]),
            "kernel_type": str(common_fields["kernel_type"]),
            "base_resolution": int(common_fields["base_resolution"]),
            "refined_resolution": int(common_fields["refined_resolution"]),
            "base_epsilon": float(common_fields["base_epsilon"]),
        },
        "refinement_hierarchy": {
            "h_definition": "h = 1 / n_side on the unit-periodic lattice; it is the lattice-spacing surrogate used to order refinements.",
            "level_parameter": "n_side",
            "construction_rule_constant": True,
            "refinement_multipliers": list(REFINEMENT_MULTIPLIERS),
            "levels": [
                {
                    "level_id": level["level_id"],
                    "refinement_multiplier": level["refinement_multiplier"],
                    "n_side": level["n_side"],
                    "h": level["h"],
                }
                for level in levels
            ],
        },
        "diagnostic_thresholds": {
            "hermiticity_tolerance": HERMITICITY_TOL,
            "semiboundedness_tolerance": SEMIBOUNDED_TOL,
            "nullspace_tolerance": NULLSPACE_TOL,
            "reordering_tolerance": REORDER_TOL,
        },
        "aggregate_diagnostics": aggregate,
        "level_summaries": levels,
        "artifacts": {
            "technical_note": relpath(FREEZE_NOTE_PATH),
            "manifest_json": relpath(MANIFEST_PATH),
            "construction_script": relpath(PHASE_ROOT / "build_operator_family.py"),
            "refinement_ledger_csv": relpath(REFINEMENT_LEDGER_PATH),
            "runs_json": relpath(RUN_LEDGER_PATH),
        },
        "success": bool(overall_success),
        "claim_boundary": CLAIM_BOUNDARY,
        "conclusion": closure_statement(overall_success),
    }


def build_runs_ledger(
    timestamp: str,
    levels: list[dict[str, Any]],
    aggregate: dict[str, Any],
    overall_success: bool,
    runtime_seconds: float,
) -> dict[str, Any]:
    return {
        "timestamp": timestamp,
        "phase": 6,
        "phase_name": "phase6-operator",
        "runner": relpath(PHASE_ROOT / "build_operator_family.py"),
        "inputs": {
            "phase4_bundle": relpath(PHASE4_BUNDLE_PATH),
            "phase4_config": relpath(PHASE4_CONFIG_PATH),
            "phase5_runs_ledger": relpath(PHASE5_LEDGER_PATH),
        },
        "selected_operator_class": OPERATOR_CLASS,
        "refinement_levels": [
            {
                "level_id": level["level_id"],
                "n_side": level["n_side"],
                "h": level["h"],
                "dimension": level["dimension"],
                "nnz": level["nnz"],
                "admissible": level["admissible"],
                "sampled_smallest_eigenvalue": level["sampled_smallest_eigenvalue"],
                "sampled_first_positive_eigenvalue": level["sampled_first_positive_eigenvalue"],
                "sampled_spectral_radius": level["sampled_spectral_radius"],
                "nullspace_estimate": level["nullspace_estimate"],
                "node_reordering_similarity_defect": level["node_reordering_similarity_defect"],
            }
            for level in levels
        ],
        "aggregate_diagnostics": aggregate,
        "runtime_seconds": round_float(runtime_seconds, 6),
        "artifacts": {
            "technical_note": relpath(FREEZE_NOTE_PATH),
            "manifest_json": relpath(MANIFEST_PATH),
            "refinement_ledger_csv": relpath(REFINEMENT_LEDGER_PATH),
            "runs_json": relpath(RUN_LEDGER_PATH),
        },
        "success": bool(overall_success),
    }


def build_freeze_note(
    timestamp: str,
    frozen_inputs: dict[str, Any],
    levels: list[dict[str, Any]],
    aggregate: dict[str, Any],
    overall_success: bool,
) -> str:
    common_fields = dict(frozen_inputs["common_fields"])
    phase5_summary = dict(frozen_inputs["phase5_summary"])
    table_header = "| level | n_side | h | dim | nnz | symmetry defect | min eigen probe | first positive | spectral radius | nullspace | reorder defect | admissible |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"
    table_rows = [
        f"| {level['level_id']} | {level['n_side']} | {level['h']:.12f} | {level['dimension']} | {level['nnz']} | {level['hermiticity_defect']:.12f} | {level['sampled_smallest_eigenvalue']:.12f} | {level['sampled_first_positive_eigenvalue']:.12f} | {level['sampled_spectral_radius']:.12f} | {level['nullspace_estimate']} | {level['node_reordering_similarity_defect']:.12f} | {level['admissible']} |"
        for level in levels
    ]
    return "\n".join(
        [
            "# Phase VI Operator Freeze",
            "",
            f"Freeze timestamp: `{timestamp}`.",
            "",
            "## Objective",
            "",
            "Freeze one branch-valid operator family across a deterministic refinement hierarchy and test whether later spectral work is well-posed on the frozen branch.",
            "",
            "## Selected Operator Class",
            "",
            f"- Selected class: `{OPERATOR_CLASS}`.",
            "- Selection rationale: it is the simplest branch-local family already exposed by the frozen DK2D complex, with direct symmetry structure, direct semiboundedness diagnostics, and no additional auxiliary machinery.",
            "- Not selected: weighted graph Laplacian, because it would discard the higher-degree cochain structure already frozen in the branch.",
            "- Not selected: discrete Dirac / Hodge-Dirac as the primary family, because the squared cochain Laplacian provides the clearer semibounded feasibility test at this stage.",
            "",
            "## Frozen Inputs",
            "",
            f"- Phase IV bundle: `{relpath(PHASE4_BUNDLE_PATH)}`.",
            f"- Phase IV config: `{relpath(PHASE4_CONFIG_PATH)}`.",
            f"- Phase V authority ledger: `{relpath(PHASE5_LEDGER_PATH)}`.",
            f"- Frozen sector identifier: `{phase5_summary['sector_identifier']}`.",
            f"- Branch operator sector: `{common_fields['operator_sector']}` with `{common_fields['boundary_type']}` `{common_fields['graph_type']}` and `{common_fields['kernel_type']}` construction discipline.",
            "",
            "## Refinement Hierarchy",
            "",
            "- Refinement parameter: `h = 1 / n_side`.",
            "- Interpretation of `h`: lattice-spacing surrogate on the unit-periodic branch-local square cell complex.",
            f"- Base resolution: `{common_fields['base_resolution']}`.",
            f"- Refinement multipliers: `{list(REFINEMENT_MULTIPLIERS)}`.",
            "- Deterministic rule: for each multiplier `m`, set `n_side = base_resolution * m`, build the periodic DK2D complex with frozen `base_epsilon`, and define `O_h = delta_h`.",
            "",
            "## Admissibility Diagnostics",
            "",
            table_header,
            *table_rows,
            "",
            "All diagnostics are computed, not asserted. Hermiticity defect is the relative Frobenius defect of `O_h - O_h^T`; semiboundedness is probed by sampled smallest eigenvalues and random Rayleigh quotients; node reordering consistency is tested by a deterministic torus translation of the full cochain basis followed by restoration to the original order.",
            "",
            "## Basic Spectral Sanity Probes",
            "",
            f"- All levels admissible: `{aggregate['all_levels_admissible']}`.",
            f"- Nullspace estimate values: `{aggregate['nullspace_estimate_values']}`.",
            f"- Spectral radius values: `{aggregate['spectral_radius_values']}`.",
            f"- First positive eigenvalue values: `{aggregate['first_positive_eigenvalue_values']}`.",
            f"- Spectral radius non-decreasing: `{aggregate['spectral_radius_non_decreasing']}`.",
            f"- Node reordering consistent: `{aggregate['node_reordering_consistent']}`.",
            "",
            "## Bounded Interpretation",
            "",
            "- The frozen operator family is branch-local only.",
            "- The diagnostics support self-adjointness, semiboundedness, stable kernel estimation, and deterministic refinement construction on the tested hierarchy.",
            "- This is an operator-feasibility result, not a continuum or geometry claim.",
            "",
            "## Explicit Non-Claims",
            "",
            "- No continuum Laplacian is claimed.",
            "- No geometric coefficients are computed.",
            "- No heat-kernel asymptotics are computed.",
            "- No spectral invariants are claimed beyond the finite diagnostic samples listed here.",
            "- No continuum or GR correspondence is asserted.",
            "",
            closure_statement(overall_success),
        ]
    )


def main() -> None:
    started_at = time.perf_counter()
    timestamp = timestamp_iso()
    frozen_inputs = load_frozen_inputs()
    base_resolution = int(frozen_inputs["common_fields"]["base_resolution"])
    epsilon = float(frozen_inputs["common_fields"]["base_epsilon"])

    levels = [
        compute_level_summary(index + 1, multiplier, base_resolution, epsilon)
        for index, multiplier in enumerate(REFINEMENT_MULTIPLIERS)
    ]
    aggregate = aggregate_level_summaries(levels)
    overall_success = (
        bool(aggregate["all_levels_admissible"])
        and bool(aggregate["nullspace_estimate_consistent"])
        and bool(aggregate["node_reordering_consistent"])
        and bool(aggregate["first_positive_eigenvalue_positive"])
    )
    runtime_seconds = time.perf_counter() - started_at

    fieldnames, csv_rows = make_csv_rows(levels)
    write_csv_rows(REFINEMENT_LEDGER_PATH, csv_rows, fieldnames)

    runs_payload = build_runs_ledger(timestamp, levels, aggregate, overall_success, runtime_seconds)
    manifest_payload = build_manifest(timestamp, frozen_inputs, levels, aggregate, overall_success)
    freeze_note = build_freeze_note(timestamp, frozen_inputs, levels, aggregate, overall_success)

    write_json(RUN_LEDGER_PATH, runs_payload)
    write_json(MANIFEST_PATH, manifest_payload)
    write_text(FREEZE_NOTE_PATH, freeze_note)

    print(f"Phase VI operator family success: {overall_success}")
    print(f"Manifest: {relpath(MANIFEST_PATH)}")
    print(f"Ledger: {relpath(REFINEMENT_LEDGER_PATH)}")
    print(closure_statement(overall_success))


if __name__ == "__main__":
    main()
