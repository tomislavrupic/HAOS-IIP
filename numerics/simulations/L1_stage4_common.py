#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Any

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import linear_sum_assignment

from L1_stage2_common import PLOTS, REPO_ROOT, save_result_payload, append_log, ensure_matplotlib
from L1_stage3_common import periodic_displacement, make_phase_plot
from L1_transverse_band_test import (
    ComplexData,
    ExactProjector,
    HarmonicProjector,
    VARIANT_LABELS,
    build_penalized_transverse_operator,
    defect_distances,
    harmonic_basis_from_reference,
    inverse_participation_ratio,
    local_support_fraction,
    project_transverse,
    radial_profile,
)

plt = ensure_matplotlib()

DEFORMATION_LABELS = {
    'none': 'undeformed',
    'radial': 'smooth radial deformation',
    'anisotropic': 'anisotropic box deformation',
    'bump': 'localized scalar bump',
}


def load_stage4_base_config(config_path: Path | None = None) -> dict[str, Any]:
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        return json.loads(path.read_text())
    return {}


def periodic_midpoint(p0: np.ndarray, p1: np.ndarray) -> np.ndarray:
    delta = periodic_displacement(np.asarray(p1, dtype=float), np.asarray(p0, dtype=float))
    return (np.asarray(p0, dtype=float) + 0.5 * delta) % 1.0


def deformation_profile(points: np.ndarray, n_side: int, family: str, eta: float) -> tuple[np.ndarray, np.ndarray]:
    center = np.full(3, 0.5 - 0.5 / n_side, dtype=float)
    delta = periodic_displacement(points, center)
    if family == 'none' or abs(eta) <= 0.0:
        return points.copy(), np.zeros(len(points), dtype=float)
    if family == 'radial':
        out = (center + (1.0 + eta) * delta) % 1.0
        profile = np.linalg.norm(delta, axis=1)
        return out, profile
    if family == 'anisotropic':
        out_delta = delta.copy()
        out_delta[:, 0] *= 1.0 + eta
        out = (center + out_delta) % 1.0
        profile = np.abs(delta[:, 0])
        return out, profile
    if family == 'bump':
        sigma = 0.18
        r2 = np.sum(delta * delta, axis=1)
        bump = np.exp(-r2 / (2.0 * sigma * sigma))
        out = (center + (1.0 + eta * bump)[:, None] * delta) % 1.0
        return out, bump
    raise ValueError(f'Unsupported deformation family: {family}')


def build_deformed_complex(
    n_side: int,
    epsilon: float,
    variant: str,
    family: str,
    eta: float,
) -> ComplexData:
    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))
    raw_points = np.array(
        [[i / n_side, j / n_side, k / n_side] for i in range(n_side) for j in range(n_side) for k in range(n_side)],
        dtype=float,
    )
    deformed_raw_points, _ = deformation_profile(raw_points, n_side, family, eta)

    active_mask = np.ones((n_side, n_side, n_side), dtype=bool)
    removed_edge_keys: set[tuple[str, int, int, int]] = set()
    if variant == 'puncture':
        ci0, ci1 = (n_side // 2 - 1, n_side // 2)
        for i in (ci0, ci1):
            for j in (ci0, ci1):
                for k in (ci0, ci1):
                    active_mask[i, j, k] = False
    elif variant == 'line_defect':
        ci = n_side // 2
        cj = n_side // 2
        for k in range(n_side):
            removed_edge_keys.add(('z', ci, cj, k))
    elif variant == 'baseline':
        pass
    else:
        raise ValueError(f'Unsupported variant: {variant}')

    active_old_indices = [int(node_index[i, j, k]) for i in range(n_side) for j in range(n_side) for k in range(n_side) if active_mask[i, j, k]]
    old_to_new = {old: new for new, old in enumerate(active_old_indices)}
    points = deformed_raw_points[active_old_indices]

    edges: list[tuple[int, int]] = []
    directions: list[np.ndarray] = []
    midpoints: list[np.ndarray] = []
    edge_weights: list[float] = []
    edge_map: dict[tuple[str, int, int, int], int] = {}
    b0_rows: list[int] = []
    b0_cols: list[int] = []
    b0_data: list[complex] = []

    def node_active(i: int, j: int, k: int) -> bool:
        return bool(active_mask[i % n_side, j % n_side, k % n_side])

    def add_edge(axis: str, i: int, j: int, k: int, ni: int, nj: int, nk: int) -> None:
        if (axis, i, j, k) in removed_edge_keys:
            return
        if not node_active(i, j, k) or not node_active(ni, nj, nk):
            return
        u_old = int(node_index[i, j, k])
        v_old = int(node_index[ni % n_side, nj % n_side, nk % n_side])
        u = old_to_new[u_old]
        v = old_to_new[v_old]
        p0 = points[u]
        p1 = points[v]
        delta = periodic_displacement(p1, p0)
        length = float(np.linalg.norm(delta))
        direction = delta / length if length > 0 else np.zeros(3, dtype=float)
        weight = math.exp(-(length * length) / (2.0 * epsilon))
        edge_idx = len(edges)
        edge_map[(axis, i, j, k)] = edge_idx
        edges.append((u, v))
        directions.append(direction)
        midpoints.append(periodic_midpoint(p0, p1))
        edge_weights.append(weight)
        b0_rows.extend([edge_idx, edge_idx])
        b0_cols.extend([u, v])
        b0_data.extend([-1.0 + 0.0j, 1.0 + 0.0j])

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_edge('x', i, j, k, (i + 1) % n_side, j, k)
                add_edge('y', i, j, k, i, (j + 1) % n_side, k)
                add_edge('z', i, j, k, i, j, (k + 1) % n_side)

    n_nodes = len(points)
    n_edges = len(edges)
    B0 = sp.coo_matrix((b0_data, (b0_rows, b0_cols)), shape=(n_edges, n_nodes), dtype=complex).tocsr()

    c_rows: list[int] = []
    c_cols: list[int] = []
    c_data: list[complex] = []
    face_weights: list[float] = []
    face_idx = 0

    def add_face(boundary_keys: list[tuple[str, int, int, int, int]]) -> None:
        nonlocal face_idx
        oriented_edges: list[tuple[int, complex]] = []
        local_weights: list[float] = []
        for axis, i, j, k, sign in boundary_keys:
            key = (axis, i % n_side, j % n_side, k % n_side)
            if key not in edge_map:
                return
            edge_id = edge_map[key]
            oriented_edges.append((edge_id, complex(sign)))
            local_weights.append(edge_weights[edge_id])
        for edge_id, coeff in oriented_edges:
            c_rows.append(face_idx)
            c_cols.append(edge_id)
            c_data.append(coeff)
        face_weights.append(float(np.mean(local_weights)))
        face_idx += 1

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_face([
                    ('x', i, j, k, +1),
                    ('y', (i + 1) % n_side, j, k, +1),
                    ('x', i, (j + 1) % n_side, k, -1),
                    ('y', i, j, k, -1),
                ])
                add_face([
                    ('x', i, j, k, +1),
                    ('z', (i + 1) % n_side, j, k, +1),
                    ('x', i, j, (k + 1) % n_side, -1),
                    ('z', i, j, k, -1),
                ])
                add_face([
                    ('y', i, j, k, +1),
                    ('z', i, (j + 1) % n_side, k, +1),
                    ('y', i, j, (k + 1) % n_side, -1),
                    ('z', i, j, k, -1),
                ])

    n_faces = face_idx
    C = sp.coo_matrix((c_data, (c_rows, c_cols)), shape=(n_faces, n_edges), dtype=complex).tocsr()
    edge_weights_arr = np.asarray(edge_weights, dtype=float)
    face_weights_arr = np.asarray(face_weights, dtype=float)
    d0 = (sp.diags(np.sqrt(edge_weights_arr)) @ B0).tocsr()
    d1 = (sp.diags(np.sqrt(face_weights_arr)) @ C).tocsr() if n_faces else sp.csr_matrix((0, n_edges), dtype=complex)
    lower = (d0 @ d0.conj().T).tocsr()
    upper = (d1.conj().T @ d1).tocsr() if n_faces else sp.csr_matrix((n_edges, n_edges), dtype=complex)
    L1 = (lower + upper).tocsr()
    node_laplacian = (d0.conj().T @ d0).tocsr()

    return ComplexData(
        points=points,
        edges=edges,
        directions=np.asarray(directions, dtype=float),
        midpoints=np.asarray(midpoints, dtype=float),
        edge_weights=edge_weights_arr,
        face_weights=face_weights_arr,
        d0=d0,
        d1=d1,
        lower=lower,
        upper=upper,
        L1=L1,
        node_laplacian=node_laplacian,
        n_side=n_side,
        variant=variant,
        removed_nodes=int(np.size(active_mask) - np.count_nonzero(active_mask)),
        removed_requested_edges=int(len(removed_edge_keys)),
    )


def low_eigensystem(operator: sp.csr_matrix | spla.LinearOperator, k: int, tol: float) -> tuple[np.ndarray, np.ndarray]:
    dim = operator.shape[0]
    k_eff = min(k, dim - 2) if dim > 2 else 1
    evals, evecs = spla.eigsh(operator, k=k_eff, which='SA', tol=tol)
    order = np.argsort(evals.real)
    return np.asarray(evals[order], dtype=float), np.asarray(evecs[:, order], dtype=complex)


def analyze_deformed_case(
    n_side: int,
    epsilon: float,
    variant: str,
    family: str,
    eta: float,
    restricted_modes: int,
    harmonic_tol: float,
    eig_tol: float,
    penalty: float,
) -> dict[str, Any]:
    analysis = build_deformed_complex(n_side=n_side, epsilon=epsilon, variant=variant, family=family, eta=eta)
    harmonic_basis = harmonic_basis_from_reference(analysis, count=max(restricted_modes, 8), harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = ExactProjector(analysis.d0, analysis.node_laplacian)
    harmonic_projector = HarmonicProjector(harmonic_basis)
    transverse_operator = build_penalized_transverse_operator(analysis, exact_projector, harmonic_projector, penalty=penalty)
    raw_evals, raw_vecs = low_eigensystem(transverse_operator, k=restricted_modes + 4, tol=eig_tol)

    restricted_spectrum: list[float] = []
    restricted_records: list[dict[str, Any]] = []
    restricted_vectors: list[np.ndarray] = []
    distances = defect_distances(analysis.midpoints, variant=variant, n_side=n_side)
    _, point_profile = deformation_profile(analysis.points, n_side, family, eta)
    midpoint_profile = edge_profile_from_node_profile(analysis, point_profile)
    for idx in range(raw_vecs.shape[1]):
        vec = project_transverse(raw_vecs[:, idx], exact_projector, harmonic_projector)
        norm = float(np.linalg.norm(vec))
        if norm <= 1e-12:
            continue
        vec = vec / norm
        lam = float(np.real(np.vdot(vec, analysis.upper @ vec)))
        if lam < 1e-10 and len(restricted_spectrum) > 0:
            continue
        restricted_spectrum.append(lam)
        restricted_vectors.append(vec)
        restricted_records.append(
            {
                'mode_index': len(restricted_spectrum) - 1,
                'eigenvalue': lam,
                'divergence_norm': float(np.linalg.norm(analysis.d0.conj().T @ vec)),
                'curl_norm': float(np.linalg.norm(analysis.d1 @ vec)) if analysis.d1.shape[0] else 0.0,
                'ipr': inverse_participation_ratio(vec),
                'near_defect_fraction': local_support_fraction(vec, distances, variant=variant, n_side=n_side),
                'radial_profile': radial_profile(vec, distances) if variant in {'puncture', 'line_defect'} else {'radius': [], 'density': []},
            }
        )
        if len(restricted_spectrum) >= restricted_modes:
            break

    return {
        'label': f'{variant}_{family}_{eta:.3f}_n{n_side}',
        'config': {
            'variant': variant,
            'variant_label': VARIANT_LABELS[variant],
            'deformation_family': family,
            'deformation_label': DEFORMATION_LABELS[family],
            'eta': float(eta),
            'n_side': int(n_side),
            'nodes': int(len(analysis.points)),
            'edges': int(len(analysis.edges)),
            'faces': int(len(analysis.face_weights)),
            'epsilon': float(epsilon),
        },
        'dimensions': {'harmonic': int(harmonic_basis.shape[1])},
        'restricted_transverse_spectrum': restricted_spectrum,
        'restricted_transverse_modes': restricted_records,
        'restricted_vectors': restricted_vectors,
        'midpoints': analysis.midpoints,
        'directions': analysis.directions,
        'analysis': analysis,
        'exact_projector': exact_projector,
        'harmonic_projector': harmonic_projector,
        'midpoint_profile': midpoint_profile,
    }


def edge_profile_from_node_profile(analysis: ComplexData, point_profile: np.ndarray) -> np.ndarray:
    values = np.zeros(len(analysis.edges), dtype=float)
    for edge_id, (u, v) in enumerate(analysis.edges):
        values[edge_id] = 0.5 * (point_profile[u] + point_profile[v])
    return values


def best_mode_overlap(baseline_vectors: list[np.ndarray], deformed_vectors: list[np.ndarray], count: int) -> dict[str, Any]:
    count_eff = min(count, len(baseline_vectors), len(deformed_vectors))
    if count_eff == 0:
        return {'matrix': [], 'matched': [], 'mean_overlap': 0.0}
    B = np.asarray(baseline_vectors[:count_eff], dtype=complex)
    D = np.asarray(deformed_vectors[:count_eff], dtype=complex)
    overlap = np.abs(B @ D.conj().T)
    rows, cols = linear_sum_assignment(-overlap)
    matched = overlap[rows, cols]
    return {
        'matrix': overlap.tolist(),
        'rows': rows.tolist(),
        'cols': cols.tolist(),
        'matched': matched.tolist(),
        'mean_overlap': float(np.mean(matched)),
    }


def restricted_operator_apply(case: dict[str, Any], vec: np.ndarray) -> np.ndarray:
    analysis = case['analysis']
    exact_projector = case['exact_projector']
    harmonic_projector = case['harmonic_projector']
    projected = project_transverse(vec, exact_projector, harmonic_projector)
    applied = np.asarray(analysis.upper @ projected, dtype=complex)
    return project_transverse(applied, exact_projector, harmonic_projector)


def coupling_matrix(baseline_case: dict[str, Any], deformed_case: dict[str, Any], k: int) -> np.ndarray:
    basis = [np.asarray(vec, dtype=complex) for vec in baseline_case['restricted_vectors'][:k]]
    matrix = np.zeros((len(basis), len(basis)), dtype=complex)
    for j, vec_j in enumerate(basis):
        delta_vec = restricted_operator_apply(deformed_case, vec_j) - restricted_operator_apply(baseline_case, vec_j)
        for i, vec_i in enumerate(basis):
            matrix[i, j] = np.vdot(vec_i, delta_vec)
    return matrix


def unique_family_representatives(lambdas: list[float], tol: float = 1e-9) -> list[int]:
    reps: list[int] = []
    current = None
    for idx, lam in enumerate(lambdas):
        if current is None or abs(lam - current) > tol:
            reps.append(idx)
            current = lam
    return reps
