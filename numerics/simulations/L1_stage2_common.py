#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from L1_transverse_band_test import (
    REPO_ROOT,
    RESULTS,
    PLOTS,
    EXPERIMENT_LOG,
    ComplexData,
    harmonic_basis_from_reference,
    plot_edge_mode,
    periodic_delta,
)
from L1_transverse_scaling_test import (
    ExactProjector as RealExactProjector,
    HarmonicProjector as RealHarmonicProjector,
    build_penalized_transverse_operator as build_penalized_transverse_operator_real,
    harmonic_basis_from_L1,
    inverse_participation_ratio,
    project_transverse as project_transverse_real,
)


def timestamp_slug() -> str:
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding='utf-8')


def stamp_plots(plot_paths: list[Path], timestamp: str) -> list[str]:
    stamped: list[str] = []
    for src in plot_paths:
        dst = PLOTS / f'{timestamp}_{src.name}'
        shutil.copy2(src, dst)
        stamped.append(str(dst.relative_to(REPO_ROOT)))
    return stamped


def save_result_payload(experiment_slug: str, result: dict[str, Any], plot_paths: list[Path]) -> tuple[Path, list[str], str]:
    timestamp = timestamp_slug()
    result_with_plots = {**result, 'plots': [str(path.relative_to(REPO_ROOT)) for path in plot_paths]}
    stamped_json = RESULTS / f'{timestamp}_{experiment_slug}.json'
    write_json(stamped_json, result_with_plots)
    stamped_plots = stamp_plots(plot_paths, timestamp)
    return stamped_json, stamped_plots, timestamp


def append_log(
    title: str,
    config_summary: str,
    result_path: Path,
    stamped_plots: list[str],
    observation: str,
    conclusion: str,
) -> None:
    with EXPERIMENT_LOG.open('a', encoding='utf-8') as handle:
        handle.write(f'\n## {title}\n')
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(f"- Config: {config_summary}\n")
        handle.write(f"- Results: `{result_path.relative_to(REPO_ROOT)}`\n")
        handle.write(f"- Plots: {', '.join(f'`{path}`' for path in stamped_plots)}\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def wrap_points(points: np.ndarray) -> np.ndarray:
    wrapped = np.mod(points, 1.0)
    wrapped[wrapped < 0.0] += 1.0
    return wrapped


def regular_points(n_side: int) -> np.ndarray:
    return np.array(
        [[i / n_side, j / n_side, k / n_side] for i in range(n_side) for j in range(n_side) for k in range(n_side)],
        dtype=float,
    )


def perturbed_points(n_side: int, sigma: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    points = regular_points(n_side)
    return wrap_points(points + sigma * rng.standard_normal(points.shape))


def build_periodic_embedded_complex(
    n_side: int,
    epsilon: float,
    points_override: np.ndarray | None = None,
) -> ComplexData:
    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))
    raw_points = regular_points(n_side) if points_override is None else np.asarray(points_override, dtype=float)
    raw_points = wrap_points(raw_points)

    edges: list[tuple[int, int]] = []
    directions: list[np.ndarray] = []
    midpoints: list[np.ndarray] = []
    edge_weights: list[float] = []
    edge_rows: list[int] = []
    edge_cols: list[int] = []
    edge_data: list[float] = []
    edge_map: dict[tuple[str, int, int, int], int] = {}

    def wrapped_delta_vec(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        delta = b - a
        return (delta + 0.5) % 1.0 - 0.5

    def add_edge(axis: str, i: int, j: int, k: int, ni: int, nj: int, nk: int) -> None:
        u = int(node_index[i, j, k])
        v = int(node_index[ni % n_side, nj % n_side, nk % n_side])
        pu = raw_points[u]
        pv = raw_points[v]
        delta = wrapped_delta_vec(pu, pv)
        length = float(np.linalg.norm(delta))
        if length <= 1e-14:
            unit = np.zeros(3, dtype=float)
        else:
            unit = delta / length
        edge_idx = len(edges)
        edge_map[(axis, i, j, k)] = edge_idx
        edges.append((u, v))
        directions.append(unit)
        midpoints.append(np.mod(pu + 0.5 * delta, 1.0))
        edge_weights.append(math.exp(-(length * length) / (2.0 * epsilon)))
        edge_rows.extend([edge_idx, edge_idx])
        edge_cols.extend([u, v])
        edge_data.extend([-1.0, 1.0])

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_edge('x', i, j, k, (i + 1) % n_side, j, k)
                add_edge('y', i, j, k, i, (j + 1) % n_side, k)
                add_edge('z', i, j, k, i, j, (k + 1) % n_side)

    n_nodes = len(raw_points)
    n_edges = len(edges)
    B0 = sp.coo_matrix((edge_data, (edge_rows, edge_cols)), shape=(n_edges, n_nodes), dtype=float).tocsr()

    face_rows: list[int] = []
    face_cols: list[int] = []
    face_data: list[float] = []
    face_weights: list[float] = []
    face_idx = 0

    def add_face(boundary_keys: list[tuple[str, int, int, int, int]]) -> None:
        nonlocal face_idx
        local_weights: list[float] = []
        for axis, i, j, k, sign in boundary_keys:
            edge_id = edge_map[(axis, i % n_side, j % n_side, k % n_side)]
            face_rows.append(face_idx)
            face_cols.append(edge_id)
            face_data.append(float(sign))
            local_weights.append(edge_weights[edge_id])
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

    C = sp.coo_matrix((face_data, (face_rows, face_cols)), shape=(face_idx, n_edges), dtype=float).tocsr()
    edge_weights_arr = np.asarray(edge_weights, dtype=float)
    face_weights_arr = np.asarray(face_weights, dtype=float)
    d0 = (sp.diags(np.sqrt(edge_weights_arr)) @ B0).tocsr()
    d1 = (sp.diags(np.sqrt(face_weights_arr)) @ C).tocsr()
    lower = (d0 @ d0.T).tocsr()
    upper = (d1.T @ d1).tocsr()
    L1 = (lower + upper).tocsr()
    node_laplacian = (d0.T @ d0).tocsr()
    return ComplexData(
        points=raw_points,
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
        variant='baseline',
        removed_nodes=0,
        removed_requested_edges=0,
    )


def low_modes(operator: sp.csr_matrix, k: int, tol: float) -> tuple[np.ndarray, np.ndarray]:
    dim = operator.shape[0]
    k_eff = min(k, dim - 2) if dim > 2 else 1
    evals, evecs = spla.eigsh(operator, k=k_eff, which='SA', tol=tol)
    order = np.argsort(evals.real)
    return np.asarray(evals[order], dtype=float), np.asarray(evecs[:, order], dtype=float)


def analyze_embedded_baseline(
    n_side: int,
    epsilon: float,
    points: np.ndarray,
    restricted_modes: int,
    phase_modes: int,
    harmonic_tol: float,
    eig_tol: float,
    penalty: float,
) -> dict[str, Any]:
    data = build_periodic_embedded_complex(n_side=n_side, epsilon=epsilon, points_override=points)
    harmonic_basis, _, _ = harmonic_basis_from_L1(data.L1, count=max(restricted_modes, phase_modes, 8), harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = RealExactProjector(data.d0, data.node_laplacian)
    harmonic_projector = RealHarmonicProjector(harmonic_basis)

    phase_evals, phase_vecs = low_modes(data.L1, phase_modes, eig_tol)
    phase_records: list[dict[str, Any]] = []
    for idx in range(len(phase_evals)):
        vec = phase_vecs[:, idx]
        div_norm = float(np.linalg.norm(data.d0.T @ vec))
        curl_norm = float(np.linalg.norm(data.d1 @ vec))
        exact_part = exact_projector.apply(vec)
        harmonic_part = harmonic_projector.apply(vec)
        norm_sq = float(np.dot(vec, vec)) or 1.0
        exact_fraction = float(np.dot(exact_part, exact_part) / norm_sq)
        harmonic_fraction = float(np.dot(harmonic_part, harmonic_part) / norm_sq)
        phase_records.append(
            {
                'mode_index': idx,
                'eigenvalue': float(phase_evals[idx]),
                'divergence_norm': div_norm,
                'curl_norm': curl_norm,
                'exact_fraction': exact_fraction,
                'harmonic_fraction': harmonic_fraction,
                'coexact_fraction': float(max(0.0, 1.0 - exact_fraction - harmonic_fraction)),
            }
        )

    transverse_operator = build_penalized_transverse_operator_real(data, exact_projector, harmonic_projector, penalty=penalty)
    raw_evals, raw_vecs = spla.eigsh(transverse_operator, k=min(restricted_modes + 4, data.upper.shape[0] - 2), which='SA', tol=eig_tol)
    order = np.argsort(raw_evals.real)
    raw_vecs = np.asarray(raw_vecs[:, order], dtype=float)

    restricted_spectrum: list[float] = []
    restricted_records: list[dict[str, Any]] = []
    restricted_vectors: list[np.ndarray] = []
    for idx in range(raw_vecs.shape[1]):
        vec = project_transverse_real(raw_vecs[:, idx], exact_projector, harmonic_projector)
        norm = float(np.linalg.norm(vec))
        if norm <= 1e-12:
            continue
        vec = vec / norm
        lam = float(np.dot(vec, data.upper @ vec))
        restricted_vectors.append(vec)
        restricted_spectrum.append(lam)
        restricted_records.append(
            {
                'mode_index': len(restricted_spectrum) - 1,
                'eigenvalue': lam,
                'divergence_norm': float(np.linalg.norm(data.d0.T @ vec)),
                'curl_norm': float(np.linalg.norm(data.d1 @ vec)),
                'ipr': inverse_participation_ratio(vec),
            }
        )
        if len(restricted_spectrum) >= restricted_modes:
            break

    return {
        'config': {
            'n_side': int(n_side),
            'epsilon': float(epsilon),
            'nodes': int(len(data.points)),
            'edges': int(len(data.edges)),
            'faces': int(len(data.face_weights)),
        },
        'dimensions': {
            'harmonic': int(harmonic_basis.shape[1]),
        },
        'phase_modes': phase_records,
        'restricted_transverse_spectrum': restricted_spectrum,
        'restricted_transverse_modes': restricted_records,
        'restricted_vectors': restricted_vectors,
        'midpoints': data.midpoints,
        'directions': data.directions,
    }


def ensure_matplotlib() -> Any:
    mplconfig = REPO_ROOT / '.mplconfig'
    mplconfig.mkdir(exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(mplconfig))
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    return plt
