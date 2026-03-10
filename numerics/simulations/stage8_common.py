#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS = REPO_ROOT / 'data'
PLOTS = REPO_ROOT / 'plots'
EXPERIMENT_LOG = REPO_ROOT / 'experiments' / 'EXPERIMENT_LOG.md'
for path in (RESULTS, PLOTS):
    path.mkdir(exist_ok=True)

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'harmonic_tol': 1.0e-8,
    'eig_tol': 1.0e-8,
    'penalty': 10.0,
    'stage8_charge_insertion': {
        'n_side': 16,
        'source_strength': 1.0,
        'radial_bins': 18,
    },
    'stage8_gauss_law_test': {
        'sizes': [12, 16, 20],
        'source_strength': 1.0,
        'flux_radii': 12,
    },
    'stage8_dispersion_test': {
        'sizes': [12, 16, 20],
        'restricted_modes': 30,
    },
    'stage8_longitudinal_decoupling': {
        'n_side': 16,
        'restricted_modes': 20,
        'alphas': [0.0, 0.25, 0.5, 1.0, 2.0],
    },
    'stage8_loop_response': {
        'n_side': 16,
        'restricted_mode_index': 0,
        'max_loop_side': 4,
    },
}


@dataclass(frozen=True)
class Stage8Complex:
    points: np.ndarray
    node_index: np.ndarray
    edges: list[tuple[int, int]]
    edge_axes: list[str]
    edge_map: dict[tuple[str, int, int, int], int]
    directions: np.ndarray
    midpoints: np.ndarray
    edge_weights: np.ndarray
    face_weights: np.ndarray
    d0: sp.csr_matrix
    d1: sp.csr_matrix
    L0: sp.csr_matrix
    L1: sp.csr_matrix
    upper: sp.csr_matrix
    n_side: int


def ensure_matplotlib() -> Any:
    mplconfig = REPO_ROOT / '.mplconfig'
    mplconfig.mkdir(exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(mplconfig))
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    return plt


def load_config(config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k not in merged or not isinstance(merged[k], dict)})
        for key, value in on_disk.items():
            if isinstance(value, dict) and isinstance(merged.get(key), dict):
                merged[key].update(value)
    return merged


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


def append_log(title: str, config_summary: str, result_path: Path, stamped_plots: list[str], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open('a', encoding='utf-8') as handle:
        handle.write(f'\n## {title}\n')
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(f"- Config: {config_summary}\n")
        handle.write(f"- Results: `{result_path.relative_to(REPO_ROOT)}`\n")
        handle.write(f"- Plots: {', '.join(f'`{path}`' for path in stamped_plots)}\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def periodic_delta(values: np.ndarray, center: np.ndarray) -> np.ndarray:
    delta = np.asarray(values, dtype=float) - np.asarray(center, dtype=float)
    return (delta + 0.5) % 1.0 - 0.5


def build_periodic_complex(n_side: int, epsilon: float) -> Stage8Complex:
    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))
    points = np.array(
        [[i / n_side, j / n_side, k / n_side] for i in range(n_side) for j in range(n_side) for k in range(n_side)],
        dtype=float,
    )
    h = 1.0 / n_side
    edge_weight = math.exp(-(h * h) / (2.0 * epsilon))
    basis = {
        'x': np.array([1.0, 0.0, 0.0]),
        'y': np.array([0.0, 1.0, 0.0]),
        'z': np.array([0.0, 0.0, 1.0]),
    }

    edges: list[tuple[int, int]] = []
    edge_axes: list[str] = []
    directions: list[np.ndarray] = []
    midpoints: list[np.ndarray] = []
    edge_weights: list[float] = []
    edge_map: dict[tuple[str, int, int, int], int] = {}
    b_rows: list[int] = []
    b_cols: list[int] = []
    b_data: list[float] = []

    def add_edge(axis: str, i: int, j: int, k: int, ni: int, nj: int, nk: int) -> None:
        u = int(node_index[i, j, k])
        v = int(node_index[ni % n_side, nj % n_side, nk % n_side])
        edge_idx = len(edges)
        edge_map[(axis, i, j, k)] = edge_idx
        edges.append((u, v))
        edge_axes.append(axis)
        directions.append(basis[axis])
        midpoints.append(
            np.array(
                [
                    ((i + 0.5) / n_side) % 1.0 if axis == 'x' else i / n_side,
                    ((j + 0.5) / n_side) % 1.0 if axis == 'y' else j / n_side,
                    ((k + 0.5) / n_side) % 1.0 if axis == 'z' else k / n_side,
                ],
                dtype=float,
            )
        )
        edge_weights.append(edge_weight)
        b_rows.extend([edge_idx, edge_idx])
        b_cols.extend([u, v])
        b_data.extend([-1.0, 1.0])

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_edge('x', i, j, k, (i + 1) % n_side, j, k)
                add_edge('y', i, j, k, i, (j + 1) % n_side, k)
                add_edge('z', i, j, k, i, j, (k + 1) % n_side)

    B0 = sp.coo_matrix((b_data, (b_rows, b_cols)), shape=(len(edges), len(points)), dtype=float).tocsr()

    c_rows: list[int] = []
    c_cols: list[int] = []
    c_data: list[float] = []
    face_weights: list[float] = []
    face_idx = 0

    def add_face(boundary_keys: list[tuple[str, int, int, int, int]]) -> None:
        nonlocal face_idx
        local_weights: list[float] = []
        for axis, i, j, k, sign in boundary_keys:
            edge_id = edge_map[(axis, i % n_side, j % n_side, k % n_side)]
            c_rows.append(face_idx)
            c_cols.append(edge_id)
            c_data.append(float(sign))
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

    C = sp.coo_matrix((c_data, (c_rows, c_cols)), shape=(face_idx, len(edges)), dtype=float).tocsr()
    edge_weights_arr = np.asarray(edge_weights, dtype=float)
    face_weights_arr = np.asarray(face_weights, dtype=float)
    d0 = (sp.diags(np.sqrt(edge_weights_arr)) @ B0).tocsr()
    d1 = (sp.diags(np.sqrt(face_weights_arr)) @ C).tocsr()
    lower = (d0 @ d0.T).tocsr()
    upper = (d1.T @ d1).tocsr()
    L1 = (lower + upper).tocsr()
    L0 = (d0.T @ d0).tocsr()

    return Stage8Complex(
        points=points,
        node_index=node_index,
        edges=edges,
        edge_axes=edge_axes,
        edge_map=edge_map,
        directions=np.asarray(directions, dtype=float),
        midpoints=np.asarray(midpoints, dtype=float),
        edge_weights=edge_weights_arr,
        face_weights=face_weights_arr,
        d0=d0,
        d1=d1,
        L0=L0,
        L1=L1,
        upper=upper,
        n_side=n_side,
    )


class ExactProjector:
    def __init__(self, d0: sp.csr_matrix, L0: sp.csr_matrix):
        self.d0 = d0.tocsr()
        reduced = L0.tocsc()[1:, 1:]
        self.solve = spla.factorized(reduced)
        self.node_dim = L0.shape[0]

    def apply(self, vec: np.ndarray) -> np.ndarray:
        rhs = np.asarray(self.d0.T @ vec, dtype=float)
        coeffs = np.zeros(self.node_dim, dtype=float)
        coeffs[1:] = self.solve(rhs[1:])
        return np.asarray(self.d0 @ coeffs, dtype=float)


class HarmonicProjector:
    def __init__(self, basis: np.ndarray):
        self.basis = np.asarray(basis, dtype=float)

    def apply(self, vec: np.ndarray) -> np.ndarray:
        if self.basis.size == 0:
            return np.zeros_like(vec)
        return self.basis @ (self.basis.T @ vec)


def low_eigensystem(operator: sp.csr_matrix | spla.LinearOperator, k: int, tol: float) -> tuple[np.ndarray, np.ndarray]:
    dim = operator.shape[0]
    k_eff = min(k, dim - 2) if dim > 2 else 1
    evals, evecs = spla.eigsh(operator, k=k_eff, which='SA', tol=tol)
    order = np.argsort(evals.real)
    return np.asarray(evals[order], dtype=float), np.asarray(evecs[:, order], dtype=float)


def harmonic_basis_from_L1(L1: sp.csr_matrix, count: int, harmonic_tol: float, tol: float) -> np.ndarray:
    evals, evecs = low_eigensystem(L1, count, tol)
    mask = evals < harmonic_tol
    basis = np.asarray(evecs[:, mask], dtype=float)
    if basis.size:
        basis, _ = np.linalg.qr(basis, mode='reduced')
    return basis


def build_penalized_transverse_operator(data: Stage8Complex, exact_projector: ExactProjector, harmonic_projector: HarmonicProjector, penalty: float) -> spla.LinearOperator:
    edge_dim = data.upper.shape[0]

    def matvec(vec: np.ndarray) -> np.ndarray:
        vec = np.asarray(vec, dtype=float)
        return np.asarray(data.upper @ vec, dtype=float) + penalty * exact_projector.apply(vec) + penalty * harmonic_projector.apply(vec)

    return spla.LinearOperator((edge_dim, edge_dim), matvec=matvec, dtype=float)


def project_transverse(vec: np.ndarray, exact_projector: ExactProjector, harmonic_projector: HarmonicProjector) -> np.ndarray:
    projected = vec - exact_projector.apply(vec)
    projected = projected - harmonic_projector.apply(projected)
    return projected


def inverse_participation_ratio(vec: np.ndarray) -> float:
    abs_sq = np.abs(vec) ** 2
    norm_sq = float(np.sum(abs_sq))
    if norm_sq == 0.0:
        return 0.0
    return float(np.sum(abs_sq * abs_sq) / (norm_sq * norm_sq))


def solve_mean_zero_poisson(L0: sp.csr_matrix, rho: np.ndarray) -> np.ndarray:
    rho = np.asarray(rho, dtype=float)
    coeffs = np.zeros(L0.shape[0], dtype=float)
    reduced = L0.tocsc()[1:, 1:]
    solve = spla.factorized(reduced)
    coeffs[1:] = solve(rho[1:])
    coeffs -= np.mean(coeffs)
    return coeffs


def continuum_transverse_q2(count: int, max_m: int = 8) -> np.ndarray:
    values: list[float] = []
    for mx in range(-max_m, max_m + 1):
        for my in range(-max_m, max_m + 1):
            for mz in range(-max_m, max_m + 1):
                if mx == 0 and my == 0 and mz == 0:
                    continue
                q2 = float(mx * mx + my * my + mz * mz)
                values.extend([q2, q2])
    values.sort()
    if len(values) < count:
        raise ValueError('Increase max_m for enough continuum modes')
    return np.asarray(values[:count], dtype=float)


def fit_linear_dispersion(k2: np.ndarray, lambdas: np.ndarray) -> dict[str, float]:
    A = np.column_stack([k2, np.ones_like(k2)])
    coeffs, *_ = np.linalg.lstsq(A, lambdas, rcond=None)
    c_eff, m_eff = coeffs
    pred = A @ coeffs
    residual = lambdas - pred
    ss_res = float(np.sum(residual * residual))
    ss_tot = float(np.sum((lambdas - np.mean(lambdas)) ** 2))
    return {
        'c_eff': float(c_eff),
        'm_eff': float(m_eff),
        'r2': float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0,
    }


def group_family_spread(lambdas: np.ndarray, q2: np.ndarray) -> float:
    spreads: list[float] = []
    start = 0
    while start < len(q2):
        stop = start + 1
        while stop < len(q2) and q2[stop] == q2[start]:
            stop += 1
        group = lambdas[start:stop]
        if len(group) > 1:
            spreads.append(float(np.std(group) / max(np.mean(np.abs(group)), 1.0e-12)))
        start = stop
    return float(np.mean(spreads)) if spreads else 0.0
