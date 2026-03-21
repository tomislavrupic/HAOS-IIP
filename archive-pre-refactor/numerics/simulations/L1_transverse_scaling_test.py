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
RESULTS.mkdir(exist_ok=True)
PLOTS = REPO_ROOT / 'plots'
PLOTS.mkdir(exist_ok=True)
EXPERIMENT_LOG = REPO_ROOT / 'experiments' / 'EXPERIMENT_LOG.md'
MPLCONFIG = REPO_ROOT / '.mplconfig'
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(MPLCONFIG))

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'L1_transverse_scaling_test': {
        'sizes': [6, 8, 10, 12, 14, 16],
        'variants': ['baseline', 'puncture', 'line_defect'],
        'phase_modes': 18,
        'restricted_modes': 6,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
    },
}

VARIANT_LABELS = {
    'baseline': 'baseline periodic torus',
    'puncture': 'single cubic puncture',
    'line_defect': 'line defect',
}

VARIANT_MARKERS = {
    'baseline': 'o',
    'puncture': 's',
    'line_defect': '^',
}


@dataclass(frozen=True)
class ComplexData:
    points: np.ndarray
    edges: list[tuple[int, int]]
    directions: np.ndarray
    midpoints: np.ndarray
    edge_weights: np.ndarray
    face_weights: np.ndarray
    d0: sp.csr_matrix
    d1: sp.csr_matrix
    lower: sp.csr_matrix
    upper: sp.csr_matrix
    L1: sp.csr_matrix
    node_laplacian: sp.csr_matrix
    n_side: int
    variant: str
    removed_nodes: int
    removed_requested_edges: int


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['L1_transverse_scaling_test'] = dict(DEFAULT_CONFIG['L1_transverse_scaling_test'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'L1_transverse_scaling_test'})
        if isinstance(on_disk.get('L1_transverse_scaling_test'), dict):
            merged['L1_transverse_scaling_test'].update(on_disk['L1_transverse_scaling_test'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'L1_transverse_scaling_test'})
        if isinstance(config.get('L1_transverse_scaling_test'), dict):
            merged['L1_transverse_scaling_test'].update(config['L1_transverse_scaling_test'])
    return merged


def central_pair(n_side: int) -> tuple[int, int]:
    return (n_side // 2 - 1, n_side // 2)


def periodic_delta(value: np.ndarray, center: float) -> np.ndarray:
    delta = value - center
    return np.minimum(np.abs(delta), 1.0 - np.abs(delta))


def build_periodic_defect_complex(n_side: int, epsilon: float, variant: str) -> ComplexData:
    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))
    raw_points = np.array(
        [[i / n_side, j / n_side, k / n_side] for i in range(n_side) for j in range(n_side) for k in range(n_side)],
        dtype=float,
    )
    active_mask = np.ones((n_side, n_side, n_side), dtype=bool)
    removed_edge_keys: set[tuple[str, int, int, int]] = set()

    if variant == 'puncture':
        ci0, ci1 = central_pair(n_side)
        for i in (ci0, ci1):
            for j in (ci0, ci1):
                for k in (ci0, ci1):
                    active_mask[i, j, k] = False
    elif variant == 'line_defect':
        ci = n_side // 2
        cj = n_side // 2
        for k in range(n_side):
            removed_edge_keys.add(('z', ci, cj, k))
    elif variant != 'baseline':
        raise ValueError(f'Unsupported variant: {variant}')

    active_old = [int(node_index[i, j, k]) for i in range(n_side) for j in range(n_side) for k in range(n_side) if active_mask[i, j, k]]
    old_to_new = {old: new for new, old in enumerate(active_old)}
    points = raw_points[active_old]

    h = 1.0 / n_side
    edge_weight = math.exp(-(h * h) / (2.0 * epsilon))
    basis = {
        'x': np.array([1.0, 0.0, 0.0]),
        'y': np.array([0.0, 1.0, 0.0]),
        'z': np.array([0.0, 0.0, 1.0]),
    }

    edges: list[tuple[int, int]] = []
    directions: list[np.ndarray] = []
    midpoints: list[np.ndarray] = []
    edge_weights: list[float] = []
    edge_map: dict[tuple[str, int, int, int], int] = {}
    b0_rows: list[int] = []
    b0_cols: list[int] = []
    b0_data: list[float] = []

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
        edge_idx = len(edges)
        edge_map[(axis, i, j, k)] = edge_idx
        edges.append((u, v))
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
        b0_rows.extend([edge_idx, edge_idx])
        b0_cols.extend([u, v])
        b0_data.extend([-1.0, 1.0])

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_edge('x', i, j, k, (i + 1) % n_side, j, k)
                add_edge('y', i, j, k, i, (j + 1) % n_side, k)
                add_edge('z', i, j, k, i, j, (k + 1) % n_side)

    n_nodes = len(points)
    n_edges = len(edges)
    B0 = sp.coo_matrix((b0_data, (b0_rows, b0_cols)), shape=(n_edges, n_nodes), dtype=float).tocsr()

    c_rows: list[int] = []
    c_cols: list[int] = []
    c_data: list[float] = []
    face_weights: list[float] = []
    face_idx = 0

    def add_face(boundary_keys: list[tuple[str, int, int, int, int]]) -> None:
        nonlocal face_idx
        oriented_edges: list[tuple[int, int]] = []
        local_weights: list[float] = []
        for axis, i, j, k, sign in boundary_keys:
            key = (axis, i % n_side, j % n_side, k % n_side)
            if key not in edge_map:
                return
            edge_id = edge_map[key]
            oriented_edges.append((edge_id, sign))
            local_weights.append(edge_weights[edge_id])
        for edge_id, sign in oriented_edges:
            c_rows.append(face_idx)
            c_cols.append(edge_id)
            c_data.append(float(sign))
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
    C = sp.coo_matrix((c_data, (c_rows, c_cols)), shape=(n_faces, n_edges), dtype=float).tocsr()

    edge_weights_arr = np.asarray(edge_weights, dtype=float)
    face_weights_arr = np.asarray(face_weights, dtype=float)
    d0 = sp.diags(np.sqrt(edge_weights_arr)) @ B0
    d1 = sp.diags(np.sqrt(face_weights_arr)) @ C if n_faces else sp.csr_matrix((0, n_edges), dtype=float)
    lower = (d0 @ d0.T).tocsr()
    upper = (d1.T @ d1).tocsr() if n_faces else sp.csr_matrix((n_edges, n_edges), dtype=float)
    L1 = (lower + upper).tocsr()
    node_laplacian = (d0.T @ d0).tocsr()

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


class ExactProjector:
    def __init__(self, d0: sp.csr_matrix, node_laplacian: sp.csr_matrix):
        self.d0 = d0.tocsr()
        reduced = node_laplacian.tocsc()[1:, 1:]
        self.solve = spla.factorized(reduced)
        self.edge_dim = d0.shape[0]
        self.node_dim = node_laplacian.shape[0]

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


def inverse_participation_ratio(vec: np.ndarray) -> float:
    abs_sq = np.abs(vec) ** 2
    norm_sq = float(np.sum(abs_sq))
    if norm_sq == 0.0:
        return 0.0
    return float(np.sum(abs_sq * abs_sq) / (norm_sq * norm_sq))


def low_eigensystem(operator: sp.csr_matrix, k: int, tol: float) -> tuple[np.ndarray, np.ndarray]:
    dim = operator.shape[0]
    k = min(k, dim - 2) if dim > 2 else 1
    evals, evecs = spla.eigsh(operator, k=k, which='SA', tol=tol)
    order = np.argsort(evals)
    return np.asarray(evals[order], dtype=float), np.asarray(evecs[:, order], dtype=float)


def harmonic_basis_from_L1(L1: sp.csr_matrix, count: int, harmonic_tol: float, tol: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    evals, evecs = low_eigensystem(L1, k=count, tol=tol)
    mask = evals < harmonic_tol
    basis = np.asarray(evecs[:, mask], dtype=float)
    if basis.size:
        basis, _ = np.linalg.qr(basis, mode='reduced')
    return basis, evals, evecs


def build_penalized_transverse_operator(
    data: ComplexData,
    exact_projector: ExactProjector,
    harmonic_projector: HarmonicProjector,
    penalty: float,
) -> spla.LinearOperator:
    edge_dim = data.upper.shape[0]

    def matvec(vec: np.ndarray) -> np.ndarray:
        vec = np.asarray(vec, dtype=float)
        return np.asarray(data.upper @ vec, dtype=float) + penalty * exact_projector.apply(vec) + penalty * harmonic_projector.apply(vec)

    return spla.LinearOperator((edge_dim, edge_dim), matvec=matvec, dtype=float)


def project_transverse(vec: np.ndarray, exact_projector: ExactProjector, harmonic_projector: HarmonicProjector) -> np.ndarray:
    projected = vec - exact_projector.apply(vec)
    projected = projected - harmonic_projector.apply(projected)
    return projected


def defect_distances(midpoints: np.ndarray, variant: str, n_side: int) -> np.ndarray:
    if variant == 'puncture':
        center = 0.5 - 0.5 / n_side
        dx = periodic_delta(midpoints[:, 0], center)
        dy = periodic_delta(midpoints[:, 1], center)
        dz = periodic_delta(midpoints[:, 2], center)
        return np.sqrt(dx * dx + dy * dy + dz * dz)
    if variant == 'line_defect':
        center = (n_side // 2) / n_side
        dx = periodic_delta(midpoints[:, 0], center)
        dy = periodic_delta(midpoints[:, 1], center)
        return np.sqrt(dx * dx + dy * dy)
    center = 0.5 - 0.5 / n_side
    dx = periodic_delta(midpoints[:, 0], center)
    dy = periodic_delta(midpoints[:, 1], center)
    dz = periodic_delta(midpoints[:, 2], center)
    return np.sqrt(dx * dx + dy * dy + dz * dz)


def local_support_fraction(mode: np.ndarray, distances: np.ndarray, variant: str, n_side: int) -> float:
    h = 1.0 / n_side
    if variant == 'puncture':
        mask = distances <= 2.5 * h
    elif variant == 'line_defect':
        mask = distances <= 1.5 * h
    else:
        mask = distances <= 2.0 * h
    norm = float(np.sum(np.abs(mode) ** 2)) or 1.0
    return float(np.sum(np.abs(mode[mask]) ** 2) / norm)


def radial_profile(mode: np.ndarray, distances: np.ndarray, bins: int = 20) -> dict[str, list[float]]:
    if len(distances) == 0:
        return {'radius': [], 'density': []}
    order = np.argsort(distances)
    distances = distances[order]
    weights = np.abs(mode[order]) ** 2
    edges = np.linspace(float(np.min(distances)), float(np.max(distances)), bins + 1)
    radii: list[float] = []
    density: list[float] = []
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (distances >= left) & (distances < right if right < edges[-1] else distances <= right)
        if not np.any(mask):
            continue
        radii.append(float(np.mean(distances[mask])))
        density.append(float(np.sum(weights[mask]) / np.sum(mask)))
    return {'radius': radii, 'density': density}


def fit_power_law(sizes: list[int], values: list[float]) -> dict[str, float]:
    x = np.log(np.asarray(sizes, dtype=float))
    y = np.log(np.asarray(values, dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    residuals = y - (slope * x + intercept)
    ss_res = float(np.sum(residuals * residuals))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        'a': float(math.exp(intercept)),
        'p': float(-slope),
        'r2': float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0,
    }


def analyze_case(
    n_side: int,
    epsilon: float,
    variant: str,
    phase_modes: int,
    restricted_modes: int,
    harmonic_tol: float,
    eig_tol: float,
    penalty: float,
) -> dict[str, Any]:
    data = build_periodic_defect_complex(n_side=n_side, epsilon=epsilon, variant=variant)
    harmonic_basis, full_evals, full_evecs = harmonic_basis_from_L1(data.L1, count=max(phase_modes, 8), harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = ExactProjector(data.d0, data.node_laplacian)
    harmonic_projector = HarmonicProjector(harmonic_basis)

    phase_count = min(phase_modes, len(full_evals))
    phase_records: list[dict[str, Any]] = []
    for idx in range(phase_count):
        vec = full_evecs[:, idx]
        div_norm = float(np.linalg.norm(data.d0.T @ vec))
        curl_norm = float(np.linalg.norm(data.d1 @ vec)) if data.d1.shape[0] else 0.0
        exact_part = exact_projector.apply(vec)
        harmonic_part = harmonic_projector.apply(vec)
        norm_sq = float(np.dot(vec, vec)) or 1.0
        exact_fraction = float(np.dot(exact_part, exact_part) / norm_sq)
        harmonic_fraction = float(np.dot(harmonic_part, harmonic_part) / norm_sq)
        coexact_fraction = float(max(0.0, 1.0 - exact_fraction - harmonic_fraction))
        phase_records.append(
            {
                'mode_index': idx,
                'eigenvalue': float(full_evals[idx]),
                'divergence_norm': div_norm,
                'curl_norm': curl_norm,
                'coexact_fraction': coexact_fraction,
            }
        )

    transverse_operator = build_penalized_transverse_operator(data, exact_projector, harmonic_projector, penalty=penalty)
    restricted_evals, restricted_vecs = spla.eigsh(transverse_operator, k=min(restricted_modes, data.upper.shape[0] - 2), which='SA', tol=eig_tol)
    order = np.argsort(restricted_evals)
    restricted_evals = np.asarray(restricted_evals[order], dtype=float)
    restricted_vecs = np.asarray(restricted_vecs[:, order], dtype=float)

    restricted_spectrum: list[float] = []
    restricted_records: list[dict[str, Any]] = []
    for idx in range(len(restricted_evals)):
        vec = project_transverse(restricted_vecs[:, idx], exact_projector, harmonic_projector)
        norm = float(np.linalg.norm(vec))
        if norm == 0.0:
            continue
        vec = vec / norm
        lam = float(np.dot(vec, data.upper @ vec))
        restricted_spectrum.append(lam)
        distances = defect_distances(data.midpoints, variant=variant, n_side=n_side)
        ipr = inverse_participation_ratio(vec)
        near_fraction = local_support_fraction(vec, distances, variant=variant, n_side=n_side)
        profile = radial_profile(vec, distances) if variant in {'puncture', 'line_defect'} else {'radius': [], 'density': []}
        div_norm = float(np.linalg.norm(data.d0.T @ vec))
        curl_norm = float(np.linalg.norm(data.d1 @ vec)) if data.d1.shape[0] else 0.0
        exact_part = exact_projector.apply(vec)
        harmonic_part = harmonic_projector.apply(vec)
        exact_fraction = float(np.dot(exact_part, exact_part))
        harmonic_fraction = float(np.dot(harmonic_part, harmonic_part))
        coexact_fraction = float(max(0.0, 1.0 - exact_fraction - harmonic_fraction))
        restricted_records.append(
            {
                'mode_index': idx,
                'eigenvalue': lam,
                'divergence_norm': div_norm,
                'curl_norm': curl_norm,
                'coexact_fraction': coexact_fraction,
                'ipr': ipr,
                'near_defect_fraction': near_fraction,
                'radial_profile': profile,
            }
        )

    return {
        'label': f'{variant}_n{n_side}',
        'config': {
            'variant': variant,
            'variant_label': VARIANT_LABELS[variant],
            'n_side': int(n_side),
            'nodes': int(len(data.points)),
            'edges': int(len(data.edges)),
            'faces': int(len(data.face_weights)),
            'removed_nodes': int(data.removed_nodes),
            'removed_requested_edges': int(data.removed_requested_edges),
            'epsilon': float(epsilon),
        },
        'dimensions': {
            'harmonic': int(harmonic_basis.shape[1]),
        },
        'phase_modes': phase_records,
        'restricted_transverse_spectrum': restricted_spectrum,
        'restricted_transverse_modes': restricted_records,
    }


def make_scaling_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> dict[str, dict[str, float]]:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    fits: dict[str, dict[str, float]] = {}
    fig, ax = plt.subplots(figsize=(7, 5))
    for variant in ['baseline', 'puncture', 'line_defect']:
        values = [cases[f'{variant}_n{n}']['restricted_transverse_spectrum'][0] for n in sizes]
        fit = fit_power_law(sizes, values)
        fits[variant] = fit
        ax.loglog(sizes, values, marker=VARIANT_MARKERS[variant], label=f"{VARIANT_LABELS[variant]} (p={fit['p']:.2f})")
        fit_curve = [fit['a'] / (n ** fit['p']) for n in sizes]
        ax.loglog(sizes, fit_curve, linestyle='--', alpha=0.6)
    ax.set_xlabel('lattice side n')
    ax.set_ylabel('lowest restricted eigenvalue')
    ax.set_title('Restricted transverse floor scaling')
    ax.grid(alpha=0.25, which='both')
    ax.legend(fontsize=8)
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    return fits


def make_ipr_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> None:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    fig, ax = plt.subplots(figsize=(7, 4))
    for variant in ['baseline', 'puncture', 'line_defect']:
        values = [cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['ipr'] for n in sizes]
        ax.plot(sizes, values, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    ax.set_xlabel('lattice side n')
    ax.set_ylabel('IPR of lowest restricted mode')
    ax.set_title('Lowest restricted mode localization')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_phase_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> None:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey=True)
    scatter = None
    for ax, n_side in zip(axes.ravel(), sizes):
        for variant in ['baseline', 'puncture', 'line_defect']:
            case = cases[f'{variant}_n{n_side}']
            x = [record['divergence_norm'] for record in case['phase_modes']]
            y = [record['curl_norm'] for record in case['phase_modes']]
            c = [record['eigenvalue'] for record in case['phase_modes']]
            scatter = ax.scatter(x, y, c=c, cmap='viridis', s=26, marker=VARIANT_MARKERS[variant], label=variant, edgecolors='none')
        ax.set_title(f'n={n_side}')
        ax.grid(alpha=0.2)
    for ax in axes[-1]:
        ax.set_xlabel(r'$||d_0^* a_k||$')
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$||d_1 a_k||$')
    handles = [plt.Line2D([0], [0], marker=VARIANT_MARKERS[v], linestyle='', color='black', label=VARIANT_LABELS[v]) for v in ['baseline', 'puncture', 'line_defect']]
    fig.legend(handles=handles, loc='upper center', ncol=3, frameon=False)
    if scatter is not None:
        cbar = fig.colorbar(scatter, ax=axes.ravel().tolist(), shrink=0.85)
        cbar.set_label(r'$\lambda_k$')
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_mode_profile_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> None:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, variant in zip(axes, ['puncture', 'line_defect']):
        for n_side in sizes:
            profile = cases[f'{variant}_n{n_side}']['restricted_transverse_modes'][0]['radial_profile']
            ax.plot(profile['radius'], profile['density'], marker='o', label=f'n={n_side}')
        ax.set_title(VARIANT_LABELS[variant])
        ax.set_xlabel('distance from defect')
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=2)
    axes[0].set_ylabel('mean mode density')
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def summarize_results(cases: dict[str, dict[str, Any]], fits: dict[str, dict[str, float]]) -> tuple[str, str, dict[str, Any]]:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    branch_metrics: dict[str, dict[str, Any]] = {}
    baseline_ipr = [cases[f'baseline_n{n}']['restricted_transverse_modes'][0]['ipr'] for n in sizes]
    baseline_extended = all(baseline_ipr[idx + 1] < baseline_ipr[idx] for idx in range(len(baseline_ipr) - 1))

    defect_localized = {}
    for variant in ['puncture', 'line_defect']:
        ipr_values = [cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['ipr'] for n in sizes]
        near_values = [cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['near_defect_fraction'] for n in sizes]
        defect_localized[variant] = {
            'ipr_decreasing': all(ipr_values[idx + 1] < ipr_values[idx] for idx in range(len(ipr_values) - 1)),
            'near_defect_decreasing': all(near_values[idx + 1] <= near_values[idx] + 1e-10 for idx in range(len(near_values) - 1)),
            'final_near_defect_fraction': float(near_values[-1]),
            'final_ipr': float(ipr_values[-1]),
        }

    if baseline_extended and all(metrics['ipr_decreasing'] and metrics['near_defect_decreasing'] for metrics in defect_localized.values()):
        observation = 'the restricted floor keeps descending while the lowest-mode IPR drops with lattice size in every branch, and the defect branches do not retain a growing fraction of norm near the defect'
        conclusion = 'the descending restricted floor is behaving more like a continuum transverse-band candidate than a defect-pinned low mode in the tested range'
    else:
        observation = 'the restricted floor descends, but one or more branches retain high localization or defect concentration as the lattice grows'
        conclusion = 'the descending restricted floor still looks more like a defect-localized low mode than a continuum transverse band in the tested range'

    for variant in ['baseline', 'puncture', 'line_defect']:
        branch_metrics[variant] = {
            'fit': fits[variant],
            'lowest_eigenvalues': {str(n): float(cases[f'{variant}_n{n}']['restricted_transverse_spectrum'][0]) for n in sizes},
            'lowest_ipr': {str(n): float(cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['ipr']) for n in sizes},
            'lowest_near_defect_fraction': {str(n): float(cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['near_defect_fraction']) for n in sizes},
        }
    branch_metrics['baseline']['extended_candidate'] = baseline_extended
    branch_metrics['puncture']['localization_flags'] = defect_localized['puncture']
    branch_metrics['line_defect']['localization_flags'] = defect_localized['line_defect']
    verdict = {'branch_metrics': branch_metrics}
    return observation, conclusion, verdict


def sanitize_case(case: dict[str, Any]) -> dict[str, Any]:
    return {
        'label': case['label'],
        'config': case['config'],
        'dimensions': case['dimensions'],
        'phase_modes': case['phase_modes'],
        'restricted_transverse_spectrum': case['restricted_transverse_spectrum'],
        'restricted_transverse_modes': case['restricted_transverse_modes'],
    }


def save_results(result: dict[str, Any], plot_paths: list[str]) -> tuple[Path, list[str], str]:
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    payload = json.dumps({**result, 'plots': plot_paths}, indent=2)
    stamped = RESULTS / f'{timestamp}_L1_transverse_scaling.json'
    latest = RESULTS / 'L1_transverse_scaling_latest.json'
    stamped.write_text(payload, encoding='utf-8')
    latest.write_text(payload, encoding='utf-8')
    stamped_plots: list[str] = []
    for rel_path in plot_paths:
        src = REPO_ROOT / rel_path
        dst = PLOTS / f'{timestamp}_{src.name}'
        shutil.copy2(src, dst)
        stamped_plots.append(str(dst.relative_to(REPO_ROOT)))
    return stamped, stamped_plots, timestamp


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_dir = REPO_ROOT / 'experiments' / 'vector_sector'
    note_dir.mkdir(parents=True, exist_ok=True)
    note_path = note_dir / 'L1_Transverse_Scaling_Test_v1.md'
    sizes = result['config']['sizes']

    def branch_table(variant: str) -> str:
        rows = []
        for n_side in sizes:
            case = result['cases'][f'{variant}_n{n_side}']
            mode = case['restricted_transverse_modes'][0]
            rows.append(
                f"| {n_side} | {case['restricted_transverse_spectrum'][0]:.6f} | {mode['ipr']:.6f} | {mode['near_defect_fraction']:.6f} | {mode['coexact_fraction']:.6f} |"
            )
        return '\n'.join(rows)

    fit_lines = '\n'.join(
        f"- {VARIANT_LABELS[variant]}: `lambda_min ~ {result['fits'][variant]['a']:.4f} / n^{result['fits'][variant]['p']:.3f}` (`R^2 = {result['fits'][variant]['r2']:.4f}`)"
        for variant in ['baseline', 'puncture', 'line_defect']
    )

    note = f"""# L1 Transverse Scaling Test

## Purpose

Determine whether the descending restricted transverse floor is organizing like a continuum transverse band or remaining as a defect-localized low mode.

## Setup

- operator unchanged:
  - `L1 = d0 d0* + d1* d1`
  - `T = d1* d1` restricted to `ker(d0*) intersect (H1)^perp`
- substrate branches:
  - baseline periodic torus
  - single cubic puncture
  - line defect
- sizes: `{sizes}`
- kernel parameter: `epsilon = {result['config']['epsilon']}`

## Lowest restricted eigenvalue scaling

{fit_lines}

## Lowest restricted mode diagnostics

### Baseline periodic torus

| `n` | lowest restricted `lambda` | IPR | near-defect fraction | coexact fraction |
| --- | ---: | ---: | ---: | ---: |
{branch_table('baseline')}

### Single cubic puncture

| `n` | lowest restricted `lambda` | IPR | near-defect fraction | coexact fraction |
| --- | ---: | ---: | ---: | ---: |
{branch_table('puncture')}

### Line defect

| `n` | lowest restricted `lambda` | IPR | near-defect fraction | coexact fraction |
| --- | ---: | ---: | ---: | ---: |
{branch_table('line_defect')}

## Direct result

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Artifacts

- results: `{result_path.relative_to(REPO_ROOT)}`
- plots: {', '.join(f'`{path}`' for path in stamped_plots)}
- timestamp: `{timestamp}`
"""
    note_path.write_text(note, encoding='utf-8')
    return note_path


def append_experiment_log(result: dict[str, Any], result_path: Path, stamped_plots: list[str]) -> None:
    with EXPERIMENT_LOG.open('a', encoding='utf-8') as handle:
        handle.write('\n## L1 transverse scaling test\n')
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            f"- Config: epsilon={result['config']['epsilon']}, sizes={result['config']['sizes']}, variants={result['config']['variants']}, phase_modes={result['config']['phase_modes']}, restricted_modes={result['config']['restricted_modes']}\n"
        )
        handle.write(f"- Results: `{result_path.relative_to(REPO_ROOT)}`\n")
        handle.write(f"- Plots: {', '.join(f'`{path}`' for path in stamped_plots)}\n")
        handle.write(f"- Observation: {result['observation']}\n")
        handle.write(f"- Conclusion: {result['conclusion']}\n")


def run_L1_transverse_scaling_test(
    config: dict[str, Any] | None = None,
    config_path: Path | None = None,
) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path=config_path)
    epsilon = float(cfg.get('epsilon', 0.2))
    experiment_cfg = cfg['L1_transverse_scaling_test']
    sizes = [int(value) for value in experiment_cfg.get('sizes', [6, 8, 10, 12, 14, 16])]
    variants = [str(value) for value in experiment_cfg.get('variants', ['baseline', 'puncture', 'line_defect'])]
    phase_modes = int(experiment_cfg.get('phase_modes', 18))
    restricted_modes = int(experiment_cfg.get('restricted_modes', 6))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))

    cases: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        for variant in variants:
            cases[f'{variant}_n{n_side}'] = analyze_case(
                n_side=n_side,
                epsilon=epsilon,
                variant=variant,
                phase_modes=phase_modes,
                restricted_modes=restricted_modes,
                harmonic_tol=harmonic_tol,
                eig_tol=eig_tol,
                penalty=penalty,
            )

    scaling_plot = PLOTS / 'transverse_floor_scaling_loglog.png'
    ipr_plot = PLOTS / 'transverse_floor_ipr_vs_n.png'
    phase_plot = PLOTS / 'divergence_curl_phase_scaling.png'
    profile_plot = PLOTS / 'transverse_mode_profiles.png'

    fits = make_scaling_plot(cases, scaling_plot)
    make_ipr_plot(cases, ipr_plot)
    make_phase_plot(cases, phase_plot)
    make_mode_profile_plot(cases, profile_plot)
    plot_paths = [
        str(scaling_plot.relative_to(REPO_ROOT)),
        str(ipr_plot.relative_to(REPO_ROOT)),
        str(profile_plot.relative_to(REPO_ROOT)),
        str(phase_plot.relative_to(REPO_ROOT)),
    ]

    observation, conclusion, verdict = summarize_results(cases, fits)
    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'phase_modes': phase_modes,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
        },
        'cases': {label: sanitize_case(case) for label, case in cases.items()},
        'fits': fits,
        'observation': observation,
        'conclusion': conclusion,
        'verdict': verdict,
    }
    result_path, stamped_plots, timestamp = save_results(result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_experiment_log(result, result_path, stamped_plots)
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_L1_transverse_scaling_test()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
