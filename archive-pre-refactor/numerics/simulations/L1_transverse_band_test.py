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
from matplotlib import cm

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'L1_transverse_band_test': {
        'sizes': [10, 12, 14, 16],
        'variants': ['baseline', 'puncture', 'line_defect', 'flux_tube'],
        'restricted_modes': 20,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': 1.5707963267948966,
    },
}

VARIANT_LABELS = {
    'baseline': 'baseline periodic torus',
    'puncture': 'single cubic puncture',
    'line_defect': 'line defect',
    'flux_tube': 'flux-tube defect',
}

VARIANT_MARKERS = {
    'baseline': 'o',
    'puncture': 's',
    'line_defect': '^',
    'flux_tube': 'D',
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
    merged['L1_transverse_band_test'] = dict(DEFAULT_CONFIG['L1_transverse_band_test'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'L1_transverse_band_test'})
        if isinstance(on_disk.get('L1_transverse_band_test'), dict):
            merged['L1_transverse_band_test'].update(on_disk['L1_transverse_band_test'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'L1_transverse_band_test'})
        if isinstance(config.get('L1_transverse_band_test'), dict):
            merged['L1_transverse_band_test'].update(config['L1_transverse_band_test'])
    return merged


def central_pair(n_side: int) -> tuple[int, int]:
    return (n_side // 2 - 1, n_side // 2)


def periodic_delta(value: np.ndarray, center: float) -> np.ndarray:
    delta = value - center
    return np.minimum(np.abs(delta), 1.0 - np.abs(delta))


def build_periodic_complex(n_side: int, epsilon: float, variant: str, flux_tube_phase: float = 0.0) -> ComplexData:
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
    elif variant in {'baseline', 'flux_tube'}:
        pass
    else:
        raise ValueError(f'Unsupported variant: {variant}')

    active_old_indices = [int(node_index[i, j, k]) for i in range(n_side) for j in range(n_side) for k in range(n_side) if active_mask[i, j, k]]
    old_to_new = {old: new for new, old in enumerate(active_old_indices)}
    points = raw_points[active_old_indices]

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
    ci, cj = n_side // 2, n_side // 2

    def add_face(face_type: str, boundary_keys: list[tuple[str, int, int, int, int]]) -> None:
        nonlocal face_idx
        oriented_edges: list[tuple[int, complex]] = []
        local_weights: list[float] = []
        phase_multiplier = 1.0 + 0.0j
        if variant == 'flux_tube' and face_type == 'xy':
            i0, j0, _ = boundary_keys[0][1:4]
            if i0 == ci and j0 == cj:
                phase_multiplier = np.exp(1j * flux_tube_phase)
        for boundary_pos, (axis, i, j, k, sign) in enumerate(boundary_keys):
            key = (axis, i % n_side, j % n_side, k % n_side)
            if key not in edge_map:
                return
            edge_id = edge_map[key]
            coeff = complex(sign)
            if boundary_pos == 0:
                coeff *= phase_multiplier
            oriented_edges.append((edge_id, coeff))
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
                add_face('xy', [
                    ('x', i, j, k, +1),
                    ('y', (i + 1) % n_side, j, k, +1),
                    ('x', i, (j + 1) % n_side, k, -1),
                    ('y', i, j, k, -1),
                ])
                add_face('xz', [
                    ('x', i, j, k, +1),
                    ('z', (i + 1) % n_side, j, k, +1),
                    ('x', i, j, (k + 1) % n_side, -1),
                    ('z', i, j, k, -1),
                ])
                add_face('yz', [
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


class ExactProjector:
    def __init__(self, d0: sp.csr_matrix, node_laplacian: sp.csr_matrix):
        self.d0 = d0.tocsr()
        reduced = node_laplacian.tocsc()[1:, 1:]
        self.solve = spla.factorized(reduced)
        self.node_dim = node_laplacian.shape[0]

    def apply(self, vec: np.ndarray) -> np.ndarray:
        rhs = np.asarray(self.d0.conj().T @ vec, dtype=complex)
        coeffs = np.zeros(self.node_dim, dtype=complex)
        coeffs[1:] = self.solve(rhs[1:])
        return np.asarray(self.d0 @ coeffs, dtype=complex)


class HarmonicProjector:
    def __init__(self, basis: np.ndarray):
        self.basis = np.asarray(basis, dtype=complex)

    def apply(self, vec: np.ndarray) -> np.ndarray:
        if self.basis.size == 0:
            return np.zeros_like(vec)
        return self.basis @ (self.basis.conj().T @ vec)


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
    order = np.argsort(evals.real)
    return np.asarray(evals[order], dtype=float), np.asarray(evecs[:, order], dtype=complex)


def harmonic_basis_from_reference(reference: ComplexData, count: int, harmonic_tol: float, tol: float) -> np.ndarray:
    evals, evecs = low_eigensystem(reference.L1, k=count, tol=tol)
    mask = evals < harmonic_tol
    basis = np.asarray(evecs[:, mask], dtype=complex)
    if basis.size:
        basis, _ = np.linalg.qr(basis, mode='reduced')
    return basis


def build_penalized_transverse_operator(
    analysis: ComplexData,
    exact_projector: ExactProjector,
    harmonic_projector: HarmonicProjector,
    penalty: float,
) -> spla.LinearOperator:
    edge_dim = analysis.upper.shape[0]

    def matvec(vec: np.ndarray) -> np.ndarray:
        vec = np.asarray(vec, dtype=complex)
        return np.asarray(analysis.upper @ vec, dtype=complex) + penalty * exact_projector.apply(vec) + penalty * harmonic_projector.apply(vec)

    return spla.LinearOperator((edge_dim, edge_dim), matvec=matvec, dtype=complex)


def project_transverse(vec: np.ndarray, exact_projector: ExactProjector, harmonic_projector: HarmonicProjector) -> np.ndarray:
    projected = vec - exact_projector.apply(vec)
    projected = projected - harmonic_projector.apply(projected)
    return projected


def defect_distances(midpoints: np.ndarray, variant: str, n_side: int) -> np.ndarray:
    center = 0.5 - 0.5 / n_side
    if variant in {'line_defect', 'flux_tube'}:
        dx = periodic_delta(midpoints[:, 0], center)
        dy = periodic_delta(midpoints[:, 1], center)
        return np.sqrt(dx * dx + dy * dy)
    dx = periodic_delta(midpoints[:, 0], center)
    dy = periodic_delta(midpoints[:, 1], center)
    dz = periodic_delta(midpoints[:, 2], center)
    return np.sqrt(dx * dx + dy * dy + dz * dz)


def local_support_fraction(mode: np.ndarray, distances: np.ndarray, variant: str, n_side: int) -> float:
    h = 1.0 / n_side
    if variant == 'puncture':
        mask = distances <= 2.5 * h
    elif variant in {'line_defect', 'flux_tube'}:
        mask = distances <= 1.5 * h
    else:
        mask = distances <= 2.0 * h
    norm = float(np.sum(np.abs(mode) ** 2)) or 1.0
    return float(np.sum(np.abs(mode[mask]) ** 2) / norm)


def radial_profile(mode: np.ndarray, distances: np.ndarray, bins: int = 20) -> dict[str, list[float]]:
    if len(distances) == 0:
        return {'radius': [], 'density': []}
    order = np.argsort(distances)
    d_sorted = distances[order]
    weights = np.abs(mode[order]) ** 2
    edges = np.linspace(float(np.min(d_sorted)), float(np.max(d_sorted)), bins + 1)
    radii: list[float] = []
    density: list[float] = []
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (d_sorted >= left) & (d_sorted < right if right < edges[-1] else d_sorted <= right)
        if not np.any(mask):
            continue
        radii.append(float(np.mean(d_sorted[mask])))
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
    restricted_modes: int,
    harmonic_tol: float,
    eig_tol: float,
    penalty: float,
    flux_tube_phase: float,
) -> dict[str, Any]:
    reference_variant = 'baseline' if variant == 'flux_tube' else variant
    reference = build_periodic_complex(n_side=n_side, epsilon=epsilon, variant=reference_variant, flux_tube_phase=0.0)
    analysis = build_periodic_complex(n_side=n_side, epsilon=epsilon, variant=variant, flux_tube_phase=flux_tube_phase)

    harmonic_basis = harmonic_basis_from_reference(reference, count=max(restricted_modes, 8), harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = ExactProjector(reference.d0, reference.node_laplacian)
    harmonic_projector = HarmonicProjector(harmonic_basis)

    transverse_operator = build_penalized_transverse_operator(analysis, exact_projector, harmonic_projector, penalty=penalty)
    k_eff = min(restricted_modes + 4, analysis.upper.shape[0] - 2) if analysis.upper.shape[0] > 2 else 1
    raw_evals, raw_vecs = spla.eigsh(transverse_operator, k=k_eff, which='SA', tol=eig_tol)
    raw_order = np.argsort(raw_evals.real)
    raw_evals = np.asarray(raw_evals[raw_order], dtype=float)
    raw_vecs = np.asarray(raw_vecs[:, raw_order], dtype=complex)

    restricted_spectrum: list[float] = []
    restricted_records: list[dict[str, Any]] = []
    restricted_vectors: list[np.ndarray] = []
    distances = defect_distances(analysis.midpoints, variant=variant, n_side=n_side)
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
                'divergence_norm': float(np.linalg.norm(reference.d0.conj().T @ vec)),
                'curl_norm': float(np.linalg.norm(analysis.d1 @ vec)) if analysis.d1.shape[0] else 0.0,
                'ipr': inverse_participation_ratio(vec),
                'near_defect_fraction': local_support_fraction(vec, distances, variant=variant, n_side=n_side),
                'radial_profile': radial_profile(vec, distances) if variant in {'puncture', 'line_defect', 'flux_tube'} else {'radius': [], 'density': []},
            }
        )
        if len(restricted_spectrum) >= restricted_modes:
            break

    return {
        'label': f'{variant}_n{n_side}',
        'config': {
            'variant': variant,
            'variant_label': VARIANT_LABELS[variant],
            'n_side': int(n_side),
            'nodes': int(len(analysis.points)),
            'edges': int(len(analysis.edges)),
            'faces': int(len(analysis.face_weights)),
            'removed_nodes': int(analysis.removed_nodes),
            'removed_requested_edges': int(analysis.removed_requested_edges),
            'epsilon': float(epsilon),
            'flux_tube_phase': float(flux_tube_phase if variant == 'flux_tube' else 0.0),
        },
        'dimensions': {
            'harmonic': int(harmonic_basis.shape[1]),
        },
        'restricted_transverse_spectrum': restricted_spectrum,
        'restricted_transverse_modes': restricted_records,
        'restricted_vectors': restricted_vectors,
        'midpoints': analysis.midpoints,
        'directions': analysis.directions,
    }


def plot_edge_mode(ax: Any, midpoints: np.ndarray, directions: np.ndarray, mode: np.ndarray, title: str) -> None:
    magnitudes = np.abs(mode)
    threshold = np.quantile(magnitudes, 0.85) if len(magnitudes) > 12 else 0.0
    select = magnitudes >= threshold
    if not np.any(select):
        select = np.ones_like(magnitudes, dtype=bool)
    max_mag = float(np.max(magnitudes[select])) or 1.0
    vectors = directions[select] * (np.real(mode[select]) / max_mag)[:, None] * 0.18
    phases = np.angle(mode[select])
    colors = cm.twilight((phases + math.pi) / (2.0 * math.pi))
    ax.scatter(
        midpoints[select, 0],
        midpoints[select, 1],
        midpoints[select, 2],
        c=phases,
        cmap='twilight',
        s=22 + 150 * magnitudes[select] / max_mag,
        vmin=-math.pi,
        vmax=math.pi,
    )
    ax.quiver(
        midpoints[select, 0],
        midpoints[select, 1],
        midpoints[select, 2],
        vectors[:, 0],
        vectors[:, 1],
        vectors[:, 2],
        colors=colors,
        linewidth=0.8,
        arrow_length_ratio=0.25,
    )
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')


def make_band_collapse_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    for ax, variant in zip(axes.ravel(), ['baseline', 'puncture', 'line_defect', 'flux_tube']):
        for n_side in sizes:
            case = cases[f'{variant}_n{n_side}']
            spectrum = np.asarray(case['restricted_transverse_spectrum'][:20], dtype=float)
            scaled = (n_side * n_side) * spectrum
            ax.plot(range(1, len(scaled) + 1), scaled, marker='o', label=f'n={n_side}')
        ax.set_title(VARIANT_LABELS[variant])
        ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel('restricted mode index k')
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$n^2 \lambda_k$')
    axes[0, 1].legend(fontsize=8, ncol=2)
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_floor_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> dict[str, dict[str, float]]:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    fits: dict[str, dict[str, float]] = {}
    fig, ax = plt.subplots(figsize=(7, 5))
    for variant in ['baseline', 'puncture', 'line_defect', 'flux_tube']:
        values = [cases[f'{variant}_n{n}']['restricted_transverse_spectrum'][0] for n in sizes]
        fit = fit_power_law(sizes, values)
        fits[variant] = fit
        ax.loglog(sizes, values, marker=VARIANT_MARKERS[variant], label=f"{VARIANT_LABELS[variant]} (p={fit['p']:.2f})")
        curve = [fit['a'] / (n ** fit['p']) for n in sizes]
        ax.loglog(sizes, curve, linestyle='--', alpha=0.6)
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
    for variant in ['baseline', 'puncture', 'line_defect', 'flux_tube']:
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
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    largest_n = max(case['config']['n_side'] for case in cases.values())
    scatter = None
    for ax, variant in zip(axes.ravel(), ['baseline', 'puncture', 'line_defect', 'flux_tube']):
        case = cases[f'{variant}_n{largest_n}']
        records = case['restricted_transverse_modes']
        x = [record['divergence_norm'] for record in records]
        y = [record['curl_norm'] for record in records]
        c = [record['eigenvalue'] for record in records]
        scatter = ax.scatter(x, y, c=c, s=30, cmap='viridis', marker=VARIANT_MARKERS[variant], edgecolors='none')
        ax.set_title(f"{VARIANT_LABELS[variant]} (n={largest_n})")
        ax.grid(alpha=0.2)
    for ax in axes[-1]:
        ax.set_xlabel(r'$||d_0^* a_k||$')
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$||d_1 a_k||$')
    if scatter is not None:
        cbar = fig.colorbar(scatter, ax=axes.ravel().tolist(), shrink=0.85)
        cbar.set_label(r'$\lambda_k$')
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_mode_profile_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> None:
    fig = plt.figure(figsize=(16, 10))
    largest_n = max(case['config']['n_side'] for case in cases.values())
    for index, variant in enumerate(['baseline', 'puncture', 'line_defect', 'flux_tube'], start=1):
        case = cases[f'{variant}_n{largest_n}']
        ax = fig.add_subplot(2, 2, index, projection='3d')
        vector = np.asarray(case['restricted_vectors'][0], dtype=complex)
        plot_edge_mode(ax, case['midpoints'], case['directions'], vector, f"{VARIANT_LABELS[variant]} (n={largest_n})")
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def summarize_results(cases: dict[str, dict[str, Any]], fits: dict[str, dict[str, float]]) -> tuple[str, str, dict[str, Any]]:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    band_metrics: dict[str, dict[str, Any]] = {}
    for variant in ['baseline', 'puncture', 'line_defect', 'flux_tube']:
        scaled_rows = []
        for n_side in sizes:
            spectrum = np.asarray(cases[f'{variant}_n{n_side}']['restricted_transverse_spectrum'][:20], dtype=float)
            scaled_rows.append((n_side * n_side) * spectrum)
        min_len = min(len(row) for row in scaled_rows)
        stack = np.vstack([row[:min_len] for row in scaled_rows])
        rel_spread = np.std(stack, axis=0) / np.maximum(np.mean(stack, axis=0), 1e-12)
        band_metrics[variant] = {
            'collapse_relative_spread_first_10': float(np.mean(rel_spread[: min(10, len(rel_spread))])),
            'fit': fits[variant],
            'lowest_ipr': {str(n): float(cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['ipr']) for n in sizes},
            'lowest_near_defect_fraction': {str(n): float(cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['near_defect_fraction']) for n in sizes},
        }

    collapse_good = all(metrics['collapse_relative_spread_first_10'] < 0.25 for metrics in band_metrics.values())
    localization_decays = True
    for variant in ['baseline', 'puncture', 'line_defect', 'flux_tube']:
        ipr_values = [cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['ipr'] for n in sizes]
        if not all(ipr_values[idx + 1] < ipr_values[idx] for idx in range(len(ipr_values) - 1)):
            localization_decays = False
    defect_fraction_decays = True
    for variant in ['puncture', 'line_defect', 'flux_tube']:
        near_values = [cases[f'{variant}_n{n}']['restricted_transverse_modes'][0]['near_defect_fraction'] for n in sizes]
        if not all(near_values[idx + 1] <= near_values[idx] + 1e-10 for idx in range(len(near_values) - 1)):
            defect_fraction_decays = False

    if collapse_good and localization_decays and defect_fraction_decays:
        observation = 'the first twenty restricted transverse eigenvalues show a partial n^2 spectral collapse while the lowest-mode IPR and defect concentration both fall with lattice size'
        conclusion = 'the restricted transverse sector is organizing as a stable low spectral band consistent with continuum scaling in the tested range'
    else:
        observation = 'the restricted floor descends, but the higher restricted modes do not collapse cleanly enough across sizes to isolate a stable band yet'
        conclusion = 'the restricted transverse sector still looks dominated by finite-size structure rather than a clean continuum band in the tested range'
    return observation, conclusion, {'band_metrics': band_metrics, 'collapse_good': collapse_good, 'localization_decays': localization_decays, 'defect_fraction_decays': defect_fraction_decays}


def sanitize_case(case: dict[str, Any]) -> dict[str, Any]:
    return {
        'label': case['label'],
        'config': case['config'],
        'dimensions': case['dimensions'],
        'restricted_transverse_spectrum': case['restricted_transverse_spectrum'],
        'restricted_transverse_modes': case['restricted_transverse_modes'],
    }


def save_results(result: dict[str, Any], plot_paths: list[str]) -> tuple[Path, list[str], str]:
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    payload = json.dumps({**result, 'plots': plot_paths}, indent=2)
    stamped = RESULTS / f'{timestamp}_L1_transverse_band_scan.json'
    latest = RESULTS / 'L1_transverse_band_scan_latest.json'
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
    note_path = note_dir / 'L1_Transverse_Band_Test_v1.md'
    sizes = result['config']['sizes']

    def branch_table(variant: str) -> str:
        rows = []
        for n_side in sizes:
            case = result['cases'][f'{variant}_n{n_side}']
            mode = case['restricted_transverse_modes'][0]
            rows.append(
                f"| {n_side} | {case['restricted_transverse_spectrum'][0]:.6f} | {mode['ipr']:.6f} | {mode['near_defect_fraction']:.6f} | {mode['divergence_norm']:.3e} | {mode['curl_norm']:.6f} |"
            )
        return '\n'.join(rows)

    fit_lines = '\n'.join(
        f"- {VARIANT_LABELS[variant]}: `lambda_min ~ {result['fits'][variant]['a']:.4f} / n^{result['fits'][variant]['p']:.3f}` (`R^2 = {result['fits'][variant]['r2']:.4f}`)"
        for variant in ['baseline', 'puncture', 'line_defect', 'flux_tube']
    )

    note = f"""# L1 Transverse Band Test

## Purpose

Determine whether the restricted transverse sector forms a stable low spectral band consistent with continuum scaling, or only an isolated descending mode.

## Setup

- operator unchanged:
  - `L1 = d0 d0* + d1* d1`
  - `T = d1* d1` restricted to `ker(d0*) intersect (H1)^perp`
- branches:
  - baseline periodic torus
  - puncture defect
  - line defect
  - flux-tube defect
- sizes: `{sizes}`
- kernel parameter: `epsilon = {result['config']['epsilon']}`
- flux-tube phase: `{result['config']['flux_tube_phase']}`

## Lowest restricted eigenvalue fits

{fit_lines}

## Lowest restricted mode diagnostics

### Baseline periodic torus

| `n` | lowest `lambda` | IPR | near-defect fraction | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: | ---: |
{branch_table('baseline')}

### Single cubic puncture

| `n` | lowest `lambda` | IPR | near-defect fraction | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: | ---: |
{branch_table('puncture')}

### Line defect

| `n` | lowest `lambda` | IPR | near-defect fraction | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: | ---: |
{branch_table('line_defect')}

### Flux-tube defect

| `n` | lowest `lambda` | IPR | near-defect fraction | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: | ---: |
{branch_table('flux_tube')}

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
        handle.write('\n## L1 transverse band test\n')
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            f"- Config: epsilon={result['config']['epsilon']}, sizes={result['config']['sizes']}, variants={result['config']['variants']}, restricted_modes={result['config']['restricted_modes']}, flux_tube_phase={result['config']['flux_tube_phase']}\n"
        )
        handle.write(f"- Results: `{result_path.relative_to(REPO_ROOT)}`\n")
        handle.write(f"- Plots: {', '.join(f'`{path}`' for path in stamped_plots)}\n")
        handle.write(f"- Observation: {result['observation']}\n")
        handle.write(f"- Conclusion: {result['conclusion']}\n")


def run_L1_transverse_band_test(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path=config_path)
    epsilon = float(cfg.get('epsilon', 0.2))
    experiment_cfg = cfg['L1_transverse_band_test']
    sizes = [int(value) for value in experiment_cfg.get('sizes', [6, 8, 10, 12, 14, 16])]
    variants = [str(value) for value in experiment_cfg.get('variants', ['baseline', 'puncture', 'line_defect', 'flux_tube'])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 20))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))

    cases: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        for variant in variants:
            cases[f'{variant}_n{n_side}'] = analyze_case(
                n_side=n_side,
                epsilon=epsilon,
                variant=variant,
                restricted_modes=restricted_modes,
                harmonic_tol=harmonic_tol,
                eig_tol=eig_tol,
                penalty=penalty,
                flux_tube_phase=flux_tube_phase,
            )

    collapse_plot = PLOTS / 'transverse_band_collapse.png'
    floor_plot = PLOTS / 'transverse_floor_vs_n.png'
    ipr_plot = PLOTS / 'transverse_ipr_vs_n.png'
    phase_plot = PLOTS / 'divergence_curl_phase_band.png'
    profile_plot = PLOTS / 'transverse_mode_profiles.png'

    make_band_collapse_plot(cases, collapse_plot)
    fits = make_floor_plot(cases, floor_plot)
    make_ipr_plot(cases, ipr_plot)
    make_phase_plot(cases, phase_plot)
    make_mode_profile_plot(cases, profile_plot)
    plot_paths = [
        str(collapse_plot.relative_to(REPO_ROOT)),
        str(floor_plot.relative_to(REPO_ROOT)),
        str(ipr_plot.relative_to(REPO_ROOT)),
        str(phase_plot.relative_to(REPO_ROOT)),
        str(profile_plot.relative_to(REPO_ROOT)),
    ]

    observation, conclusion, verdict = summarize_results(cases, fits)
    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'flux_tube_phase': flux_tube_phase,
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
    result, result_path, _, _, note_path = run_L1_transverse_band_test()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
