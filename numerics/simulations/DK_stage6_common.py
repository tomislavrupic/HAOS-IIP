#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
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
    'DK_stage6': {
        'cochain_sizes': [8, 12, 16],
        'spectrum_n': 12,
        'flux_n': 12,
        'cycle_phases': [0.0, 0.2617993877991494, 0.5235987755982988, 0.7853981633974483],
        'square_modes': 40,
        'spectrum_modes': 60,
        'eig_tol': 1e-9,
    },
}


@dataclass(frozen=True)
class DK2DComplex:
    n_side: int
    epsilon: float
    cycle_phase_x: float
    cycle_phase_y: float
    points: np.ndarray
    edge_midpoints: np.ndarray
    face_centers: np.ndarray
    edge_axes: np.ndarray
    edge_weights: np.ndarray
    face_weights: np.ndarray
    d0: sp.csr_matrix
    d1: sp.csr_matrix
    delta1: sp.csr_matrix
    delta2: sp.csr_matrix
    delta_h: sp.csr_matrix
    dirac_kahler: sp.csr_matrix
    block_sizes: tuple[int, int, int]


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['DK_stage6'] = dict(DEFAULT_CONFIG['DK_stage6'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'DK_stage6'})
        if isinstance(on_disk.get('DK_stage6'), dict):
            merged['DK_stage6'].update(on_disk['DK_stage6'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'DK_stage6'})
        if isinstance(config.get('DK_stage6'), dict):
            merged['DK_stage6'].update(config['DK_stage6'])
    return merged


def timestamp_slug() -> str:
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding='utf-8')


def append_log(title: str, config_summary: str, result_path: Path, plot_paths: list[Path], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open('a', encoding='utf-8') as handle:
        handle.write(f'\n## {title}\n')
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(f'- Config: {config_summary}\n')
        handle.write(f"- Results: `{result_path.relative_to(REPO_ROOT)}`\n")
        handle.write(f"- Plots: {', '.join(f'`{path.relative_to(REPO_ROOT)}`' for path in plot_paths)}\n")
        handle.write(f'- Observation: {observation}\n')
        handle.write(f'- Conclusion: {conclusion}\n')


def sparse_frobenius_norm(matrix: sp.spmatrix) -> float:
    matrix = matrix.tocsr()
    return float(np.sqrt(np.sum(np.abs(matrix.data) ** 2)))


def low_eigs(operator: sp.csr_matrix, k: int, which: str = 'SA', sigma: float | None = None, tol: float = 1e-9) -> tuple[np.ndarray, np.ndarray]:
    dim = operator.shape[0]
    k_eff = max(1, min(k, dim - 2)) if dim > 2 else 1
    kwargs: dict[str, Any] = {'k': k_eff, 'tol': tol}
    if sigma is None:
        kwargs['which'] = which
    else:
        kwargs['sigma'] = sigma
        kwargs['which'] = 'LM'
    evals, evecs = spla.eigsh(operator, **kwargs)
    order = np.argsort(evals.real)
    return np.asarray(evals[order], dtype=float), np.asarray(evecs[:, order], dtype=complex)


def inverse_participation_ratio(vec: np.ndarray) -> float:
    abs_sq = np.abs(vec) ** 2
    norm_sq = float(np.sum(abs_sq)) or 1.0
    return float(np.sum(abs_sq * abs_sq) / (norm_sq * norm_sq))


def pairing_error(evals: np.ndarray) -> float:
    vals = np.sort(np.asarray(evals, dtype=float))
    negatives = np.abs(vals[vals < -1e-12])[::-1]
    positives = vals[vals > 1e-12]
    if negatives.size == 0 or positives.size == 0:
        return 0.0
    count = min(len(negatives), len(positives))
    ref = np.maximum(np.abs(positives[:count]), 1e-12)
    return float(np.mean(np.abs(negatives[:count] - positives[:count]) / ref))


def degeneracy_summary(evals: np.ndarray, tol: float = 1e-8) -> list[dict[str, Any]]:
    vals = np.sort(np.abs(np.asarray(evals, dtype=float)))
    groups: list[list[float]] = []
    for val in vals:
        if not groups or abs(val - groups[-1][-1]) > tol:
            groups.append([float(val)])
        else:
            groups[-1].append(float(val))
    return [{'abs_eigenvalue': float(np.mean(group)), 'degeneracy': len(group)} for group in groups[:16]]


def sector_weights(vec: np.ndarray, block_sizes: tuple[int, int, int]) -> dict[str, float]:
    n0, n1, n2 = block_sizes
    split0 = n0
    split1 = n0 + n1
    parts = [vec[:split0], vec[split0:split1], vec[split1:split1 + n2]]
    norms = [float(np.sum(np.abs(part) ** 2)) for part in parts]
    total = sum(norms) or 1.0
    return {
        'omega0_fraction': norms[0] / total,
        'omega1_fraction': norms[1] / total,
        'omega2_fraction': norms[2] / total,
    }


def edge_phase_x(i: int, n_side: int, cycle_phase_x: float) -> complex:
    if cycle_phase_x == 0.0:
        return 1.0 + 0.0j
    return np.exp(1j * cycle_phase_x) if i == n_side - 1 else 1.0 + 0.0j


def edge_phase_y(j: int, n_side: int, cycle_phase_y: float) -> complex:
    if cycle_phase_y == 0.0:
        return 1.0 + 0.0j
    return np.exp(1j * cycle_phase_y) if j == n_side - 1 else 1.0 + 0.0j


def build_dk2d_complex(n_side: int, epsilon: float, cycle_phase_x: float = 0.0, cycle_phase_y: float = 0.0) -> DK2DComplex:
    h = 1.0 / n_side
    edge_weight = math.exp(-(h * h) / (2.0 * epsilon))
    face_weight = edge_weight
    root_edge = math.sqrt(edge_weight)
    root_face = math.sqrt(face_weight)

    node_index = np.arange(n_side * n_side).reshape((n_side, n_side))
    points = np.array([[i / n_side, j / n_side] for i in range(n_side) for j in range(n_side)], dtype=float)

    edge_rows: list[int] = []
    edge_cols: list[int] = []
    edge_data: list[complex] = []
    edge_midpoints: list[np.ndarray] = []
    edge_axes: list[int] = []
    edge_map: dict[tuple[str, int, int], int] = {}

    edge_id = 0
    for i in range(n_side):
        for j in range(n_side):
            u = int(node_index[i, j])
            v = int(node_index[(i + 1) % n_side, j])
            phase = edge_phase_x(i, n_side, cycle_phase_x)
            edge_map[('x', i, j)] = edge_id
            edge_rows.extend([edge_id, edge_id])
            edge_cols.extend([u, v])
            edge_data.extend([-root_edge + 0.0j, root_edge * phase])
            edge_midpoints.append(np.array([((i + 0.5) / n_side) % 1.0, j / n_side], dtype=float))
            edge_axes.append(0)
            edge_id += 1
    for i in range(n_side):
        for j in range(n_side):
            u = int(node_index[i, j])
            v = int(node_index[i, (j + 1) % n_side])
            phase = edge_phase_y(j, n_side, cycle_phase_y)
            edge_map[('y', i, j)] = edge_id
            edge_rows.extend([edge_id, edge_id])
            edge_cols.extend([u, v])
            edge_data.extend([-root_edge + 0.0j, root_edge * phase])
            edge_midpoints.append(np.array([i / n_side, ((j + 0.5) / n_side) % 1.0], dtype=float))
            edge_axes.append(1)
            edge_id += 1

    n_nodes = n_side * n_side
    n_edges = edge_id
    d0 = sp.coo_matrix((edge_data, (edge_rows, edge_cols)), shape=(n_edges, n_nodes), dtype=complex).tocsr()

    face_rows: list[int] = []
    face_cols: list[int] = []
    face_data: list[complex] = []
    face_centers: list[np.ndarray] = []
    face_id = 0
    for i in range(n_side):
        for j in range(n_side):
            ex0 = edge_map[('x', i, j)]
            ey1 = edge_map[('y', (i + 1) % n_side, j)]
            ex1 = edge_map[('x', i, (j + 1) % n_side)]
            ey0 = edge_map[('y', i, j)]
            ux = edge_phase_x(i, n_side, cycle_phase_x)
            uy = edge_phase_y(j, n_side, cycle_phase_y)
            entries = [
                (ex0, +root_face),
                (ey1, +root_face * ux),
                (ex1, -root_face * uy),
                (ey0, -root_face),
            ]
            for edge_idx, value in entries:
                face_rows.append(face_id)
                face_cols.append(edge_idx)
                face_data.append(value)
            face_centers.append(np.array([((i + 0.5) / n_side) % 1.0, ((j + 0.5) / n_side) % 1.0], dtype=float))
            face_id += 1

    n_faces = face_id
    d1 = sp.coo_matrix((face_data, (face_rows, face_cols)), shape=(n_faces, n_edges), dtype=complex).tocsr()

    delta1 = d0.getH().tocsr()
    delta2 = d1.getH().tocsr()

    delta0_op = (delta1 @ d0).tocsr()
    delta1_op = (d0 @ delta1 + delta2 @ d1).tocsr()
    delta2_op = (d1 @ delta2).tocsr()
    delta_h = sp.block_diag((delta0_op, delta1_op, delta2_op), format='csr')

    zero00 = sp.csr_matrix((n_nodes, n_nodes), dtype=complex)
    zero01 = sp.csr_matrix((n_nodes, n_faces), dtype=complex)
    zero10 = sp.csr_matrix((n_edges, n_edges), dtype=complex)
    zero20 = sp.csr_matrix((n_faces, n_nodes), dtype=complex)
    zero22 = sp.csr_matrix((n_faces, n_faces), dtype=complex)
    dirac_kahler = sp.bmat(
        [
            [zero00, delta1, zero01],
            [d0, zero10, delta2],
            [zero20, d1, zero22],
        ],
        format='csr',
    )

    return DK2DComplex(
        n_side=n_side,
        epsilon=epsilon,
        cycle_phase_x=cycle_phase_x,
        cycle_phase_y=cycle_phase_y,
        points=points,
        edge_midpoints=np.asarray(edge_midpoints, dtype=float),
        face_centers=np.asarray(face_centers, dtype=float),
        edge_axes=np.asarray(edge_axes, dtype=int),
        edge_weights=np.full(n_edges, edge_weight, dtype=float),
        face_weights=np.full(n_faces, face_weight, dtype=float),
        d0=d0,
        d1=d1,
        delta1=delta1,
        delta2=delta2,
        delta_h=delta_h,
        dirac_kahler=dirac_kahler,
        block_sizes=(n_nodes, n_edges, n_faces),
    )


def cochain_identity_errors(complex_data: DK2DComplex) -> dict[str, float]:
    d1d0 = (complex_data.d1 @ complex_data.d0).tocsr()
    delta1delta2 = (complex_data.delta1 @ complex_data.delta2).tocsr()
    dd = (complex_data.dirac_kahler @ complex_data.dirac_kahler - complex_data.delta_h).tocsr()
    return {
        'd1_d0_error': sparse_frobenius_norm(d1d0),
        'delta1_delta2_error': sparse_frobenius_norm(delta1delta2),
        'square_residual': sparse_frobenius_norm(dd),
        'square_relative_residual': sparse_frobenius_norm(dd) / (sparse_frobenius_norm(complex_data.delta_h) or 1.0),
    }


def graded_low_spectrum(complex_data: DK2DComplex, modes: int, tol: float) -> dict[str, Any]:
    evals, evecs = low_eigs(complex_data.dirac_kahler, k=modes, sigma=0.0, tol=tol)
    ipr = np.array([inverse_participation_ratio(evecs[:, idx]) for idx in range(evecs.shape[1])], dtype=float)
    weights = [sector_weights(evecs[:, idx], complex_data.block_sizes) for idx in range(evecs.shape[1])]
    return {
        'eigenvalues': evals.tolist(),
        'ipr': ipr.tolist(),
        'pairing_error': pairing_error(evals),
        'lowest_abs_eigenvalue': float(np.min(np.abs(evals))),
        'mean_ipr': float(np.mean(ipr)),
        'max_ipr': float(np.max(ipr)),
        'degeneracy_summary': degeneracy_summary(evals),
        'sector_weights': weights,
    }


def plot_square_vs_hodge(records: dict[str, dict[str, Any]], path: Path) -> None:
    plt.figure(figsize=(6.8, 4.6))
    for label, record in records.items():
        hodge = np.asarray(record['hodge_eigenvalues'], dtype=float)
        square = np.asarray(record['dk_square_eigenvalues'], dtype=float)
        count = min(len(hodge), len(square))
        plt.plot(hodge[:count], square[:count], marker='o', linestyle='-', label=label)
    mins, maxs = plt.xlim(), plt.ylim()
    lo = min(mins[0], maxs[0])
    hi = max(mins[1], maxs[1])
    plt.plot([lo, hi], [lo, hi], linestyle='--', color='black', linewidth=1.0)
    plt.xlabel('Hodge eigenvalue')
    plt.ylabel('Dirac-Kaehler square eigenvalue')
    plt.title('Dirac-Kaehler square vs Hodge spectrum')
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def plot_dk_spectrum(summary: dict[str, Any], spectrum_path: Path, grading_path: Path, ipr_path: Path) -> None:
    evals = np.asarray(summary['eigenvalues'], dtype=float)
    ipr = np.asarray(summary['ipr'], dtype=float)
    weights = summary['sector_weights']
    mode_index = np.arange(len(evals))

    plt.figure(figsize=(6.6, 4.2))
    plt.scatter(mode_index, evals, c=np.abs(evals), cmap='viridis', s=28)
    plt.axhline(0.0, color='black', linestyle='--', linewidth=1.0)
    plt.xlabel('mode index')
    plt.ylabel('eigenvalue')
    plt.title('Low Dirac-Kaehler spectrum')
    plt.tight_layout()
    plt.savefig(spectrum_path, dpi=220)
    plt.close()

    plt.figure(figsize=(6.6, 4.2))
    omega0 = [record['omega0_fraction'] for record in weights]
    omega1 = [record['omega1_fraction'] for record in weights]
    omega2 = [record['omega2_fraction'] for record in weights]
    plt.stackplot(mode_index, omega0, omega1, omega2, labels=[r'$\Omega^0$', r'$\Omega^1$', r'$\Omega^2$'], alpha=0.85)
    plt.xlabel('mode index')
    plt.ylabel('mode weight fraction')
    plt.title('Graded weight distribution of low modes')
    plt.legend(loc='upper right', fontsize=8)
    plt.tight_layout()
    plt.savefig(grading_path, dpi=220)
    plt.close()

    plt.figure(figsize=(6.6, 4.2))
    plt.plot(mode_index, ipr, marker='o')
    plt.xlabel('mode index')
    plt.ylabel('IPR')
    plt.title('Dirac-Kaehler mode localization')
    plt.tight_layout()
    plt.savefig(ipr_path, dpi=220)
    plt.close()


def plot_flux_flow(cycle_phases: list[float], branches: list[np.ndarray], flow_path: Path, residuals: list[float], residual_path: Path) -> None:
    plt.figure(figsize=(6.6, 4.2))
    branch_arr = np.asarray(branches, dtype=float)
    for idx in range(branch_arr.shape[1]):
        plt.plot(cycle_phases, branch_arr[:, idx], marker='o', linewidth=1.2)
    plt.xlabel('cycle phase')
    plt.ylabel('positive eigenvalue')
    plt.title('Dirac-Kaehler low spectral flow')
    plt.tight_layout()
    plt.savefig(flow_path, dpi=220)
    plt.close()

    plt.figure(figsize=(6.4, 4.0))
    plt.plot(cycle_phases, residuals, marker='o')
    plt.xlabel('cycle phase')
    plt.ylabel('relative square residual')
    plt.title('Square identity residual under cycle holonomy')
    plt.tight_layout()
    plt.savefig(residual_path, dpi=220)
    plt.close()
