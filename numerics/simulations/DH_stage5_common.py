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
    'DH_stage5': {
        'two_d_n': 6,
        'three_d_n': 8,
        'flux_quanta': [0, 1, 2, 3],
        'spectrum_modes': 50,
        'square_root_modes': 48,
        'spectral_flow_modes': 12,
        'eig_tol': 1e-9,
    },
}

PAULI = [
    np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex),
    np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
]


@dataclass(frozen=True)
class DiracSystem:
    dim: int
    n_side: int
    epsilon: float
    flux_quanta: int
    points: np.ndarray
    diff_ops: list[sp.csr_matrix]
    L0: sp.csr_matrix
    Q: sp.csr_matrix
    DH: sp.csr_matrix
    DH2: sp.csr_matrix
    spin_dim: int
    total_dim: int
    plaquette_angle: float


@dataclass(frozen=True)
class SpectrumSummary:
    evals: np.ndarray
    evecs: np.ndarray
    ipr: np.ndarray


def ensure_matplotlib():
    return plt


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['DH_stage5'] = dict(DEFAULT_CONFIG['DH_stage5'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'DH_stage5'})
        if isinstance(on_disk.get('DH_stage5'), dict):
            merged['DH_stage5'].update(on_disk['DH_stage5'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'DH_stage5'})
        if isinstance(config.get('DH_stage5'), dict):
            merged['DH_stage5'].update(config['DH_stage5'])
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


def inverse_participation_ratio(vec: np.ndarray) -> float:
    abs_sq = np.abs(vec) ** 2
    norm_sq = float(np.sum(abs_sq)) or 1.0
    return float(np.sum(abs_sq * abs_sq) / (norm_sq * norm_sq))


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


def repeated_scalar_spectrum(l0_evals: np.ndarray, multiplicity: int) -> np.ndarray:
    return np.repeat(np.asarray(l0_evals, dtype=float), multiplicity)


def sparse_frobenius_norm(matrix: sp.csr_matrix) -> float:
    matrix = matrix.tocsr()
    return float(np.sqrt(np.sum(np.abs(matrix.data) ** 2)))


def link_phase(dim: int, axis: int, coord: tuple[int, ...], n_side: int, flux_quanta: int) -> complex:
    if dim < 2 or flux_quanta == 0:
        return 1.0 + 0.0j
    i = coord[0]
    j = coord[1]
    phi = 2.0 * math.pi * flux_quanta / max(n_side * n_side, 1)
    if axis == 0:
        return np.exp(-1j * phi * n_side * j) if i == n_side - 1 else 1.0 + 0.0j
    if axis == 1:
        return np.exp(1j * phi * i)
    return 1.0 + 0.0j


def regular_points(dim: int, n_side: int) -> np.ndarray:
    if dim == 2:
        return np.array([[i / n_side, j / n_side] for i in range(n_side) for j in range(n_side)], dtype=float)
    if dim == 3:
        return np.array([[i / n_side, j / n_side, k / n_side] for i in range(n_side) for j in range(n_side) for k in range(n_side)], dtype=float)
    raise ValueError(f'Unsupported dimension: {dim}')


def build_covariant_difference_ops(dim: int, n_side: int, epsilon: float, flux_quanta: int = 0) -> tuple[np.ndarray, list[sp.csr_matrix]]:
    points = regular_points(dim, n_side)
    node_count = points.shape[0]
    h = 1.0 / n_side
    edge_weight = math.exp(-(h * h) / (2.0 * epsilon))
    rootw = math.sqrt(edge_weight)

    if dim == 2:
        node_index = np.arange(n_side * n_side).reshape((n_side, n_side))
    else:
        node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))

    diff_ops: list[sp.csr_matrix] = []
    for axis in range(dim):
        rows: list[int] = []
        cols: list[int] = []
        data: list[complex] = []
        if dim == 2:
            for i in range(n_side):
                for j in range(n_side):
                    coord = (i, j)
                    u = int(node_index[i, j])
                    nxt = [i, j]
                    nxt[axis] = (nxt[axis] + 1) % n_side
                    v = int(node_index[tuple(nxt)])
                    phase = link_phase(dim, axis, coord, n_side, flux_quanta)
                    rows.extend([u, u])
                    cols.extend([u, v])
                    data.extend([-rootw + 0.0j, rootw * phase])
        else:
            for i in range(n_side):
                for j in range(n_side):
                    for k in range(n_side):
                        coord = (i, j, k)
                        u = int(node_index[i, j, k])
                        nxt = [i, j, k]
                        nxt[axis] = (nxt[axis] + 1) % n_side
                        v = int(node_index[tuple(nxt)])
                        phase = link_phase(dim, axis, coord, n_side, flux_quanta)
                        rows.extend([u, u])
                        cols.extend([u, v])
                        data.extend([-rootw + 0.0j, rootw * phase])
        diff = sp.coo_matrix((data, (rows, cols)), shape=(node_count, node_count), dtype=complex).tocsr()
        diff_ops.append(diff)
    return points, diff_ops


def build_dirac_system(dim: int, n_side: int, epsilon: float, flux_quanta: int = 0) -> DiracSystem:
    points, diff_ops = build_covariant_difference_ops(dim=dim, n_side=n_side, epsilon=epsilon, flux_quanta=flux_quanta)
    L0 = sum(diff.getH() @ diff for diff in diff_ops).tocsr()
    spin_dim = 2
    q_operator = None
    for sigma, diff in zip(PAULI[:dim], diff_ops):
        term = sp.kron(sp.csr_matrix(sigma), diff, format='csr')
        q_operator = term if q_operator is None else (q_operator + term)
    assert q_operator is not None
    zero = sp.csr_matrix(q_operator.shape, dtype=complex)
    DH = sp.bmat([[zero, q_operator.getH()], [q_operator, zero]], format='csr')
    DH2 = (DH.getH() @ DH).tocsr()
    plaquette_angle = 2.0 * math.pi * flux_quanta / max(n_side * n_side, 1) if dim >= 2 else 0.0
    return DiracSystem(
        dim=dim,
        n_side=n_side,
        epsilon=epsilon,
        flux_quanta=flux_quanta,
        points=points,
        diff_ops=diff_ops,
        L0=L0,
        Q=q_operator.tocsr(),
        DH=DH,
        DH2=DH2,
        spin_dim=spin_dim,
        total_dim=DH.shape[0],
        plaquette_angle=plaquette_angle,
    )


def compare_square_root(system: DiracSystem, modes: int, tol: float) -> dict[str, Any]:
    l0_evals, _ = low_eigs(system.L0, k=max(4, min(modes, system.L0.shape[0] - 2)), which='SA', tol=tol)
    qhq = (system.Q.getH() @ system.Q).tocsr()
    q_evals, _ = low_eigs(qhq, k=max(4, min(2 * len(l0_evals), qhq.shape[0] - 2)), which='SA', tol=tol)
    target = repeated_scalar_spectrum(l0_evals, multiplicity=2)[: len(q_evals)]
    valid = np.maximum(np.abs(target), 1e-14)
    relative_error = np.abs(q_evals - target) / valid
    correlation = float(np.corrcoef(q_evals, target)[0, 1]) if len(q_evals) > 1 else 1.0
    spin_scalar = sp.kron(sp.eye(system.spin_dim, format='csr', dtype=complex), system.L0, format='csr')
    diff = (qhq - spin_scalar).tocsr()
    diff_norm = sparse_frobenius_norm(diff)
    ref_norm = sparse_frobenius_norm(spin_scalar) or 1.0
    return {
        'l0_eigenvalues': l0_evals.tolist(),
        'qhq_eigenvalues': q_evals.tolist(),
        'target_repeated_l0': target.tolist(),
        'mean_relative_error': float(np.mean(relative_error)),
        'max_relative_error': float(np.max(relative_error)),
        'eigenvalue_correlation': correlation,
        'operator_relative_frobenius_error': float(diff_norm / ref_norm),
    }


def low_dirac_spectrum(system: DiracSystem, modes: int, tol: float) -> SpectrumSummary:
    evals, evecs = low_eigs(system.DH, k=modes, sigma=0.0, tol=tol)
    ipr = np.asarray([inverse_participation_ratio(vec) for vec in evecs.T], dtype=float)
    return SpectrumSummary(evals=evals, evecs=evecs, ipr=ipr)


def pairing_error(evals: np.ndarray) -> float:
    vals = np.asarray(evals, dtype=float)
    order = np.argsort(vals)
    vals = vals[order]
    pairs = min(np.sum(vals < 0.0), np.sum(vals > 0.0))
    if pairs == 0:
        return 0.0
    neg = vals[:pairs]
    pos = vals[-pairs:][::-1]
    return float(np.max(np.abs(neg + pos)))


def degeneracy_summary(evals: np.ndarray, tol: float = 1e-8) -> list[dict[str, Any]]:
    abs_vals = np.sort(np.abs(np.asarray(evals, dtype=float)))
    groups: list[list[float]] = []
    for value in abs_vals:
        if not groups or abs(value - groups[-1][0]) > tol:
            groups.append([value])
        else:
            groups[-1].append(value)
    return [
        {'abs_eigenvalue': float(group[0]), 'degeneracy': int(len(group))}
        for group in groups[:12]
    ]


def plot_square_root_comparison(system_records: dict[str, dict[str, Any]], path: Path) -> None:
    fig, axes = plt.subplots(1, len(system_records), figsize=(5 * len(system_records), 4), squeeze=False)
    for ax, (label, record) in zip(axes[0], system_records.items()):
        qvals = np.asarray(record['qhq_eigenvalues'], dtype=float)
        target = np.asarray(record['target_repeated_l0'], dtype=float)
        count = min(len(qvals), len(target))
        ax.plot(range(count), target[:count], marker='o', label='repeated $L_0$')
        ax.plot(range(count), qvals[:count], marker='s', label='$Q^\dagger Q$')
        ax.set_title(label)
        ax.set_xlabel('mode index')
        ax.set_ylabel('eigenvalue')
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_dirac_spectrum(summary_records: dict[str, SpectrumSummary], path_spectrum: Path, path_ipr: Path) -> None:
    fig, axes = plt.subplots(1, len(summary_records), figsize=(5 * len(summary_records), 4), squeeze=False)
    for ax, (label, summary) in zip(axes[0], summary_records.items()):
        ax.axhline(0.0, color='black', linewidth=0.8, alpha=0.5)
        ax.scatter(range(len(summary.evals)), summary.evals, s=18)
        ax.set_title(label)
        ax.set_xlabel('mode index')
        ax.set_ylabel(r'$\lambda$')
        ax.grid(alpha=0.25)
    fig.savefig(path_spectrum, dpi=180, bbox_inches='tight')
    plt.close(fig)

    fig, axes = plt.subplots(1, len(summary_records), figsize=(5 * len(summary_records), 4), squeeze=False)
    for ax, (label, summary) in zip(axes[0], summary_records.items()):
        ax.scatter(range(len(summary.ipr)), summary.ipr, s=18)
        ax.set_title(label)
        ax.set_xlabel('mode index')
        ax.set_ylabel('IPR')
        ax.grid(alpha=0.25)
    fig.savefig(path_ipr, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_flux_flow(flux_angles: list[float], positive_branches: list[np.ndarray], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    mode_count = len(positive_branches[0]) if positive_branches else 0
    for mode_index in range(mode_count):
        ax.plot(flux_angles, [branch[mode_index] for branch in positive_branches], marker='o', label=f'mode {mode_index + 1}' if mode_index < 6 else None)
    ax.set_xlabel('plaquette angle')
    ax.set_ylabel('positive near-zero $D_H$ eigenvalue')
    ax.set_title('DH flux spectral flow')
    ax.grid(alpha=0.25)
    if mode_count > 0:
        ax.legend(fontsize=8, ncol=2)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)
