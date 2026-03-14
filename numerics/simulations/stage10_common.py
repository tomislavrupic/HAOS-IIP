#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

import numpy as np
import scipy.sparse as sp

from stage8_common import (
    REPO_ROOT,
    RESULTS,
    PLOTS,
    ExactProjector,
    HarmonicProjector,
    append_log,
    build_periodic_complex,
    ensure_matplotlib,
    harmonic_basis_from_L1,
    project_transverse,
)
from stage9_common import estimate_lambda_max

plt = ensure_matplotlib()

ATLAS_NOTES = REPO_ROOT / 'experiments' / 'pre_geometry_atlas'
ATLAS_NOTES.mkdir(exist_ok=True)

DEFAULT_STAGE10: dict[str, Any] = {
    'n_side': 16,
    'epsilon': 0.2,
    'harmonic_tol': 1.0e-8,
    'packet_center': [0.25, 0.5, 0.5],
    'kick_axis': 0,
    'kick_strength_scalar': 1.0,
    'kick_strength_transverse': 1.0,
    't_final': 0.9,
    'dt_safety': 0.28,
    'slice_axis': 2,
    'random_seed': 10,
}


@dataclass(frozen=True)
class Stage10Complex:
    points: np.ndarray
    node_index: np.ndarray
    edges: list[tuple[int, int]]
    edge_axes: list[str]
    edge_map: dict[tuple[str, int, int, int], int]
    midpoints: np.ndarray
    edge_weights: np.ndarray
    face_weights: np.ndarray
    d0: sp.csr_matrix
    d1: sp.csr_matrix
    L0: sp.csr_matrix
    L1: sp.csr_matrix
    upper: sp.csr_matrix
    n_side: int
    boundary_type: str


@dataclass(frozen=True)
class AtlasRun:
    run_id: str
    atlas_phase: str
    graph_type: str
    kernel_type: str
    operator_sector: str
    boundary_type: str
    initial_seed_type: str
    central_k: float
    bandwidth: float
    amplitude: float
    phase_pattern: str
    packet_count: int
    random_seed: int
    notes: str


def load_stage10_defaults(config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_STAGE10))
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        stage10 = on_disk.get('stage10', {})
        if isinstance(stage10, dict):
            merged.update(stage10)
        for key in ('epsilon', 'harmonic_tol'):
            if key in on_disk:
                merged[key] = on_disk[key]
    return merged


def timestamp_slug() -> str:
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def load_run_sheet(path: Path) -> list[AtlasRun]:
    payload = json.loads(path.read_text())
    return [AtlasRun(**item) for item in payload['runs']]


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_open_complex(n_side: int, epsilon: float) -> Stage10Complex:
    if n_side < 2:
        raise ValueError('open boundary requires n_side >= 2')

    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))
    scale = n_side - 1
    points = np.array(
        [[i / scale, j / scale, k / scale] for i in range(n_side) for j in range(n_side) for k in range(n_side)],
        dtype=float,
    )
    h = 1.0 / scale
    edge_weight = math.exp(-(h * h) / (2.0 * epsilon))

    edges: list[tuple[int, int]] = []
    edge_axes: list[str] = []
    midpoints: list[np.ndarray] = []
    edge_weights: list[float] = []
    edge_map: dict[tuple[str, int, int, int], int] = {}
    b_rows: list[int] = []
    b_cols: list[int] = []
    b_data: list[float] = []

    def add_edge(axis: str, i: int, j: int, k: int, ni: int, nj: int, nk: int) -> None:
        u = int(node_index[i, j, k])
        v = int(node_index[ni, nj, nk])
        edge_idx = len(edges)
        edge_map[(axis, i, j, k)] = edge_idx
        edges.append((u, v))
        edge_axes.append(axis)
        midpoints.append(
            np.array(
                [
                    (i + 0.5) / scale if axis == 'x' else i / scale,
                    (j + 0.5) / scale if axis == 'y' else j / scale,
                    (k + 0.5) / scale if axis == 'z' else k / scale,
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
                if i + 1 < n_side:
                    add_edge('x', i, j, k, i + 1, j, k)
                if j + 1 < n_side:
                    add_edge('y', i, j, k, i, j + 1, k)
                if k + 1 < n_side:
                    add_edge('z', i, j, k, i, j, k + 1)

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
            edge_id = edge_map[(axis, i, j, k)]
            c_rows.append(face_idx)
            c_cols.append(edge_id)
            c_data.append(float(sign))
            local_weights.append(edge_weights[edge_id])
        face_weights.append(float(np.mean(local_weights)))
        face_idx += 1

    for i in range(n_side - 1):
        for j in range(n_side - 1):
            for k in range(n_side):
                add_face([
                    ('x', i, j, k, +1),
                    ('y', i + 1, j, k, +1),
                    ('x', i, j + 1, k, -1),
                    ('y', i, j, k, -1),
                ])
    for i in range(n_side - 1):
        for j in range(n_side):
            for k in range(n_side - 1):
                add_face([
                    ('x', i, j, k, +1),
                    ('z', i + 1, j, k, +1),
                    ('x', i, j, k + 1, -1),
                    ('z', i, j, k, -1),
                ])
    for i in range(n_side):
        for j in range(n_side - 1):
            for k in range(n_side - 1):
                add_face([
                    ('y', i, j, k, +1),
                    ('z', i, j + 1, k, +1),
                    ('y', i, j, k + 1, -1),
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

    return Stage10Complex(
        points=points,
        node_index=node_index,
        edges=edges,
        edge_axes=edge_axes,
        edge_map=edge_map,
        midpoints=np.asarray(midpoints, dtype=float),
        edge_weights=edge_weights_arr,
        face_weights=face_weights_arr,
        d0=d0,
        d1=d1,
        L0=L0,
        L1=L1,
        upper=upper,
        n_side=n_side,
        boundary_type='open',
    )


def build_regular_complex(n_side: int, epsilon: float, boundary_type: str) -> Stage10Complex:
    if boundary_type == 'periodic':
        data = build_periodic_complex(n_side=n_side, epsilon=epsilon)
        return Stage10Complex(
            points=data.points,
            node_index=data.node_index,
            edges=data.edges,
            edge_axes=data.edge_axes,
            edge_map=data.edge_map,
            midpoints=data.midpoints,
            edge_weights=data.edge_weights,
            face_weights=data.face_weights,
            d0=data.d0,
            d1=data.d1,
            L0=data.L0,
            L1=data.L1,
            upper=data.upper,
            n_side=data.n_side,
            boundary_type='periodic',
        )
    if boundary_type == 'open':
        return build_open_complex(n_side=n_side, epsilon=epsilon)
    raise ValueError(f'unsupported boundary_type: {boundary_type}')


def displacement(points: np.ndarray, anchor: np.ndarray, boundary_type: str) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    anchor = np.asarray(anchor, dtype=float)
    delta = points - anchor
    if boundary_type == 'periodic':
        return (delta + 0.5) % 1.0 - 0.5
    return delta


def gaussian_packet(
    points: np.ndarray,
    center: np.ndarray,
    sigma: float,
    boundary_type: str,
    central_k: float,
    phase_pattern: str,
    kick_axis: int,
) -> tuple[np.ndarray, np.ndarray]:
    delta = displacement(points, center, boundary_type)
    r2 = np.sum(delta * delta, axis=1)
    envelope = np.exp(-r2 / (2.0 * sigma * sigma))
    phase_coord = delta[:, kick_axis]
    if phase_pattern == 'cosine-carrier':
        carrier = np.cos(2.0 * math.pi * central_k * phase_coord)
    elif phase_pattern == 'sine-carrier':
        carrier = np.sin(2.0 * math.pi * central_k * phase_coord)
    else:
        carrier = np.ones_like(phase_coord)
    return envelope * carrier, delta


def weighted_center(points: np.ndarray, weights: np.ndarray, boundary_type: str) -> np.ndarray:
    weights = np.asarray(weights, dtype=float)
    total = float(np.sum(weights))
    if total <= 0.0:
        return np.zeros(3, dtype=float)
    if boundary_type == 'periodic':
        center = np.zeros(3, dtype=float)
        for axis in range(3):
            angles = 2.0 * math.pi * np.asarray(points[:, axis], dtype=float)
            s = float(np.sum(weights * np.sin(angles)))
            c = float(np.sum(weights * np.cos(angles)))
            center[axis] = ((math.atan2(s, c) / (2.0 * math.pi)) + 1.0) % 1.0
        return center
    return np.sum(points * weights[:, None], axis=0) / total


def packet_width(points: np.ndarray, weights: np.ndarray, center: np.ndarray, boundary_type: str) -> float:
    total = float(np.sum(weights))
    if total <= 0.0:
        return 0.0
    delta = displacement(points, center, boundary_type)
    r2 = np.sum(delta * delta, axis=1)
    return float(np.sqrt(np.sum(weights * r2) / total))


def anisotropy_ratio(points: np.ndarray, weights: np.ndarray, center: np.ndarray, boundary_type: str) -> float:
    total = float(np.sum(weights))
    if total <= 0.0:
        return 0.0
    delta = displacement(points, center, boundary_type)
    cov = (delta.T * weights) @ delta / total
    evals = np.linalg.eigvalsh(cov)
    smallest = float(max(np.min(evals), 1.0e-12))
    largest = float(np.max(evals))
    return largest / smallest


def state_energy(operator: sp.csr_matrix, q: np.ndarray, v: np.ndarray) -> float:
    kinetic = float(np.vdot(v, v).real)
    potential = float(np.vdot(q, operator @ q).real)
    return 0.5 * (kinetic + potential)


def spectral_moments(operator: sp.csr_matrix, q: np.ndarray) -> tuple[float, float]:
    norm_sq = float(np.vdot(q, q).real)
    if norm_sq <= 0.0:
        return 0.0, 0.0
    lq = np.asarray(operator @ q, dtype=float)
    centroid = float(np.vdot(q, lq).real / norm_sq)
    second = float(np.vdot(lq, lq).real / norm_sq)
    spread_sq = max(second - centroid * centroid, 0.0)
    return centroid, math.sqrt(spread_sq)


def coherence_score(weights: np.ndarray) -> float:
    total = float(np.sum(weights))
    if total <= 0.0:
        return 0.0
    probs = np.asarray(weights, dtype=float) / total
    return float(np.sum(probs * probs))


def suggested_dt(operator: sp.csr_matrix, safety: float) -> float:
    lam_max = estimate_lambda_max(operator)
    return float(safety * 2.0 / math.sqrt(max(lam_max, 1.0e-12)))


def center_shift(centers: np.ndarray, boundary_type: str) -> float:
    if len(centers) < 2:
        return 0.0
    delta = displacement(centers[-1][None, :], centers[0], boundary_type)[0]
    return float(np.linalg.norm(delta))


def path_length(centers: np.ndarray, boundary_type: str) -> float:
    if len(centers) < 2:
        return 0.0
    total = 0.0
    for prev, cur in zip(centers[:-1], centers[1:]):
        total += float(np.linalg.norm(displacement(cur[None, :], prev, boundary_type)[0]))
    return total


def build_transverse_setup(n_side: int, epsilon: float, harmonic_tol: float, boundary_type: str) -> tuple[Stage10Complex, Callable[[np.ndarray], np.ndarray]]:
    data = build_regular_complex(n_side=n_side, epsilon=epsilon, boundary_type=boundary_type)
    harmonic_basis = harmonic_basis_from_L1(data.L1, count=8, harmonic_tol=harmonic_tol, tol=1.0e-8)
    exact_projector = ExactProjector(data.d0, data.L0)
    harmonic_projector = HarmonicProjector(harmonic_basis)

    def projector(vec: np.ndarray) -> np.ndarray:
        return project_transverse(np.asarray(vec, dtype=float), exact_projector, harmonic_projector)

    return data, projector


def build_scalar_initial_state(run: AtlasRun, points: np.ndarray, center: np.ndarray, kick_axis: int, boundary_type: str) -> tuple[np.ndarray, np.ndarray]:
    profile, delta = gaussian_packet(
        points=points,
        center=center,
        sigma=run.bandwidth,
        boundary_type=boundary_type,
        central_k=run.central_k,
        phase_pattern=run.phase_pattern,
        kick_axis=kick_axis,
    )
    q0 = profile / max(np.linalg.norm(profile), 1.0e-12)
    base_v = (delta[:, kick_axis] / max(run.bandwidth * run.bandwidth, 1.0e-12)) * profile
    v0 = base_v / max(np.linalg.norm(base_v), 1.0e-12)
    return run.amplitude * q0, run.amplitude * v0


def build_transverse_initial_state(
    run: AtlasRun,
    midpoints: np.ndarray,
    edge_axes: list[str],
    center: np.ndarray,
    kick_axis: int,
    boundary_type: str,
    projector: Callable[[np.ndarray], np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    profile, delta = gaussian_packet(
        points=midpoints,
        center=center,
        sigma=run.bandwidth,
        boundary_type=boundary_type,
        central_k=run.central_k,
        phase_pattern=run.phase_pattern,
        kick_axis=kick_axis,
    )
    axis_mask = np.asarray([1.0 if axis == 'y' else 0.0 for axis in edge_axes], dtype=float)
    raw = profile * axis_mask
    q0 = projector(raw)
    q0 = q0 / max(np.linalg.norm(q0), 1.0e-12)

    raw_v = (delta[:, kick_axis] / max(run.bandwidth * run.bandwidth, 1.0e-12)) * raw
    v0 = projector(raw_v)
    v0 = v0 / max(np.linalg.norm(v0), 1.0e-12)
    return run.amplitude * q0, run.amplitude * v0


def simulate_atlas_run(
    operator: sp.csr_matrix,
    q0: np.ndarray,
    v0: np.ndarray,
    positions: np.ndarray,
    boundary_type: str,
    dt: float,
    steps: int,
    project: Callable[[np.ndarray], np.ndarray] | None = None,
    divergence_op: sp.csr_matrix | None = None,
    leakage_fn: Callable[[np.ndarray], float] | None = None,
) -> dict[str, Any]:
    q = np.asarray(q0, dtype=float).copy()
    v = np.asarray(v0, dtype=float).copy()
    if project is not None:
        q = np.asarray(project(q), dtype=float)
        v = np.asarray(project(v), dtype=float)

    times: list[float] = []
    centers: list[list[float]] = []
    widths: list[float] = []
    norms: list[float] = []
    energies: list[float] = []
    anisotropies: list[float] = []
    spectral_centroids: list[float] = []
    spectral_spreads: list[float] = []
    sector_leakages: list[float] = []
    constraint_norms: list[float] = []
    coherences: list[float] = []

    def record(t: float) -> None:
        weights = np.abs(q) ** 2
        center = weighted_center(positions, weights, boundary_type)
        centroid, spread = spectral_moments(operator, q)
        times.append(float(t))
        centers.append(center.tolist())
        widths.append(packet_width(positions, weights, center, boundary_type))
        norms.append(float(np.linalg.norm(q)))
        energies.append(state_energy(operator, q, v))
        anisotropies.append(anisotropy_ratio(positions, weights, center, boundary_type))
        spectral_centroids.append(centroid)
        spectral_spreads.append(spread)
        sector_leakages.append(float(leakage_fn(q)) if leakage_fn is not None else 0.0)
        constraint_norms.append(float(np.linalg.norm(divergence_op @ q)) if divergence_op is not None else 0.0)
        coherences.append(coherence_score(weights))

    record(0.0)
    for step in range(steps):
        acc = -np.asarray(operator @ q, dtype=float)
        if project is not None:
            acc = np.asarray(project(acc), dtype=float)
        v_half = v + 0.5 * dt * acc
        q = q + dt * v_half
        if project is not None:
            q = np.asarray(project(q), dtype=float)
        acc_new = -np.asarray(operator @ q, dtype=float)
        if project is not None:
            acc_new = np.asarray(project(acc_new), dtype=float)
        v = v_half + 0.5 * dt * acc_new
        if project is not None:
            v = np.asarray(project(v), dtype=float)
        record((step + 1) * dt)

    centers_arr = np.asarray(centers, dtype=float)
    widths_arr = np.asarray(widths, dtype=float)
    energies_arr = np.asarray(energies, dtype=float)
    path = path_length(centers_arr, boundary_type)
    shift = center_shift(centers_arr, boundary_type)
    return {
        'times': times,
        'centers': centers,
        'widths': widths,
        'norms': norms,
        'energies': energies,
        'anisotropies': anisotropies,
        'spectral_centroids': spectral_centroids,
        'spectral_spreads': spectral_spreads,
        'sector_leakages': sector_leakages,
        'constraint_norms': constraint_norms,
        'coherences': coherences,
        'summary': {
            'center_shift': shift,
            'path_length': path,
            'path_efficiency': float(shift / max(path, 1.0e-12)),
            'initial_width': float(widths_arr[0]),
            'final_width': float(widths_arr[-1]),
            'width_change': float(widths_arr[-1] - widths_arr[0]),
            'width_ratio': float(widths_arr[-1] / max(widths_arr[0], 1.0e-12)),
            'relative_norm_change': float((norms[-1] - norms[0]) / max(norms[0], 1.0e-12)),
            'relative_energy_drift': float((energies_arr[-1] - energies_arr[0]) / max(abs(energies_arr[0]), 1.0e-12)),
            'max_anisotropy': float(np.max(anisotropies)),
            'max_constraint_norm': float(np.max(constraint_norms)),
            'max_sector_leakage': float(np.max(sector_leakages)),
            'coherence_initial': float(coherences[0]),
            'coherence_final': float(coherences[-1]),
            'coherence_ratio': float(coherences[-1] / max(coherences[0], 1.0e-12)),
            'spectral_centroid_initial': float(spectral_centroids[0]),
            'spectral_centroid_final': float(spectral_centroids[-1]),
            'spectral_spread_initial': float(spectral_spreads[0]),
            'spectral_spread_final': float(spectral_spreads[-1]),
        },
    }


def create_field_snapshot(path: Path, run_id: str, positions: np.ndarray, values: np.ndarray, boundary_type: str, slice_axis: int) -> None:
    weights = np.abs(values) ** 2
    coords = np.asarray(positions, dtype=float)
    slice_value = 0.5 if boundary_type == 'periodic' else 0.5
    distances = np.abs(coords[:, slice_axis] - slice_value)
    threshold = np.min(distances) + 1.0e-9
    mask = distances <= threshold
    plane_axes = [axis for axis in range(3) if axis != slice_axis]
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    scatter = ax.scatter(coords[mask, plane_axes[0]], coords[mask, plane_axes[1]], c=weights[mask], s=20, cmap='viridis')
    ax.set_xlabel(f'axis {plane_axes[0]}')
    ax.set_ylabel(f'axis {plane_axes[1]}')
    ax.set_title(f'Field snapshot: {run_id}')
    fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_spectral_trace(path: Path, run_id: str, times: list[float], centroids: list[float], spreads: list[float]) -> None:
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    ax.plot(times, centroids, label='centroid')
    ax.plot(times, spreads, label='spread')
    ax.set_xlabel('time')
    ax.set_ylabel('spectral moment')
    ax.set_title(f'Spectral morphology: {run_id}')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_trace_plot(path: Path, run_id: str, times: list[float], values: list[float], ylabel: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    ax.plot(times, values)
    ax.set_xlabel('time')
    ax.set_ylabel(ylabel)
    ax.set_title(f'{title}: {run_id}')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_summary_scatter(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    colors = {
        'Ballistic coherent': 'tab:blue',
        'Ballistic dispersive': 'tab:green',
        'Diffusive': 'tab:orange',
        'Localized': 'tab:purple',
        'Oscillatory trapped': 'tab:brown',
        'Fragmenting': 'tab:red',
        'Interfering coherent': 'tab:pink',
        'Synchronizing': 'tab:cyan',
        'Chaotic or irregular': 'tab:gray',
        'Metastable structured': 'tab:olive',
    }
    for row in rows:
        ax.scatter(
            float(row['center_shift']),
            float(row['width_ratio']),
            color=colors.get(str(row['regime_label']), 'black'),
            label=str(row['regime_label']),
            alpha=0.8,
        )
    handles, labels = ax.get_legend_handles_labels()
    unique: dict[str, Any] = {}
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    ax.legend(unique.values(), unique.keys(), fontsize=7)
    ax.set_xlabel('center shift')
    ax.set_ylabel('width ratio')
    ax.set_title('Atlas-0 morphology map')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def save_atlas_payload(
    experiment_slug: str,
    result: dict[str, Any],
    csv_rows: list[dict[str, Any]],
    csv_fieldnames: list[str],
    plot_paths: list[Path],
) -> tuple[Path, Path, list[str], str]:
    timestamp = timestamp_slug()
    json_path = RESULTS / f'{timestamp}_{experiment_slug}.json'
    csv_path = RESULTS / f'{timestamp}_{experiment_slug}.csv'
    json_path.write_text(json.dumps(result, indent=2), encoding='utf-8')
    write_csv(csv_path, csv_rows, csv_fieldnames)
    stamped_plots: list[str] = []
    for src in plot_paths:
        dst = PLOTS / f'{timestamp}_{src.name}'
        dst.write_bytes(src.read_bytes())
        stamped_plots.append(str(dst.relative_to(REPO_ROOT)))
    return json_path, csv_path, stamped_plots, timestamp


__all__ = [
    'ATLAS_NOTES',
    'AtlasRun',
    'PLOTS',
    'REPO_ROOT',
    'RESULTS',
    'append_log',
    'build_regular_complex',
    'build_scalar_initial_state',
    'build_transverse_initial_state',
    'build_transverse_setup',
    'create_field_snapshot',
    'create_spectral_trace',
    'create_summary_scatter',
    'create_trace_plot',
    'load_run_sheet',
    'load_stage10_defaults',
    'plt',
    'save_atlas_payload',
    'simulate_atlas_run',
    'suggested_dt',
]
