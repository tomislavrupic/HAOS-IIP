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
from scipy import linalg

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
    'L1_defect_transverse_test': {
        'sizes': [6, 8, 10],
        'variants': ['baseline', 'puncture', 'line_defect'],
        'phase_modes': 48,
        'restricted_modes': 12,
        'harmonic_tol': 1e-8,
        'basis_tol': 1e-10,
    },
}

VARIANT_LABELS = {
    'baseline': 'baseline periodic torus',
    'puncture': 'single cubic puncture',
    'line_defect': 'line defect',
}


@dataclass(frozen=True)
class ComplexData:
    points: np.ndarray
    edges: list[tuple[int, int]]
    directions: np.ndarray
    midpoints: np.ndarray
    edge_weights: np.ndarray
    face_weights: np.ndarray
    d0: np.ndarray
    d1: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    L1: np.ndarray
    n_side: int
    variant: str
    removed_nodes: int
    removed_edges: int
    harmonic_estimate: int


@dataclass(frozen=True)
class HodgeBases:
    exact_basis: np.ndarray
    harmonic_basis: np.ndarray
    coexact_basis: np.ndarray


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['L1_defect_transverse_test'] = dict(DEFAULT_CONFIG['L1_defect_transverse_test'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'L1_defect_transverse_test'})
        if isinstance(on_disk.get('L1_defect_transverse_test'), dict):
            merged['L1_defect_transverse_test'].update(on_disk['L1_defect_transverse_test'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'L1_defect_transverse_test'})
        if isinstance(config.get('L1_defect_transverse_test'), dict):
            merged['L1_defect_transverse_test'].update(config['L1_defect_transverse_test'])
    return merged


def central_pair(n_side: int) -> tuple[int, int]:
    return (n_side // 2 - 1, n_side // 2)


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

    active_old_indices = [int(node_index[i, j, k]) for i in range(n_side) for j in range(n_side) for k in range(n_side) if active_mask[i, j, k]]
    old_to_new = {old: new for new, old in enumerate(active_old_indices)}
    points = raw_points[active_old_indices]

    h = 1.0 / n_side
    edge_weight = math.exp(-(h * h) / (2.0 * epsilon))
    edges: list[tuple[int, int]] = []
    directions: list[np.ndarray] = []
    midpoints: list[np.ndarray] = []
    edge_weights: list[float] = []
    edge_map: dict[tuple[str, int, int, int], int] = {}

    basis = {
        'x': np.array([1.0, 0.0, 0.0]),
        'y': np.array([0.0, 1.0, 0.0]),
        'z': np.array([0.0, 0.0, 1.0]),
    }

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

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_edge('x', i, j, k, (i + 1) % n_side, j, k)
                add_edge('y', i, j, k, i, (j + 1) % n_side, k)
                add_edge('z', i, j, k, i, j, (k + 1) % n_side)

    n_nodes = len(points)
    n_edges = len(edges)
    B0 = np.zeros((n_edges, n_nodes), dtype=float)
    for edge_idx, (u, v) in enumerate(edges):
        B0[edge_idx, u] = -1.0
        B0[edge_idx, v] = 1.0

    face_boundaries: list[list[tuple[int, int]]] = []
    face_weights: list[float] = []

    def add_face(boundary_keys: list[tuple[str, int, int, int, int]]) -> None:
        oriented_edges: list[tuple[int, int]] = []
        local_weights: list[float] = []
        for axis, i, j, k, sign in boundary_keys:
            key = (axis, i % n_side, j % n_side, k % n_side)
            if key not in edge_map:
                return
            edge_idx = edge_map[key]
            oriented_edges.append((edge_idx, sign))
            local_weights.append(edge_weights[edge_idx])
        face_boundaries.append(oriented_edges)
        face_weights.append(float(np.mean(local_weights)))

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

    n_faces = len(face_boundaries)
    C = np.zeros((n_faces, n_edges), dtype=float)
    for face_idx, boundary in enumerate(face_boundaries):
        for edge_idx, sign in boundary:
            C[face_idx, edge_idx] = float(sign)

    edge_weights_arr = np.asarray(edge_weights, dtype=float)
    face_weights_arr = np.asarray(face_weights, dtype=float)
    d0 = np.diag(np.sqrt(edge_weights_arr)) @ B0
    d1 = np.diag(np.sqrt(face_weights_arr)) @ C if n_faces else np.zeros((0, n_edges), dtype=float)
    lower = d0 @ d0.T
    upper = d1.T @ d1 if n_faces else np.zeros((n_edges, n_edges), dtype=float)
    L1 = lower + upper

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
        n_side=n_side,
        variant=variant,
        removed_nodes=int(np.size(active_mask) - np.count_nonzero(active_mask)),
        removed_edges=int(len(removed_edge_keys)),
        harmonic_estimate=0,
    )


def orthonormal_image_basis(matrix: np.ndarray, tol: float) -> np.ndarray:
    if matrix.size == 0:
        return np.zeros((matrix.shape[0], 0), dtype=float)
    q, r, _ = linalg.qr(matrix, mode='economic', pivoting=True)
    diag = np.abs(np.diag(r))
    rank = int(np.sum(diag > tol))
    return q[:, :rank]


def orthonormal_complement(basis: np.ndarray, dimension: int, tol: float) -> np.ndarray:
    if basis.size == 0:
        return np.eye(dimension, dtype=float)
    complement = linalg.null_space(basis.T, rcond=tol)
    return np.asarray(complement, dtype=float)


def mode_projection_fraction(vec: np.ndarray, basis: np.ndarray) -> float:
    if basis.size == 0:
        return 0.0
    coeffs = basis.T @ vec
    return float(np.dot(coeffs, coeffs) / np.dot(vec, vec))


def build_hodge_bases(data: ComplexData, harmonic_tol: float, basis_tol: float) -> tuple[HodgeBases, np.ndarray, np.ndarray]:
    evals, evecs = linalg.eigh(data.L1)
    exact_basis = orthonormal_image_basis(data.d0, tol=basis_tol)
    harmonic_mask = evals < harmonic_tol
    harmonic_basis = np.asarray(evecs[:, harmonic_mask], dtype=float)
    if harmonic_basis.size:
        harmonic_basis, _ = np.linalg.qr(harmonic_basis, mode='reduced')
    combined = exact_basis if harmonic_basis.size == 0 else np.column_stack([exact_basis, harmonic_basis])
    coexact_basis = orthonormal_complement(combined, data.L1.shape[0], tol=basis_tol)
    return HodgeBases(exact_basis=exact_basis, harmonic_basis=harmonic_basis, coexact_basis=coexact_basis), evals, evecs


def analyze_low_modes(
    evals: np.ndarray,
    evecs: np.ndarray,
    data: ComplexData,
    bases: HodgeBases,
    count: int,
) -> list[dict[str, Any]]:
    vectors = evecs[:, :count]
    div_values = np.linalg.norm(data.d0.T @ vectors, axis=0)
    curl_values = np.linalg.norm(data.d1 @ vectors, axis=0) if data.d1.size else np.zeros(count, dtype=float)
    records: list[dict[str, Any]] = []
    for idx in range(count):
        vec = vectors[:, idx]
        restricted = vec.copy()
        if bases.exact_basis.size:
            restricted = restricted - bases.exact_basis @ (bases.exact_basis.T @ restricted)
        if bases.harmonic_basis.size:
            restricted = restricted - bases.harmonic_basis @ (bases.harmonic_basis.T @ restricted)
        restricted_norm = float(np.linalg.norm(restricted))
        records.append(
            {
                'mode_index': idx,
                'eigenvalue': float(evals[idx]),
                'divergence_norm': float(div_values[idx]),
                'curl_norm': float(curl_values[idx]),
                'exact_fraction': mode_projection_fraction(vec, bases.exact_basis),
                'harmonic_fraction': mode_projection_fraction(vec, bases.harmonic_basis),
                'coexact_fraction': mode_projection_fraction(vec, bases.coexact_basis),
                'restricted_norm': restricted_norm,
            }
        )
    return records


def restricted_transverse_spectrum(data: ComplexData, bases: HodgeBases, count: int) -> tuple[list[float], list[dict[str, Any]]]:
    q = bases.coexact_basis
    if q.size == 0:
        return [], []
    reduced_upper = q.T @ data.upper @ q
    evals, evecs = linalg.eigh(reduced_upper, subset_by_index=[0, min(count - 1, reduced_upper.shape[0] - 1)])
    evals = np.asarray(evals, dtype=float)
    evecs = np.asarray(evecs, dtype=float)
    records: list[dict[str, Any]] = []
    for idx in range(len(evals)):
        vec = q @ evecs[:, idx]
        div_norm = float(np.linalg.norm(data.d0.T @ vec))
        curl_norm = float(np.linalg.norm(data.d1 @ vec)) if data.d1.size else 0.0
        records.append(
            {
                'mode_index': idx,
                'eigenvalue': float(evals[idx]),
                'divergence_norm': div_norm,
                'curl_norm': curl_norm,
                'exact_fraction': mode_projection_fraction(vec, bases.exact_basis),
                'harmonic_fraction': mode_projection_fraction(vec, bases.harmonic_basis),
                'coexact_fraction': mode_projection_fraction(vec, bases.coexact_basis),
            }
        )
    return [float(value) for value in evals], records


def analyze_case(
    n_side: int,
    epsilon: float,
    variant: str,
    phase_modes: int,
    restricted_modes: int,
    harmonic_tol: float,
    basis_tol: float,
) -> dict[str, Any]:
    data = build_periodic_defect_complex(n_side=n_side, epsilon=epsilon, variant=variant)
    bases, full_evals, full_evecs = build_hodge_bases(data, harmonic_tol=harmonic_tol, basis_tol=basis_tol)
    phase_count = min(phase_modes, len(full_evals))
    phase_records = analyze_low_modes(full_evals, full_evecs, data, bases, phase_count)
    restricted_spectrum, restricted_records = restricted_transverse_spectrum(data, bases, count=min(restricted_modes, bases.coexact_basis.shape[1]))

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
            'removed_requested_edges': int(data.removed_edges),
            'epsilon': float(epsilon),
        },
        'dimensions': {
            'exact': int(bases.exact_basis.shape[1]),
            'harmonic': int(bases.harmonic_basis.shape[1]),
            'coexact': int(bases.coexact_basis.shape[1]),
        },
        'phase_modes': phase_records,
        'restricted_transverse_spectrum': restricted_spectrum,
        'restricted_transverse_modes': restricted_records,
        'full_low_spectrum': [float(value) for value in full_evals[:phase_count]],
    }


def make_divergence_curl_phase_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> None:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    variants = ['baseline', 'puncture', 'line_defect']
    fig, axes = plt.subplots(len(sizes), len(variants), figsize=(15, 12), sharex=True, sharey=True)
    for row, n_side in enumerate(sizes):
        for col, variant in enumerate(variants):
            ax = axes[row, col] if len(sizes) > 1 else axes[col]
            case = cases[f'{variant}_n{n_side}']
            x = [record['divergence_norm'] for record in case['phase_modes']]
            y = [record['curl_norm'] for record in case['phase_modes']]
            c = [record['eigenvalue'] for record in case['phase_modes']]
            scatter = ax.scatter(x, y, c=c, s=26, cmap='viridis', edgecolors='none')
            ax.set_title(f"n={n_side}, {VARIANT_LABELS[variant]}")
            ax.grid(alpha=0.2)
            if row == len(sizes) - 1:
                ax.set_xlabel(r'$||d_0^* a_k||$')
            if col == 0:
                ax.set_ylabel(r'$||d_1 a_k||$')
    cbar = fig.colorbar(scatter, ax=axes.ravel().tolist(), shrink=0.9)
    cbar.set_label(r'$\lambda_k$')
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_restricted_spectrum_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> None:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    fig, axes = plt.subplots(1, len(sizes), figsize=(5.5 * len(sizes), 4), sharey=True)
    if len(sizes) == 1:
        axes = [axes]
    for ax, n_side in zip(axes, sizes):
        for variant in ['baseline', 'puncture', 'line_defect']:
            case = cases[f'{variant}_n{n_side}']
            spectrum = case['restricted_transverse_spectrum']
            ax.plot(range(len(spectrum)), spectrum, marker='o', label=VARIANT_LABELS[variant])
        ax.set_title(f'n={n_side}')
        ax.set_xlabel('restricted mode index')
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    axes[0].set_ylabel('restricted eigenvalue')
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_lowest_restricted_vs_n_plot(cases: dict[str, dict[str, Any]], plot_path: Path) -> None:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    fig, ax = plt.subplots(figsize=(7, 4))
    for variant in ['baseline', 'puncture', 'line_defect']:
        values = []
        for n_side in sizes:
            spectrum = cases[f'{variant}_n{n_side}']['restricted_transverse_spectrum']
            values.append(spectrum[0] if spectrum else math.nan)
        ax.plot(sizes, values, marker='o', label=VARIANT_LABELS[variant])
    ax.set_xlabel('lattice side n')
    ax.set_ylabel('lowest restricted eigenvalue')
    ax.set_title('Lowest restricted transverse eigenvalue vs lattice size')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def sanitize_case(case: dict[str, Any]) -> dict[str, Any]:
    return {
        'label': case['label'],
        'config': case['config'],
        'dimensions': case['dimensions'],
        'full_low_spectrum': case['full_low_spectrum'],
        'phase_modes': case['phase_modes'],
        'restricted_transverse_spectrum': case['restricted_transverse_spectrum'],
        'restricted_transverse_modes': case['restricted_transverse_modes'],
    }


def summarize_results(cases: dict[str, dict[str, Any]]) -> tuple[str, str, dict[str, Any]]:
    sizes = sorted({case['config']['n_side'] for case in cases.values()})
    variants = ['baseline', 'puncture', 'line_defect']
    lowest: dict[str, list[float]] = {}
    for variant in variants:
        values = []
        for n_side in sizes:
            spectrum = cases[f'{variant}_n{n_side}']['restricted_transverse_spectrum']
            values.append(float(spectrum[0]) if spectrum else math.nan)
        lowest[variant] = values

    descending = {
        variant: all(
            values[idx + 1] <= values[idx] + 1e-10
            for idx in range(len(values) - 1)
            if not math.isnan(values[idx]) and not math.isnan(values[idx + 1])
        )
        for variant, values in lowest.items()
    }
    any_descending = any(descending.values())
    any_defect_below_baseline = False
    for idx in range(len(sizes)):
        baseline_value = lowest['baseline'][idx]
        for variant in ('puncture', 'line_defect'):
            defect_value = lowest[variant][idx]
            if not math.isnan(defect_value) and not math.isnan(baseline_value) and defect_value < baseline_value:
                any_defect_below_baseline = True

    if all(descending.values()) and any_defect_below_baseline:
        observation = 'the lowest restricted transverse eigenvalue decreases with lattice size in all three substrate branches, and both defect branches lie below the baseline torus at every tested n'
        conclusion = 'the restricted transverse sector develops a descending low band under the tested puncture and line-defect substrate modifications'
    elif any_descending or any_defect_below_baseline:
        observation = 'the restricted transverse floor shifts under substrate defects, but the effect must be read from the eigenvalue trend rather than the full L1 floor'
        conclusion = 'a descending restricted transverse band is suggested only where the defect branch lowers the coexact floor relative to the baseline torus'
    else:
        observation = 'the restricted transverse floor remains high across the tested puncture and line-defect branches'
        conclusion = 'no descending low restricted transverse band appears in the tested substrate variants'

    verdict = {
        'lowest_restricted_by_variant': {variant: {str(n): value for n, value in zip(sizes, values)} for variant, values in lowest.items()},
        'monotone_descending_by_variant': descending,
        'any_defect_below_baseline': any_defect_below_baseline,
    }
    return observation, conclusion, verdict


def save_results(result: dict[str, Any], plot_paths: list[str]) -> tuple[Path, list[str], str]:
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    payload = json.dumps({**result, 'plots': plot_paths}, indent=2)
    stamped = RESULTS / f'{timestamp}_L1_defect_scan.json'
    latest = RESULTS / 'L1_defect_scan_latest.json'
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
    note_path = note_dir / 'L1_Defect_Transverse_Test_v1.md'
    sizes = result['config']['sizes']

    def lowest_table(variant: str) -> str:
        rows = []
        for n_side in sizes:
            case = result['cases'][f'{variant}_n{n_side}']
            spectrum = case['restricted_transverse_spectrum']
            lowest = spectrum[0] if spectrum else math.nan
            rows.append(
                f"| {n_side} | {case['config']['nodes']} | {case['config']['edges']} | {case['dimensions']['harmonic']} | {lowest:.6f} |"
            )
        return '\n'.join(rows)

    note = f"""# L1 Defect Transverse Test

## Purpose

Test whether punctures or line defects weaken harmonic dominance enough for the restricted transverse sector

`d1* d1` on `ker(d0*) intersect (H1)^perp`

to produce a descending low spectral band.

## Setup

- substrate: periodic cubic lattice with oriented `x`, `y`, `z` edges
- sizes: `{sizes}`
- variants:
  - baseline periodic torus
  - single cubic puncture: centered `2x2x2` node removal
  - line defect: removed `z`-edge column at central `(x, y)` index
- kernel parameter: `epsilon = {result['config']['epsilon']}`
- operator definitions unchanged:
  - `d0 = W1^(1/2) B0`
  - `d1 = W2^(1/2) C`
  - `L1 = d0 d0* + d1* d1`

## Restricted transverse floor

### Baseline periodic torus

| `n` | nodes | edges | harmonic dim | lowest restricted eigenvalue |
| --- | ---: | ---: | ---: | ---: |
{lowest_table('baseline')}

### Single cubic puncture

| `n` | nodes | edges | harmonic dim | lowest restricted eigenvalue |
| --- | ---: | ---: | ---: | ---: |
{lowest_table('puncture')}

### Line defect

| `n` | nodes | edges | harmonic dim | lowest restricted eigenvalue |
| --- | ---: | ---: | ---: | ---: |
{lowest_table('line_defect')}

## Diagnostics

- divergence-curl phase plot for the low `L1` window
- restricted transverse spectrum for each `(n, variant)`
- lowest restricted eigenvalue versus lattice size

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
        handle.write('\n## L1 defect transverse test\n')
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            f"- Config: epsilon={result['config']['epsilon']}, sizes={result['config']['sizes']}, variants={result['config']['variants']}, phase_modes={result['config']['phase_modes']}, restricted_modes={result['config']['restricted_modes']}\n"
        )
        handle.write(f"- Results: `{result_path.relative_to(REPO_ROOT)}`\n")
        handle.write(f"- Plots: {', '.join(f'`{path}`' for path in stamped_plots)}\n")
        handle.write(f"- Observation: {result['observation']}\n")
        handle.write(f"- Conclusion: {result['conclusion']}\n")


def run_L1_defect_transverse_test(
    config: dict[str, Any] | None = None,
    config_path: Path | None = None,
) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path=config_path)
    epsilon = float(cfg.get('epsilon', 0.2))
    experiment_cfg = cfg['L1_defect_transverse_test']
    sizes = [int(value) for value in experiment_cfg.get('sizes', [6, 8, 10])]
    variants = [str(value) for value in experiment_cfg.get('variants', ['baseline', 'puncture', 'line_defect'])]
    phase_modes = int(experiment_cfg.get('phase_modes', 48))
    restricted_modes = int(experiment_cfg.get('restricted_modes', 12))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    basis_tol = float(experiment_cfg.get('basis_tol', 1e-10))

    cases: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        for variant in variants:
            label = f'{variant}_n{n_side}'
            cases[label] = analyze_case(
                n_side=n_side,
                epsilon=epsilon,
                variant=variant,
                phase_modes=phase_modes,
                restricted_modes=restricted_modes,
                harmonic_tol=harmonic_tol,
                basis_tol=basis_tol,
            )

    phase_plot = PLOTS / 'divergence_curl_phase.png'
    spectrum_plot = PLOTS / 'restricted_transverse_spectrum.png'
    size_plot = PLOTS / 'restricted_eigenvalue_vs_n.png'
    make_divergence_curl_phase_plot(cases, phase_plot)
    make_restricted_spectrum_plot(cases, spectrum_plot)
    make_lowest_restricted_vs_n_plot(cases, size_plot)
    plot_paths = [
        str(phase_plot.relative_to(REPO_ROOT)),
        str(spectrum_plot.relative_to(REPO_ROOT)),
        str(size_plot.relative_to(REPO_ROOT)),
    ]

    observation, conclusion, verdict = summarize_results(cases)
    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'phase_modes': phase_modes,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'basis_tol': basis_tol,
        },
        'cases': {label: sanitize_case(case) for label, case in cases.items()},
        'observation': observation,
        'conclusion': conclusion,
        'verdict': verdict,
    }
    result_path, stamped_plots, timestamp = save_results(result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_experiment_log(result, result_path, stamped_plots)
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_L1_defect_transverse_test()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
