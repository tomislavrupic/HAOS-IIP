#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment

from L1_stage3_common import (
    PLOTS,
    REPO_ROOT,
    VARIANT_LABELS,
    VARIANT_MARKERS,
    append_log,
    save_result_payload,
    plt,
)
from L1_transverse_band_test import analyze_case

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_active_sector_transport': {
        'sizes': [12, 16, 20, 24, 28],
        'variants': ['baseline', 'puncture', 'line_defect', 'flux_tube', 'mild_disorder'],
        'restricted_modes': 6,
        'transport_modes': 4,
        'probe_n_side': 6,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': math.pi / 2.0,
        'disorder_strength': 0.12,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_active_sector_transport'] = dict(DEFAULT_CONFIG['transverse_active_sector_transport'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_active_sector_transport'})
        if isinstance(on_disk.get('transverse_active_sector_transport'), dict):
            merged['transverse_active_sector_transport'].update(on_disk['transverse_active_sector_transport'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_active_sector_transport'})
        if isinstance(config.get('transverse_active_sector_transport'), dict):
            merged['transverse_active_sector_transport'].update(config['transverse_active_sector_transport'])
    return merged


def probe_cell_ids(midpoints: np.ndarray, probe_n_side: int) -> np.ndarray:
    coords = np.floor(np.mod(np.asarray(midpoints, dtype=float), 1.0) * probe_n_side).astype(int)
    coords = np.clip(coords, 0, probe_n_side - 1)
    return coords[:, 0] + probe_n_side * (coords[:, 1] + probe_n_side * coords[:, 2])


def transport_mode_to_probe(midpoints: np.ndarray, directions: np.ndarray, mode: np.ndarray, probe_n_side: int) -> np.ndarray:
    cell_ids = probe_cell_ids(midpoints, probe_n_side)
    axis_ids = np.argmax(np.abs(np.asarray(directions, dtype=float)), axis=1)
    probe_dim = 3 * (probe_n_side**3)
    transported = np.zeros(probe_dim, dtype=complex)
    counts = np.zeros(probe_dim, dtype=float)
    for edge_idx, cell_id in enumerate(cell_ids):
        channel = 3 * int(cell_id) + int(axis_ids[edge_idx])
        transported[channel] += mode[edge_idx]
        counts[channel] += 1.0
    mask = counts > 0.0
    transported[mask] /= np.sqrt(counts[mask])
    norm = float(np.linalg.norm(transported))
    if norm <= 1.0e-12:
        return transported
    return transported / norm


def transported_modes(case: dict[str, Any], probe_n_side: int, transport_modes_count: int) -> list[np.ndarray]:
    transported: list[np.ndarray] = []
    vectors = list(case.get('restricted_vectors', []))[:transport_modes_count]
    for vec in vectors:
        moved = transport_mode_to_probe(case['midpoints'], case['directions'], np.asarray(vec, dtype=complex), probe_n_side)
        if float(np.linalg.norm(moved)) > 1.0e-12:
            transported.append(moved)
    return transported


def normalized_matrix(vectors: list[np.ndarray]) -> np.ndarray:
    if not vectors:
        return np.zeros((0, 0), dtype=complex)
    return np.column_stack(vectors)


def subspace_cosines(vectors_a: list[np.ndarray], vectors_b: list[np.ndarray]) -> list[float]:
    if not vectors_a or not vectors_b:
        return []
    matrix_a = normalized_matrix(vectors_a)
    matrix_b = normalized_matrix(vectors_b)
    qa, _ = np.linalg.qr(matrix_a, mode='reduced')
    qb, _ = np.linalg.qr(matrix_b, mode='reduced')
    singular_values = np.linalg.svd(qa.conj().T @ qb, compute_uv=False)
    return [float(value) for value in singular_values]


def pair_transport_metrics(
    variant: str,
    n_from: int,
    n_to: int,
    case_from: dict[str, Any],
    case_to: dict[str, Any],
    probe_n_side: int,
    transport_modes_count: int,
) -> dict[str, Any]:
    modes_from = transported_modes(case_from, probe_n_side, transport_modes_count)
    modes_to = transported_modes(case_to, probe_n_side, transport_modes_count)
    if not modes_from or not modes_to:
        return {
            'variant': variant,
            'n_from': int(n_from),
            'n_to': int(n_to),
            'probe_n_side': int(probe_n_side),
            'transport_modes': 0,
            'matched_pairs': [],
            'mean_overlap': math.nan,
            'min_overlap': math.nan,
            'leading_mode_overlap': math.nan,
            'principal_cosines': [],
            'mean_principal_cosine': math.nan,
            'min_principal_cosine': math.nan,
            'mean_scaled_eigen_drift': math.nan,
            'max_scaled_eigen_drift': math.nan,
        }

    matrix_from = normalized_matrix(modes_from)
    matrix_to = normalized_matrix(modes_to)
    overlap = np.abs(matrix_from.conj().T @ matrix_to)
    row_ind, col_ind = linear_sum_assignment(-overlap)
    order = np.argsort(row_ind)
    row_ind = row_ind[order]
    col_ind = col_ind[order]

    scaled_from = [(n_from**2) * float(value) for value in case_from['restricted_transverse_spectrum'][: len(modes_from)]]
    scaled_to = [(n_to**2) * float(value) for value in case_to['restricted_transverse_spectrum'][: len(modes_to)]]
    matched_pairs: list[dict[str, Any]] = []
    exact_matches: list[float] = []
    mode_shifts: list[float] = []
    assignment_margins: list[float] = []
    for left, right in zip(row_ind, col_ind):
        overlap_value = float(overlap[left, right])
        scaled_drift = abs(scaled_from[left] - scaled_to[right]) / max(abs(scaled_from[left]), abs(scaled_to[right]), 1.0e-12)
        row_alternatives = np.delete(overlap[left], right)
        col_alternatives = np.delete(overlap[:, right], left)
        row_runner_up = float(np.max(row_alternatives)) if row_alternatives.size else 0.0
        col_runner_up = float(np.max(col_alternatives)) if col_alternatives.size else 0.0
        assignment_margin = float(max(0.0, min(overlap_value - row_runner_up, overlap_value - col_runner_up)))
        exact_match = float(left == right)
        mode_shift = float(abs(left - right))
        matched_pairs.append(
            {
                'from_mode': int(left),
                'to_mode': int(right),
                'overlap': overlap_value,
                'scaled_eigen_from': float(scaled_from[left]),
                'scaled_eigen_to': float(scaled_to[right]),
                'scaled_eigen_drift': float(scaled_drift),
                'assignment_margin': assignment_margin,
                'exact_match': exact_match,
                'mode_shift': mode_shift,
            }
        )
        exact_matches.append(exact_match)
        mode_shifts.append(mode_shift)
        assignment_margins.append(assignment_margin)

    principal = subspace_cosines(modes_from, modes_to)
    leading_pair = next((pair for pair in matched_pairs if pair['from_mode'] == 0), None)
    overlaps = [pair['overlap'] for pair in matched_pairs]
    drifts = [pair['scaled_eigen_drift'] for pair in matched_pairs]
    return {
        'variant': variant,
        'n_from': int(n_from),
        'n_to': int(n_to),
        'probe_n_side': int(probe_n_side),
        'transport_modes': int(len(matched_pairs)),
        'matched_pairs': matched_pairs,
        'mean_overlap': float(np.mean(overlaps)),
        'min_overlap': float(np.min(overlaps)),
        'leading_mode_overlap': float(leading_pair['overlap']) if leading_pair is not None else math.nan,
        'principal_cosines': principal,
        'mean_principal_cosine': float(np.mean(principal)) if principal else math.nan,
        'min_principal_cosine': float(np.min(principal)) if principal else math.nan,
        'mean_scaled_eigen_drift': float(np.mean(drifts)),
        'max_scaled_eigen_drift': float(np.max(drifts)),
        'exact_match_fraction': float(np.mean(exact_matches)),
        'mean_mode_shift': float(np.mean(mode_shifts)),
        'max_mode_shift': float(np.max(mode_shifts)),
        'mean_assignment_margin': float(np.mean(assignment_margins)),
        'min_assignment_margin': float(np.min(assignment_margins)),
    }


def make_overlap_plot(metrics: dict[str, dict[str, Any]], sizes: list[int], variants: list[str], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    x_values = sizes[:-1]
    for variant in variants:
        mean_overlap = []
        min_overlap = []
        for n_from, n_to in zip(sizes[:-1], sizes[1:]):
            item = metrics[f'{variant}_n{n_from}_to_n{n_to}']
            mean_overlap.append(float(item['mean_overlap']))
            min_overlap.append(float(item['min_overlap']))
        axes[0].plot(x_values, mean_overlap, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
        axes[1].plot(x_values, min_overlap, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    axes[0].set_xlabel('coarse size n')
    axes[0].set_ylabel('mean matched overlap')
    axes[0].set_title('Active-sector mode transport overlap')
    axes[0].set_ylim(0.0, 1.02)
    axes[0].grid(alpha=0.25)
    axes[1].set_xlabel('coarse size n')
    axes[1].set_ylabel('minimum matched overlap')
    axes[1].set_title('Worst matched mode overlap')
    axes[1].set_ylim(0.0, 1.02)
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_stability_plot(metrics: dict[str, dict[str, Any]], sizes: list[int], variants: list[str], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    x_values = sizes[:-1]
    for variant in variants:
        principal = []
        drifts = []
        for n_from, n_to in zip(sizes[:-1], sizes[1:]):
            item = metrics[f'{variant}_n{n_from}_to_n{n_to}']
            principal.append(float(item['mean_principal_cosine']))
            drifts.append(float(item['max_scaled_eigen_drift']))
        axes[0].plot(x_values, principal, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
        axes[1].plot(x_values, drifts, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    axes[0].set_xlabel('coarse size n')
    axes[0].set_ylabel('mean principal cosine')
    axes[0].set_title('Transported active-subspace alignment')
    axes[0].set_ylim(0.0, 1.02)
    axes[0].grid(alpha=0.25)
    axes[1].set_xlabel('coarse size n')
    axes[1].set_ylabel('max relative scaled-eigen drift')
    axes[1].set_title('Matched active-branch eigen drift')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Active_Sector_Transport_v1.md'
    rows = []
    for key, item in sorted(result['transport_metrics'].items()):
        pair_text = ', '.join(f"{pair['from_mode']}->{pair['to_mode']}" for pair in item['matched_pairs'])
        rows.append(
            f"| {item['variant']} | {item['n_from']} -> {item['n_to']} | {item['mean_overlap']:.3f} | {item['min_overlap']:.3f} | {item['mean_principal_cosine']:.3f} | {item['max_scaled_eigen_drift']:.3f} | {pair_text} |"
        )
    note = f"""# Transverse Active Sector Transport

## Purpose

Execute `N1` of the continuum-closure program by making the `T8` comparison maps `J_n` explicit on the restricted coexact sector, then use those maps for a first low-window branch-identity test.

## Transport map

For each refinement level, the active restricted edge modes are transported into one common probe space:

- the unit torus is partitioned into a fixed `{result['config']['probe_n_side']}^3` probe lattice
- each edge amplitude is assigned to the probe cell containing its midpoint
- amplitudes are separated into `x`, `y`, and `z` channels
- each channel is normalized by the square root of the number of contributing edges

This defines the explicit comparison map `J_n` used in the present test.

## Setup

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- restricted modes: `{result['config']['restricted_modes']}`
- transported low-window modes: `{result['config']['transport_modes']}`
- probe lattice side: `{result['config']['probe_n_side']}`
- mild-disorder strength: `{result['config']['disorder_strength']}`

## Transport metrics

| branch | refinement | mean overlap | min overlap | mean principal cosine | max scaled-eigen drift | mode matching |
| --- | --- | ---: | ---: | ---: | ---: | --- |
{chr(10).join(rows)}

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


def run_transverse_active_sector_transport(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_active_sector_transport']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20, 24, 28])]
    variants = [str(v) for v in experiment_cfg.get('variants', ['baseline', 'puncture', 'line_defect', 'flux_tube', 'mild_disorder'])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 6))
    transport_modes_count = int(experiment_cfg.get('transport_modes', 4))
    probe_n_side = int(experiment_cfg.get('probe_n_side', 6))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))
    disorder_strength = float(experiment_cfg.get('disorder_strength', 0.12))

    cases: dict[str, Any] = {}
    for n_side in sizes:
        for variant in variants:
            key = f'{variant}_n{n_side}'
            cases[key] = analyze_case(
                n_side=n_side,
                epsilon=epsilon,
                variant=variant,
                restricted_modes=restricted_modes,
                harmonic_tol=harmonic_tol,
                eig_tol=eig_tol,
                penalty=penalty,
                flux_tube_phase=flux_tube_phase,
                disorder_strength=disorder_strength,
            )

    transport_metrics: dict[str, dict[str, Any]] = {}
    for variant in variants:
        for n_from, n_to in zip(sizes[:-1], sizes[1:]):
            key = f'{variant}_n{n_from}_to_n{n_to}'
            transport_metrics[key] = pair_transport_metrics(
                variant=variant,
                n_from=n_from,
                n_to=n_to,
                case_from=cases[f'{variant}_n{n_from}'],
                case_to=cases[f'{variant}_n{n_to}'],
                probe_n_side=probe_n_side,
                transport_modes_count=transport_modes_count,
            )

    overlap_plot = PLOTS / 'transverse_active_sector_transport_overlap.png'
    stability_plot = PLOTS / 'transverse_active_sector_transport_stability.png'
    make_overlap_plot(transport_metrics, sizes, variants, overlap_plot)
    make_stability_plot(transport_metrics, sizes, variants, stability_plot)
    plot_paths = [overlap_plot, stability_plot]

    all_pairs = list(transport_metrics.values())
    strongest_pair = max(all_pairs, key=lambda item: item['mean_overlap'])
    weakest_pair = min(all_pairs, key=lambda item: item['mean_overlap'])
    most_stable_pair = min(all_pairs, key=lambda item: item['max_scaled_eigen_drift'])
    observation = (
        'explicit probe-space comparison maps now carry the restricted coexact branch across refinement, so active-sector identity can be measured directly by transported overlaps, principal-angle alignment, and scaled-eigenvalue drift rather than inferred only from separate spectra'
    )
    conclusion = (
        f"within the tested low transported window, the strongest refinement pair is {strongest_pair['variant']} {strongest_pair['n_from']}->{strongest_pair['n_to']} with mean overlap {strongest_pair['mean_overlap']:.3f}, the weakest pair is {weakest_pair['variant']} {weakest_pair['n_from']}->{weakest_pair['n_to']} with mean overlap {weakest_pair['mean_overlap']:.3f}, and the smallest max scaled-eigen drift is {most_stable_pair['variant']} {most_stable_pair['n_from']}->{most_stable_pair['n_to']} at {most_stable_pair['max_scaled_eigen_drift']:.3f}; this gives `T8` a real `J_n` mechanism and an initial branch-identity test, but the continuum question still depends on whether these transported overlaps stay strong as the window and background family are enlarged"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'restricted_modes': restricted_modes,
            'transport_modes': transport_modes_count,
            'probe_n_side': probe_n_side,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'flux_tube_phase': flux_tube_phase,
            'disorder_strength': disorder_strength,
        },
        'transport_metrics': transport_metrics,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_active_sector_transport', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse active-sector transport',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variants={variants}, restricted_modes={restricted_modes}, transport_modes={transport_modes_count}, probe_n_side={probe_n_side}, disorder_strength={disorder_strength}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_active_sector_transport()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
