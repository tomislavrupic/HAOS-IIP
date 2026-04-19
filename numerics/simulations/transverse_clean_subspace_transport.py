#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage3_common import PLOTS, REPO_ROOT, VARIANT_LABELS, VARIANT_MARKERS, append_log, save_result_payload, plt
from L1_transverse_band_test import analyze_case
from transverse_active_sector_transport import subspace_cosines, transport_mode_to_probe

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_clean_subspace_transport': {
        'sizes': [12, 16, 20],
        'variants': ['baseline', 'mild_disorder'],
        'restricted_modes': 8,
        'window_size': 4,
        'scan_depth': 8,
        'probe_n_side': 6,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': math.pi / 2.0,
        'disorder_strength': 0.12,
        'thresholds': {
            'min_best_projector_affinity': 0.80,
            'min_best_mean_principal_cosine': 0.80,
            'max_best_scaled_window_drift': 0.05,
            'min_affinity_gain': 0.20,
        },
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_clean_subspace_transport'] = dict(DEFAULT_CONFIG['transverse_clean_subspace_transport'])
    merged['transverse_clean_subspace_transport']['thresholds'] = dict(DEFAULT_CONFIG['transverse_clean_subspace_transport']['thresholds'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_clean_subspace_transport'})
        if isinstance(on_disk.get('transverse_clean_subspace_transport'), dict):
            block = dict(on_disk['transverse_clean_subspace_transport'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_clean_subspace_transport'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_clean_subspace_transport']['thresholds'].update(thresholds)
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_clean_subspace_transport'})
        if isinstance(config.get('transverse_clean_subspace_transport'), dict):
            block = dict(config['transverse_clean_subspace_transport'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_clean_subspace_transport'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_clean_subspace_transport']['thresholds'].update(thresholds)
    return merged


def transported_window(case: dict[str, Any], probe_n_side: int, start: int, window_size: int) -> list[np.ndarray]:
    vectors = list(case.get('restricted_vectors', []))[start:start + window_size]
    transported: list[np.ndarray] = []
    for vec in vectors:
        moved = transport_mode_to_probe(case['midpoints'], case['directions'], np.asarray(vec, dtype=complex), probe_n_side)
        if float(np.linalg.norm(moved)) > 1.0e-12:
            transported.append(moved)
    return transported


def orthonormal_basis(vectors: list[np.ndarray]) -> np.ndarray:
    if not vectors:
        return np.zeros((0, 0), dtype=complex)
    matrix = np.column_stack(vectors)
    q, _ = np.linalg.qr(matrix, mode='reduced')
    return q


def projector_affinity(vectors_a: list[np.ndarray], vectors_b: list[np.ndarray]) -> float:
    if not vectors_a or not vectors_b:
        return math.nan
    qa = orthonormal_basis(vectors_a)
    qb = orthonormal_basis(vectors_b)
    overlap = qa.conj().T @ qb
    singular_values = np.linalg.svd(overlap, compute_uv=False)
    if singular_values.size == 0:
        return math.nan
    return float(np.mean(np.square(np.clip(np.real(singular_values), 0.0, 1.0))))


def scaled_window_drift(
    case_from: dict[str, Any],
    case_to: dict[str, Any],
    n_from: int,
    n_to: int,
    start_from: int,
    start_to: int,
    window_size: int,
) -> float:
    spec_from = np.asarray(case_from['restricted_transverse_spectrum'][start_from:start_from + window_size], dtype=float)
    spec_to = np.asarray(case_to['restricted_transverse_spectrum'][start_to:start_to + window_size], dtype=float)
    count = min(len(spec_from), len(spec_to))
    if count == 0:
        return math.nan
    lhs = np.sort((n_from**2) * spec_from[:count])
    rhs = np.sort((n_to**2) * spec_to[:count])
    denom = np.maximum(np.maximum(np.abs(lhs), np.abs(rhs)), 1.0e-12)
    return float(np.max(np.abs(lhs - rhs) / denom))


def window_metrics(
    case_from: dict[str, Any],
    case_to: dict[str, Any],
    n_from: int,
    n_to: int,
    start_from: int,
    start_to: int,
    window_size: int,
    probe_n_side: int,
) -> dict[str, Any]:
    modes_from = transported_window(case_from, probe_n_side, start_from, window_size)
    modes_to = transported_window(case_to, probe_n_side, start_to, window_size)
    principal = subspace_cosines(modes_from, modes_to)
    return {
        'start_from': int(start_from),
        'start_to': int(start_to),
        'projector_affinity': projector_affinity(modes_from, modes_to),
        'mean_principal_cosine': float(np.mean(principal)) if principal else math.nan,
        'min_principal_cosine': float(np.min(principal)) if principal else math.nan,
        'scaled_window_drift': scaled_window_drift(case_from, case_to, n_from, n_to, start_from, start_to, window_size),
    }


def best_window_metrics(
    case_from: dict[str, Any],
    case_to: dict[str, Any],
    n_from: int,
    n_to: int,
    scan_depth: int,
    window_size: int,
    probe_n_side: int,
) -> dict[str, Any]:
    available_from = max(0, min(scan_depth, len(case_from.get('restricted_vectors', []))) - window_size + 1)
    available_to = max(0, min(scan_depth, len(case_to.get('restricted_vectors', []))) - window_size + 1)
    candidates: list[dict[str, Any]] = []
    for start_from in range(available_from):
        for start_to in range(available_to):
            metrics = window_metrics(case_from, case_to, n_from, n_to, start_from, start_to, window_size, probe_n_side)
            candidates.append(metrics)
    if not candidates:
        return {
            'start_from': 0,
            'start_to': 0,
            'projector_affinity': math.nan,
            'mean_principal_cosine': math.nan,
            'min_principal_cosine': math.nan,
            'scaled_window_drift': math.nan,
        }
    candidates.sort(
        key=lambda item: (
            float(item['projector_affinity']),
            float(item['mean_principal_cosine']),
            -float(item['scaled_window_drift']),
        ),
        reverse=True,
    )
    return candidates[0]


def evaluate_pair(best_metrics: dict[str, Any], fixed_metrics: dict[str, Any], thresholds: dict[str, float]) -> dict[str, Any]:
    fixed_checks = {
        'projector_affinity': float(fixed_metrics['projector_affinity']) >= float(thresholds['min_best_projector_affinity']),
        'mean_principal_cosine': float(fixed_metrics['mean_principal_cosine']) >= float(thresholds['min_best_mean_principal_cosine']),
        'scaled_window_drift': float(fixed_metrics['scaled_window_drift']) <= float(thresholds['max_best_scaled_window_drift']),
    }
    best_checks = {
        'projector_affinity': float(best_metrics['projector_affinity']) >= float(thresholds['min_best_projector_affinity']),
        'mean_principal_cosine': float(best_metrics['mean_principal_cosine']) >= float(thresholds['min_best_mean_principal_cosine']),
        'scaled_window_drift': float(best_metrics['scaled_window_drift']) <= float(thresholds['max_best_scaled_window_drift']),
    }
    affinity_gain = float(best_metrics['projector_affinity']) - float(fixed_metrics['projector_affinity'])
    rescued = bool(all(best_checks.values()) and affinity_gain >= float(thresholds['min_affinity_gain']))
    fixed_stable = bool(all(fixed_checks.values()))
    if fixed_stable:
        status = 'stable_fixed'
    elif rescued:
        status = 'rescued_by_window'
    else:
        status = 'open'
    return {
        'fixed_checks': fixed_checks,
        'best_checks': best_checks,
        'affinity_gain': affinity_gain,
        'status': status,
        'resolved': bool(status != 'open'),
    }


def make_affinity_plot(pair_rows_by_variant: dict[str, list[dict[str, Any]]], sizes: list[int], path: Path) -> None:
    x_values = sizes[:-1]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant, rows in pair_rows_by_variant.items():
        axes[0].plot(
            x_values,
            [row['fixed_window']['projector_affinity'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=f"{VARIANT_LABELS[variant]} fixed",
        )
        axes[0].plot(
            x_values,
            [row['best_window']['projector_affinity'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            linestyle='--',
            label=f"{VARIANT_LABELS[variant]} best",
        )
        axes[1].plot(
            x_values,
            [row['best_window']['mean_principal_cosine'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=VARIANT_LABELS[variant],
        )
    axes[0].axhline(0.80, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[0].set_xlabel('coarse size n')
    axes[0].set_ylabel('projector affinity')
    axes[0].set_title('Fixed vs best low-window subspace affinity')
    axes[0].set_ylim(0.0, 1.02)
    axes[0].grid(alpha=0.25)
    axes[1].axhline(0.80, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel('coarse size n')
    axes[1].set_ylabel('best mean principal cosine')
    axes[1].set_title('Best-window subspace alignment')
    axes[1].set_ylim(0.0, 1.02)
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_shift_plot(pair_rows_by_variant: dict[str, list[dict[str, Any]]], sizes: list[int], path: Path) -> None:
    x_values = sizes[:-1]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant, rows in pair_rows_by_variant.items():
        axes[0].plot(
            x_values,
            [row['best_window']['start_from'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=f"{VARIANT_LABELS[variant]} coarse",
        )
        axes[0].plot(
            x_values,
            [row['best_window']['start_to'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            linestyle='--',
            label=f"{VARIANT_LABELS[variant]} fine",
        )
        axes[1].plot(
            x_values,
            [row['best_window']['scaled_window_drift'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=VARIANT_LABELS[variant],
        )
    axes[0].set_xlabel('coarse size n')
    axes[0].set_ylabel('best window start')
    axes[0].set_title('Best window shift under refinement')
    axes[0].grid(alpha=0.25)
    axes[1].axhline(0.05, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel('coarse size n')
    axes[1].set_ylabel('max scaled-window drift')
    axes[1].set_title('Best-window spectral drift')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Clean_Subspace_Transport_v1.md'
    threshold_lines = [
        f"- best projector affinity >= `{result['thresholds']['min_best_projector_affinity']:.2f}`",
        f"- best mean principal cosine >= `{result['thresholds']['min_best_mean_principal_cosine']:.2f}`",
        f"- best scaled-window drift <= `{result['thresholds']['max_best_scaled_window_drift']:.2f}`",
        f"- best minus fixed affinity gain >= `{result['thresholds']['min_affinity_gain']:.2f}`",
    ]
    rows = []
    for item in result['pair_rows']:
        fixed = item['fixed_window']
        best = item['best_window']
        rows.append(
            f"| {item['variant']} | {item['n_from']} -> {item['n_to']} | {fixed['projector_affinity']:.3f} | {best['projector_affinity']:.3f} | {best['mean_principal_cosine']:.3f} | {best['scaled_window_drift']:.3f} | {best['start_from']}->{best['start_to']} | {item['evaluation']['status']} |"
        )
    summary_rows = []
    for item in result['variant_summary']:
        summary_rows.append(
            f"| {item['variant']} | {item['resolved_pair_count']}/{item['total_pairs']} | {item['worst_fixed_affinity']:.3f} | {item['worst_best_affinity']:.3f} | {item['largest_affinity_gain']:.3f} | {item['worst_best_mean_principal_cosine']:.3f} | {item['worst_best_scaled_window_drift']:.3f} | {'RESOLVED' if item['all_pairs_resolved'] else 'OPEN'} |"
        )
    note = f"""# Transverse Clean Subspace Transport

## Purpose

Fast-track the clean-baseline identity question by treating the low active coexact branch as a potentially rotating low-dimensional subspace rather than as a fixed mode-by-mode labeling problem.

This bounded `CB1-lite` scan keeps the same probe-space transport map, the same active coexact sector, and the same refinement pairs. The only added freedom is that the low transported window may slide inside the first `{result['config']['scan_depth']}` restricted coexact modes.

## Setup

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- restricted modes analyzed: `{result['config']['restricted_modes']}`
- window size: `{result['config']['window_size']}`
- scan depth: `{result['config']['scan_depth']}`
- probe lattice side: `{result['config']['probe_n_side']}`

## CB1-lite criteria

{chr(10).join(threshold_lines)}

If the fixed window already meets them, the branch is treated as already stable. If only the shifted best window meets them with a substantial affinity gain, the pair is treated as rescued by windowed subspace transport. If neither happens, the clean-baseline ambiguity remains open.

## Pair results

| branch | refinement | fixed affinity | best affinity | best mean principal cosine | best scaled-window drift | best window shift | status |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
{chr(10).join(rows)}

## Variant summary

| branch | resolved pairs | worst fixed affinity | worst best affinity | largest affinity gain | worst best mean principal cosine | worst best drift | CB1-lite status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
{chr(10).join(summary_rows)}

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


def run_transverse_clean_subspace_transport(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_clean_subspace_transport']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20])]
    variants = [str(v) for v in experiment_cfg.get('variants', ['baseline', 'mild_disorder'])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 8))
    window_size = int(experiment_cfg.get('window_size', 4))
    scan_depth = int(experiment_cfg.get('scan_depth', 8))
    probe_n_side = int(experiment_cfg.get('probe_n_side', 6))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))
    disorder_strength = float(experiment_cfg.get('disorder_strength', 0.12))
    thresholds = {k: float(v) for k, v in experiment_cfg.get('thresholds', {}).items()}

    cases: dict[str, Any] = {}
    for n_side in sizes:
        for variant in variants:
            print(f'[cb1-lite] analyze_case variant={variant} n={n_side}', flush=True)
            kwargs: dict[str, Any] = {
                'n_side': n_side,
                'epsilon': epsilon,
                'variant': variant,
                'restricted_modes': restricted_modes,
                'harmonic_tol': harmonic_tol,
                'eig_tol': eig_tol,
                'penalty': penalty,
                'flux_tube_phase': flux_tube_phase,
            }
            if variant == 'mild_disorder':
                kwargs['disorder_strength'] = disorder_strength
            cases[f'{variant}_n{n_side}'] = analyze_case(**kwargs)

    pair_rows: list[dict[str, Any]] = []
    pair_rows_by_variant: dict[str, list[dict[str, Any]]] = {variant: [] for variant in variants}
    for variant in variants:
        for n_from, n_to in zip(sizes[:-1], sizes[1:]):
            case_from = cases[f'{variant}_n{n_from}']
            case_to = cases[f'{variant}_n{n_to}']
            fixed_metrics = window_metrics(case_from, case_to, n_from, n_to, 0, 0, window_size, probe_n_side)
            best_metrics = best_window_metrics(case_from, case_to, n_from, n_to, scan_depth, window_size, probe_n_side)
            evaluation = evaluate_pair(best_metrics, fixed_metrics, thresholds)
            row = {
                'variant': variant,
                'n_from': int(n_from),
                'n_to': int(n_to),
                'fixed_window': fixed_metrics,
                'best_window': best_metrics,
                'evaluation': evaluation,
            }
            pair_rows.append(row)
            pair_rows_by_variant[variant].append(row)
            print(
                f"[cb1-lite] variant={variant} {n_from}->{n_to} "
                f"fixed_affinity={fixed_metrics['projector_affinity']:.3f} "
                f"best_affinity={best_metrics['projector_affinity']:.3f} "
                f"best_shift={best_metrics['start_from']}->{best_metrics['start_to']}",
                flush=True,
            )

    variant_summary: list[dict[str, Any]] = []
    for variant, rows in pair_rows_by_variant.items():
        variant_summary.append(
            {
                'variant': variant,
                'resolved_pair_count': int(sum(1 for row in rows if row['evaluation']['resolved'])),
                'total_pairs': int(len(rows)),
                'all_pairs_resolved': bool(all(row['evaluation']['resolved'] for row in rows)),
                'worst_fixed_affinity': float(min(row['fixed_window']['projector_affinity'] for row in rows)),
                'worst_best_affinity': float(min(row['best_window']['projector_affinity'] for row in rows)),
                'largest_affinity_gain': float(max(row['best_window']['projector_affinity'] - row['fixed_window']['projector_affinity'] for row in rows)),
                'worst_best_mean_principal_cosine': float(min(row['best_window']['mean_principal_cosine'] for row in rows)),
                'worst_best_scaled_window_drift': float(max(row['best_window']['scaled_window_drift'] for row in rows)),
            }
        )

    strongest_gain = max(pair_rows, key=lambda row: row['best_window']['projector_affinity'] - row['fixed_window']['projector_affinity'])
    weakest_best = min(pair_rows, key=lambda row: row['best_window']['projector_affinity'])
    observation = (
        "the clean-baseline identity question can be fast-tracked by scanning low contiguous coexact windows as transported subspaces; "
        "this keeps the same active-sector transport map while allowing the clean basis and the low window to rotate or slide inside a bounded scan depth"
    )
    if any(row['variant'] == 'baseline' and row['evaluation']['status'] == 'rescued_by_window' for row in pair_rows):
        conclusion = (
            f"the fast-track CB1-lite scan supports a degeneracy reading: the strongest gain is on "
            f"{strongest_gain['variant']} {strongest_gain['n_from']}->{strongest_gain['n_to']}, where projector affinity rises from "
            f"{strongest_gain['fixed_window']['projector_affinity']:.3f} to {strongest_gain['best_window']['projector_affinity']:.3f}; "
            f"this means at least part of the clean-baseline weakness is consistent with rotating or shifted low-window identity rather than branch death"
        )
    else:
        conclusion = (
            f"the fast-track CB1-lite scan improves low-window identification but does not yet close the clean baseline: "
            f"the weakest best-window case is {weakest_best['variant']} {weakest_best['n_from']}->{weakest_best['n_to']} with best projector affinity "
            f"{weakest_best['best_window']['projector_affinity']:.3f}; this means the next shortest move is still CB2-style infinitesimal symmetry pinning, "
            f"but we now know whether simple window sliding already resolves the issue"
        )

    result = {
        'config': {
            'sizes': sizes,
            'variants': variants,
            'restricted_modes': restricted_modes,
            'window_size': window_size,
            'scan_depth': scan_depth,
            'probe_n_side': probe_n_side,
            'penalty': penalty,
            'disorder_strength': disorder_strength,
        },
        'thresholds': thresholds,
        'pair_rows': pair_rows,
        'variant_summary': variant_summary,
        'observation': observation,
        'conclusion': conclusion,
    }

    affinity_plot = PLOTS / 'transverse_clean_subspace_transport_affinity.png'
    shift_plot = PLOTS / 'transverse_clean_subspace_transport_shift.png'
    make_affinity_plot(pair_rows_by_variant, sizes, affinity_plot)
    make_shift_plot(pair_rows_by_variant, sizes, shift_plot)
    plot_paths = [affinity_plot, shift_plot]
    result_path, stamped_plots, timestamp = save_result_payload('transverse_clean_subspace_transport', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse clean subspace transport',
        config_summary=(
            f"epsilon={epsilon}, sizes={sizes}, variants={variants}, restricted_modes={restricted_modes}, "
            f"window_size={window_size}, scan_depth={scan_depth}, probe_n_side={probe_n_side}, "
            f"disorder_strength={disorder_strength}, thresholds={thresholds}"
        ),
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


if __name__ == '__main__':
    result, result_path, _, _, note_path = run_transverse_clean_subspace_transport()
    print(json.dumps(result, indent=2))
    print(f'results: {result_path.relative_to(REPO_ROOT)}')
    print(f'note: {note_path.relative_to(REPO_ROOT)}')
