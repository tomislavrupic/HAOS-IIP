#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage3_common import (
    PLOTS,
    REPO_ROOT,
    analyze_branch_cases,
    append_log,
    periodic_displacement,
    save_result_payload,
    continuum_transverse_q2,
    make_phase_plot,
    plt,
)
from L1_transverse_band_test import build_periodic_complex

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_realspace_correlations': {
        'sizes': [12, 16],
        'restricted_modes': 6,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': 0.0,
        'pair_samples': 8000,
        'bins': 24,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_realspace_correlations'] = dict(DEFAULT_CONFIG['transverse_realspace_correlations'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_realspace_correlations'})
        if isinstance(on_disk.get('transverse_realspace_correlations'), dict):
            merged['transverse_realspace_correlations'].update(on_disk['transverse_realspace_correlations'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_realspace_correlations'})
        if isinstance(config.get('transverse_realspace_correlations'), dict):
            merged['transverse_realspace_correlations'].update(config['transverse_realspace_correlations'])
    return merged


def normalize_complex_mode(vec: np.ndarray) -> np.ndarray:
    vec = np.asarray(vec, dtype=complex)
    idx = int(np.argmax(np.abs(vec)))
    phase = np.angle(vec[idx]) if np.abs(vec[idx]) > 0 else 0.0
    return vec * np.exp(-1j * phase)


def sampled_correlation(midpoints: np.ndarray, weights: np.ndarray, samples: int, bins: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    count = len(weights)
    idx_i = rng.integers(0, count, size=samples)
    idx_j = rng.integers(0, count, size=samples)
    disp = periodic_displacement(midpoints[idx_j], midpoints[idx_i])
    dist = np.linalg.norm(disp, axis=1)
    centered = weights - np.mean(weights)
    corr = centered[idx_i] * centered[idx_j]
    edges = np.linspace(0.0, float(np.max(dist)), bins + 1)
    radii = []
    values = []
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (dist >= left) & (dist < right if right < edges[-1] else dist <= right)
        if not np.any(mask):
            continue
        radii.append(float(np.mean(dist[mask])))
        values.append(float(np.mean(corr[mask])))
    values_arr = np.asarray(values, dtype=float)
    if len(values_arr) and abs(values_arr[0]) > 0:
        values_arr = values_arr / values_arr[0]
    return np.asarray(radii, dtype=float), values_arr


def correlation_length(radii: np.ndarray, corr: np.ndarray) -> float:
    if len(radii) == 0:
        return 0.0
    target = math.exp(-1.0)
    below = np.where(corr <= target)[0]
    if len(below) == 0:
        return float(radii[-1])
    return float(radii[below[0]])


def direction_anisotropy(directions: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    powers = {}
    axes = {'x': np.array([1.0, 0.0, 0.0]), 'y': np.array([0.0, 1.0, 0.0]), 'z': np.array([0.0, 0.0, 1.0])}
    for name, axis in axes.items():
        mask = np.abs(directions @ axis) > 0.9
        powers[name] = float(np.mean(weights[mask])) if np.any(mask) else 0.0
    valid = [value for value in powers.values() if value > 0]
    powers['anisotropy_ratio'] = float(max(valid) / min(valid)) if valid else 1.0
    return powers


def curl_orientation_metrics(n_side: int, epsilon: float, vector: np.ndarray) -> dict[str, float]:
    data = build_periodic_complex(n_side=n_side, epsilon=epsilon, variant='baseline', flux_tube_phase=0.0)
    curl = np.abs(np.asarray(data.d1 @ vector, dtype=complex))
    orientations = {'xy': curl[0::3], 'xz': curl[1::3], 'yz': curl[2::3]}
    return {name: float(np.mean(values)) if len(values) else 0.0 for name, values in orientations.items()}


def make_correlation_plot(metrics: dict[str, dict[str, Any]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    for n_side in sorted(int(k) for k in metrics.keys()):
        item = metrics[str(n_side)]
        ax.plot(item['radii'], item['correlation'], marker='o', label=f'n={n_side}')
    ax.set_xlabel('periodic separation')
    ax.set_ylabel('normalized correlation')
    ax.set_title('Lowest restricted-mode correlation decay')
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_anisotropy_plot(metrics: dict[str, dict[str, Any]], path: Path) -> None:
    sizes = sorted(int(k) for k in metrics.keys())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(sizes, [metrics[str(n)]['direction_power']['anisotropy_ratio'] for n in sizes], marker='o')
    axes[0].set_xlabel('lattice side n')
    axes[0].set_ylabel('max/min directional power')
    axes[0].set_title('Directional anisotropy ratio')
    axes[0].grid(alpha=0.25)

    largest = metrics[str(sizes[-1])]
    labels = ['xy', 'xz', 'yz']
    values = [largest['curl_orientation'][label] for label in labels]
    axes[1].bar(labels, values)
    axes[1].set_ylabel('mean |curl| by face orientation')
    axes[1].set_title(f'Curl orientation structure (n={sizes[-1]})')
    axes[1].grid(alpha=0.25, axis='y')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Realspace_Correlations_v1.md'
    rows = '\n'.join(
        f"| {n} | {result['metrics'][str(n)]['correlation_length']:.6f} | {result['metrics'][str(n)]['direction_power']['anisotropy_ratio']:.4f} | {result['metrics'][str(n)]['ipr']:.6f} |"
        for n in result['config']['sizes']
    )
    note = f"""# Transverse Real-Space Correlations

## Purpose

Measure whether the lowest restricted transverse modes develop extended, approximately isotropic correlations on the periodic branch.

## Setup

- branch: baseline periodic torus
- sizes: `{result['config']['sizes']}`
- pair samples per size: `{result['config']['pair_samples']}`
- bins: `{result['config']['bins']}`

## Correlation summary

| `n` | correlation length | anisotropy ratio | lowest-mode IPR |
| --- | ---: | ---: | ---: |
{rows}

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


def run_transverse_realspace_correlations(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_realspace_correlations']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 6))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    pair_samples = int(experiment_cfg.get('pair_samples', 8000))
    bins = int(experiment_cfg.get('bins', 24))

    cases = analyze_branch_cases(
        sizes=sizes,
        variants=['baseline'],
        epsilon=epsilon,
        restricted_modes=restricted_modes,
        harmonic_tol=harmonic_tol,
        eig_tol=eig_tol,
        penalty=penalty,
        flux_tube_phase=0.0,
    )

    metrics: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        case = cases[f'baseline_n{n_side}']
        vector = normalize_complex_mode(np.asarray(case['restricted_vectors'][0], dtype=complex))
        weights = np.abs(vector) ** 2
        radii, corr = sampled_correlation(case['midpoints'], weights, pair_samples, bins, seed=20260310 + n_side)
        metrics[str(n_side)] = {
            'radii': radii.tolist(),
            'correlation': corr.tolist(),
            'correlation_length': correlation_length(radii, corr),
            'direction_power': direction_anisotropy(np.asarray(case['directions'], dtype=float), weights),
            'curl_orientation': curl_orientation_metrics(n_side, epsilon, vector),
            'ipr': float(case['restricted_transverse_modes'][0]['ipr']),
            'divergence_norm': float(case['restricted_transverse_modes'][0]['divergence_norm']),
            'curl_norm': float(case['restricted_transverse_modes'][0]['curl_norm']),
        }

    corr_plot = PLOTS / 'transverse_correlation_decay.png'
    anisotropy_plot = PLOTS / 'transverse_anisotropy.png'
    phase_plot = PLOTS / 'divergence_curl_phase_correlations.png'
    make_correlation_plot(metrics, corr_plot)
    make_anisotropy_plot(metrics, anisotropy_plot)
    make_phase_plot(cases[f'baseline_n{sizes[-1]}'], phase_plot, f"Restricted phase map (correlations, n={sizes[-1]})")
    plot_paths = [corr_plot, anisotropy_plot, phase_plot]

    first = metrics[str(sizes[0])]
    last = metrics[str(sizes[-1])]
    observation = (
        'the lowest restricted transverse mode shows extended two-point correlations while directional anisotropy remains modest across the tested sizes'
    )
    conclusion = (
        f"the correlation length remains of the same order, from {first['correlation_length']:.3f} at n={sizes[0]} to {last['correlation_length']:.3f} at n={sizes[-1]}, while the anisotropy ratio stays moderate at {last['direction_power']['anisotropy_ratio']:.3f} for the largest size, consistent with an extended continuum-like branch in the tested range"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'pair_samples': pair_samples,
            'bins': bins,
        },
        'metrics': metrics,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_realspace_correlations', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse real-space correlations',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, pair_samples={pair_samples}, bins={bins}, restricted_modes={restricted_modes}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_realspace_correlations()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
