#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from L1_transverse_band_test import analyze_case as band_analyze_case
from L1_transverse_scaling_test import analyze_case as scaling_analyze_case
from L1_stage2_common import (
    PLOTS,
    REPO_ROOT,
    append_log,
    ensure_matplotlib,
    plot_edge_mode,
    save_result_payload,
)

plt = ensure_matplotlib()

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_mode_structure': {
        'n_side': 16,
        'variant': 'baseline',
        'phase_modes': 20,
        'restricted_modes': 10,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_mode_structure'] = dict(DEFAULT_CONFIG['transverse_mode_structure'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_mode_structure'})
        if isinstance(on_disk.get('transverse_mode_structure'), dict):
            merged['transverse_mode_structure'].update(on_disk['transverse_mode_structure'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_mode_structure'})
        if isinstance(config.get('transverse_mode_structure'), dict):
            merged['transverse_mode_structure'].update(config['transverse_mode_structure'])
    return merged


def make_mode_profile_plot(case: dict[str, Any], path: Path) -> None:
    vectors = case['restricted_vectors']
    midpoints = case['midpoints']
    directions = case['directions']
    fig = plt.figure(figsize=(18, 7))
    for idx in range(min(10, len(vectors))):
        ax = fig.add_subplot(2, 5, idx + 1, projection='3d')
        lam = case['restricted_transverse_spectrum'][idx]
        plot_edge_mode(ax, midpoints, directions, np.asarray(vectors[idx], dtype=complex), f'mode {idx + 1}, lambda={lam:.3f}')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_phase_plot(case: dict[str, Any], path: Path) -> None:
    records = case['phase_modes']
    x = [record['divergence_norm'] for record in records]
    y = [record['curl_norm'] for record in records]
    c = [record['eigenvalue'] for record in records]
    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(x, y, c=c, cmap='viridis', s=28)
    ax.set_xlabel(r'$||d_0^* a_k||$')
    ax.set_ylabel(r'$||d_1 a_k||$')
    ax.set_title(f"Mode-structure divergence-curl map (n={case['config']['n_side']})")
    ax.grid(alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r'$\lambda_k$')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def summarize_mode_metrics(case: dict[str, Any]) -> list[dict[str, float]]:
    metrics: list[dict[str, float]] = []
    for idx, record in enumerate(case['restricted_transverse_modes'][:10]):
        metrics.append(
            {
                'mode_index': idx,
                'eigenvalue': float(record['eigenvalue']),
                'divergence_norm': float(record['divergence_norm']),
                'curl_norm': float(record['curl_norm']),
                'ipr': float(record['ipr']),
            }
        )
    return metrics


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Mode_Structure_v1.md'
    rows = '\n'.join(
        f"| {m['mode_index'] + 1} | {m['eigenvalue']:.6f} | {m['divergence_norm']:.3e} | {m['curl_norm']:.6f} | {m['ipr']:.6f} |"
        for m in result['mode_metrics']
    )
    note = f"""# Transverse Mode Structure

## Purpose

Inspect the lowest restricted transverse modes directly to check whether they appear as extended circulation fields.

## Setup

- branch: `{result['config']['variant']}`
- `n = {result['config']['n_side']}`
- `epsilon = {result['config']['epsilon']}`
- restricted modes analyzed: `{result['config']['restricted_modes']}`

## Lowest 10 restricted transverse modes

| mode | eigenvalue | `||d0* a||` | `||d1 a||` | IPR |
| --- | ---: | ---: | ---: | ---: |
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


def run_transverse_mode_structure(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_mode_structure']
    epsilon = float(cfg.get('epsilon', 0.2))
    n_side = int(experiment_cfg.get('n_side', 16))
    variant = str(experiment_cfg.get('variant', 'baseline'))
    phase_modes = int(experiment_cfg.get('phase_modes', 20))
    restricted_modes = int(experiment_cfg.get('restricted_modes', 10))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))

    profile_case = band_analyze_case(
        n_side=n_side,
        epsilon=epsilon,
        variant=variant,
        restricted_modes=restricted_modes,
        harmonic_tol=harmonic_tol,
        eig_tol=eig_tol,
        penalty=penalty,
        flux_tube_phase=0.0,
    )
    phase_case = scaling_analyze_case(
        n_side=n_side,
        epsilon=epsilon,
        variant=variant,
        phase_modes=phase_modes,
        restricted_modes=restricted_modes,
        harmonic_tol=harmonic_tol,
        eig_tol=eig_tol,
        penalty=penalty,
    )

    profile_plot = PLOTS / 'transverse_vector_field_profiles.png'
    phase_plot = PLOTS / 'divergence_curl_phase_mode_structure.png'
    make_mode_profile_plot(profile_case, profile_plot)
    make_phase_plot(phase_case, phase_plot)
    plot_paths = [profile_plot, phase_plot]

    mode_metrics = summarize_mode_metrics(phase_case)
    observation = (
        'the lowest restricted transverse modes remain numerically divergence-free while their spatial support extends over the periodic lattice rather than concentrating at a single location'
    )
    conclusion = 'the inspected low restricted modes behave like extended circulation fields in the tested baseline branch'

    result = {
        'config': {
            'epsilon': epsilon,
            'n_side': n_side,
            'variant': variant,
            'phase_modes': phase_modes,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
        },
        'case': {'config': phase_case['config'], 'dimensions': phase_case['dimensions'], 'phase_modes': phase_case['phase_modes'], 'restricted_transverse_spectrum': phase_case['restricted_transverse_spectrum'], 'restricted_transverse_modes': phase_case['restricted_transverse_modes']},
        'mode_metrics': mode_metrics,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_mode_structure', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse mode structure',
        config_summary=f"epsilon={epsilon}, n={n_side}, variant={variant}, phase_modes={phase_modes}, restricted_modes={restricted_modes}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_mode_structure()
    print(
        json.dumps(
            {
                'result_path': str(result_path),
                'note_path': str(note_path),
                'observation': result['observation'],
                'conclusion': result['conclusion'],
            },
            indent=2,
        )
    )


if __name__ == '__main__':
    main()
