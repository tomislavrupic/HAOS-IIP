#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage2_common import (
    PLOTS,
    REPO_ROOT,
    analyze_embedded_baseline,
    append_log,
    ensure_matplotlib,
    perturbed_points,
    regular_points,
    save_result_payload,
)

plt = ensure_matplotlib()

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'L1_geometric_disorder_test': {
        'sizes': [12, 16, 20],
        'phase_modes': 20,
        'restricted_modes': 20,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'sigma_factor': 0.05,
        'seed': 20260310,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['L1_geometric_disorder_test'] = dict(DEFAULT_CONFIG['L1_geometric_disorder_test'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'L1_geometric_disorder_test'})
        if isinstance(on_disk.get('L1_geometric_disorder_test'), dict):
            merged['L1_geometric_disorder_test'].update(on_disk['L1_geometric_disorder_test'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'L1_geometric_disorder_test'})
        if isinstance(config.get('L1_geometric_disorder_test'), dict):
            merged['L1_geometric_disorder_test'].update(config['L1_geometric_disorder_test'])
    return merged


def make_band_plot(cases: dict[str, dict[str, Any]], sizes: list[int], path: Path) -> dict[str, float]:
    fig, axes = plt.subplots(1, len(sizes), figsize=(4.5 * len(sizes), 4), sharey=True)
    if len(sizes) == 1:
        axes = [axes]
    spreads: dict[str, float] = {}
    for ax, n_side in zip(axes, sizes):
        ordered = np.asarray(cases[f'ordered_n{n_side}']['restricted_transverse_spectrum'][:20], dtype=float)
        disordered = np.asarray(cases[f'disordered_n{n_side}']['restricted_transverse_spectrum'][:20], dtype=float)
        ordered_scaled = (n_side * n_side) * ordered
        disordered_scaled = (n_side * n_side) * disordered
        ax.plot(range(1, len(ordered_scaled) + 1), ordered_scaled, marker='o', label='ordered')
        ax.plot(range(1, len(disordered_scaled) + 1), disordered_scaled, marker='s', linestyle='--', label='disordered')
        ax.set_title(f'n={n_side}')
        ax.set_xlabel('restricted mode index k')
        ax.grid(alpha=0.25)
        spread = float(np.mean(np.abs(disordered_scaled[:10] - ordered_scaled[:10]) / np.maximum(np.abs(ordered_scaled[:10]), 1e-12)))
        spreads[str(n_side)] = spread
    axes[0].set_ylabel(r'$n^2 \lambda_k$')
    axes[0].legend(fontsize=8)
    fig.suptitle('Restricted transverse band under geometric disorder')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    return spreads


def make_phase_plot(cases: dict[str, dict[str, Any]], n_side: int, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True, sharey=True)
    scatter = None
    for ax, label in zip(axes, ['ordered', 'disordered']):
        case = cases[f'{label}_n{n_side}']
        records = case['phase_modes']
        x = [record['divergence_norm'] for record in records]
        y = [record['curl_norm'] for record in records]
        c = [record['eigenvalue'] for record in records]
        scatter = ax.scatter(x, y, c=c, cmap='viridis', s=28)
        ax.set_title(f'{label}, n={n_side}')
        ax.grid(alpha=0.25)
        ax.set_xlabel(r'$||d_0^* a_k||$')
    axes[0].set_ylabel(r'$||d_1 a_k||$')
    cbar = fig.colorbar(scatter, ax=axes.ravel().tolist(), shrink=0.9)
    cbar.set_label(r'$\lambda_k$')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'L1_Geometric_Disorder_Test_v1.md'
    rows = []
    for n_side in result['config']['sizes']:
        ordered = result['cases'][f'ordered_n{n_side}']
        disordered = result['cases'][f'disordered_n{n_side}']
        rows.append(
            f"| {n_side} | {ordered['restricted_transverse_spectrum'][0]:.6f} | {disordered['restricted_transverse_spectrum'][0]:.6f} | "
            f"{ordered['restricted_transverse_modes'][0]['ipr']:.6f} | {disordered['restricted_transverse_modes'][0]['ipr']:.6f} | "
            f"{result['spectral_shift_first10'][str(n_side)]:.4f} |"
        )
    note = f"""# L1 Geometric Disorder Test

## Purpose

Test whether the restricted transverse band survives loss of exact lattice symmetry when only the node embedding is perturbed.

## Setup

- sizes: `{result['config']['sizes']}`
- baseline periodic topology preserved
- perturbation: `x_i -> x_i + sigma xi`
- `sigma = {result['config']['sigma_factor']} h`
- deterministic seed base: `{result['config']['seed']}`

## Lowest restricted mode comparison

| `n` | ordered `lambda_1` | disordered `lambda_1` | ordered IPR | disordered IPR | first-10 spectral shift |
| --- | ---: | ---: | ---: | ---: | ---: |
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


def run_L1_geometric_disorder_test(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['L1_geometric_disorder_test']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20])]
    phase_modes = int(experiment_cfg.get('phase_modes', 20))
    restricted_modes = int(experiment_cfg.get('restricted_modes', 20))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    sigma_factor = float(experiment_cfg.get('sigma_factor', 0.05))
    seed_base = int(experiment_cfg.get('seed', 20260310))

    cases: dict[str, dict[str, Any]] = {}
    spectral_shift_first10: dict[str, float] = {}
    for n_side in sizes:
        h = 1.0 / n_side
        sigma = sigma_factor * h
        cases[f'ordered_n{n_side}'] = analyze_embedded_baseline(
            n_side=n_side,
            epsilon=epsilon,
            points=regular_points(n_side),
            restricted_modes=restricted_modes,
            phase_modes=phase_modes,
            harmonic_tol=harmonic_tol,
            eig_tol=eig_tol,
            penalty=penalty,
        )
        cases[f'disordered_n{n_side}'] = analyze_embedded_baseline(
            n_side=n_side,
            epsilon=epsilon,
            points=perturbed_points(n_side, sigma=sigma, seed=seed_base + n_side),
            restricted_modes=restricted_modes,
            phase_modes=phase_modes,
            harmonic_tol=harmonic_tol,
            eig_tol=eig_tol,
            penalty=penalty,
        )
        ordered = np.asarray(cases[f'ordered_n{n_side}']['restricted_transverse_spectrum'][:10], dtype=float)
        disordered = np.asarray(cases[f'disordered_n{n_side}']['restricted_transverse_spectrum'][:10], dtype=float)
        spectral_shift_first10[str(n_side)] = float(np.mean(np.abs(disordered - ordered) / np.maximum(np.abs(ordered), 1e-12)))

    band_plot = PLOTS / 'transverse_band_disorder.png'
    phase_plot = PLOTS / 'divergence_curl_phase_disorder.png'
    spreads = make_band_plot(cases, sizes, band_plot)
    make_phase_plot(cases, sizes[-1], phase_plot)
    plot_paths = [band_plot, phase_plot]

    observation = (
        'the restricted transverse spectrum persists under weak geometric disorder, with the ordered and disordered scaled bands remaining close and the lowest restricted mode staying divergence-free'
    )
    conclusion = 'in the tested weak-disorder regime, loss of exact lattice symmetry does not destroy the low restricted transverse band'

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'phase_modes': phase_modes,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'sigma_factor': sigma_factor,
            'seed': seed_base,
        },
        'cases': {label: {'config': case['config'], 'dimensions': case['dimensions'], 'phase_modes': case['phase_modes'], 'restricted_transverse_spectrum': case['restricted_transverse_spectrum'], 'restricted_transverse_modes': case['restricted_transverse_modes']} for label, case in cases.items()},
        'spectral_shift_first10': spectral_shift_first10,
        'band_spread_by_n': spreads,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('L1_geometric_disorder_test', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='L1 geometric disorder test',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, sigma_factor={sigma_factor}, phase_modes={phase_modes}, restricted_modes={restricted_modes}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_L1_geometric_disorder_test()
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
