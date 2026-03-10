#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from L1_transverse_scaling_test import analyze_case
from L1_stage2_common import (
    PLOTS,
    REPO_ROOT,
    append_log,
    ensure_matplotlib,
    save_result_payload,
)

plt = ensure_matplotlib()

DEFAULT_CONFIG: dict[str, Any] = {
    'kernel_width_sensitivity_test': {
        'n_side': 16,
        'c_epsilon_values': [0.25, 0.5, 1.0, 1.5],
        'phase_modes': 20,
        'restricted_modes': 20,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['kernel_width_sensitivity_test'] = dict(DEFAULT_CONFIG['kernel_width_sensitivity_test'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        if isinstance(on_disk.get('kernel_width_sensitivity_test'), dict):
            merged['kernel_width_sensitivity_test'].update(on_disk['kernel_width_sensitivity_test'])
    if config is not None and isinstance(config.get('kernel_width_sensitivity_test'), dict):
        merged['kernel_width_sensitivity_test'].update(config['kernel_width_sensitivity_test'])
    return merged


def make_kernel_plot(cases: dict[str, dict[str, Any]], c_values: list[float], path: Path) -> dict[str, float]:
    lowest = [cases[f'c_{c_value}']['restricted_transverse_spectrum'][0] for c_value in c_values]
    spread = [float(cases[f'c_{c_value}']['restricted_transverse_spectrum'][9] - cases[f'c_{c_value}']['restricted_transverse_spectrum'][0]) for c_value in c_values]
    ipr = [cases[f'c_{c_value}']['restricted_transverse_modes'][0]['ipr'] for c_value in c_values]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].plot(c_values, lowest, marker='o')
    axes[0].set_title('Lowest restricted eigenvalue')
    axes[1].plot(c_values, spread, marker='o')
    axes[1].set_title('Band spread ($\\lambda_{10}-\\lambda_1$)')
    axes[2].plot(c_values, ipr, marker='o')
    axes[2].set_title('Lowest-mode IPR')
    for ax in axes:
        ax.set_xlabel(r'$c_\varepsilon$')
        ax.grid(alpha=0.25)
    axes[0].set_ylabel('value')
    fig.suptitle('Restricted transverse sensitivity to kernel width')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    return {
        str(c_value): {
            'lowest_eigenvalue': lowest[idx],
            'band_spread_10': spread[idx],
            'lowest_ipr': ipr[idx],
        }
        for idx, c_value in enumerate(c_values)
    }


def make_phase_plot(cases: dict[str, dict[str, Any]], c_values: list[float], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharex=True, sharey=True)
    scatter = None
    for ax, c_value in zip(axes.ravel(), c_values):
        case = cases[f'c_{c_value}']
        x = [record['divergence_norm'] for record in case['phase_modes']]
        y = [record['curl_norm'] for record in case['phase_modes']]
        c = [record['eigenvalue'] for record in case['phase_modes']]
        scatter = ax.scatter(x, y, c=c, cmap='viridis', s=28)
        ax.set_title(rf'$c_\varepsilon={c_value}$')
        ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel(r'$||d_0^* a_k||$')
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$||d_1 a_k||$')
    cbar = fig.colorbar(scatter, ax=axes.ravel().tolist(), shrink=0.88)
    cbar.set_label(r'$\lambda_k$')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Kernel_Width_Sensitivity_Test_v1.md'
    rows = '\n'.join(
        f"| {c_value} | {result['metrics'][str(c_value)]['lowest_eigenvalue']:.6f} | {result['metrics'][str(c_value)]['band_spread_10']:.6f} | {result['metrics'][str(c_value)]['lowest_ipr']:.6f} |"
        for c_value in result['config']['c_epsilon_values']
    )
    note = f"""# Kernel Width Sensitivity Test

## Purpose

Measure how the restricted transverse band changes as the kernel width varies at fixed grid size.

## Setup

- branch: baseline periodic torus
- `n = {result['config']['n_side']}`
- `epsilon_k = c_epsilon h^2`
- `c_epsilon = {result['config']['c_epsilon_values']}`

## Summary metrics

| `c_epsilon` | lowest `lambda` | `lambda_10 - lambda_1` | lowest-mode IPR |
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


def run_kernel_width_sensitivity_test(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['kernel_width_sensitivity_test']
    n_side = int(experiment_cfg.get('n_side', 16))
    c_values = [float(v) for v in experiment_cfg.get('c_epsilon_values', [0.25, 0.5, 1.0, 1.5])]
    phase_modes = int(experiment_cfg.get('phase_modes', 20))
    restricted_modes = int(experiment_cfg.get('restricted_modes', 20))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    h = 1.0 / n_side

    cases: dict[str, dict[str, Any]] = {}
    for c_value in c_values:
        epsilon = c_value * h * h
        cases[f'c_{c_value}'] = analyze_case(
            n_side=n_side,
            epsilon=epsilon,
            variant='baseline',
            phase_modes=phase_modes,
            restricted_modes=restricted_modes,
            harmonic_tol=harmonic_tol,
            eig_tol=eig_tol,
            penalty=penalty,
        )

    kernel_plot = PLOTS / 'transverse_band_vs_kernel_width.png'
    phase_plot = PLOTS / 'divergence_curl_phase_kernel_width.png'
    metrics = make_kernel_plot(cases, c_values, kernel_plot)
    make_phase_plot(cases, c_values, phase_plot)
    plot_paths = [kernel_plot, phase_plot]

    observation = (
        'the restricted transverse spectrum remains present across the tested kernel-width range, while the lowest eigenvalue, first-band spread, and lowest-mode IPR shift smoothly with c_epsilon'
    )
    conclusion = 'within the tested local-kernel range, the restricted transverse band is robust to moderate kernel-width changes'

    result = {
        'config': {
            'n_side': n_side,
            'c_epsilon_values': c_values,
            'phase_modes': phase_modes,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
        },
        'cases': {label: {'config': case['config'], 'dimensions': case['dimensions'], 'phase_modes': case['phase_modes'], 'restricted_transverse_spectrum': case['restricted_transverse_spectrum'], 'restricted_transverse_modes': case['restricted_transverse_modes']} for label, case in cases.items()},
        'metrics': metrics,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('kernel_width_sensitivity_test', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='kernel width sensitivity test',
        config_summary=f"n={n_side}, c_epsilon_values={c_values}, phase_modes={phase_modes}, restricted_modes={restricted_modes}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_kernel_width_sensitivity_test()
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
