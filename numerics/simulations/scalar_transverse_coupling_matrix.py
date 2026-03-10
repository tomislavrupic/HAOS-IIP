#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage2_common import PLOTS, REPO_ROOT, append_log, save_result_payload
from L1_stage4_common import analyze_deformed_case, coupling_matrix, load_stage4_base_config, plt

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'scalar_transverse_coupling_matrix': {
        'n_side': 12,
        'variant': 'baseline',
        'family': 'anisotropic',
        'etas': [0.02, 0.05, 0.10],
        'restricted_modes': 12,
        'matrix_modes': 10,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['scalar_transverse_coupling_matrix'] = dict(DEFAULT_CONFIG['scalar_transverse_coupling_matrix'])
    on_disk = load_stage4_base_config(config_path)
    if on_disk:
        merged.update({k: v for k, v in on_disk.items() if k != 'scalar_transverse_coupling_matrix'})
        if isinstance(on_disk.get('scalar_transverse_coupling_matrix'), dict):
            merged['scalar_transverse_coupling_matrix'].update(on_disk['scalar_transverse_coupling_matrix'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'scalar_transverse_coupling_matrix'})
        if isinstance(config.get('scalar_transverse_coupling_matrix'), dict):
            merged['scalar_transverse_coupling_matrix'].update(config['scalar_transverse_coupling_matrix'])
    return merged


def matrix_metrics(matrix: np.ndarray) -> dict[str, float]:
    abs_m = np.abs(matrix)
    diag = np.sum(np.diag(abs_m))
    total = np.sum(abs_m) or 1.0
    offdiag = total - diag
    nn = float(np.mean([abs_m[i, i + 1] for i in range(abs_m.shape[0] - 1)])) if abs_m.shape[0] > 1 else 0.0
    return {
        'diagonal_fraction': float(diag / total),
        'offdiagonal_fraction': float(offdiag / total),
        'nearest_neighbor_mean': nn,
        'frobenius_norm': float(np.linalg.norm(matrix)),
    }


def make_heatmap(matrix: np.ndarray, eta: float, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 5))
    image = ax.imshow(np.abs(matrix), cmap='magma')
    ax.set_title(f'|C_ij| heatmap (eta={eta:.2f})')
    ax.set_xlabel('j')
    ax.set_ylabel('i')
    fig.colorbar(image, ax=ax, label='|C_ij|')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_mixing_plot(etas: list[float], metrics: dict[str, dict[str, float]], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
    axes[0].plot(etas, [metrics[f'{eta:.2f}']['diagonal_fraction'] for eta in etas], marker='o', label='diagonal fraction')
    axes[0].plot(etas, [metrics[f'{eta:.2f}']['offdiagonal_fraction'] for eta in etas], marker='o', label='off-diagonal fraction')
    axes[0].set_xlabel('deformation amplitude eta')
    axes[0].set_ylabel('matrix weight fraction')
    axes[0].set_title('Coupling-matrix weight distribution')
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    axes[1].plot(etas, [metrics[f'{eta:.2f}']['nearest_neighbor_mean'] for eta in etas], marker='o', label='nearest-neighbor mean')
    axes[1].plot(etas, [metrics[f'{eta:.2f}']['frobenius_norm'] for eta in etas], marker='o', label='Frobenius norm')
    axes[1].set_xlabel('deformation amplitude eta')
    axes[1].set_ylabel('magnitude')
    axes[1].set_title('Mode mixing versus deformation')
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Scalar_Transverse_Coupling_Matrix_v1.md'
    rows = '\n'.join(
        f"| {eta} | {result['metrics'][f'{eta:.2f}']['diagonal_fraction']:.4f} | {result['metrics'][f'{eta:.2f}']['offdiagonal_fraction']:.4f} | {result['metrics'][f'{eta:.2f}']['nearest_neighbor_mean']:.6f} |"
        for eta in result['config']['etas']
    )
    note = f"""# Scalar-Transverse Coupling Matrix

## Purpose

Measure the low-dimensional response matrix `C_ij(eta) = <a_i, (T_eta - T_0) a_j>` for the low restricted transverse basis.

## Setup

- branch: `{result['config']['variant']}`
- deformation family: `{result['config']['family']}`
- lattice size: `{result['config']['n_side']}`
- matrix modes: `{result['config']['matrix_modes']}`

## Matrix summary

| eta | diagonal fraction | off-diagonal fraction | nearest-neighbor mean |
| ---: | ---: | ---: | ---: |
{rows}

## Direct result

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Artifacts

- results: `{result_path.relative_to(REPO_ROOT)}`
- plots: {', '.join(f'`{p}`' for p in stamped_plots)}
- timestamp: `{timestamp}`
"""
    note_path.write_text(note, encoding='utf-8')
    return note_path


def run_scalar_transverse_coupling_matrix(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['scalar_transverse_coupling_matrix']
    epsilon = float(cfg.get('epsilon', 0.2))
    n_side = int(experiment_cfg.get('n_side', 12))
    variant = str(experiment_cfg.get('variant', 'baseline'))
    family = str(experiment_cfg.get('family', 'anisotropic'))
    etas = [float(v) for v in experiment_cfg.get('etas', [0.02, 0.05, 0.10])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 12))
    matrix_modes = int(experiment_cfg.get('matrix_modes', 10))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))

    baseline = analyze_deformed_case(
        n_side=n_side,
        epsilon=epsilon,
        variant=variant,
        family='none',
        eta=0.0,
        restricted_modes=restricted_modes,
        harmonic_tol=harmonic_tol,
        eig_tol=eig_tol,
        penalty=penalty,
    )

    matrices: dict[str, Any] = {}
    metrics: dict[str, Any] = {}
    for eta in etas:
        deformed = analyze_deformed_case(
            n_side=n_side,
            epsilon=epsilon,
            variant=variant,
            family=family,
            eta=eta,
            restricted_modes=restricted_modes,
            harmonic_tol=harmonic_tol,
            eig_tol=eig_tol,
            penalty=penalty,
        )
        matrix = coupling_matrix(baseline, deformed, matrix_modes)
        matrices[f'{eta:.2f}'] = {
            'real': np.real(matrix).tolist(),
            'imag': np.imag(matrix).tolist(),
            'abs': np.abs(matrix).tolist(),
        }
        metrics[f'{eta:.2f}'] = matrix_metrics(matrix)

    heatmap_plot = PLOTS / 'transverse_coupling_matrix_heatmap.png'
    mixing_plot = PLOTS / 'transverse_mode_mixing_vs_eta.png'
    make_heatmap(np.asarray(matrices[f'{etas[-1]:.2f}']['abs'], dtype=float), etas[-1], heatmap_plot)
    make_mixing_plot(etas, metrics, mixing_plot)
    plot_paths = [heatmap_plot, mixing_plot]

    observation = (
        'the deformation-induced coupling matrix remains structured, but its weight is broadly distributed across the low-mode block rather than sharply diagonal-dominant'
    )
    conclusion = (
        f"at eta={etas[-1]:.2f}, the diagonal fraction is {metrics[f'{etas[-1]:.2f}']['diagonal_fraction']:.4f} and the off-diagonal fraction is {metrics[f'{etas[-1]:.2f}']['offdiagonal_fraction']:.4f}, indicating structured but non-diagonal mode mixing in the low transverse sector"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'n_side': n_side,
            'variant': variant,
            'family': family,
            'etas': etas,
            'restricted_modes': restricted_modes,
            'matrix_modes': matrix_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
        },
        'metrics': metrics,
        'matrices': matrices,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('scalar_transverse_coupling_matrix', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='scalar-transverse coupling matrix',
        config_summary=f"epsilon={epsilon}, n_side={n_side}, variant={variant}, family={family}, etas={etas}, matrix_modes={matrix_modes}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_scalar_transverse_coupling_matrix()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
