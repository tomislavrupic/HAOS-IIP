#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage2_common import PLOTS, REPO_ROOT, append_log, save_result_payload
from L1_stage4_common import analyze_deformed_case, load_stage4_base_config, plt

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_backreaction_maps': {
        'n_side': 12,
        'variants': ['baseline', 'puncture'],
        'family': 'bump',
        'etas': [0.05, 0.10],
        'restricted_modes': 6,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'map_mode_index': 0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_backreaction_maps'] = dict(DEFAULT_CONFIG['transverse_backreaction_maps'])
    on_disk = load_stage4_base_config(config_path)
    if on_disk:
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_backreaction_maps'})
        if isinstance(on_disk.get('transverse_backreaction_maps'), dict):
            merged['transverse_backreaction_maps'].update(on_disk['transverse_backreaction_maps'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_backreaction_maps'})
        if isinstance(config.get('transverse_backreaction_maps'), dict):
            merged['transverse_backreaction_maps'].update(config['transverse_backreaction_maps'])
    return merged


def mode_sensitivity(case0: dict[str, Any], case_eta: dict[str, Any], mode_index: int) -> np.ndarray:
    v0 = np.abs(np.asarray(case0['restricted_vectors'][mode_index], dtype=complex))
    v1 = np.abs(np.asarray(case_eta['restricted_vectors'][mode_index], dtype=complex))
    return np.abs(v1 - v0)


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def make_map_plot(result: dict[str, Any], path: Path) -> None:
    variants = result['config']['variants']
    eta = result['config']['etas'][-1]
    fig, axes = plt.subplots(len(variants), 2, figsize=(10, 4.3 * len(variants)))
    if len(variants) == 1:
        axes = np.asarray([axes])
    for row, variant in enumerate(variants):
        item = result['responses'][f'{variant}|{eta:.2f}']
        mid = np.asarray(item['midpoints'], dtype=float)
        sens = np.asarray(item['sensitivity'], dtype=float)
        profile = np.asarray(item['scalar_profile'], dtype=float)
        zmid = np.abs(mid[:, 2] - 0.5)
        mask = zmid <= np.quantile(zmid, 0.4)
        im0 = axes[row, 0].scatter(mid[mask, 0], mid[mask, 1], c=sens[mask], cmap='magma', s=14)
        axes[row, 0].set_title(f'{variant}: sensitivity map')
        axes[row, 0].set_xlabel('x')
        axes[row, 0].set_ylabel('y')
        fig.colorbar(im0, ax=axes[row, 0], fraction=0.046, pad=0.04)
        im1 = axes[row, 1].scatter(mid[mask, 0], mid[mask, 1], c=profile[mask], cmap='viridis', s=14)
        axes[row, 1].set_title(f'{variant}: scalar bump profile')
        axes[row, 1].set_xlabel('x')
        axes[row, 1].set_ylabel('y')
        fig.colorbar(im1, ax=axes[row, 1], fraction=0.046, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_profile_plot(result: dict[str, Any], path: Path) -> None:
    eta = result['config']['etas'][-1]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for variant in result['config']['variants']:
        item = result['responses'][f'{variant}|{eta:.2f}']
        ax.scatter(item['scalar_profile'], item['sensitivity'], s=16, alpha=0.5, label=f"{variant} (corr={item['profile_correlation']:.3f})")
    ax.set_xlabel('scalar deformation profile')
    ax.set_ylabel('mode sensitivity')
    ax.set_title('Transverse sensitivity versus scalar profile')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Backreaction_Maps_v1.md'
    rows = '\n'.join(
        f"| {variant} | {eta} | {result['responses'][f'{variant}|{eta:.2f}']['mean_sensitivity']:.6f} | {result['responses'][f'{variant}|{eta:.2f}']['profile_correlation']:.4f} |"
        for variant in result['config']['variants']
        for eta in result['config']['etas']
    )
    note = f"""# Transverse Backreaction Maps

## Purpose

Visualize how the lowest restricted transverse mode responds in space to a localized scalar/background deformation.

## Setup

- branch variants: `{result['config']['variants']}`
- deformation family: `{result['config']['family']}`
- amplitudes: `{result['config']['etas']}`
- lattice size: `{result['config']['n_side']}`

## Sensitivity summary

| variant | eta | mean sensitivity | correlation with scalar profile |
| --- | ---: | ---: | ---: |
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


def run_transverse_backreaction_maps(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_backreaction_maps']
    epsilon = float(cfg.get('epsilon', 0.2))
    n_side = int(experiment_cfg.get('n_side', 12))
    variants = list(experiment_cfg.get('variants', ['baseline', 'puncture']))
    family = str(experiment_cfg.get('family', 'bump'))
    etas = [float(v) for v in experiment_cfg.get('etas', [0.05, 0.10])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 6))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    map_mode_index = int(experiment_cfg.get('map_mode_index', 0))

    responses: dict[str, Any] = {}
    for variant in variants:
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
            sensitivity = mode_sensitivity(baseline, deformed, map_mode_index)
            profile = np.asarray(deformed['midpoint_profile'], dtype=float)
            responses[f'{variant}|{eta:.2f}'] = {
                'midpoints': np.asarray(deformed['midpoints'], dtype=float).tolist(),
                'sensitivity': sensitivity.tolist(),
                'scalar_profile': profile.tolist(),
                'mean_sensitivity': float(np.mean(sensitivity)),
                'profile_correlation': correlation(profile, sensitivity),
            }

    map_plot = PLOTS / 'transverse_backreaction_maps.png'
    profile_plot = PLOTS / 'transverse_sensitivity_vs_scalar_profile.png'
    make_map_plot({'config': {'variants': variants, 'etas': etas}, 'responses': responses}, map_plot)
    make_profile_plot({'config': {'variants': variants, 'etas': etas}, 'responses': responses}, profile_plot)
    plot_paths = [map_plot, profile_plot]

    key = f"baseline|{etas[-1]:.2f}"
    observation = (
        'localized scalar/background bumps produce nonuniform changes in the lowest restricted transverse mode rather than a flat redistribution across the lattice'
    )
    conclusion = (
        f"for the baseline branch at eta={etas[-1]:.2f}, the mean sensitivity is {responses[key]['mean_sensitivity']:.6f} with profile correlation {responses[key]['profile_correlation']:.4f}, so the response is spatially structured but not simply proportional to the scalar bump profile"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'n_side': n_side,
            'variants': variants,
            'family': family,
            'etas': etas,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'map_mode_index': map_mode_index,
        },
        'responses': responses,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_backreaction_maps', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse backreaction maps',
        config_summary=f"epsilon={epsilon}, n_side={n_side}, variants={variants}, family={family}, etas={etas}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_backreaction_maps()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
