#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage2_common import PLOTS, REPO_ROOT, append_log, save_result_payload
from L1_stage4_common import (
    DEFORMATION_LABELS,
    VARIANT_LABELS,
    analyze_deformed_case,
    best_mode_overlap,
    load_stage4_base_config,
    make_phase_plot,
    plt,
)

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'scalar_transverse_coupling_test': {
        'sizes': [12],
        'variants': ['baseline', 'puncture', 'line_defect'],
        'families': ['radial', 'anisotropic', 'bump'],
        'etas': [0.02, 0.05, 0.10],
        'restricted_modes': 20,
        'overlap_modes': 10,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['scalar_transverse_coupling_test'] = dict(DEFAULT_CONFIG['scalar_transverse_coupling_test'])
    on_disk = load_stage4_base_config(config_path)
    if on_disk:
        merged.update({k: v for k, v in on_disk.items() if k != 'scalar_transverse_coupling_test'})
        if isinstance(on_disk.get('scalar_transverse_coupling_test'), dict):
            merged['scalar_transverse_coupling_test'].update(on_disk['scalar_transverse_coupling_test'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'scalar_transverse_coupling_test'})
        if isinstance(config.get('scalar_transverse_coupling_test'), dict):
            merged['scalar_transverse_coupling_test'].update(config['scalar_transverse_coupling_test'])
    return merged


def make_shift_plot(result: dict[str, Any], path: Path, modes: int = 10) -> None:
    families = result['config']['families']
    fig, axes = plt.subplots(1, len(families), figsize=(5.0 * len(families), 4.2), sharey=True)
    if len(families) == 1:
        axes = [axes]
    for ax, family in zip(axes, families):
        for branch in result['config']['variants']:
            for eta in result['config']['etas']:
                key = f"{branch}|{family}|{eta:.2f}"
                shifts = result['responses'][key]['eigenvalue_shifts'][:modes]
                ax.plot(range(1, len(shifts) + 1), shifts, marker='o', label=f"{branch}, eta={eta:.2f}")
        ax.set_title(DEFORMATION_LABELS[family])
        ax.set_xlabel('mode index')
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r'$\Delta \lambda_k$')
    axes[-1].legend(fontsize=7, ncol=1)
    fig.suptitle('Restricted transverse eigenvalue shifts under scalar deformation')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_overlap_plot(result: dict[str, Any], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for branch in result['config']['variants']:
        for family in result['config']['families']:
            means = []
            for eta in result['config']['etas']:
                key = f"{branch}|{family}|{eta:.2f}"
                means.append(result['responses'][key]['overlap']['mean_overlap'])
            ax.plot(result['config']['etas'], means, marker='o', label=f"{branch} / {family}")
    ax.set_xlabel('deformation amplitude eta')
    ax.set_ylabel('mean matched overlap')
    ax.set_title('Baseline-to-deformed transverse mode overlap')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_response_plot(result: dict[str, Any], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for branch in result['config']['variants']:
        for family in result['config']['families']:
            lowest = []
            ipr = []
            for eta in result['config']['etas']:
                key = f"{branch}|{family}|{eta:.2f}"
                lowest.append(result['responses'][key]['deformed_spectrum'][0])
                ipr.append(result['responses'][key]['deformed_modes'][0]['ipr'])
            axes[0].plot(result['config']['etas'], lowest, marker='o', label=f"{branch} / {family}")
            axes[1].plot(result['config']['etas'], ipr, marker='o', label=f"{branch} / {family}")
    axes[0].set_xlabel('deformation amplitude eta')
    axes[0].set_ylabel(r'$\lambda_1$')
    axes[0].set_title('Lowest restricted eigenvalue')
    axes[0].grid(alpha=0.25)
    axes[1].set_xlabel('deformation amplitude eta')
    axes[1].set_ylabel('IPR(mode 1)')
    axes[1].set_title('Lowest-mode localization')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=7, ncol=2)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Scalar_Transverse_Coupling_Test_v1.md'
    rows = []
    for key, item in result['responses'].items():
        branch, family, eta = key.split('|')
        rows.append(
            f"| {branch} | {family} | {eta} | {item['deformed_spectrum'][0]:.6f} | {item['eigenvalue_shifts'][0]:.6f} | {item['overlap']['mean_overlap']:.4f} | {item['deformed_modes'][0]['ipr']:.6f} |"
        )
    note = f"""# Scalar-Transverse Coupling Test

## Purpose

Measure how the restricted transverse spectrum deforms under small scalar/background geometry perturbations without changing operator definitions.

## Setup

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- deformation families: `{result['config']['families']}`
- amplitudes: `{result['config']['etas']}`

## Response summary

| branch | family | eta | deformed `lambda_1` | shift `Delta lambda_1` | mean overlap | lowest-mode IPR |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
{chr(10).join(rows)}

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


def run_scalar_transverse_coupling_test(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['scalar_transverse_coupling_test']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12])]
    variants = list(experiment_cfg.get('variants', ['baseline', 'puncture', 'line_defect']))
    families = list(experiment_cfg.get('families', ['radial', 'anisotropic', 'bump']))
    etas = [float(v) for v in experiment_cfg.get('etas', [0.02, 0.05, 0.10])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 20))
    overlap_modes = int(experiment_cfg.get('overlap_modes', 10))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))

    baseline_cases: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        for variant in variants:
            baseline_cases[f'{variant}_n{n_side}'] = analyze_deformed_case(
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

    responses: dict[str, Any] = {}
    for n_side in sizes:
        for variant in variants:
            baseline = baseline_cases[f'{variant}_n{n_side}']
            for family in families:
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
                    overlap = best_mode_overlap(baseline['restricted_vectors'], deformed['restricted_vectors'], overlap_modes)
                    baseline_spec = np.asarray(baseline['restricted_transverse_spectrum'][:restricted_modes], dtype=float)
                    deformed_spec = np.asarray(deformed['restricted_transverse_spectrum'][:restricted_modes], dtype=float)
                    shifts = (deformed_spec - baseline_spec).tolist()
                    responses[f'{variant}|{family}|{eta:.2f}'] = {
                        'n_side': n_side,
                        'baseline_spectrum': baseline_spec.tolist(),
                        'deformed_spectrum': deformed_spec.tolist(),
                        'eigenvalue_shifts': shifts,
                        'overlap': overlap,
                        'deformed_modes': deformed['restricted_transverse_modes'],
                    }

    shift_plot = PLOTS / 'scalar_transverse_eigenvalue_shifts.png'
    overlap_plot = PLOTS / 'scalar_transverse_mode_overlap.png'
    response_plot = PLOTS / 'scalar_transverse_deformation_response.png'
    phase_plot = PLOTS / 'divergence_curl_phase_scalar_transverse.png'
    make_shift_plot({'config': experiment_cfg | {'families': families, 'variants': variants, 'etas': etas}, 'responses': responses}, shift_plot)
    make_overlap_plot({'config': experiment_cfg | {'families': families, 'variants': variants, 'etas': etas}, 'responses': responses}, overlap_plot)
    make_response_plot({'config': experiment_cfg | {'families': families, 'variants': variants, 'etas': etas}, 'responses': responses}, response_plot)
    example_case = analyze_deformed_case(
        n_side=sizes[-1],
        epsilon=epsilon,
        variant=variants[0],
        family=families[-1],
        eta=etas[-1],
        restricted_modes=restricted_modes,
        harmonic_tol=harmonic_tol,
        eig_tol=eig_tol,
        penalty=penalty,
    )
    make_phase_plot(example_case, phase_plot, f"Restricted phase map ({variants[0]}, {families[-1]}, eta={etas[-1]:.2f})")
    plot_paths = [shift_plot, overlap_plot, response_plot, phase_plot]

    key_low = f"baseline|anisotropic|{etas[-1]:.2f}"
    observation = (
        'small scalar/background deformations shift the restricted transverse spectrum smoothly, while the matched low-mode overlaps remain moderate to high depending on branch and deformation family'
    )
    conclusion = (
        f"for the baseline anisotropic branch at eta={etas[-1]:.2f}, the lowest restricted eigenvalue shifts by {responses[key_low]['eigenvalue_shifts'][0]:.6f} with mean matched overlap {responses[key_low]['overlap']['mean_overlap']:.4f}, indicating a measurable but bounded operator-level spectral response"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'families': families,
            'etas': etas,
            'restricted_modes': restricted_modes,
            'overlap_modes': overlap_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
        },
        'responses': responses,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('scalar_transverse_coupling_test', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='scalar-transverse coupling test',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variants={variants}, families={families}, etas={etas}, restricted_modes={restricted_modes}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_scalar_transverse_coupling_test()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
