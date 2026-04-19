#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage3_common import PLOTS, REPO_ROOT, append_log, save_result_payload, plt
from L1_transverse_band_test import analyze_case
from transverse_active_sector_transport import pair_transport_metrics

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_active_sector_disorder_sweep': {
        'sizes': [12, 16, 20, 24, 28],
        'disorder_strengths': [0.0, 0.015, 0.03, 0.06, 0.09, 0.12],
        'restricted_modes': 6,
        'transport_modes': 4,
        'probe_n_side': 6,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': math.pi / 2.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_active_sector_disorder_sweep'] = dict(DEFAULT_CONFIG['transverse_active_sector_disorder_sweep'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_active_sector_disorder_sweep'})
        if isinstance(on_disk.get('transverse_active_sector_disorder_sweep'), dict):
            merged['transverse_active_sector_disorder_sweep'].update(on_disk['transverse_active_sector_disorder_sweep'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_active_sector_disorder_sweep'})
        if isinstance(config.get('transverse_active_sector_disorder_sweep'), dict):
            merged['transverse_active_sector_disorder_sweep'].update(config['transverse_active_sector_disorder_sweep'])
    return merged


def summarize_strength(strength: float, pair_metrics: dict[str, dict[str, Any]], sizes: list[int]) -> dict[str, Any]:
    pairs = [pair_metrics[f'n{n_from}_to_n{n_to}'] for n_from, n_to in zip(sizes[:-1], sizes[1:])]
    last_pair = pair_metrics[f'n{sizes[-2]}_to_n{sizes[-1]}']
    return {
        'disorder_strength': float(strength),
        'mean_pair_overlap': float(np.mean([item['mean_overlap'] for item in pairs])),
        'worst_pair_mean_overlap': float(np.min([item['mean_overlap'] for item in pairs])),
        'worst_pair_min_overlap': float(np.min([item['min_overlap'] for item in pairs])),
        'mean_pair_principal_cosine': float(np.mean([item['mean_principal_cosine'] for item in pairs])),
        'max_pair_scaled_eigen_drift': float(np.max([item['max_scaled_eigen_drift'] for item in pairs])),
        'last_pair_mean_overlap': float(last_pair['mean_overlap']),
        'last_pair_min_overlap': float(last_pair['min_overlap']),
        'last_pair_mean_principal_cosine': float(last_pair['mean_principal_cosine']),
        'last_pair_max_scaled_eigen_drift': float(last_pair['max_scaled_eigen_drift']),
    }


def make_overlap_plot(summary_rows: list[dict[str, Any]], path: Path) -> None:
    strengths = [row['disorder_strength'] for row in summary_rows]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    axes[0].plot(strengths, [row['mean_pair_overlap'] for row in summary_rows], marker='o', label='mean pair overlap')
    axes[0].plot(strengths, [row['worst_pair_mean_overlap'] for row in summary_rows], marker='s', label='worst pair mean overlap')
    axes[0].plot(strengths, [row['worst_pair_min_overlap'] for row in summary_rows], marker='^', label='worst matched overlap')
    axes[0].set_xlabel('disorder strength')
    axes[0].set_ylabel('transport overlap')
    axes[0].set_title('Active-sector overlap versus disorder strength')
    axes[0].set_ylim(0.0, 1.02)
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    axes[1].plot(strengths, [row['last_pair_mean_overlap'] for row in summary_rows], marker='o', label='24->28 mean overlap')
    axes[1].plot(strengths, [row['last_pair_min_overlap'] for row in summary_rows], marker='s', label='24->28 min overlap')
    axes[1].plot(strengths, [row['last_pair_mean_principal_cosine'] for row in summary_rows], marker='D', label='24->28 mean principal cosine')
    axes[1].set_xlabel('disorder strength')
    axes[1].set_ylabel('last-pair transport stability')
    axes[1].set_title('Largest refinement pair versus disorder strength')
    axes[1].set_ylim(0.0, 1.02)
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_drift_plot(summary_rows: list[dict[str, Any]], path: Path) -> None:
    strengths = [row['disorder_strength'] for row in summary_rows]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(strengths, [row['max_pair_scaled_eigen_drift'] for row in summary_rows], marker='o', label='max pair drift')
    ax.plot(strengths, [row['last_pair_max_scaled_eigen_drift'] for row in summary_rows], marker='s', label='24->28 max drift')
    ax.set_xlabel('disorder strength')
    ax.set_ylabel('max relative scaled-eigen drift')
    ax.set_title('Active-sector eigen drift versus disorder strength')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Active_Sector_Disorder_Sweep_v1.md'
    rows = []
    for item in result['summary_rows']:
        rows.append(
            f"| {item['disorder_strength']:.3f} | {item['mean_pair_overlap']:.3f} | {item['worst_pair_mean_overlap']:.3f} | {item['worst_pair_min_overlap']:.3f} | {item['last_pair_mean_overlap']:.3f} | {item['last_pair_mean_principal_cosine']:.3f} | {item['last_pair_max_scaled_eigen_drift']:.3f} |"
        )
    note = f"""# Transverse Active Sector Disorder Sweep

## Purpose

Keep the same active-sector transport map `J_n` fixed and ask whether the clean-torus branch-mixing weakness decays continuously once bounded smooth disorder is introduced.

## Setup

- sizes: `{result['config']['sizes']}`
- disorder strengths: `{result['config']['disorder_strengths']}`
- restricted modes: `{result['config']['restricted_modes']}`
- transported low-window modes: `{result['config']['transport_modes']}`
- probe lattice side: `{result['config']['probe_n_side']}`

## Strength summary

| disorder strength | mean pair overlap | worst pair mean overlap | worst matched overlap | `24->28` mean overlap | `24->28` mean principal cosine | `24->28` max drift |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
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


def run_transverse_active_sector_disorder_sweep(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_active_sector_disorder_sweep']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20, 24, 28])]
    disorder_strengths = [float(v) for v in experiment_cfg.get('disorder_strengths', [0.0, 0.015, 0.03, 0.06, 0.09, 0.12])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 6))
    transport_modes_count = int(experiment_cfg.get('transport_modes', 4))
    probe_n_side = int(experiment_cfg.get('probe_n_side', 6))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))

    all_metrics: dict[str, dict[str, Any]] = {}
    summary_rows: list[dict[str, Any]] = []

    for strength in disorder_strengths:
        print(f'[disorder-sweep] starting strength={strength:.3f}', flush=True)
        cases: dict[str, Any] = {}
        for n_side in sizes:
            cases[f'n{n_side}'] = analyze_case(
                n_side=n_side,
                epsilon=epsilon,
                variant='mild_disorder',
                restricted_modes=restricted_modes,
                harmonic_tol=harmonic_tol,
                eig_tol=eig_tol,
                penalty=penalty,
                flux_tube_phase=flux_tube_phase,
                disorder_strength=strength,
            )
        pair_metrics: dict[str, Any] = {}
        for n_from, n_to in zip(sizes[:-1], sizes[1:]):
            pair_key = f'n{n_from}_to_n{n_to}'
            pair_metrics[pair_key] = pair_transport_metrics(
                variant='mild_disorder',
                n_from=n_from,
                n_to=n_to,
                case_from=cases[f'n{n_from}'],
                case_to=cases[f'n{n_to}'],
                probe_n_side=probe_n_side,
                transport_modes_count=transport_modes_count,
            )
        strength_key = f"{strength:.3f}"
        all_metrics[strength_key] = pair_metrics
        summary = summarize_strength(strength, pair_metrics, sizes)
        summary_rows.append(summary)
        print(
            '[disorder-sweep] completed '
            f"strength={strength:.3f} "
            f"last_pair_overlap={summary['last_pair_mean_overlap']:.3f} "
            f"last_pair_principal_cosine={summary['last_pair_mean_principal_cosine']:.3f} "
            f"last_pair_max_drift={summary['last_pair_max_scaled_eigen_drift']:.3f}",
            flush=True,
        )

    overlap_plot = PLOTS / 'transverse_active_sector_disorder_sweep_overlap.png'
    drift_plot = PLOTS / 'transverse_active_sector_disorder_sweep_drift.png'
    make_overlap_plot(summary_rows, overlap_plot)
    make_drift_plot(summary_rows, drift_plot)
    plot_paths = [overlap_plot, drift_plot]

    stabilization = next((row for row in summary_rows if row['last_pair_mean_overlap'] >= 0.9), None)
    weakest = min(summary_rows, key=lambda row: row['last_pair_mean_overlap'])
    strongest = max(summary_rows, key=lambda row: row['last_pair_mean_overlap'])
    observation = (
        'with the active-sector comparison map held fixed, the disorder sweep isolates whether branch identity improves continuously under bounded symmetry breaking instead of only under topological defects'
    )
    if stabilization is None:
        conclusion = (
            f"across the tested strengths, the largest refinement pair improves from {weakest['last_pair_mean_overlap']:.3f} at strength {weakest['disorder_strength']:.3f} to {strongest['last_pair_mean_overlap']:.3f} at strength {strongest['disorder_strength']:.3f}, but it never crosses a 0.9 mean-overlap stabilization threshold in this window; that would mean the clean-baseline weakness shrinks under disorder but is not yet fully resolved"
        )
    else:
        conclusion = (
            f"the largest refinement pair rises from mean overlap {weakest['last_pair_mean_overlap']:.3f} at disorder strength {weakest['disorder_strength']:.3f} to {strongest['last_pair_mean_overlap']:.3f} at strength {strongest['disorder_strength']:.3f}, and it first crosses a 0.9 mean-overlap stabilization threshold at strength {stabilization['disorder_strength']:.3f}; because the same `J_n` map is used throughout, this points to a clean-background degeneracy / branch-mixing problem rather than a broken transport construction"
        )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'disorder_strengths': disorder_strengths,
            'restricted_modes': restricted_modes,
            'transport_modes': transport_modes_count,
            'probe_n_side': probe_n_side,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'flux_tube_phase': flux_tube_phase,
        },
        'transport_metrics_by_strength': all_metrics,
        'summary_rows': summary_rows,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_active_sector_disorder_sweep', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse active-sector disorder sweep',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, disorder_strengths={disorder_strengths}, restricted_modes={restricted_modes}, transport_modes={transport_modes_count}, probe_n_side={probe_n_side}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_active_sector_disorder_sweep()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
