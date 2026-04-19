#!/usr/bin/env python3

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage3_common import PLOTS, REPO_ROOT, append_log, save_result_payload, plt
from L1_transverse_band_test import analyze_case
from transverse_pde_reconstruction import cell_center_average, curl, edge_fields_from_mode, project_centered_transverse

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_clean_holonomy_symmetry_labels': {
        'sizes': [12, 16, 20],
        'holonomy_phase': 0.05,
        'restricted_modes': 8,
        'active_window_start': 2,
        'window_size': 4,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_clean_holonomy_symmetry_labels'] = dict(DEFAULT_CONFIG['transverse_clean_holonomy_symmetry_labels'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_clean_holonomy_symmetry_labels'})
        if isinstance(on_disk.get('transverse_clean_holonomy_symmetry_labels'), dict):
            merged['transverse_clean_holonomy_symmetry_labels'].update(on_disk['transverse_clean_holonomy_symmetry_labels'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_clean_holonomy_symmetry_labels'})
        if isinstance(config.get('transverse_clean_holonomy_symmetry_labels'), dict):
            merged['transverse_clean_holonomy_symmetry_labels'].update(config['transverse_clean_holonomy_symmetry_labels'])
    return merged


def signed_index(index: int, n_side: int) -> int:
    return index if index <= n_side // 2 else index - n_side


def real_projected_mode_fields(case: dict[str, Any], mode_index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ex, ey, ez = edge_fields_from_mode(case, mode_index)
    Ax, Ay, Az = cell_center_average(ex, ey, ez)
    return project_centered_transverse(Ax, Ay, Az)


def fft_energy(Ax: np.ndarray, Ay: np.ndarray, Az: np.ndarray) -> np.ndarray:
    return np.abs(np.fft.fftn(Ax)) ** 2 + np.abs(np.fft.fftn(Ay)) ** 2 + np.abs(np.fft.fftn(Az)) ** 2


def energy_fractions_by_k_axis(energy: np.ndarray) -> dict[str, float]:
    n_side = int(energy.shape[0])
    weights = {'x': 0.0, 'y': 0.0, 'z': 0.0}
    total = float(np.sum(energy)) or 1.0
    for index, value in np.ndenumerate(energy):
        if index == (0, 0, 0):
            continue
        k_vec = tuple(signed_index(idx, n_side) for idx in index)
        axis = 'xyz'[int(np.argmax(np.abs(k_vec)))]
        weights[axis] += float(value) / total
    return weights


def first_shell_fraction(energy: np.ndarray) -> float:
    n_side = int(energy.shape[0])
    total = float(np.sum(energy)) or 1.0
    shell_weight = 0.0
    for index, value in np.ndenumerate(energy):
        if index == (0, 0, 0):
            continue
        k_vec = tuple(signed_index(idx, n_side) for idx in index)
        shell = tuple(sorted(abs(v) for v in k_vec))
        if shell == (0, 0, 1):
            shell_weight += float(value)
    return float(shell_weight / total)


def dominant_momentum_metrics(energy: np.ndarray) -> dict[str, Any]:
    energy = np.array(energy, copy=True)
    n_side = int(energy.shape[0])
    energy[0, 0, 0] = 0.0
    dominant_index = np.unravel_index(int(np.argmax(energy)), energy.shape)
    k_vec = tuple(int(signed_index(idx, n_side)) for idx in dominant_index)
    shell = tuple(sorted(abs(v) for v in k_vec))
    axis = 'xyz'[int(np.argmax(np.abs(k_vec)))]
    total = float(np.sum(energy)) or 1.0
    return {
        'dominant_momentum_vector': k_vec,
        'dominant_momentum_axis': axis,
        'dominant_shell': shell,
        'dominant_fraction': float(energy[dominant_index] / total),
    }


def dominant_axis_from_field(Ax: np.ndarray, Ay: np.ndarray, Az: np.ndarray) -> str:
    energy = {
        'x': float(np.sum(np.abs(Ax) ** 2)),
        'y': float(np.sum(np.abs(Ay) ** 2)),
        'z': float(np.sum(np.abs(Az) ** 2)),
    }
    return max(energy, key=energy.get)


def mode_label(mode_metrics: dict[str, Any]) -> str:
    return f"k{mode_metrics['dominant_momentum_axis']}-A{mode_metrics['dominant_field_axis']}-curl{mode_metrics['dominant_curl_axis']}"


def analyze_window(case: dict[str, Any], start: int, window_size: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    mode_rows: list[dict[str, Any]] = []
    aggregate_fft = None
    aggregate_field = {'x': 0.0, 'y': 0.0, 'z': 0.0}
    aggregate_curl = {'x': 0.0, 'y': 0.0, 'z': 0.0}
    for mode_index in range(start, start + window_size):
        Ax, Ay, Az = real_projected_mode_fields(case, mode_index)
        Cx, Cy, Cz = curl(Ax, Ay, Az)
        energy = fft_energy(Ax, Ay, Az)
        aggregate_fft = energy if aggregate_fft is None else aggregate_fft + energy
        aggregate_field['x'] += float(np.sum(np.abs(Ax) ** 2))
        aggregate_field['y'] += float(np.sum(np.abs(Ay) ** 2))
        aggregate_field['z'] += float(np.sum(np.abs(Az) ** 2))
        aggregate_curl['x'] += float(np.sum(np.abs(Cx) ** 2))
        aggregate_curl['y'] += float(np.sum(np.abs(Cy) ** 2))
        aggregate_curl['z'] += float(np.sum(np.abs(Cz) ** 2))
        momentum = dominant_momentum_metrics(energy)
        row = {
            'mode_index': int(mode_index),
            'eigenvalue': float(case['restricted_transverse_spectrum'][mode_index]),
            'coexact_fraction': float(case['restricted_transverse_modes'][mode_index]['coexact_fraction']),
            'harmonic_fraction': float(case['restricted_transverse_modes'][mode_index]['harmonic_fraction']),
            'divergence_norm': float(case['restricted_transverse_modes'][mode_index]['divergence_norm']),
            'curl_norm': float(case['restricted_transverse_modes'][mode_index]['curl_norm']),
            **momentum,
            'dominant_field_axis': dominant_axis_from_field(Ax, Ay, Az),
            'dominant_curl_axis': dominant_axis_from_field(Cx, Cy, Cz),
        }
        row['family_label'] = mode_label(row)
        mode_rows.append(row)

    total_field = sum(aggregate_field.values()) or 1.0
    total_curl = sum(aggregate_curl.values()) or 1.0
    axis_weights = energy_fractions_by_k_axis(aggregate_fft if aggregate_fft is not None else np.zeros((1, 1, 1), dtype=float))
    summary = {
        'first_shell_fraction': first_shell_fraction(aggregate_fft if aggregate_fft is not None else np.zeros((1, 1, 1), dtype=float)),
        'k_axis_fractions': axis_weights,
        'field_axis_fractions': {axis: float(value / total_field) for axis, value in aggregate_field.items()},
        'curl_axis_fractions': {axis: float(value / total_curl) for axis, value in aggregate_curl.items()},
        'family_counts': dict(Counter(row['family_label'] for row in mode_rows)),
        'family_sequence': [row['family_label'] for row in mode_rows],
        'scaled_eigenvalues': [float((case['config']['n_side'] ** 2) * row['eigenvalue']) for row in mode_rows],
    }
    return mode_rows, summary


def make_fraction_plot(rows: list[dict[str, Any]], path: Path) -> None:
    n_values = [row['n_side'] for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    axes[0].plot(n_values, [row['summary']['first_shell_fraction'] for row in rows], marker='o', label='first-shell fraction')
    axes[0].plot(n_values, [row['summary']['field_axis_fractions']['x'] for row in rows], marker='s', label='x-field fraction')
    axes[0].set_xlabel('lattice side n')
    axes[0].set_ylabel('window fraction')
    axes[0].set_ylim(0.0, 1.02)
    axes[0].set_title('Resolved-window shell and field support')
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    axes[1].plot(n_values, [row['summary']['k_axis_fractions']['y'] for row in rows], marker='o', label='k-axis y fraction')
    axes[1].plot(n_values, [row['summary']['k_axis_fractions']['z'] for row in rows], marker='s', label='k-axis z fraction')
    axes[1].plot(n_values, [row['summary']['curl_axis_fractions']['y'] for row in rows], marker='^', label='curl-axis y fraction')
    axes[1].plot(n_values, [row['summary']['curl_axis_fractions']['z'] for row in rows], marker='d', label='curl-axis z fraction')
    axes[1].set_xlabel('lattice side n')
    axes[1].set_ylabel('window fraction')
    axes[1].set_ylim(0.0, 1.02)
    axes[1].set_title('Resolved-window family balance')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_scaled_eigen_plot(rows: list[dict[str, Any]], path: Path) -> None:
    family_colors = {
        'ky-Ax-curlz': '#1f77b4',
        'kz-Ax-curly': '#d62728',
    }
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for row in rows:
        n_side = row['n_side']
        for mode in row['mode_rows']:
            ax.scatter(
                n_side,
                (n_side ** 2) * mode['eigenvalue'],
                s=60,
                color=family_colors.get(mode['family_label'], '#444444'),
                alpha=0.85,
            )
    handles = [
        plt.Line2D([0], [0], marker='o', linestyle='', color=color, label=label)
        for label, color in family_colors.items()
    ]
    ax.set_xlabel('lattice side n')
    ax.set_ylabel(r'scaled eigenvalue $n^2 \lambda$')
    ax.set_title('Resolved holonomy-window scaled eigenvalues by family label')
    ax.grid(alpha=0.25)
    ax.legend(handles=handles, fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Clean_Holonomy_Symmetry_Labels_v1.md'
    summary_rows = []
    sequence_rows = []
    for row in result['rows']:
        summary_rows.append(
            f"| {row['n_side']} | {row['summary']['first_shell_fraction']:.3f} | {row['summary']['field_axis_fractions']['x']:.3f} | {row['summary']['k_axis_fractions']['y']:.3f} | {row['summary']['k_axis_fractions']['z']:.3f} | {row['summary']['curl_axis_fractions']['y']:.3f} | {row['summary']['curl_axis_fractions']['z']:.3f} | `{row['summary']['family_counts']}` |"
        )
        sequence_rows.append(f"- `n = {row['n_side']}`: `{row['summary']['family_sequence']}`")
    note = f"""# Transverse Clean Holonomy Symmetry Labels

## Purpose

Execute `CB3` in bounded form on the already-resolved clean holonomy family. The goal is not a full point-group irrep package, but a stable momentum-family / circulation-axis label on the rescued active coexact window.

## Setup

- branch: `baseline`
- sizes: `{result['config']['sizes']}`
- holonomy phase: `{result['config']['holonomy_phase']}`
- restricted modes: `{result['config']['restricted_modes']}`
- resolved active window: start `{result['config']['active_window_start']}`, size `{result['config']['window_size']}`

## Window summary

| `n` | first-shell fraction | `A_x` fraction | `k_y` fraction | `k_z` fraction | `curl_y` fraction | `curl_z` fraction | family counts |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
{chr(10).join(summary_rows)}

## Family sequences

{chr(10).join(sequence_rows)}

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


def run_transverse_clean_holonomy_symmetry_labels(
    config: dict[str, Any] | None = None,
    config_path: Path | None = None,
) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_clean_holonomy_symmetry_labels']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20])]
    holonomy_phase = float(experiment_cfg.get('holonomy_phase', 0.05))
    restricted_modes = int(experiment_cfg.get('restricted_modes', 8))
    active_window_start = int(experiment_cfg.get('active_window_start', 2))
    window_size = int(experiment_cfg.get('window_size', 4))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))

    rows: list[dict[str, Any]] = []
    for n_side in sizes:
        print(f'[cb3] analyze_case holonomy={holonomy_phase:.3f} n={n_side}', flush=True)
        case = analyze_case(
            n_side=n_side,
            epsilon=epsilon,
            variant='baseline',
            restricted_modes=restricted_modes,
            harmonic_tol=harmonic_tol,
            eig_tol=eig_tol,
            penalty=penalty,
            flux_tube_phase=0.0,
            cycle_holonomy_phase=holonomy_phase,
        )
        mode_rows, summary = analyze_window(case, active_window_start, window_size)
        rows.append(
            {
                'n_side': int(n_side),
                'mode_rows': mode_rows,
                'summary': summary,
            }
        )
        print(
            f"[cb3] n={n_side} shell={summary['first_shell_fraction']:.3f} "
            f"Ax={summary['field_axis_fractions']['x']:.3f} "
            f"ky={summary['k_axis_fractions']['y']:.3f} "
            f"kz={summary['k_axis_fractions']['z']:.3f}",
            flush=True,
        )

    first_shell_min = min(row['summary']['first_shell_fraction'] for row in rows)
    x_field_min = min(row['summary']['field_axis_fractions']['x'] for row in rows)
    ky_values = [row['summary']['k_axis_fractions']['y'] for row in rows]
    kz_values = [row['summary']['k_axis_fractions']['z'] for row in rows]
    family_sets_stable = all(set(row['summary']['family_counts'].keys()) == set(rows[0]['summary']['family_counts'].keys()) for row in rows[1:])
    raw_sequences_shuffle = len({tuple(row['summary']['family_sequence']) for row in rows}) > 1

    fractions_plot = PLOTS / 'transverse_clean_holonomy_symmetry_labels_fractions.png'
    eigen_plot = PLOTS / 'transverse_clean_holonomy_symmetry_labels_scaled_eigenvalues.png'
    make_fraction_plot(rows, fractions_plot)
    make_scaled_eigen_plot(rows, eigen_plot)
    plot_paths = [fractions_plot, eigen_plot]

    observation = (
        'the already-resolved holonomy window can be labeled at the family level by combining first-shell momentum support with coarse field-axis and curl-axis readout, rather than demanding a rigid per-mode identity at exact clean symmetry'
    )
    conclusion = (
        f"the holonomy-resolved clean window carries the same bounded first-shell family across n={sizes}: the first-shell support stays at least {first_shell_min:.3f}, "
        f"the coarse field stays fully x-polarized at {x_field_min:.3f}, and the window-level ky/kz split only ranges from {min(ky_values):.3f}-{max(ky_values):.3f} versus {min(kz_values):.3f}-{max(kz_values):.3f}; "
        f"the raw mode ordering still shuffles={raw_sequences_shuffle}, but the family label set is stable={family_sets_stable}, so the clean-baseline identity problem closes in the bounded sense as a symmetry-degenerate first-shell family rather than a missing branch"
    )

    result = {
        'config': {
            'sizes': sizes,
            'holonomy_phase': holonomy_phase,
            'restricted_modes': restricted_modes,
            'active_window_start': active_window_start,
            'window_size': window_size,
            'penalty': penalty,
        },
        'rows': rows,
        'observation': observation,
        'conclusion': conclusion,
        'family_sets_stable': family_sets_stable,
        'raw_sequences_shuffle': raw_sequences_shuffle,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_clean_holonomy_symmetry_labels', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse clean holonomy symmetry labels',
        config_summary=(
            f"epsilon={epsilon}, sizes={sizes}, holonomy_phase={holonomy_phase}, restricted_modes={restricted_modes}, "
            f"active_window_start={active_window_start}, window_size={window_size}"
        ),
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_clean_holonomy_symmetry_labels()
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
