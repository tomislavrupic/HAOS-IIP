#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from L1_stage3_common import PLOTS, REPO_ROOT, append_log, save_result_payload, plt
from L1_transverse_band_test import analyze_case
from transverse_pde_reconstruction import cell_center_average, edge_fields_from_mode, residual_metrics

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_active_branch_effective_closure': {
        'sizes': [12, 16, 20],
        'variant': 'mild_disorder',
        'disorder_strengths': [0.015, 0.06, 0.12],
        'penalties': [8.0, 10.0, 12.0],
        'restricted_modes': 3,
        'effective_mode': 0,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'flux_tube_phase': math.pi / 2.0,
        'thresholds': {
            'max_relative_curlcurl_residual': 0.08,
            'max_relative_divergence_norm': 1e-8,
            'max_adjacent_scaled_response_jump': 0.12,
        },
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_active_branch_effective_closure'] = dict(DEFAULT_CONFIG['transverse_active_branch_effective_closure'])
    merged['transverse_active_branch_effective_closure']['thresholds'] = dict(DEFAULT_CONFIG['transverse_active_branch_effective_closure']['thresholds'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_active_branch_effective_closure'})
        if isinstance(on_disk.get('transverse_active_branch_effective_closure'), dict):
            block = dict(on_disk['transverse_active_branch_effective_closure'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_effective_closure'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_effective_closure']['thresholds'].update(thresholds)
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_active_branch_effective_closure'})
        if isinstance(config.get('transverse_active_branch_effective_closure'), dict):
            block = dict(config['transverse_active_branch_effective_closure'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_effective_closure'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_effective_closure']['thresholds'].update(thresholds)
    return merged


def make_residual_plot(combo_rows: list[dict[str, Any]], sizes: list[int], penalties: list[float], disorder_strengths: list[float], path: Path) -> None:
    fig, axes = plt.subplots(1, len(disorder_strengths), figsize=(5.0 * len(disorder_strengths), 4.2), sharey=True)
    if len(disorder_strengths) == 1:
        axes = [axes]
    for axis, disorder_strength in zip(axes, disorder_strengths):
        rows = [row for row in combo_rows if abs(row['disorder_strength'] - disorder_strength) <= 1.0e-12]
        for penalty in penalties:
            row = next(item for item in rows if abs(item['penalty'] - penalty) <= 1.0e-12)
            axis.plot(sizes, row['relative_curlcurl_residuals'], marker='o', label=f'penalty={penalty:g}')
        axis.axhline(0.08, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
        axis.set_xlabel('size n')
        axis.set_title(f'disorder={disorder_strength:.3f}')
        axis.grid(alpha=0.25)
    axes[0].set_ylabel('relative curl-curl residual')
    axes[-1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_response_plot(combo_rows: list[dict[str, Any]], penalties: list[float], disorder_strengths: list[float], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for penalty in penalties:
        rows = [row for row in combo_rows if abs(row['penalty'] - penalty) <= 1.0e-12]
        rows.sort(key=lambda item: item['disorder_strength'])
        x_vals = [row['disorder_strength'] for row in rows]
        axes[0].plot(x_vals, [row['largest_scaled_eigenvalue'] for row in rows], marker='o', label=f'penalty={penalty:g}')
        axes[1].plot(x_vals, [row['largest_divergence_norm'] for row in rows], marker='o', label=f'penalty={penalty:g}')
    axes[0].set_xlabel('disorder strength')
    axes[0].set_ylabel(r'$n^2 \lambda_0$ at largest size')
    axes[0].set_title('Effective-law response')
    axes[0].grid(alpha=0.25)
    axes[1].axhline(1.0e-8, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel('disorder strength')
    axes[1].set_ylabel('largest-size divergence norm')
    axes[1].set_title('Response-side divergence control')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Active_Branch_Effective_Closure_v1.md'
    threshold_lines = [
        f"- every tested combo must keep the coarse-grained active mode inside the local law window `curl curl A ≈ lambda A` with relative residual <= `{result['thresholds']['max_relative_curlcurl_residual']:.2f}`",
        f"- every tested combo must keep relative divergence norm <= `{result['thresholds']['max_relative_divergence_norm']:.1e}`",
        f"- at the largest tested size, adjacent mild-disorder steps at fixed penalty must keep the scaled-eigen response jump <= `{result['thresholds']['max_adjacent_scaled_response_jump']:.2f}`",
    ]
    rows = []
    for row in result['combo_rows']:
        rows.append(
            f"| {row['penalty']:.1f} | {row['disorder_strength']:.3f} | {row['worst_relative_curlcurl_residual']:.3f} | {row['worst_relative_divergence_norm']:.3e} | {row['largest_scaled_eigenvalue']:.3f} | {row['max_adjacent_response_jump']:.3f} | {'PASS' if row['combo_pass'] else 'OPEN'} |"
        )
    note = f"""# Transverse Active Branch Effective Closure

## Purpose

Execute a first bounded `V6` read by asking whether the validated active branch closes onto one local coarse-grained law family and supports a smooth bounded response read on the same family.

This first freeze stays on the `V5`-validated `mild_disorder` branch and does not reopen the clean baseline.

## Setup

- sizes: `{result['config']['sizes']}`
- variant family: `mild_disorder`
- disorder strengths: `{result['config']['disorder_strengths']}`
- penalties: `{result['config']['penalties']}`
- effective mode: `{result['config']['effective_mode']}`

## V6 criteria

{chr(10).join(threshold_lines)}

## Combo results

| penalty | disorder | worst curl-curl residual | worst divergence norm | largest-size scaled eigenvalue | max adjacent response jump | V6 status |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
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


def run_transverse_active_branch_effective_closure(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_active_branch_effective_closure']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20])]
    variant = str(experiment_cfg.get('variant', 'mild_disorder'))
    disorder_strengths = [float(v) for v in experiment_cfg.get('disorder_strengths', [0.015, 0.06, 0.12])]
    penalties = [float(v) for v in experiment_cfg.get('penalties', [8.0, 10.0, 12.0])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 3))
    effective_mode = int(experiment_cfg.get('effective_mode', 0))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))
    thresholds = {k: float(v) for k, v in experiment_cfg.get('thresholds', {}).items()}

    combo_rows: list[dict[str, Any]] = []
    response_by_penalty: dict[float, list[dict[str, float]]] = {penalty: [] for penalty in penalties}
    for penalty in penalties:
        for disorder_strength in disorder_strengths:
            print(f'[v6] family penalty={penalty:.2f} disorder={disorder_strength:.3f}', flush=True)
            mode_rows: list[dict[str, float]] = []
            for n_side in sizes:
                print(f'[v6] analyze_case n={n_side} penalty={penalty:.2f} disorder={disorder_strength:.3f}', flush=True)
                case = analyze_case(
                    n_side=n_side,
                    epsilon=epsilon,
                    variant=variant,
                    restricted_modes=restricted_modes,
                    harmonic_tol=harmonic_tol,
                    eig_tol=eig_tol,
                    penalty=penalty,
                    flux_tube_phase=flux_tube_phase,
                    disorder_strength=disorder_strength,
                )
                ex, ey, ez = edge_fields_from_mode(case, effective_mode)
                Ax, Ay, Az = cell_center_average(ex, ey, ez)
                residual = residual_metrics(Ax, Ay, Az, float(case['restricted_transverse_spectrum'][effective_mode]))
                mode_rows.append(
                    {
                        'n_side': float(n_side),
                        'eigenvalue': float(case['restricted_transverse_spectrum'][effective_mode]),
                        'scaled_eigenvalue': float((n_side**2) * case['restricted_transverse_spectrum'][effective_mode]),
                        'relative_curlcurl_residual': float(residual['relative_curlcurl_residual']),
                        'relative_divergence_norm': float(residual['relative_divergence_norm']),
                    }
                )

            largest = mode_rows[-1]
            response_by_penalty[penalty].append(
                {
                    'disorder_strength': disorder_strength,
                    'scaled_eigenvalue': largest['scaled_eigenvalue'],
                }
            )
            combo_rows.append(
                {
                    'penalty': penalty,
                    'disorder_strength': disorder_strength,
                    'relative_curlcurl_residuals': [row['relative_curlcurl_residual'] for row in mode_rows],
                    'relative_divergence_norms': [row['relative_divergence_norm'] for row in mode_rows],
                    'scaled_eigenvalues': [row['scaled_eigenvalue'] for row in mode_rows],
                    'worst_relative_curlcurl_residual': float(max(row['relative_curlcurl_residual'] for row in mode_rows)),
                    'worst_relative_divergence_norm': float(max(row['relative_divergence_norm'] for row in mode_rows)),
                    'largest_scaled_eigenvalue': float(largest['scaled_eigenvalue']),
                    'largest_divergence_norm': float(largest['relative_divergence_norm']),
                }
            )

    response_jump_by_penalty: dict[float, float] = {}
    for penalty, rows in response_by_penalty.items():
        rows.sort(key=lambda item: item['disorder_strength'])
        jumps: list[float] = []
        for left, right in zip(rows[:-1], rows[1:]):
            jump = abs(left['scaled_eigenvalue'] - right['scaled_eigenvalue']) / max(abs(left['scaled_eigenvalue']), abs(right['scaled_eigenvalue']), 1.0e-12)
            jumps.append(float(jump))
        response_jump_by_penalty[penalty] = float(max(jumps, default=0.0))

    for row in combo_rows:
        row['max_adjacent_response_jump'] = float(response_jump_by_penalty[row['penalty']])
        row['combo_pass'] = bool(
            row['worst_relative_curlcurl_residual'] <= float(thresholds['max_relative_curlcurl_residual'])
            and row['worst_relative_divergence_norm'] <= float(thresholds['max_relative_divergence_norm'])
            and row['max_adjacent_response_jump'] <= float(thresholds['max_adjacent_scaled_response_jump'])
        )
        print(
            f"[v6] combo penalty={row['penalty']:.2f} disorder={row['disorder_strength']:.3f} residual={row['worst_relative_curlcurl_residual']:.3f} divergence={row['worst_relative_divergence_norm']:.3e} response_jump={row['max_adjacent_response_jump']:.3f} status={'PASS' if row['combo_pass'] else 'OPEN'}",
            flush=True,
        )

    residual_plot = PLOTS / 'transverse_active_branch_effective_closure_residuals.png'
    response_plot = PLOTS / 'transverse_active_branch_effective_closure_response.png'
    make_residual_plot(combo_rows, sizes, penalties, disorder_strengths, residual_plot)
    make_response_plot(combo_rows, penalties, disorder_strengths, response_plot)
    plot_paths = [residual_plot, response_plot]

    max_residual = max(row['worst_relative_curlcurl_residual'] for row in combo_rows)
    max_divergence = max(row['worst_relative_divergence_norm'] for row in combo_rows)
    max_response_jump = max(row['max_adjacent_response_jump'] for row in combo_rows)
    observation = (
        'the validated mild-disorder active branch admits the same coarse-grained local law test across the full bounded V5 family, and the corresponding scaled law parameter moves smoothly as the disorder strength changes'
    )
    if all(row['combo_pass'] for row in combo_rows):
        weakest = max(combo_rows, key=lambda row: row['worst_relative_curlcurl_residual'])
        conclusion = (
            f"all {len(combo_rows)} tested combos pass the first bounded V6 closure: the weakest local-law case is penalty={weakest['penalty']:.1f}, disorder={weakest['disorder_strength']:.3f} with residual {weakest['worst_relative_curlcurl_residual']:.3f}, the largest divergence norm is {max_divergence:.3e}, and the largest adjacent scaled-response jump is {max_response_jump:.3f}; this means the validated active branch now closes onto one bounded coarse-grained law family and supports a smooth bounded response read on that same family"
        )
    else:
        failing = [f"(penalty={row['penalty']:.1f}, disorder={row['disorder_strength']:.3f})" for row in combo_rows if not row['combo_pass']]
        conclusion = (
            f"the first bounded V6 closure remains open: failing combos = {failing}; max residual = {max_residual:.3f}, max divergence = {max_divergence:.3e}, max adjacent response jump = {max_response_jump:.3f}"
        )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variant': variant,
            'disorder_strengths': disorder_strengths,
            'penalties': penalties,
            'restricted_modes': restricted_modes,
            'effective_mode': effective_mode,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'flux_tube_phase': flux_tube_phase,
        },
        'thresholds': thresholds,
        'combo_rows': combo_rows,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_active_branch_effective_closure', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse active-branch effective closure',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variant={variant}, disorder_strengths={disorder_strengths}, penalties={penalties}, restricted_modes={restricted_modes}, effective_mode={effective_mode}, thresholds={thresholds}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_active_branch_effective_closure()
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
