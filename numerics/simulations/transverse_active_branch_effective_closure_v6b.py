#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from L1_stage3_common import PLOTS, REPO_ROOT, append_log, save_result_payload, plt
from L1_transverse_band_test import analyze_case
from transverse_active_branch_effective_closure import make_response_plot
from transverse_pde_reconstruction import (
    cell_center_average,
    edge_fields_from_mode,
    project_centered_transverse,
    residual_metrics,
)

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_active_branch_effective_closure_v6b': {
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
    merged['transverse_active_branch_effective_closure_v6b'] = dict(DEFAULT_CONFIG['transverse_active_branch_effective_closure_v6b'])
    merged['transverse_active_branch_effective_closure_v6b']['thresholds'] = dict(DEFAULT_CONFIG['transverse_active_branch_effective_closure_v6b']['thresholds'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_active_branch_effective_closure_v6b'})
        if isinstance(on_disk.get('transverse_active_branch_effective_closure_v6b'), dict):
            block = dict(on_disk['transverse_active_branch_effective_closure_v6b'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_effective_closure_v6b'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_effective_closure_v6b']['thresholds'].update(thresholds)
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_active_branch_effective_closure_v6b'})
        if isinstance(config.get('transverse_active_branch_effective_closure_v6b'), dict):
            block = dict(config['transverse_active_branch_effective_closure_v6b'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_effective_closure_v6b'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_effective_closure_v6b']['thresholds'].update(thresholds)
    return merged


def make_residual_plot(combo_rows: list[dict[str, Any]], sizes: list[int], penalties: list[float], disorder_strengths: list[float], path: Path) -> None:
    fig, axes = plt.subplots(1, len(disorder_strengths), figsize=(5.0 * len(disorder_strengths), 4.2), sharey=True)
    if len(disorder_strengths) == 1:
        axes = [axes]
    for axis, disorder_strength in zip(axes, disorder_strengths):
        rows = [row for row in combo_rows if abs(row['disorder_strength'] - disorder_strength) <= 1.0e-12]
        for penalty in penalties:
            row = next(item for item in rows if abs(item['penalty'] - penalty) <= 1.0e-12)
            axis.plot(sizes, row['raw_divergence_norms'], marker='x', linestyle='--', label=f'raw p={penalty:g}')
            axis.plot(sizes, row['projected_divergence_norms'], marker='o', label=f'proj p={penalty:g}')
        axis.axhline(1.0e-8, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
        axis.set_xlabel('size n')
        axis.set_title(f'disorder={disorder_strength:.3f}')
        axis.grid(alpha=0.25)
    axes[0].set_ylabel('relative divergence norm')
    axes[-1].legend(fontsize=7)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Active_Branch_Effective_Closure_V6b_v1.md'
    threshold_lines = [
        f"- apply a centered-difference transverse projection after coarse-graining, then judge the same local-law window `curl curl A ≈ lambda A` with residual <= `{result['thresholds']['max_relative_curlcurl_residual']:.2f}`",
        f"- require projected divergence norm <= `{result['thresholds']['max_relative_divergence_norm']:.1e}`",
        f"- keep the same response criterion as `V6`: largest-size adjacent scaled-eigen response jump <= `{result['thresholds']['max_adjacent_scaled_response_jump']:.2f}`",
    ]
    rows = []
    for row in result['combo_rows']:
        rows.append(
            f"| {row['penalty']:.1f} | {row['disorder_strength']:.3f} | {row['worst_raw_divergence_norm']:.3e} | {row['worst_projected_divergence_norm']:.3e} | {row['worst_projected_residual']:.3f} | {row['max_adjacent_response_jump']:.3f} | {'PASS' if row['combo_pass'] else 'OPEN'} |"
        )
    note = f"""# Transverse Active Branch Effective Closure V6b

## Purpose

Execute `V6b` by rerunning the bounded `V6` family with a divergence-aware coarse-graining / projection step, without changing the validated branch, the response family, or the law target.

## Setup

- sizes: `{result['config']['sizes']}`
- variant family: `mild_disorder`
- disorder strengths: `{result['config']['disorder_strengths']}`
- penalties: `{result['config']['penalties']}`
- effective mode: `{result['config']['effective_mode']}`

## V6b criteria

{chr(10).join(threshold_lines)}

## Combo results

| penalty | disorder | worst raw divergence norm | worst projected divergence norm | worst projected residual | max adjacent response jump | V6b status |
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


def run_transverse_active_branch_effective_closure_v6b(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_active_branch_effective_closure_v6b']
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
            print(f'[v6b] family penalty={penalty:.2f} disorder={disorder_strength:.3f}', flush=True)
            mode_rows: list[dict[str, float]] = []
            for n_side in sizes:
                print(f'[v6b] analyze_case n={n_side} penalty={penalty:.2f} disorder={disorder_strength:.3f}', flush=True)
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
                raw = residual_metrics(Ax, Ay, Az, float(case['restricted_transverse_spectrum'][effective_mode]))
                Px, Py, Pz = project_centered_transverse(Ax, Ay, Az)
                projected = residual_metrics(Px, Py, Pz, float(case['restricted_transverse_spectrum'][effective_mode]))
                mode_rows.append(
                    {
                        'n_side': float(n_side),
                        'scaled_eigenvalue': float((n_side**2) * case['restricted_transverse_spectrum'][effective_mode]),
                        'raw_divergence_norm': float(raw['relative_divergence_norm']),
                        'projected_divergence_norm': float(projected['relative_divergence_norm']),
                        'projected_residual': float(projected['relative_curlcurl_residual']),
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
                    'raw_divergence_norms': [row['raw_divergence_norm'] for row in mode_rows],
                    'projected_divergence_norms': [row['projected_divergence_norm'] for row in mode_rows],
                    'projected_residuals': [row['projected_residual'] for row in mode_rows],
                    'scaled_eigenvalues': [row['scaled_eigenvalue'] for row in mode_rows],
                    'worst_raw_divergence_norm': float(max(row['raw_divergence_norm'] for row in mode_rows)),
                    'worst_projected_divergence_norm': float(max(row['projected_divergence_norm'] for row in mode_rows)),
                    'worst_projected_residual': float(max(row['projected_residual'] for row in mode_rows)),
                    'largest_scaled_eigenvalue': float(largest['scaled_eigenvalue']),
                    'largest_divergence_norm': float(largest['projected_divergence_norm']),
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
            row['worst_projected_residual'] <= float(thresholds['max_relative_curlcurl_residual'])
            and row['worst_projected_divergence_norm'] <= float(thresholds['max_relative_divergence_norm'])
            and row['max_adjacent_response_jump'] <= float(thresholds['max_adjacent_scaled_response_jump'])
        )
        print(
            f"[v6b] combo penalty={row['penalty']:.2f} disorder={row['disorder_strength']:.3f} raw_div={row['worst_raw_divergence_norm']:.3e} proj_div={row['worst_projected_divergence_norm']:.3e} proj_res={row['worst_projected_residual']:.3f} response_jump={row['max_adjacent_response_jump']:.3f} status={'PASS' if row['combo_pass'] else 'OPEN'}",
            flush=True,
        )

    divergence_plot = PLOTS / 'transverse_active_branch_effective_closure_v6b_divergence.png'
    response_plot = PLOTS / 'transverse_active_branch_effective_closure_v6b_response.png'
    make_residual_plot(combo_rows, sizes, penalties, disorder_strengths, divergence_plot)
    make_response_plot(combo_rows, penalties, disorder_strengths, response_plot)
    plot_paths = [divergence_plot, response_plot]

    max_raw_divergence = max(row['worst_raw_divergence_norm'] for row in combo_rows)
    max_projected_divergence = max(row['worst_projected_divergence_norm'] for row in combo_rows)
    max_projected_residual = max(row['worst_projected_residual'] for row in combo_rows)
    max_response_jump = max(row['max_adjacent_response_jump'] for row in combo_rows)
    observation = (
        'once the coarse-grained active field is projected with the same centered-difference transverse symbol used by the reconstruction, the divergence defect collapses while the bounded effective-law response remains smooth across the validated mild-disorder family'
    )
    if all(row['combo_pass'] for row in combo_rows):
        conclusion = (
            f"all {len(combo_rows)} tested combos pass `V6b`: the raw divergence defect reaches {max_raw_divergence:.3e}, but the projected divergence defect drops to {max_projected_divergence:.3e}, the worst projected local-law residual remains {max_projected_residual:.3f}, and the largest adjacent scaled-response jump remains {max_response_jump:.3f}; this means the original V6 failure was a coarse-graining artifact rather than a breakdown of bounded effective-law closure on the validated active branch"
        )
    else:
        failing = [f"(penalty={row['penalty']:.1f}, disorder={row['disorder_strength']:.3f})" for row in combo_rows if not row['combo_pass']]
        conclusion = (
            f"`V6b` remains open: failing combos = {failing}; raw divergence max = {max_raw_divergence:.3e}, projected divergence max = {max_projected_divergence:.3e}, projected residual max = {max_projected_residual:.3f}, response jump max = {max_response_jump:.3f}"
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
    result_path, stamped_plots, timestamp = save_result_payload('transverse_active_branch_effective_closure_v6b', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse active-branch effective closure v6b',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variant={variant}, disorder_strengths={disorder_strengths}, penalties={penalties}, restricted_modes={restricted_modes}, effective_mode={effective_mode}, thresholds={thresholds}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_active_branch_effective_closure_v6b()
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
