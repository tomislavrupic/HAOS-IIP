#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from L1_stage3_common import PLOTS, REPO_ROOT, append_log, save_result_payload, plt
from L1_transverse_band_test import analyze_case
from transverse_active_branch_identity import evaluate_pair
from transverse_active_branch_purity import evaluate_case as evaluate_purity_case
from transverse_active_branch_scaling import V1_THRESHOLDS, scaling_metrics
from transverse_active_sector_transport import pair_transport_metrics

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_active_branch_universality': {
        'sizes': [12, 16, 20],
        'variant': 'mild_disorder',
        'disorder_strengths': [0.015, 0.06, 0.12],
        'penalties': [8.0, 10.0, 12.0],
        'restricted_modes': 6,
        'purity_modes': 4,
        'transport_modes': 4,
        'probe_n_side': 6,
        'scaling_exponent': 2.0,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'flux_tube_phase': math.pi / 2.0,
        'thresholds': {
            'max_scaling_drift': 0.05,
            'max_exact_fraction': 1e-8,
            'max_harmonic_fraction': 1e-8,
            'min_coexact_fraction': 0.999999,
        },
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_active_branch_universality'] = dict(DEFAULT_CONFIG['transverse_active_branch_universality'])
    merged['transverse_active_branch_universality']['thresholds'] = dict(DEFAULT_CONFIG['transverse_active_branch_universality']['thresholds'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_active_branch_universality'})
        if isinstance(on_disk.get('transverse_active_branch_universality'), dict):
            block = dict(on_disk['transverse_active_branch_universality'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_universality'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_universality']['thresholds'].update(thresholds)
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_active_branch_universality'})
        if isinstance(config.get('transverse_active_branch_universality'), dict):
            block = dict(config['transverse_active_branch_universality'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_universality'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_universality']['thresholds'].update(thresholds)
    return merged


def make_universality_plot(combo_rows: list[dict[str, Any]], disorder_strengths: list[float], penalties: list[float], path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), sharex=True)
    for penalty in penalties:
        rows = [row for row in combo_rows if abs(row['penalty'] - penalty) <= 1.0e-12]
        rows.sort(key=lambda item: item['disorder_strength'])
        x_vals = [row['disorder_strength'] for row in rows]
        axes[0].plot(x_vals, [row['worst_mean_overlap'] for row in rows], marker='o', label=f'penalty={penalty:g}')
        axes[1].plot(x_vals, [row['worst_scaling_drift'] for row in rows], marker='o', label=f'penalty={penalty:g}')
        axes[2].plot(x_vals, [row['worst_harmonic_fraction'] for row in rows], marker='o', label=f'penalty={penalty:g}')
    axes[0].axhline(V1_THRESHOLDS['min_mean_overlap'], color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[0].set_xlabel('disorder strength')
    axes[0].set_ylabel('worst mean overlap')
    axes[0].set_title('V5 identity floor')
    axes[0].grid(alpha=0.25)
    axes[1].axhline(0.05, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel('disorder strength')
    axes[1].set_ylabel('worst scaling drift')
    axes[1].set_title('V5 scaling floor')
    axes[1].grid(alpha=0.25)
    axes[2].axhline(1.0e-8, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[2].set_xlabel('disorder strength')
    axes[2].set_ylabel('worst harmonic fraction')
    axes[2].set_title('V5 purity floor')
    axes[2].grid(alpha=0.25)
    axes[2].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Active_Branch_Universality_v1.md'
    threshold_lines = [
        f"- keep the same active-window contract fixed: `V1` identity thresholds, `V3` purity thresholds, and `V4` scaling exponent `alpha = {result['config']['scaling_exponent']:.2f}`",
        f"- test bounded operator-family variation by changing penalty across `{result['config']['penalties']}`",
        f"- test bounded substrate-family variation by changing mild-disorder strength across `{result['config']['disorder_strengths']}`",
        f"- a combo counts as `PASS` only if every adjacent refinement pair passes `V1`, every tested size passes `V3`, and the fixed-`alpha` scaling drift stays <= `{result['thresholds']['max_scaling_drift']:.2f}`",
    ]
    rows = []
    for row in result['combo_rows']:
        rows.append(
            f"| {row['penalty']:.1f} | {row['disorder_strength']:.3f} | {row['identity_pass_count']}/{row['identity_total_pairs']} | {row['worst_mean_overlap']:.3f} | {row['worst_scaling_drift']:.3f} | {row['worst_harmonic_fraction']:.3e} | {row['worst_coexact_fraction']:.9f} | {'PASS' if row['combo_pass'] else 'OPEN'} |"
        )
    note = f"""# Transverse Active Branch Universality

## Purpose

Execute a first bounded `V5` read by asking whether the current active-branch contract survives a small operator-family and substrate-family variation without changing the contract itself.

This first freeze stays on the already-validated `mild_disorder` branch. It does not reopen the clean baseline or enlarge the physical-sector definition.

## Setup

- sizes: `{result['config']['sizes']}`
- variant family: `mild_disorder`
- disorder strengths: `{result['config']['disorder_strengths']}`
- penalties: `{result['config']['penalties']}`
- restricted modes: `{result['config']['restricted_modes']}`
- transported active window: `{result['config']['transport_modes']}`
- fixed scaling exponent: `{result['config']['scaling_exponent']:.2f}`

## V5 criteria

{chr(10).join(threshold_lines)}

## Combo results

| penalty | disorder | identity passes | worst mean overlap | worst scaling drift | worst harmonic fraction | worst coexact fraction | V5 status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
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


def run_transverse_active_branch_universality(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_active_branch_universality']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20])]
    variant = str(experiment_cfg.get('variant', 'mild_disorder'))
    disorder_strengths = [float(v) for v in experiment_cfg.get('disorder_strengths', [0.015, 0.06, 0.12])]
    penalties = [float(v) for v in experiment_cfg.get('penalties', [8.0, 10.0, 12.0])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 6))
    purity_modes = int(experiment_cfg.get('purity_modes', 4))
    transport_modes = int(experiment_cfg.get('transport_modes', 4))
    probe_n_side = int(experiment_cfg.get('probe_n_side', 6))
    scaling_exponent = float(experiment_cfg.get('scaling_exponent', 2.0))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))
    thresholds = {k: float(v) for k, v in experiment_cfg.get('thresholds', {}).items()}

    combo_rows: list[dict[str, Any]] = []
    for penalty in penalties:
        for disorder_strength in disorder_strengths:
            print(f'[v5] family penalty={penalty:.2f} disorder={disorder_strength:.3f}', flush=True)
            cases: dict[int, dict[str, Any]] = {}
            for n_side in sizes:
                print(f'[v5] analyze_case n={n_side} penalty={penalty:.2f} disorder={disorder_strength:.3f}', flush=True)
                cases[n_side] = analyze_case(
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

            identity_passes = 0
            identity_total = 0
            worst_mean_overlap = math.inf
            worst_scaling_drift = 0.0
            for n_from, n_to in zip(sizes[:-1], sizes[1:]):
                transport = pair_transport_metrics(
                    variant=variant,
                    n_from=n_from,
                    n_to=n_to,
                    case_from=cases[n_from],
                    case_to=cases[n_to],
                    probe_n_side=probe_n_side,
                    transport_modes_count=transport_modes,
                )
                if evaluate_pair(transport, V1_THRESHOLDS)['passed']:
                    identity_passes += 1
                identity_total += 1
                worst_mean_overlap = min(worst_mean_overlap, float(transport['mean_overlap']))
                scaling = scaling_metrics(
                    cases[n_from],
                    cases[n_to],
                    n_from=n_from,
                    n_to=n_to,
                    probe_n_side=probe_n_side,
                    transport_modes_count=transport_modes,
                    exponent=scaling_exponent,
                )
                worst_scaling_drift = max(worst_scaling_drift, float(scaling['max_scaled_eigen_drift']))

            worst_exact_fraction = 0.0
            worst_harmonic_fraction = 0.0
            worst_coexact_fraction = 1.0
            purity_sizes_pass = 0
            for n_side in sizes:
                purity = evaluate_purity_case(cases[n_side]['restricted_transverse_modes'], thresholds, purity_modes)
                if purity['passed']:
                    purity_sizes_pass += 1
                worst_exact_fraction = max(worst_exact_fraction, float(purity['worst_exact_fraction']))
                worst_harmonic_fraction = max(worst_harmonic_fraction, float(purity['worst_harmonic_fraction']))
                worst_coexact_fraction = min(worst_coexact_fraction, float(purity['worst_coexact_fraction']))

            combo_pass = bool(
                identity_passes == identity_total
                and purity_sizes_pass == len(sizes)
                and worst_scaling_drift <= float(thresholds['max_scaling_drift'])
            )
            combo_rows.append(
                {
                    'penalty': penalty,
                    'disorder_strength': disorder_strength,
                    'identity_pass_count': identity_passes,
                    'identity_total_pairs': identity_total,
                    'purity_pass_count': purity_sizes_pass,
                    'purity_total_sizes': len(sizes),
                    'worst_mean_overlap': float(worst_mean_overlap),
                    'worst_scaling_drift': float(worst_scaling_drift),
                    'worst_exact_fraction': float(worst_exact_fraction),
                    'worst_harmonic_fraction': float(worst_harmonic_fraction),
                    'worst_coexact_fraction': float(worst_coexact_fraction),
                    'combo_pass': combo_pass,
                }
            )
            print(
                f"[v5] combo penalty={penalty:.2f} disorder={disorder_strength:.3f} identity={identity_passes}/{identity_total} scaling={worst_scaling_drift:.3f} purity_harmonic={worst_harmonic_fraction:.3e} status={'PASS' if combo_pass else 'OPEN'}",
                flush=True,
            )

    plot_path = PLOTS / 'transverse_active_branch_universality.png'
    make_universality_plot(combo_rows, disorder_strengths, penalties, plot_path)
    plot_paths = [plot_path]

    passing = [row for row in combo_rows if row['combo_pass']]
    strongest = max(combo_rows, key=lambda row: row['worst_mean_overlap'])
    weakest = min(combo_rows, key=lambda row: row['worst_mean_overlap'])
    observation = (
        'the first bounded V5 scan keeps the active-branch contract fixed and asks whether it survives a small family of penalty choices and mild-disorder substrates instead of being re-earned after each parameter change'
    )
    if len(passing) == len(combo_rows):
        conclusion = (
            f"all {len(combo_rows)} tested operator/substrate combos pass the first bounded V5 universality check: the weakest combo is penalty={weakest['penalty']:.1f}, disorder={weakest['disorder_strength']:.3f} with worst mean overlap {weakest['worst_mean_overlap']:.3f}, while the largest fixed-alpha drift over the whole family is {max(row['worst_scaling_drift'] for row in combo_rows):.3f}; this means the current active-branch contract survives a real bounded mild-disorder and penalty family without changing the physical-sector rule, purity rule, or scaling exponent"
        )
    else:
        failing = [f"(penalty={row['penalty']:.1f}, disorder={row['disorder_strength']:.3f})" for row in combo_rows if not row['combo_pass']]
        conclusion = (
            f"the first bounded V5 scan remains open: {len(passing)}/{len(combo_rows)} tested combos pass, while failing combos = {failing}; this means the current active-branch contract is not yet stable across even the first bounded operator/substrate family"
        )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variant': variant,
            'disorder_strengths': disorder_strengths,
            'penalties': penalties,
            'restricted_modes': restricted_modes,
            'purity_modes': purity_modes,
            'transport_modes': transport_modes,
            'probe_n_side': probe_n_side,
            'scaling_exponent': scaling_exponent,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'flux_tube_phase': flux_tube_phase,
        },
        'thresholds': thresholds,
        'combo_rows': combo_rows,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_active_branch_universality', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse active-branch universality',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variant={variant}, disorder_strengths={disorder_strengths}, penalties={penalties}, restricted_modes={restricted_modes}, purity_modes={purity_modes}, transport_modes={transport_modes}, probe_n_side={probe_n_side}, scaling_exponent={scaling_exponent}, thresholds={thresholds}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_active_branch_universality()
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
