#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import (
    CSV_FIELDS as BASE_FIELDS,
    WORK_PLOT_DIR,
    load_runsheet,
    lookup_by_id,
    selected_runs,
    simulate_run,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_2_dk_grade_exchange_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_2_DK_Grade_Exchange_Texture_Scan_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'grade_transfer_onset_time',
    'grade_transfer_duty_cycle',
    'omega1_peak_weight',
    'omega2_peak_weight',
    'grade_exchange_signature',
    'texture_stable_under_refinement',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.2 DK grade-exchange and collision-texture scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def onset_time(times: list[float], grade_hist: np.ndarray, amplitude: float) -> float:
    if len(times) != len(grade_hist):
        return 0.0
    threshold = max(0.02, 0.25 * amplitude)
    baseline = grade_hist[0]
    deviations = np.linalg.norm(grade_hist - baseline[None, :], axis=1)
    idx = np.where(deviations >= threshold)[0]
    return float(times[int(idx[0])]) if idx.size else 0.0


def duty_cycle(grade_hist: np.ndarray, amplitude: float) -> float:
    threshold = max(0.02, 0.25 * amplitude)
    baseline = grade_hist[0]
    deviations = np.linalg.norm(grade_hist - baseline[None, :], axis=1)
    return float(np.mean(deviations >= threshold)) if deviations.size else 0.0


def signature_for(row: dict[str, Any]) -> str:
    amp = float(row['grade_transfer_amplitude'])
    label = str(row['collision_label'])
    if label == 'pass-through with grade exchange' and amp >= 0.15:
        return 'strong_pass_through_exchange'
    if label == 'pass-through with grade exchange':
        return 'weak_pass_through_exchange'
    if label == 'unresolved / mixed' and float(row['binding_persistence']) >= 0.5:
        return 'clustered_mixed_encounter'
    if label == 'pass-through dispersive' and amp >= 0.06:
        return 'dispersive_exchange_tail'
    return 'texture_null'



def plot_grade_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_2_', '') for row in rows]
    amps = [float(row['grade_transfer_amplitude']) for row in rows]
    onsets = [float(row['grade_transfer_onset_time']) for row in rows]
    fig, axes = plt.subplots(2, 1, figsize=(11.0, 7.0), sharex=True)
    axes[0].bar(range(len(rows)), amps, color='tab:green')
    axes[0].set_ylabel('grade transfer amplitude')
    axes[0].set_title('Stage 23.2 grade-exchange summary')
    axes[0].grid(axis='y', alpha=0.25)
    axes[1].bar(range(len(rows)), onsets, color='tab:blue')
    axes[1].set_ylabel('onset time')
    axes[1].set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    axes[1].grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)



def plot_phase_sensitivity(path: Path, rows: list[dict[str, Any]]) -> None:
    grouped: dict[tuple[str, int], dict[str, float]] = defaultdict(dict)
    for row in rows:
        grouped[(str(row['geometry_id']), int(row['resolution']))][str(row['phase_id'])] = float(row['grade_transfer_amplitude'])
    items = sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1]))
    labels = [f"{geom.replace('_pair','').replace('_collision','')}\n n={res}" for (geom, res), _ in items]
    deltas = [abs(values.get('in_phase', 0.0) - values.get('out_of_phase', 0.0)) for _, values in items]
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.bar(range(len(items)), deltas, color='tab:orange')
    ax.set_ylabel('|in - out| grade transfer')
    ax.set_xticks(range(len(items)), labels)
    ax.set_title('Phase sensitivity by geometry and resolution')
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)



def plot_texture_stability(path: Path, rows: list[dict[str, Any]]) -> None:
    geometries = ['counter_propagating_corridor_pair', 'offset_glancing_collision', 'tight_clustered_pair']
    phases = ['in_phase', 'out_of_phase']
    resolutions = [12, 24]
    label_to_value = {
        'texture_null': 0,
        'weak_pass_through_exchange': 1,
        'strong_pass_through_exchange': 2,
        'clustered_mixed_encounter': 3,
        'dispersive_exchange_tail': 4,
    }
    grid = np.full((len(geometries), len(phases) * len(resolutions)), -1, dtype=int)
    text_labels: dict[tuple[int, int], str] = {}
    for row in rows:
        i = geometries.index(str(row['geometry_id']))
        phase_idx = phases.index(str(row['phase_id']))
        res_idx = resolutions.index(int(row['resolution']))
        j = phase_idx * len(resolutions) + res_idx
        sig = str(row['grade_exchange_signature'])
        grid[i, j] = label_to_value[sig]
        text_labels[(i, j)] = sig.replace('_', ' ')
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    im = ax.imshow(grid, cmap='viridis', aspect='auto', vmin=0, vmax=max(label_to_value.values()))
    ax.set_yticks(range(len(geometries)), [g.replace('_pair', '').replace('_collision', '') for g in geometries])
    xlabels = [f"{phase}\n n={res}" for phase in phases for res in resolutions]
    ax.set_xticks(range(len(xlabels)), xlabels)
    ax.set_title('DK collision-texture signatures')
    for (i, j), label in text_labels.items():
        ax.text(j, i, label.replace('pass through', 'pass-thru'), ha='center', va='center', color='white', fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)



def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], summary: dict[str, Any]) -> None:
    lines = [
        '# Stage 23.2 DK Grade-Exchange Texture Scan v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f"Total runs: {len(rows)}",
        '',
        '## Collision labels',
    ]
    for label, count in sorted(summary['collision_labels'].items()):
        lines.append(f'- {label}: {count}')
    lines.extend(['', '## Grade-exchange signatures'])
    for label, count in sorted(summary['grade_exchange_signatures'].items()):
        lines.append(f'- {label}: {count}')
    lines.extend([
        '',
        '## Texture families stable under refinement',
        f"- stable families: `{summary['stable_family_count']}`",
        f"- stable keys: `{summary['stable_keys']}`",
        '',
        '## Plots',
    ])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')
    NOTE_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')



def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    runs = selected_runs(runsheet['runs'], args.run_ids)
    geometry_by_id = lookup_by_id(runsheet['geometries'], 'geometry_id')
    fixed = runsheet['fixed_parameters']

    rows: list[dict[str, Any]] = []
    details_by_run: dict[str, dict[str, Any]] = {}
    plot_paths: list[Path] = []

    for run in runs:
        geometry = geometry_by_id[run['geometry_id']]
        row, detail, run_plots = simulate_run(run, geometry, fixed)
        grade_hist = np.asarray(detail['total_grade_weights'], dtype=float)
        times = [float(t) for t in detail['times']]
        amp = float(row['grade_transfer_amplitude'])
        row['grade_transfer_onset_time'] = onset_time(times, grade_hist, amp)
        row['grade_transfer_duty_cycle'] = duty_cycle(grade_hist, amp)
        row['omega1_peak_weight'] = float(np.max(grade_hist[:, 1])) if grade_hist.size else 0.0
        row['omega2_peak_weight'] = float(np.max(grade_hist[:, 2])) if grade_hist.size else 0.0
        row['grade_exchange_signature'] = signature_for(row)
        row['texture_stable_under_refinement'] = 0
        detail['metrics'] = row
        rows.append(row)
        details_by_run[str(row['run_id'])] = detail
        plot_paths.extend(run_plots)

    rows.sort(key=lambda item: item['run_id'])

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row['geometry_id']), str(row['phase_id']))].append(row)
    stable_keys: list[str] = []
    for key, group in grouped.items():
        if len(group) != 2:
            continue
        group = sorted(group, key=lambda item: int(item['resolution']))
        stable = (
            group[0]['grade_exchange_signature'] == group[1]['grade_exchange_signature']
            and group[0]['grade_exchange_signature'] != 'texture_null'
        )
        if stable:
            stable_keys.append(f"{key[0]}::{key[1]}")
        for item in group:
            item['texture_stable_under_refinement'] = int(stable)
    details = []
    for row in rows:
        detail = details_by_run[str(row['run_id'])]
        detail['metrics'] = row
        details.append(detail)

    grade_matrix = WORK_PLOT_DIR / 'stage23_2_grade_exchange_summary.png'
    phase_plot = WORK_PLOT_DIR / 'stage23_2_phase_sensitivity_summary.png'
    texture_plot = WORK_PLOT_DIR / 'stage23_2_texture_stability_matrix.png'
    plot_grade_matrix(grade_matrix, rows)
    plot_phase_sensitivity(phase_plot, rows)
    plot_texture_stability(texture_plot, rows)
    plot_paths.extend([grade_matrix, phase_plot, texture_plot])

    summary = {
        'collision_labels': dict(Counter(row['collision_label'] for row in rows)),
        'persistence_labels': dict(Counter(row['persistence_label'] for row in rows)),
        'grade_exchange_signatures': dict(Counter(row['grade_exchange_signature'] for row in rows)),
        'stable_family_count': len(stable_keys),
        'stable_keys': stable_keys,
        'max_grade_transfer_amplitude': float(max(float(row['grade_transfer_amplitude']) for row in rows)) if rows else 0.0,
    }
    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'fixed_parameters': fixed,
        'runs': details,
        'summary': summary,
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        runsheet['experiment_slug'],
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, summary)


if __name__ == '__main__':
    main()
