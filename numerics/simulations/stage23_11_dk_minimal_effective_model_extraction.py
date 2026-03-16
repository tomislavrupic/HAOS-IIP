#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_8_dk_topology_protection_scan import simulate_run as simulate_observed_run

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_11_dk_minimal_effective_model_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_III_C2_Minimal_Effective_Model_Extraction_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage23_11_minimal_effective_model'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'run_order',
    'phase_offset_fraction_of_pi',
    'detuning_delta',
    'grade0_scale',
    'grade1_scale',
    'observed_topology_label',
    'predicted_topology_label',
    'topology_match_flag',
    'observed_y_peak',
    'predicted_y_peak',
    'observed_x_final',
    'predicted_x_final',
    'x_trace_rmse',
    'y_trace_rmse',
    'observed_braid_survival_proxy',
    'predicted_braid_survival_proxy',
    'survival_time_abs_error',
    'flow_concentration_index',
    'grade_exchange_coherence',
    'notes',
]

TOPOLOGY_CODE = {
    'braid_like_exchange': 0,
    'transfer_smeared': 1,
    'unresolved': 2,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 23.11 DK minimal effective model extraction.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def grade_scales_for_detuning(detuning_delta: float, max_span: float, max_detuning: float) -> tuple[float, float]:
    alpha = float(detuning_delta) / max(max_detuning, 1.0e-12)
    return 1.0 + float(max_span) * alpha, 1.0 - float(max_span) * alpha


def compact_topology(label: str) -> str:
    mapping = {
        'braid_like_exchange': 'braid_like_exchange',
        'transfer_smeared': 'transfer_smeared',
        'fragmented_exchange': 'unresolved',
        'unresolved_mixed_topology': 'unresolved',
        'dispersive_pass': 'unresolved',
    }
    return mapping.get(label, label)


def signed_separation_trace(peak_tracks: list[list[list[float]]]) -> np.ndarray:
    track_a = np.asarray(peak_tracks[0], dtype=float)
    track_b = np.asarray(peak_tracks[1], dtype=float)
    return np.asarray(track_b[:, 0] - track_a[:, 0], dtype=float)


def build_observed_run(run: dict[str, Any], common: dict[str, Any]) -> dict[str, Any]:
    max_detuning = max(float(value) for value in common['detuning_values'])
    grade0_scale, grade1_scale = grade_scales_for_detuning(
        float(run['detuning_delta']),
        float(common['max_grade_detuning_span']),
        max_detuning,
    )
    observed_run = {
        'run_id': run['run_id'],
        'branch_id': 'A_grade_balance',
        'branch_label': 'Stage III-C2 effective extraction',
        'branch_order': 0,
        'variant_order': int(run['run_order']),
        'role': 'calibration',
        'variant_label': run['run_id'],
        'phase_offset_fraction_of_pi': float(run['phase_offset_fraction_of_pi']),
        'grade0_scale': grade0_scale,
        'grade1_scale': grade1_scale,
        'width_ratio_a_to_b': 1.0,
        'anisotropy_x': 1.0,
        'anisotropy_y': 1.0,
        'notes': run['notes'],
    }
    row, detail, _plot_paths = simulate_observed_run(observed_run, common)
    x_trace = signed_separation_trace(detail['peak_tracks'])
    y_trace = np.asarray(detail['grade_exchange_trace'], dtype=float)
    times = np.asarray(detail['times'], dtype=float)
    dt = float(times[1] - times[0]) if len(times) >= 2 else 0.0
    return {
        'run': run,
        'observed_row': row,
        'observed_detail': detail,
        'grade0_scale': grade0_scale,
        'grade1_scale': grade1_scale,
        'x_trace': x_trace,
        'y_trace': y_trace,
        'times': times,
        'dt': dt,
        'phase_centered': float(run['phase_offset_fraction_of_pi']) - 0.5,
        'observed_topology_label': compact_topology(str(row['topology_label'])),
    }


def fit_effective_maps(runs: list[dict[str, Any]]) -> tuple[np.ndarray, np.ndarray]:
    x_design: list[list[float]] = []
    x_target: list[float] = []
    y_design: list[list[float]] = []
    y_target: list[float] = []
    for item in runs:
        x_trace = np.asarray(item['x_trace'], dtype=float)
        y_trace = np.asarray(item['y_trace'], dtype=float)
        phase = float(item['phase_centered'])
        detuning = float(item['run']['detuning_delta'])
        for idx in range(len(x_trace) - 1):
            x_design.append([1.0, y_trace[idx], y_trace[idx] ** 2, phase, phase * y_trace[idx]])
            x_target.append(float(x_trace[idx + 1] - x_trace[idx]))
            y_design.append([1.0, x_trace[idx], x_trace[idx] ** 2, detuning, detuning * x_trace[idx]])
            y_target.append(float(y_trace[idx + 1] - y_trace[idx]))

    x_coeffs = np.linalg.lstsq(np.asarray(x_design, dtype=float), np.asarray(x_target, dtype=float), rcond=None)[0]
    y_coeffs = np.linalg.lstsq(np.asarray(y_design, dtype=float), np.asarray(y_target, dtype=float), rcond=None)[0]
    return x_coeffs, y_coeffs


def simulate_surrogate(
    x0: float,
    y0: float,
    steps: int,
    phase_centered: float,
    detuning_delta: float,
    x_coeffs: np.ndarray,
    y_coeffs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    x_values = [float(x0)]
    y_values = [float(y0)]
    for _ in range(steps - 1):
        x_n = float(x_values[-1])
        y_n = float(y_values[-1])
        dx = float(np.dot(x_coeffs, np.asarray([1.0, y_n, y_n * y_n, phase_centered, phase_centered * y_n], dtype=float)))
        dy = float(np.dot(y_coeffs, np.asarray([1.0, x_n, x_n * x_n, detuning_delta, detuning_delta * x_n], dtype=float)))
        x_values.append(float(np.clip(x_n + dx, -0.2, 0.2)))
        y_values.append(float(np.clip(y_n + dy, 0.0, 1.0)))
    return np.asarray(x_values, dtype=float), np.asarray(y_values, dtype=float)


def midpoint_phase_halfwidth(runs: list[dict[str, Any]]) -> float:
    braid_offsets = sorted(
        abs(float(item['phase_centered']))
        for item in runs
        if item['observed_topology_label'] == 'braid_like_exchange'
    )
    return 0.5 * braid_offsets[0] if braid_offsets else 0.0625


def unresolved_detuning_threshold(runs: list[dict[str, Any]]) -> float:
    unresolved = sorted(
        float(item['run']['detuning_delta'])
        for item in runs
        if item['observed_topology_label'] == 'unresolved'
    )
    low_non_unresolved = sorted(
        float(item['run']['detuning_delta'])
        for item in runs
        if item['observed_topology_label'] != 'unresolved'
    )
    if not unresolved or not low_non_unresolved:
        return 0.25
    lower = max(value for value in low_non_unresolved if value < unresolved[0]) if any(value < unresolved[0] for value in low_non_unresolved) else 0.0
    return 0.5 * (lower + unresolved[0])


def smear_peak_threshold(surrogate_runs: list[dict[str, Any]], unresolved_threshold: float) -> float:
    high_smear = [
        float(item['predicted_y_peak'])
        for item in surrogate_runs
        if item['observed_topology_label'] == 'transfer_smeared'
        and float(item['run']['detuning_delta']) >= unresolved_threshold
    ]
    other = [
        float(item['predicted_y_peak'])
        for item in surrogate_runs
        if not (
            item['observed_topology_label'] == 'transfer_smeared'
            and float(item['run']['detuning_delta']) >= unresolved_threshold
        )
    ]
    if not high_smear or not other:
        return 0.18
    return 0.5 * (min(high_smear) + max(other))


def classify_surrogate_topology(
    phase_centered: float,
    detuning_delta: float,
    y_peak: float,
    phase_halfwidth: float,
    unresolved_threshold: float,
    smear_threshold: float,
) -> str:
    if y_peak >= smear_threshold:
        return 'transfer_smeared'
    if detuning_delta >= unresolved_threshold:
        return 'unresolved'
    if abs(phase_centered) <= phase_halfwidth:
        return 'transfer_smeared'
    return 'braid_like_exchange'


def braid_survival_proxy(
    y_trace: np.ndarray,
    dt: float,
    phase_centered: float,
    detuning_delta: float,
    phase_halfwidth: float,
    unresolved_threshold: float,
    smear_threshold: float,
) -> float:
    if detuning_delta >= unresolved_threshold or abs(phase_centered) <= phase_halfwidth:
        return 0.0
    return float(np.sum(y_trace < smear_threshold) * dt)


def detuning_threshold(labels: dict[tuple[float, float], str], outer_phases: list[float], detuning_values: list[float]) -> float:
    thresholds: list[float] = []
    for phase in outer_phases:
        for detuning in sorted(detuning_values):
            if labels[(phase, detuning)] != 'braid_like_exchange':
                thresholds.append(float(detuning))
                break
    if not thresholds:
        return float(detuning_values[-1])
    return float(np.mean(thresholds))


def phase_band_width(labels: dict[tuple[float, float], str], phase_values: list[float], detuning: float) -> tuple[float, bool]:
    smear_phases = [float(phase) for phase in phase_values if labels[(phase, detuning)] == 'transfer_smeared']
    if not smear_phases:
        return 0.0, False
    width = float(max(smear_phases) - min(smear_phases))
    sorted_positions = sorted(phase_values.index(phase) for phase in smear_phases)
    contiguous = all(b - a == 1 for a, b in zip(sorted_positions, sorted_positions[1:])) if len(sorted_positions) > 1 else True
    return width, contiguous


def plot_phase_plane(path: Path, row: dict[str, Any]) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    ax.plot(row['observed_x_trace'], row['observed_y_trace'], color='tab:blue', linewidth=2.0, label='observed')
    ax.plot(row['predicted_x_trace'], row['predicted_y_trace'], color='tab:orange', linestyle='--', linewidth=2.0, label='surrogate')
    ax.scatter([row['observed_x_trace'][0]], [row['observed_y_trace'][0]], color='tab:blue', s=20)
    ax.scatter([row['predicted_x_trace'][0]], [row['predicted_y_trace'][0]], color='tab:orange', s=20)
    ax.set_xlabel('X (signed separation)')
    ax.set_ylabel('Y (grade exchange)')
    ax.set_title(
        f"{row['run_id']}\nobs={row['observed_topology_label']}, pred={row['predicted_topology_label']}"
    )
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_phase_diagram(path: Path, rows: list[dict[str, Any]], phase_values: list[float], detuning_values: list[float]) -> None:
    observed = np.zeros((len(detuning_values), len(phase_values)), dtype=float)
    predicted = np.zeros((len(detuning_values), len(phase_values)), dtype=float)
    for row in rows:
        i = detuning_values.index(float(row['detuning_delta']))
        j = phase_values.index(float(row['phase_offset_fraction_of_pi']))
        observed[i, j] = TOPOLOGY_CODE[str(row['observed_topology_label'])]
        predicted[i, j] = TOPOLOGY_CODE[str(row['predicted_topology_label'])]

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.8), sharey=True)
    for ax, grid, title in zip(axes, [observed, predicted], ['Observed DK phase diagram', 'Reduced surrogate phase diagram']):
        im = ax.imshow(grid, aspect='auto', cmap='viridis', vmin=0.0, vmax=2.0)
        ax.set_xticks(range(len(phase_values)), [f'{phase:.3f}' for phase in phase_values])
        ax.set_yticks(range(len(detuning_values)), [f'{delta:.2f}' for delta in detuning_values])
        ax.set_xlabel('phase offset (fraction of pi)')
        ax.set_title(title)
    axes[0].set_ylabel('detuning delta')
    fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_failure_map(path: Path, rows: list[dict[str, Any]], phase_values: list[float], detuning_values: list[float]) -> None:
    grid = np.zeros((len(detuning_values), len(phase_values)), dtype=float)
    for row in rows:
        i = detuning_values.index(float(row['detuning_delta']))
        j = phase_values.index(float(row['phase_offset_fraction_of_pi']))
        grid[i, j] = 0.0 if int(row['topology_match_flag']) == 1 else 1.0

    fig, ax = plt.subplots(figsize=(4.8, 4.4))
    im = ax.imshow(grid, aspect='auto', cmap='Reds', vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(phase_values)), [f'{phase:.3f}' for phase in phase_values])
    ax.set_yticks(range(len(detuning_values)), [f'{delta:.2f}' for delta in detuning_values])
    ax.set_xlabel('phase offset (fraction of pi)')
    ax.set_ylabel('detuning delta')
    ax.set_title('Stage III-C2 surrogate breakdown map')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    json_path: Path,
    csv_path: Path,
    rows: list[dict[str, Any]],
    stamped_plots: list[str],
    x_coeffs: np.ndarray,
    y_coeffs: np.ndarray,
    metrics: dict[str, Any],
    thresholds: dict[str, float],
    valid: bool,
) -> None:
    lines = [
        '# Stage III-C2 Minimal Effective Model Extraction v1',
        '',
        'This reduction stays inside the frozen HAOS (Harmonic Address Operating System) / IIP (Interaction Invariance Physics) Dirac-Kahler clustered braid window.',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f"effective_model_valid = `{'TRUE' if valid else 'FALSE'}`",
        '',
        '## Calibration note',
        'The source brief listed more phase and detuning values than a 9-run matrix can contain. This implementation therefore uses the resolved 3 x 3 anchor lattice `phase in {0.375, 0.500, 0.625}` crossed with bounded clustered-seed grade detuning `delta in {0.0, 0.5, 0.75}`.',
        'The extracted `X` coordinate is the signed separation envelope of the ordered peak tracks, because the unsigned clustered peak distance remains nearly constant on this branch.',
        '',
        '## Fitted map',
        'Use the reduced state',
        '- `X_n`: signed separation envelope',
        '- `Y_n`: net grade-exchange amplitude',
        '',
        '`X_{n+1} = X_n + c0 + c1 Y_n + c2 Y_n^2 + c3 phi + c4 phi Y_n`',
        '`Y_{n+1} = Y_n + d0 + d1 X_n + d2 X_n^2 + d3 delta + d4 delta X_n`',
        '',
        '### X-map coefficients',
    ]
    for idx, value in enumerate(x_coeffs):
        lines.append(f'- `c{idx}` = `{float(value):.6f}`')
    lines.extend([
        '',
        '### Y-map coefficients',
    ])
    for idx, value in enumerate(y_coeffs):
        lines.append(f'- `d{idx}` = `{float(value):.6f}`')

    lines.extend([
        '',
        '## Learned classifier thresholds',
        f"- midpoint phase halfwidth: `{thresholds['midpoint_phase_halfwidth']:.4f}`",
        f"- unresolved detuning threshold: `{thresholds['unresolved_detuning_threshold']:.4f}`",
        f"- smear peak threshold: `{thresholds['smear_peak_threshold']:.4f}`",
        '',
        '## Validation metrics',
        f"- topology agreement fraction: `{metrics['topology_agreement_fraction']:.4f}`",
        f"- detuning threshold error: `{metrics['detuning_threshold_error']:.4f}`",
        f"- phase-band width error: `{metrics['phase_band_width_error']:.4f}`",
        f"- single smeared band reproduced: `{metrics['single_smeared_band_reproduced']}`",
        f"- mean survival-time envelope mismatch: `{metrics['mean_survival_time_envelope_mismatch']:.4f}`",
        '',
        '## Per-run summary',
    ])
    for row in sorted(rows, key=lambda item: int(item['run_order'])):
        lines.append(
            f"- `{row['run_id']}`: obs=`{row['observed_topology_label']}`, pred=`{row['predicted_topology_label']}`, "
            f"`x_rmse={row['x_trace_rmse']:.4f}`, `y_rmse={row['y_trace_rmse']:.4f}`"
        )

    lines.extend([
        '',
        '## Interpretation',
    ])
    if valid:
        lines.append(
            'The clustered DK braid / smear sector admits a compact two-state surrogate on this anchor lattice. The reduction is descriptive rather than ontic: it compresses the clustered texture without overturning the Stage 23.10 result that the braid does not generalize into a family-wide intrinsic mechanism.'
        )
    else:
        lines.append(
            'The surrogate does not reproduce the clustered braid / smear atlas well enough on the anchor lattice. At this resolution the texture remains effectively higher-dimensional, so further progress should come from sector change or stronger structural constraints rather than more surrogate forcing.'
        )

    lines.extend([
        '',
        '## Plots',
    ])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')

    NOTE_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    common = runsheet['common_fields']

    observed_runs = [build_observed_run(run, common) for run in sorted(runsheet['runs'], key=lambda item: int(item['run_order']))]
    x_coeffs, y_coeffs = fit_effective_maps(observed_runs)

    phase_halfwidth = midpoint_phase_halfwidth(observed_runs)
    unresolved_threshold = unresolved_detuning_threshold(observed_runs)

    surrogate_runs: list[dict[str, Any]] = []
    for item in observed_runs:
        x_pred, y_pred = simulate_surrogate(
            x0=float(item['x_trace'][0]),
            y0=float(item['y_trace'][0]),
            steps=len(item['x_trace']),
            phase_centered=float(item['phase_centered']),
            detuning_delta=float(item['run']['detuning_delta']),
            x_coeffs=x_coeffs,
            y_coeffs=y_coeffs,
        )
        surrogate_runs.append(
            {
                **item,
                'predicted_x_trace': x_pred,
                'predicted_y_trace': y_pred,
                'predicted_y_peak': float(np.max(y_pred)),
                'predicted_x_final': float(x_pred[-1]),
            }
        )

    peak_threshold = smear_peak_threshold(surrogate_runs, unresolved_threshold)

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []

    for item in surrogate_runs:
        predicted_topology = classify_surrogate_topology(
            phase_centered=float(item['phase_centered']),
            detuning_delta=float(item['run']['detuning_delta']),
            y_peak=float(item['predicted_y_peak']),
            phase_halfwidth=phase_halfwidth,
            unresolved_threshold=unresolved_threshold,
            smear_threshold=peak_threshold,
        )
        observed_survival = braid_survival_proxy(
            y_trace=np.asarray(item['y_trace'], dtype=float),
            dt=float(item['dt']),
            phase_centered=float(item['phase_centered']),
            detuning_delta=float(item['run']['detuning_delta']),
            phase_halfwidth=phase_halfwidth,
            unresolved_threshold=unresolved_threshold,
            smear_threshold=peak_threshold,
        )
        predicted_survival = braid_survival_proxy(
            y_trace=np.asarray(item['predicted_y_trace'], dtype=float),
            dt=float(item['dt']),
            phase_centered=float(item['phase_centered']),
            detuning_delta=float(item['run']['detuning_delta']),
            phase_halfwidth=phase_halfwidth,
            unresolved_threshold=unresolved_threshold,
            smear_threshold=peak_threshold,
        )
        phase_plane_path = WORK_PLOT_DIR / f"{item['run']['run_id']}_phase_plane.png"
        plot_phase_plane(
            phase_plane_path,
            {
                'run_id': item['run']['run_id'],
                'observed_topology_label': item['observed_topology_label'],
                'predicted_topology_label': predicted_topology,
                'observed_x_trace': np.asarray(item['x_trace'], dtype=float),
                'observed_y_trace': np.asarray(item['y_trace'], dtype=float),
                'predicted_x_trace': np.asarray(item['predicted_x_trace'], dtype=float),
                'predicted_y_trace': np.asarray(item['predicted_y_trace'], dtype=float),
            },
        )
        plot_paths.append(phase_plane_path)

        row = {
            'run_id': item['run']['run_id'],
            'run_order': int(item['run']['run_order']),
            'phase_offset_fraction_of_pi': float(item['run']['phase_offset_fraction_of_pi']),
            'detuning_delta': float(item['run']['detuning_delta']),
            'grade0_scale': float(item['grade0_scale']),
            'grade1_scale': float(item['grade1_scale']),
            'observed_topology_label': item['observed_topology_label'],
            'predicted_topology_label': predicted_topology,
            'topology_match_flag': int(predicted_topology == item['observed_topology_label']),
            'observed_y_peak': float(np.max(item['y_trace'])),
            'predicted_y_peak': float(item['predicted_y_peak']),
            'observed_x_final': float(item['x_trace'][-1]),
            'predicted_x_final': float(item['predicted_x_final']),
            'x_trace_rmse': float(np.sqrt(np.mean((np.asarray(item['predicted_x_trace']) - np.asarray(item['x_trace'])) ** 2))),
            'y_trace_rmse': float(np.sqrt(np.mean((np.asarray(item['predicted_y_trace']) - np.asarray(item['y_trace'])) ** 2))),
            'observed_braid_survival_proxy': observed_survival,
            'predicted_braid_survival_proxy': predicted_survival,
            'survival_time_abs_error': abs(predicted_survival - observed_survival),
            'flow_concentration_index': float(item['observed_row']['flow_concentration_index']),
            'grade_exchange_coherence': float(item['observed_row']['grade_exchange_coherence']),
            'notes': str(item['run']['notes']),
        }
        detail = {
            'run_id': item['run']['run_id'],
            'phase_offset_fraction_of_pi': float(item['run']['phase_offset_fraction_of_pi']),
            'detuning_delta': float(item['run']['detuning_delta']),
            'grade0_scale': float(item['grade0_scale']),
            'grade1_scale': float(item['grade1_scale']),
            'observed_topology_label': item['observed_topology_label'],
            'predicted_topology_label': predicted_topology,
            'times': item['times'].tolist(),
            'observed_x_trace': np.asarray(item['x_trace']).tolist(),
            'observed_y_trace': np.asarray(item['y_trace']).tolist(),
            'predicted_x_trace': np.asarray(item['predicted_x_trace']).tolist(),
            'predicted_y_trace': np.asarray(item['predicted_y_trace']).tolist(),
            'observed_flow_concentration_trace': np.asarray(item['observed_detail']['flow_concentration_trace']).tolist(),
        }
        rows.append(row)
        details.append(detail)

    phase_values = [float(value) for value in common['phase_anchor_values_fraction_of_pi']]
    detuning_values = [float(value) for value in common['detuning_values']]

    observed_labels = {
        (float(row['phase_offset_fraction_of_pi']), float(row['detuning_delta'])): str(row['observed_topology_label'])
        for row in rows
    }
    predicted_labels = {
        (float(row['phase_offset_fraction_of_pi']), float(row['detuning_delta'])): str(row['predicted_topology_label'])
        for row in rows
    }

    observed_threshold = detuning_threshold(observed_labels, [phase_values[0], phase_values[-1]], detuning_values)
    predicted_threshold = detuning_threshold(predicted_labels, [phase_values[0], phase_values[-1]], detuning_values)
    observed_width, observed_single_band = phase_band_width(observed_labels, phase_values, detuning_values[0])
    predicted_width, predicted_single_band = phase_band_width(predicted_labels, phase_values, detuning_values[0])

    phase_diagram_path = WORK_PLOT_DIR / 'stage23_11_effective_phase_diagram.png'
    failure_map_path = WORK_PLOT_DIR / 'stage23_11_failure_map.png'
    plot_phase_diagram(phase_diagram_path, rows, phase_values, detuning_values)
    plot_failure_map(failure_map_path, rows, phase_values, detuning_values)
    plot_paths.extend([phase_diagram_path, failure_map_path])

    topology_agreement = float(np.mean([float(row['topology_match_flag']) for row in rows])) if rows else 0.0
    threshold_error = abs(predicted_threshold - observed_threshold)
    width_error = abs(predicted_width - observed_width)
    mean_survival_mismatch = float(np.mean([float(row['survival_time_abs_error']) for row in rows])) if rows else 0.0
    valid = topology_agreement >= 0.75 and threshold_error <= 0.1 and bool(predicted_single_band)

    metrics = {
        'topology_agreement_fraction': topology_agreement,
        'detuning_threshold_error': threshold_error,
        'phase_band_width_error': width_error,
        'single_smeared_band_reproduced': bool(predicted_single_band and observed_single_band),
        'mean_survival_time_envelope_mismatch': mean_survival_mismatch,
        'observed_detuning_threshold': observed_threshold,
        'predicted_detuning_threshold': predicted_threshold,
    }
    thresholds = {
        'midpoint_phase_halfwidth': phase_halfwidth,
        'unresolved_detuning_threshold': unresolved_threshold,
        'smear_peak_threshold': peak_threshold,
    }

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'fitted_x_coefficients': x_coeffs.tolist(),
        'fitted_y_coefficients': y_coeffs.tolist(),
        'classifier_thresholds': thresholds,
        'runs': details,
        'summary': {
            'effective_model_valid': bool(valid),
            'metrics': metrics,
            'observed_topology_counts': dict(Counter(str(row['observed_topology_label']) for row in rows)),
            'predicted_topology_counts': dict(Counter(str(row['predicted_topology_label']) for row in rows)),
        },
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        'stage23_11_dk_minimal_effective_model_extraction',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, x_coeffs, y_coeffs, metrics, thresholds, valid)


if __name__ == '__main__':
    main()
