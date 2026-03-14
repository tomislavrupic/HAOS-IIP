#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import (
    ATLAS_NOTES,
    PLOTS,
    REPO_ROOT,
    append_log,
    anisotropy_ratio,
    build_transverse_setup,
    coherence_score,
    create_field_snapshot,
    displacement,
    gaussian_packet,
    load_stage10_defaults,
    packet_width,
    plt,
    save_atlas_payload,
    spectral_moments,
    state_energy,
    suggested_dt,
    weighted_center,
)

NOTE_PATH = ATLAS_NOTES / 'Stage_11_Collective_Wave_Interaction_Atlas_v1.md'
RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage11_collective_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage11_collective_plots'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'stage',
    'operator_sector',
    'boundary_type',
    'packet_count',
    'interaction_regime',
    'collective_regime_label',
    'center_shift',
    'width_ratio',
    'coherence_score',
    'anisotropy_score',
    'constraint_max',
    'sector_leakage',
    'initial_mean_pair_distance',
    'final_mean_pair_distance',
    'min_mean_pair_distance',
    'max_mean_overlap',
    'phase_lock_indicator',
    'phase_lock_event_count',
    'composite_lifetime',
    'oscillation_count',
    'exchange_or_merger_flag',
    'notes',
]

LABELS = [
    'dispersive independent regime',
    'transient binding regime',
    'metastable composite regime',
    'oscillatory exchange regime',
    'large-scale drift field regime',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 11 collective wave interaction pilot.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    parser.add_argument('--max-runs', type=int, default=0)
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_runs(all_runs: list[dict[str, Any]], run_ids: list[str], max_runs: int) -> list[dict[str, Any]]:
    runs = all_runs
    if run_ids:
        wanted = set(run_ids)
        runs = [run for run in runs if run['run_id'] in wanted]
    if max_runs > 0:
        runs = runs[:max_runs]
    return runs


def make_leakage_fn(projector):
    def leakage(vec: np.ndarray) -> float:
        denom = max(float(np.linalg.norm(vec)), 1.0e-12)
        return float(np.linalg.norm(np.asarray(projector(vec), dtype=float) - vec) / denom)

    return leakage


def build_component_packet(
    center: np.ndarray,
    amplitude: float,
    phase_offset: float,
    bandwidth: float,
    central_k: float,
    phase_pattern: str,
    kick_axis: int,
    boundary_type: str,
    midpoints: np.ndarray,
    edge_axes: list[str],
    projector,
) -> tuple[np.ndarray, np.ndarray]:
    profile, delta = gaussian_packet(
        points=midpoints,
        center=center,
        sigma=bandwidth,
        boundary_type=boundary_type,
        central_k=0.0,
        phase_pattern='flat',
        kick_axis=kick_axis,
    )
    phase_coord = delta[:, kick_axis]
    if phase_pattern == 'cosine-carrier':
        carrier = np.cos(2.0 * math.pi * central_k * phase_coord + phase_offset)
    elif phase_pattern == 'sine-carrier':
        carrier = np.sin(2.0 * math.pi * central_k * phase_coord + phase_offset)
    else:
        carrier = np.cos(phase_offset) * np.ones_like(phase_coord)

    axis_mask = np.asarray([1.0 if axis == 'y' else 0.0 for axis in edge_axes], dtype=float)
    raw = profile * carrier * axis_mask
    q0 = np.asarray(projector(raw), dtype=float)
    q0 = q0 / max(float(np.linalg.norm(q0)), 1.0e-12)

    raw_v = (delta[:, kick_axis] / max(bandwidth * bandwidth, 1.0e-12)) * raw
    v0 = np.asarray(projector(raw_v), dtype=float)
    v0 = v0 / max(float(np.linalg.norm(v0)), 1.0e-12)
    return amplitude * q0, amplitude * v0


def pairwise_distance(centers: list[np.ndarray], boundary_type: str) -> list[float]:
    values: list[float] = []
    for i, j in combinations(range(len(centers)), 2):
        delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float), boundary_type)[0]
        values.append(float(np.linalg.norm(delta)))
    return values


def pairwise_overlap(component_qs: list[np.ndarray]) -> list[float]:
    values: list[float] = []
    weights = [np.abs(q) ** 2 for q in component_qs]
    for i, j in combinations(range(len(weights)), 2):
        wi = weights[i]
        wj = weights[j]
        denom = math.sqrt(max(float(np.sum(wi)) * float(np.sum(wj)), 1.0e-12))
        values.append(float(np.sum(np.sqrt(wi * wj)) / denom))
    return values


def pairwise_phase_alignment(component_qs: list[np.ndarray]) -> list[float]:
    values: list[float] = []
    for i, j in combinations(range(len(component_qs)), 2):
        qi = component_qs[i]
        qj = component_qs[j]
        denom = max(float(np.linalg.norm(qi)) * float(np.linalg.norm(qj)), 1.0e-12)
        values.append(float(np.vdot(qi, qj).real / denom))
    return values


def count_phase_lock_events(series: list[float], threshold: float = 0.85) -> int:
    if not series:
        return 0
    active = np.asarray(series, dtype=float) >= threshold
    transitions = np.diff(active.astype(int))
    return int(np.sum(transitions == 1) + int(active[0]))


def count_oscillations(series: list[float], eps: float = 1.0e-4) -> int:
    values = np.asarray(series, dtype=float)
    if values.size < 3:
        return 0
    deriv = np.diff(values)
    deriv[np.abs(deriv) < eps] = 0.0
    signs = np.sign(deriv)
    filtered = [s for s in signs if s != 0.0]
    if len(filtered) < 2:
        return 0
    return int(sum(1 for a, b in zip(filtered[:-1], filtered[1:]) if a != b))


def classify_collective(summary: dict[str, float | int | str]) -> str:
    packet_count = int(summary['packet_count'])
    center_shift = float(summary['center_shift'])
    width_ratio = float(summary['width_ratio'])
    max_overlap = float(summary['max_mean_overlap'])
    phase_lock = float(summary['phase_lock_indicator'])
    composite_lifetime = float(summary['composite_lifetime'])
    oscillations = int(summary['oscillation_count'])
    initial_distance = max(float(summary['initial_mean_pair_distance']), 1.0e-12)
    min_distance = float(summary['min_mean_pair_distance'])
    exchange_flag = str(summary['exchange_or_merger_flag'])
    t_final = float(summary['t_final'])

    if packet_count >= 3 and center_shift > 0.075 and phase_lock > 0.55:
        return 'large-scale drift field regime'
    if exchange_flag == 'exchange' or (oscillations >= 2 and max_overlap > 0.18):
        return 'oscillatory exchange regime'
    if exchange_flag == 'merger' or (composite_lifetime > 0.45 * t_final and max_overlap > 0.5 and width_ratio <= 1.2):
        return 'metastable composite regime'
    if max_overlap > 0.2 or min_distance < 0.8 * initial_distance:
        return 'transient binding regime'
    return 'dispersive independent regime'


def run_case(run: dict[str, Any], defaults: dict[str, Any], base: dict[str, Any]) -> tuple[dict[str, Any], list[Path]]:
    data, projector = build_transverse_setup(
        n_side=int(defaults['n_side']),
        epsilon=float(defaults['epsilon']),
        harmonic_tol=float(defaults['harmonic_tol']),
        boundary_type=str(base['boundary_type']),
    )

    kick_axis = int(defaults['kick_axis'])
    packet_qs: list[np.ndarray] = []
    packet_vs: list[np.ndarray] = []
    for center, amp, phase in zip(run['packet_centers'], run['packet_amplitudes'], run['phase_offsets_rad']):
        q0, v0 = build_component_packet(
            center=np.asarray(center, dtype=float),
            amplitude=float(amp),
            phase_offset=float(phase),
            bandwidth=float(base['bandwidth']),
            central_k=float(base['central_k']),
            phase_pattern=str(base['phase_pattern']),
            kick_axis=kick_axis,
            boundary_type=str(base['boundary_type']),
            midpoints=data.midpoints,
            edge_axes=data.edge_axes,
            projector=projector,
        )
        packet_qs.append(q0)
        packet_vs.append(v0)

    operator = data.L1
    dt = suggested_dt(operator, float(defaults['dt_safety']))
    steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
    times: list[float] = []
    total_centers: list[list[float]] = []
    widths: list[float] = []
    norms: list[float] = []
    energies: list[float] = []
    anisotropies: list[float] = []
    spectral_centroids: list[float] = []
    spectral_spreads: list[float] = []
    sector_leakages: list[float] = []
    constraint_norms: list[float] = []
    coherences: list[float] = []
    pair_distance_series: list[float] = []
    overlap_series: list[float] = []
    phase_alignment_series: list[float] = []
    component_center_traces: list[list[list[float]]] = [[] for _ in range(len(packet_qs))]

    leakage_fn = make_leakage_fn(projector)

    def record(t: float) -> None:
        total_q = np.sum(packet_qs, axis=0)
        total_v = np.sum(packet_vs, axis=0)
        total_weights = np.abs(total_q) ** 2
        total_center = weighted_center(data.midpoints, total_weights, str(base['boundary_type']))
        centroid, spread = spectral_moments(operator, total_q)
        times.append(float(t))
        total_centers.append(total_center.tolist())
        widths.append(packet_width(data.midpoints, total_weights, total_center, str(base['boundary_type'])))
        norms.append(float(np.linalg.norm(total_q)))
        energies.append(state_energy(operator, total_q, total_v))
        anisotropies.append(anisotropy_ratio(data.midpoints, total_weights, total_center, str(base['boundary_type'])))
        spectral_centroids.append(centroid)
        spectral_spreads.append(spread)
        sector_leakages.append(float(leakage_fn(total_q)))
        constraint_norms.append(float(np.linalg.norm(data.d0.T @ total_q)))
        coherences.append(coherence_score(total_weights))

        component_centers: list[np.ndarray] = []
        for idx, comp_q in enumerate(packet_qs):
            comp_weights = np.abs(comp_q) ** 2
            comp_center = weighted_center(data.midpoints, comp_weights, str(base['boundary_type']))
            component_centers.append(comp_center)
            component_center_traces[idx].append(comp_center.tolist())

        distances = pairwise_distance(component_centers, str(base['boundary_type']))
        overlaps = pairwise_overlap(packet_qs)
        alignments = [abs(v) for v in pairwise_phase_alignment(packet_qs)]
        pair_distance_series.append(float(np.mean(distances)) if distances else 0.0)
        overlap_series.append(float(np.mean(overlaps)) if overlaps else 0.0)
        phase_alignment_series.append(float(np.mean(alignments)) if alignments else 0.0)

    record(0.0)
    for _ in range(steps):
        next_qs: list[np.ndarray] = []
        next_vs: list[np.ndarray] = []
        for q, v in zip(packet_qs, packet_vs):
            acc = np.asarray(projector(-(operator @ q)), dtype=float)
            v_half = v + 0.5 * dt * acc
            q_new = np.asarray(projector(q + dt * v_half), dtype=float)
            acc_new = np.asarray(projector(-(operator @ q_new)), dtype=float)
            v_new = np.asarray(projector(v_half + 0.5 * dt * acc_new), dtype=float)
            next_qs.append(q_new)
            next_vs.append(v_new)
        packet_qs = next_qs
        packet_vs = next_vs
        record(times[-1] + dt)

    centers_arr = np.asarray(total_centers, dtype=float)
    path = 0.0
    for prev, cur in zip(centers_arr[:-1], centers_arr[1:]):
        path += float(np.linalg.norm(displacement(cur[None, :], prev, str(base['boundary_type']))[0]))
    shift = float(np.linalg.norm(displacement(centers_arr[-1][None, :], centers_arr[0], str(base['boundary_type']))[0]))

    phase_lock_indicator = float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0
    phase_lock_events = count_phase_lock_events(phase_alignment_series)
    composite_lifetime = float(np.sum(np.asarray(overlap_series) >= 0.35) * dt)
    oscillation_count = count_oscillations(pair_distance_series)
    max_mean_overlap = float(np.max(overlap_series)) if overlap_series else 0.0
    initial_mean_pair_distance = float(pair_distance_series[0]) if pair_distance_series else 0.0
    final_mean_pair_distance = float(pair_distance_series[-1]) if pair_distance_series else 0.0
    min_mean_pair_distance = float(np.min(pair_distance_series)) if pair_distance_series else 0.0

    if max_mean_overlap > 0.65 and composite_lifetime > 0.18:
        exchange_flag = 'merger'
    elif oscillation_count >= 2 and max_mean_overlap > 0.18:
        exchange_flag = 'exchange'
    else:
        exchange_flag = 'none'

    summary = {
        'packet_count': int(run['packet_count']),
        'center_shift': shift,
        'path_length': path,
        'path_efficiency': float(shift / max(path, 1.0e-12)),
        'initial_width': float(widths[0]),
        'final_width': float(widths[-1]),
        'width_change': float(widths[-1] - widths[0]),
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'relative_norm_change': float((norms[-1] - norms[0]) / max(norms[0], 1.0e-12)),
        'relative_energy_drift': float((energies[-1] - energies[0]) / max(abs(energies[0]), 1.0e-12)),
        'max_anisotropy': float(np.max(anisotropies)),
        'max_constraint_norm': float(np.max(constraint_norms)),
        'max_sector_leakage': float(np.max(sector_leakages)),
        'coherence_initial': float(coherences[0]),
        'coherence_final': float(coherences[-1]),
        'coherence_ratio': float(coherences[-1] / max(coherences[0], 1.0e-12)),
        'spectral_centroid_initial': float(spectral_centroids[0]),
        'spectral_centroid_final': float(spectral_centroids[-1]),
        'spectral_spread_initial': float(spectral_spreads[0]),
        'spectral_spread_final': float(spectral_spreads[-1]),
        'initial_mean_pair_distance': initial_mean_pair_distance,
        'final_mean_pair_distance': final_mean_pair_distance,
        'min_mean_pair_distance': min_mean_pair_distance,
        'max_mean_overlap': max_mean_overlap,
        'phase_lock_indicator': phase_lock_indicator,
        'phase_lock_event_count': phase_lock_events,
        'composite_lifetime': composite_lifetime,
        'oscillation_count': oscillation_count,
        'exchange_or_merger_flag': exchange_flag,
        't_final': float(times[-1]),
    }
    label = classify_collective(summary)

    row = {
        'run_id': run['run_id'],
        'stage': 'Stage 11',
        'operator_sector': base['operator_sector'],
        'boundary_type': base['boundary_type'],
        'packet_count': int(run['packet_count']),
        'interaction_regime': run['interaction_regime'],
        'collective_regime_label': label,
        'center_shift': float(summary['center_shift']),
        'width_ratio': float(summary['width_ratio']),
        'coherence_score': float(summary['coherence_ratio']),
        'anisotropy_score': float(summary['max_anisotropy']),
        'constraint_max': float(summary['max_constraint_norm']),
        'sector_leakage': float(summary['max_sector_leakage']),
        'initial_mean_pair_distance': initial_mean_pair_distance,
        'final_mean_pair_distance': final_mean_pair_distance,
        'min_mean_pair_distance': min_mean_pair_distance,
        'max_mean_overlap': max_mean_overlap,
        'phase_lock_indicator': phase_lock_indicator,
        'phase_lock_event_count': phase_lock_events,
        'composite_lifetime': composite_lifetime,
        'oscillation_count': oscillation_count,
        'exchange_or_merger_flag': exchange_flag,
        'notes': run['selection_note'],
    }

    payload = {
        'run': run,
        'summary': summary,
        'collective_regime_label': label,
        'times': times,
        'total_centers': total_centers,
        'component_center_traces': component_center_traces,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'phase_alignment_series': phase_alignment_series,
        'coherences': coherences,
        'constraint_norms': constraint_norms,
        'sector_leakages': sector_leakages,
    }

    plot_paths: list[Path] = []
    total_q = np.sum(packet_qs, axis=0)
    field_path = WORK_PLOT_DIR / f"stage11_{run['run_id']}_field_snapshot.png"
    create_field_snapshot(field_path, run['run_id'], data.midpoints, total_q, str(base['boundary_type']), int(defaults['slice_axis']))
    plot_paths.append(field_path)

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    for idx, trace in enumerate(component_center_traces):
        xs = [item[0] for item in trace]
        ax.plot(times, xs, label=f'packet {idx+1}')
    ax.set_xlabel('time')
    ax.set_ylabel('x-centroid')
    ax.set_title(f'Centroid traces: {run["run_id"]}')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    centroid_path = WORK_PLOT_DIR / f"stage11_{run['run_id']}_centroid_traces.png"
    fig.savefig(centroid_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(centroid_path)

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.plot(times, pair_distance_series)
    ax.set_xlabel('time')
    ax.set_ylabel('mean pair distance')
    ax.set_title(f'Pair distance: {run["run_id"]}')
    ax.grid(alpha=0.25)
    distance_path = WORK_PLOT_DIR / f"stage11_{run['run_id']}_pair_distance.png"
    fig.savefig(distance_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(distance_path)

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.plot(times, overlap_series, label='mean overlap')
    ax.plot(times, phase_alignment_series, label='phase lock')
    ax.set_xlabel('time')
    ax.set_ylabel('interaction observable')
    ax.set_title(f'Overlap and phase lock: {run["run_id"]}')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    overlap_path = WORK_PLOT_DIR / f"stage11_{run['run_id']}_overlap_phase.png"
    fig.savefig(overlap_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(overlap_path)

    return payload, row, plot_paths


def create_summary_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    colors = {
        'dispersive independent regime': 'tab:blue',
        'transient binding regime': 'tab:green',
        'metastable composite regime': 'tab:orange',
        'oscillatory exchange regime': 'tab:red',
        'large-scale drift field regime': 'tab:purple',
    }
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for row in rows:
        ax.scatter(
            float(row['max_mean_overlap']),
            float(row['composite_lifetime']),
            color=colors.get(str(row['collective_regime_label']), 'black'),
            s=80,
            alpha=0.85,
            label=str(row['collective_regime_label']),
        )
        ax.text(float(row['max_mean_overlap']) + 0.01, float(row['composite_lifetime']), str(row['packet_count']), fontsize=8)
    handles, labels = ax.get_legend_handles_labels()
    unique: dict[str, Any] = {}
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    ax.legend(unique.values(), unique.keys(), fontsize=8)
    ax.set_xlabel('max mean overlap')
    ax.set_ylabel('composite lifetime')
    ax.set_title('Stage 11 collective morphology map')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(path: Path, json_rel: str, csv_rel: str, rows: list[dict[str, Any]], regime_counts: dict[str, int]) -> None:
    lines = [
        '# Stage 11 Collective Wave Interaction Atlas v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        '',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Regime counts: {regime_counts}',
        '',
        'This pilot remains classification-first and interpretation-minimal.',
        'Because the current frozen dynamics is linear, the collective observables should be read as superposition and interference morphology rather than as proof of nonlinear effective forces.',
        '',
        'Per-run summary:',
        '',
    ]
    for row in rows:
        lines.extend([
            f"- `{row['run_id']}`",
            f"  - label: `{row['collective_regime_label']}`",
            f"  - max overlap: `{float(row['max_mean_overlap']):.4f}`",
            f"  - composite lifetime: `{float(row['composite_lifetime']):.4f}`",
            f"  - phase-lock indicator: `{float(row['phase_lock_indicator']):.4f}`",
            f"  - exchange/merger: `{row['exchange_or_merger_flag']}`",
        ])
    lines.extend([
        '',
        'Classification labels remain descriptive only:',
        '- dispersive independent regime',
        '- transient binding regime',
        '- metastable composite regime',
        '- oscillatory exchange regime',
        '- large-scale drift field regime',
        '',
        'Do not read these runs as gravity, curvature, or force claims. They are interaction-morphology diagnostics only.',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    runs = selected_runs(runsheet['runs'], args.run_ids, args.max_runs)
    base = runsheet['base_seed_reference']

    payload_runs: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        payload, row, per_run_plots = run_case(run, defaults, base)
        payload_runs.append(payload)
        rows.append(row)
        plot_paths.extend(per_run_plots)

    regime_counts = dict(Counter(row['collective_regime_label'] for row in rows))
    summary_plot = WORK_PLOT_DIR / 'stage11_collective_morphology_map.png'
    create_summary_plot(summary_plot, rows)
    plot_paths.append(summary_plot)

    result = {
        'stage': 'Stage 11',
        'description': 'Collective wave interaction atlas pilot on the frozen projected transverse architecture.',
        'base_seed_reference': base,
        'runs': payload_runs,
        'regime_counts': regime_counts,
        'labels': LABELS,
        'conclusion': 'Stage 11 pilot maps interaction morphology for five canonical multi-packet transverse cases without introducing any new operators.',
    }

    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        experiment_slug='stage11_collective_wave_interaction',
        result=result,
        csv_rows=rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )

    write_note(
        NOTE_PATH,
        json_rel=str(json_path.relative_to(REPO_ROOT)),
        csv_rel=str(csv_path.relative_to(REPO_ROOT)),
        rows=rows,
        regime_counts=regime_counts,
    )

    append_log(
        title=f'Stage 11 Collective Wave Interaction ({json_path.stem})',
        config_summary=f"runs={len(runs)}, sector={base['operator_sector']}, boundary={base['boundary_type']}, bandwidth={base['bandwidth']}, k={base['central_k']}",
        result_path=json_path,
        stamped_plots=stamped_plots,
        observation=f'regime_counts={regime_counts}',
        conclusion='the Stage 11 pilot adds a first interaction-morphology layer on the frozen projected transverse architecture without extending the interpretation boundary',
    )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(f'regime_counts={regime_counts}')


if __name__ == '__main__':
    main()
