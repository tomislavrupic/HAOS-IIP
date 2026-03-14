#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import REPO_ROOT, build_transverse_setup, load_stage10_defaults, suggested_dt
from stage11_collective_wave_interaction import build_component_packet
from stage14_effective_binding_test import edge_feedback_profile

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage14_binding_runs.json'
OUTPUT_DIR = REPO_ROOT / 'data'


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Probe Stage 14A feedback-term survival through the transverse projector.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--representative-id', default='tight_clustered_pair')
    parser.add_argument('--condition-id', default='weak_feedback')
    parser.add_argument('--resolution', type=int, default=0)
    parser.add_argument('--sample-count', type=int, default=5)
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sample_indices(length: int, count: int) -> list[int]:
    if length <= count:
        return list(range(length))
    grid = np.linspace(0, length - 1, count)
    return sorted({int(round(x)) for x in grid})


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    base = runsheet['base_seed_reference']
    rep = next(item for item in runsheet['representatives'] if item['representative_id'] == args.representative_id)
    cond = next(item for item in runsheet['conditions'] if item['condition_id'] == args.condition_id)
    resolution = int(args.resolution or rep['base_resolution'])
    feedback_eps = float(cond['feedback_eps'])
    feedback_sigma = float(runsheet['feedback_sigma'])

    data, projector = build_transverse_setup(
        n_side=resolution,
        epsilon=float(defaults['epsilon']),
        harmonic_tol=float(defaults['harmonic_tol']),
        boundary_type=str(base['boundary_type']),
    )

    packet_qs: list[np.ndarray] = []
    packet_vs: list[np.ndarray] = []
    kick_axis = int(defaults['kick_axis'])
    for center, amp, phase, kick_sign in zip(rep['packet_centers'], rep['packet_amplitudes'], rep['phase_offsets_rad'], rep['kick_signs']):
        q0, v0 = build_component_packet(
            center=np.asarray(center, dtype=float),
            amplitude=float(amp),
            phase_offset=float(phase),
            kick_sign=float(kick_sign),
            bandwidth=float(rep['bandwidth']),
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

    trace_rows: list[dict[str, float]] = []
    for step_idx in range(steps + 1):
        total_q = np.sum(packet_qs, axis=0)
        raw_feedback_profile, _env = edge_feedback_profile(data.midpoints, total_q, resolution, feedback_sigma)
        raw_feedback = feedback_eps * raw_feedback_profile * total_q
        projected_feedback = np.asarray(projector(raw_feedback), dtype=float)
        base_force = np.asarray(projector(-(operator @ total_q)), dtype=float)

        raw_norm = float(np.linalg.norm(raw_feedback))
        projected_norm = float(np.linalg.norm(projected_feedback))
        base_force_norm = float(np.linalg.norm(base_force))
        trace_rows.append({
            'time': float(step_idx * dt),
            'raw_feedback_norm': raw_norm,
            'projected_feedback_norm': projected_norm,
            'projection_survival_ratio': float(projected_norm / max(raw_norm, 1.0e-12)),
            'projected_vs_base_force_ratio': float(projected_norm / max(base_force_norm, 1.0e-12)),
            'raw_feedback_max': float(np.max(np.abs(raw_feedback))),
            'projected_feedback_max': float(np.max(np.abs(projected_feedback))),
            'base_force_norm': base_force_norm,
        })

        if step_idx == steps:
            break

        next_qs: list[np.ndarray] = []
        next_vs: list[np.ndarray] = []
        for q, v in zip(packet_qs, packet_vs):
            feedback_term = feedback_eps * raw_feedback_profile * q
            acc = np.asarray(projector(-(operator @ q) + feedback_term), dtype=float)
            v_half = v + 0.5 * dt * acc
            q_new = np.asarray(projector(q + dt * v_half), dtype=float)
            acc_new = np.asarray(projector(-(operator @ q_new) + feedback_eps * raw_feedback_profile * q_new), dtype=float)
            v_new = np.asarray(projector(v_half + 0.5 * dt * acc_new), dtype=float)
            next_qs.append(q_new)
            next_vs.append(v_new)
        packet_qs = next_qs
        packet_vs = next_vs

    sample_ids = sample_indices(len(trace_rows), int(args.sample_count))
    samples = [dict(index=i, **trace_rows[i]) for i in sample_ids]
    summary = {
        'representative_id': rep['representative_id'],
        'condition_id': cond['condition_id'],
        'feedback_eps': feedback_eps,
        'resolution': resolution,
        'sample_count': len(samples),
        'mean_projection_survival_ratio': float(np.mean([row['projection_survival_ratio'] for row in trace_rows])),
        'min_projection_survival_ratio': float(np.min([row['projection_survival_ratio'] for row in trace_rows])),
        'max_projection_survival_ratio': float(np.max([row['projection_survival_ratio'] for row in trace_rows])),
        'mean_projected_vs_base_force_ratio': float(np.mean([row['projected_vs_base_force_ratio'] for row in trace_rows])),
        'max_projected_vs_base_force_ratio': float(np.max([row['projected_vs_base_force_ratio'] for row in trace_rows])),
        'samples': samples,
    }

    stem = f"stage14_feedback_probe_{rep['representative_id']}_{cond['condition_id']}_n{resolution}"
    json_path = OUTPUT_DIR / f"{stem}.json"
    csv_path = OUTPUT_DIR / f"{stem}.csv"
    json_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_csv(csv_path, trace_rows, fieldnames=list(trace_rows[0].keys()))

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    print(json.dumps({k: v for k, v in summary.items() if k != 'samples'}, indent=2))


if __name__ == '__main__':
    main()
