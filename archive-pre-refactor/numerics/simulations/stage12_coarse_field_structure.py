#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, deque
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import (
    ATLAS_NOTES,
    REPO_ROOT,
    append_log,
    anisotropy_ratio,
    build_transverse_setup,
    coherence_score,
    create_field_snapshot,
    displacement,
    load_stage10_defaults,
    packet_width,
    plt,
    save_atlas_payload,
    spectral_moments,
    state_energy,
    suggested_dt,
    weighted_center,
)
from stage11_collective_wave_interaction import build_component_packet, make_leakage_fn

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage12_coarse_runs.json'
CSV_FIELDS = [
    'run_id',
    'stage',
    'operator_sector',
    'boundary_type',
    'packet_count',
    'cluster_pattern',
    'coarse_regime_label',
    'center_shift',
    'width_ratio',
    'coherence_score',
    'anisotropy_score',
    'constraint_max',
    'sector_leakage',
    'mean_basin_count_sigma2',
    'mean_basin_count_sigma4',
    'mean_dominant_area_fraction_sigma2',
    'mean_dominant_area_fraction_sigma4',
    'basin_lifetime_sigma2',
    'basin_lifetime_sigma4',
    'coarse_persistence_sigma2',
    'coarse_persistence_sigma4',
    'mean_envelope_variance_sigma2',
    'mean_envelope_variance_sigma4',
    'notes',
]
LABELS = [
    'basin-dominated regime',
    'multi-basin fluctuating regime',
    'diffuse coarse field regime',
]
SIGMAS = (2.0, 4.0)
ALPHA = 0.6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 12 coarse-grained field structure pilot.')
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


def edge_field_to_grid(midpoints: np.ndarray, values: np.ndarray, n_side: int) -> np.ndarray:
    grid = np.zeros((n_side, n_side, n_side), dtype=float)
    coords = np.floor(midpoints * n_side).astype(int) % n_side
    for idx, amp in zip(coords, values):
        grid[idx[0], idx[1], idx[2]] += float(amp)
    return grid


def gaussian_transfer(n_side: int, sigma: float) -> np.ndarray:
    freqs = np.fft.fftfreq(n_side) * (2.0 * math.pi)
    kx, ky, kz = np.meshgrid(freqs, freqs, freqs, indexing='ij')
    return np.exp(-0.5 * sigma * sigma * (kx * kx + ky * ky + kz * kz))


def gaussian_smooth_periodic(grid: np.ndarray, sigma: float) -> np.ndarray:
    transfer = gaussian_transfer(grid.shape[0], sigma)
    return np.fft.ifftn(np.fft.fftn(grid) * transfer).real


def connected_components_periodic(mask: np.ndarray) -> list[list[tuple[int, int, int]]]:
    n_side = mask.shape[0]
    visited = np.zeros_like(mask, dtype=bool)
    comps: list[list[tuple[int, int, int]]] = []
    neighbors = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    for start in zip(*np.where(mask)):
        if visited[start]:
            continue
        comp: list[tuple[int, int, int]] = []
        queue: deque[tuple[int, int, int]] = deque([start])
        visited[start] = True
        while queue:
            i, j, k = queue.popleft()
            comp.append((i, j, k))
            for di, dj, dk in neighbors:
                ni = (i + di) % n_side
                nj = (j + dj) % n_side
                nk = (k + dk) % n_side
                if mask[ni, nj, nk] and not visited[ni, nj, nk]:
                    visited[ni, nj, nk] = True
                    queue.append((ni, nj, nk))
        comps.append(comp)
    return comps


def component_mask(component: list[tuple[int, int, int]], n_side: int) -> np.ndarray:
    mask = np.zeros((n_side, n_side, n_side), dtype=bool)
    if component:
        idx = np.asarray(component, dtype=int)
        mask[idx[:, 0], idx[:, 1], idx[:, 2]] = True
    return mask


def jaccard(a: np.ndarray | None, b: np.ndarray | None) -> float:
    if a is None or b is None:
        return 0.0
    inter = np.logical_and(a, b).sum()
    union = np.logical_or(a, b).sum()
    if union == 0:
        return 0.0
    return float(inter / union)


def classify_coarse(summary: dict[str, float]) -> str:
    mean_basin_count_2 = float(summary['mean_basin_count_sigma2'])
    mean_dom_frac_2 = float(summary['mean_dominant_area_fraction_sigma2'])
    mean_dom_frac_4 = float(summary['mean_dominant_area_fraction_sigma4'])
    basin_lifetime = float(summary['basin_lifetime_sigma4'])
    persistence_4 = float(summary['coarse_persistence_sigma4'])
    env_var_4 = float(summary['mean_envelope_variance_sigma4'])
    t_final = float(summary['t_final'])

    if mean_basin_count_2 > 1.05:
        return 'multi-basin fluctuating regime'
    if basin_lifetime >= 0.45 * t_final and mean_dom_frac_2 <= 0.06 and mean_dom_frac_4 <= 0.35 and persistence_4 >= 0.9 and env_var_4 >= 0.02:
        return 'basin-dominated regime'
    if mean_dom_frac_4 >= 0.35 or env_var_4 < 0.02:
        return 'diffuse coarse field regime'
    return 'basin-dominated regime'


def render_grid_slice(ax, grid: np.ndarray, title: str, cmap: str = 'viridis') -> None:
    n_side = grid.shape[0]
    slice_idx = n_side // 2
    img = grid[:, :, slice_idx].T
    im = ax.imshow(img, origin='lower', cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)


def run_case(run: dict[str, Any], defaults: dict[str, Any], base: dict[str, Any], work_plot_dir: Path) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    data, projector = build_transverse_setup(
        n_side=int(defaults['n_side']),
        epsilon=float(defaults['epsilon']),
        harmonic_tol=float(defaults['harmonic_tol']),
        boundary_type=str(base['boundary_type']),
    )
    kick_axis = int(defaults['kick_axis'])
    packet_qs: list[np.ndarray] = []
    packet_vs: list[np.ndarray] = []
    kick_signs = run.get('kick_signs', [1.0] * int(run['packet_count']))
    for center, amp, phase, kick_sign in zip(run['packet_centers'], run['packet_amplitudes'], run['phase_offsets_rad'], kick_signs):
        q0, v0 = build_component_packet(
            center=np.asarray(center, dtype=float),
            amplitude=float(amp),
            phase_offset=float(phase),
            kick_sign=float(kick_sign),
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
    leakage_fn = make_leakage_fn(projector)

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

    coarse: dict[float, dict[str, Any]] = {
        sigma: {
            'envelopes': [],
            'basin_counts': [],
            'dominant_area_fractions': [],
            'envelope_variances': [],
            'stable_flags': [],
            'dominant_masks': [],
            'basin_lifetime_flags': [],
        }
        for sigma in SIGMAS
    }

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

        point_grid = edge_field_to_grid(data.midpoints, np.abs(total_q), int(defaults['n_side']))
        for sigma in SIGMAS:
            env = gaussian_smooth_periodic(point_grid, sigma)
            env = np.maximum(env, 0.0)
            coarse[sigma]['envelopes'].append(env)
            max_env = float(np.max(env))
            norm_env = env / max(max_env, 1.0e-12)
            coarse[sigma]['envelope_variances'].append(float(np.var(norm_env)))
            mask = env > (ALPHA * max_env) if max_env > 0.0 else np.zeros_like(env, dtype=bool)
            comps = connected_components_periodic(mask)
            coarse[sigma]['basin_counts'].append(float(len(comps)))
            if comps:
                dominant = max(comps, key=len)
                dominant_mask = component_mask(dominant, env.shape[0])
                dom_frac = float(len(dominant) / env.size)
            else:
                dominant_mask = None
                dom_frac = 0.0
            coarse[sigma]['dominant_area_fractions'].append(dom_frac)
            prev_mask = coarse[sigma]['dominant_masks'][-1] if coarse[sigma]['dominant_masks'] else None
            if dominant_mask is None:
                stable = 0.0
            elif prev_mask is None:
                stable = 1.0
            else:
                stable = 1.0 if jaccard(prev_mask, dominant_mask) >= 0.5 else 0.0
            coarse[sigma]['stable_flags'].append(stable)
            coarse[sigma]['basin_lifetime_flags'].append(1.0 if dom_frac >= 0.05 else 0.0)
            coarse[sigma]['dominant_masks'].append(dominant_mask)

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
    shift = float(np.linalg.norm(displacement(centers_arr[-1][None, :], centers_arr[0], str(base['boundary_type']))[0]))
    summary: dict[str, float] = {
        'center_shift': shift,
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'relative_norm_change': float((norms[-1] - norms[0]) / max(norms[0], 1.0e-12)),
        'relative_energy_drift': float((energies[-1] - energies[0]) / max(abs(energies[0]), 1.0e-12)),
        'max_anisotropy': float(np.max(anisotropies)),
        'max_constraint_norm': float(np.max(constraint_norms)),
        'max_sector_leakage': float(np.max(sector_leakages)),
        'coherence_ratio': float(coherences[-1] / max(coherences[0], 1.0e-12)),
        't_final': float(times[-1]),
    }

    for sigma in SIGMAS:
        stable_flags = np.asarray(coarse[sigma]['stable_flags'], dtype=float)
        lifetime_flags = np.asarray(coarse[sigma]['basin_lifetime_flags'], dtype=float)
        summary[f'mean_basin_count_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['basin_counts']))
        summary[f'mean_dominant_area_fraction_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['dominant_area_fractions']))
        summary[f'coarse_persistence_sigma{int(sigma)}'] = float(np.mean(stable_flags)) if stable_flags.size else 0.0
        summary[f'basin_lifetime_sigma{int(sigma)}'] = float(np.sum(lifetime_flags) * dt)
        summary[f'mean_envelope_variance_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['envelope_variances']))

    label = classify_coarse(summary)
    row = {
        'run_id': run['run_id'],
        'stage': 'Stage 12',
        'operator_sector': base['operator_sector'],
        'boundary_type': base['boundary_type'],
        'packet_count': int(run['packet_count']),
        'cluster_pattern': run['cluster_pattern'],
        'coarse_regime_label': label,
        'center_shift': summary['center_shift'],
        'width_ratio': summary['width_ratio'],
        'coherence_score': summary['coherence_ratio'],
        'anisotropy_score': summary['max_anisotropy'],
        'constraint_max': summary['max_constraint_norm'],
        'sector_leakage': summary['max_sector_leakage'],
        'mean_basin_count_sigma2': summary['mean_basin_count_sigma2'],
        'mean_basin_count_sigma4': summary['mean_basin_count_sigma4'],
        'mean_dominant_area_fraction_sigma2': summary['mean_dominant_area_fraction_sigma2'],
        'mean_dominant_area_fraction_sigma4': summary['mean_dominant_area_fraction_sigma4'],
        'basin_lifetime_sigma2': summary['basin_lifetime_sigma2'],
        'basin_lifetime_sigma4': summary['basin_lifetime_sigma4'],
        'coarse_persistence_sigma2': summary['coarse_persistence_sigma2'],
        'coarse_persistence_sigma4': summary['coarse_persistence_sigma4'],
        'mean_envelope_variance_sigma2': summary['mean_envelope_variance_sigma2'],
        'mean_envelope_variance_sigma4': summary['mean_envelope_variance_sigma4'],
        'notes': run['selection_note'],
    }

    payload = {
        'run': run,
        'summary': summary,
        'coarse_regime_label': label,
        'times': times,
        'total_centers': total_centers,
        'widths': widths,
        'coherences': coherences,
        'constraint_norms': constraint_norms,
        'sector_leakages': sector_leakages,
        'coarse_metrics': {
            f'sigma{int(sigma)}': {
                'basin_counts': coarse[sigma]['basin_counts'],
                'dominant_area_fractions': coarse[sigma]['dominant_area_fractions'],
                'envelope_variances': coarse[sigma]['envelope_variances'],
                'stable_flags': coarse[sigma]['stable_flags'],
            }
            for sigma in SIGMAS
        },
    }

    plot_paths: list[Path] = []
    total_q = np.sum(packet_qs, axis=0)
    field_path = work_plot_dir / f"stage12_{run['run_id']}_field_snapshot.png"
    create_field_snapshot(field_path, run['run_id'], data.midpoints, total_q, str(base['boundary_type']), int(defaults['slice_axis']))
    plot_paths.append(field_path)

    env = coarse[4.0]['envelopes'][-1]
    env_path = work_plot_dir / f"stage12_{run['run_id']}_envelope_snapshot.png"
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    render_grid_slice(ax, env, f'Envelope sigma=4: {run["run_id"]}')
    fig.savefig(env_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(env_path)

    basin_path = work_plot_dir / f"stage12_{run['run_id']}_basin_segmentation.png"
    max_env = float(np.max(env))
    basin_mask = env > (ALPHA * max_env) if max_env > 0.0 else np.zeros_like(env, dtype=bool)
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    render_grid_slice(ax, basin_mask.astype(float), f'Basin segmentation sigma=4: {run["run_id"]}', cmap='magma')
    fig.savefig(basin_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(basin_path)

    area_path = work_plot_dir / f"stage12_{run['run_id']}_dominant_basin_area.png"
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    ax.plot(times, coarse[2.0]['dominant_area_fractions'], label='sigma=2')
    ax.plot(times, coarse[4.0]['dominant_area_fractions'], label='sigma=4')
    ax.set_xlabel('time')
    ax.set_ylabel('dominant basin area fraction')
    ax.set_title(f'Dominant basin area: {run["run_id"]}')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(area_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(area_path)

    persistence_path = work_plot_dir / f"stage12_{run['run_id']}_persistence_trace.png"
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    for sigma in SIGMAS:
        stable = np.asarray(coarse[sigma]['stable_flags'], dtype=float)
        trace = np.cumsum(stable) / np.arange(1, stable.size + 1)
        ax.plot(times, trace, label=f'sigma={int(sigma)}')
    ax.set_xlabel('time')
    ax.set_ylabel('coarse persistence')
    ax.set_ylim(0.0, 1.05)
    ax.set_title(f'Persistence trace: {run["run_id"]}')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(persistence_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(persistence_path)

    return payload, row, plot_paths


def create_summary_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    colors = {
        'basin-dominated regime': 'tab:green',
        'multi-basin fluctuating regime': 'tab:orange',
        'diffuse coarse field regime': 'tab:blue',
    }
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for row in rows:
        ax.scatter(
            float(row['mean_dominant_area_fraction_sigma4']),
            float(row['coarse_persistence_sigma4']),
            color=colors.get(str(row['coarse_regime_label']), 'black'),
            s=90,
            alpha=0.85,
            label=str(row['coarse_regime_label']),
        )
        ax.text(
            float(row['mean_dominant_area_fraction_sigma4']) + 0.01,
            float(row['coarse_persistence_sigma4']),
            str(row['packet_count']),
            fontsize=8,
        )
    handles, labels = ax.get_legend_handles_labels()
    unique: dict[str, Any] = {}
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    ax.legend(unique.values(), unique.keys(), fontsize=8)
    ax.set_xlabel('mean dominant basin area fraction (sigma=4)')
    ax.set_ylabel('coarse persistence (sigma=4)')
    ax.set_title('Stage 12 coarse regime map')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(path: Path, json_rel: str, csv_rel: str, rows: list[dict[str, Any]], regime_counts: dict[str, int]) -> None:
    lines = [
        '# Stage 12 Coarse-Grained Field Structure v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        '',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Regime counts: {regime_counts}',
        '',
        'This pilot remains classification-first and interpretation-minimal.',
        'Stage 12 identifies coarse-grained organizational morphology and does not imply nonlinear effective interactions.',
        '',
        'Per-run summary:',
        '',
    ]
    for row in rows:
        lines.extend([
            f"- `{row['run_id']}`",
            f"  - label: `{row['coarse_regime_label']}`",
            f"  - mean dominant area fraction (sigma=4): `{float(row['mean_dominant_area_fraction_sigma4']):.4f}`",
            f"  - coarse persistence (sigma=4): `{float(row['coarse_persistence_sigma4']):.4f}`",
            f"  - basin lifetime (sigma=4): `{float(row['basin_lifetime_sigma4']):.4f}`",
            f"  - mean basin count (sigma=4): `{float(row['mean_basin_count_sigma4']):.4f}`",
        ])
    lines.extend([
        '',
        'Classification labels remain descriptive only:',
        '- basin-dominated regime',
        '- multi-basin fluctuating regime',
        '- diffuse coarse field regime',
        '',
        'Do not read these runs as force, curvature, or nonlinear interaction claims. They are coarse-grained superposition diagnostics only.',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    runs = selected_runs(runsheet['runs'], args.run_ids, args.max_runs)
    base = runsheet['base_seed_reference']
    stem = args.runsheet.stem
    experiment_slug = stem.replace('_runs', '')
    note_name = runsheet.get('note_name', 'Stage_12_Coarse_Grained_Field_Structure_v1.md')
    note_path = ATLAS_NOTES / note_name
    work_plot_dir = Path('/tmp') / f'haos_iip_{experiment_slug}_plots'
    work_plot_dir.mkdir(parents=True, exist_ok=True)

    payload_runs: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        payload, row, per_run_plots = run_case(run, defaults, base, work_plot_dir)
        payload_runs.append(payload)
        rows.append(row)
        plot_paths.extend(per_run_plots)

    regime_counts = dict(Counter(row['coarse_regime_label'] for row in rows))
    summary_plot = work_plot_dir / f'{experiment_slug}_coarse_regime_map.png'
    create_summary_plot(summary_plot, rows)
    plot_paths.append(summary_plot)

    result = {
        'stage': 'Stage 12',
        'description': runsheet['description'],
        'base_seed_reference': base,
        'labels': LABELS,
        'coarse_scales': list(runsheet.get('coarse_scales', [2, 4])),
        'alpha': float(runsheet.get('alpha', ALPHA)),
        'regime_counts': regime_counts,
        'runs': payload_runs,
        'conclusion': 'Stage 12 identifies coarse-grained organizational morphology on top of frozen linear superposition and does not imply nonlinear effective interactions.',
    }

    json_path, csv_path, stamped_plots, timestamp = save_atlas_payload(
        experiment_slug=experiment_slug,
        result=result,
        csv_rows=rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), rows, regime_counts)
    append_log(
        title=f'Stage 12 Coarse-Grained Field Structure ({json_path.stem})',
        config_summary=f"runs={len(rows)}, sector={base['operator_sector']}, boundary={base['boundary_type']}, sigmas={runsheet.get('coarse_scales', [2, 4])}, alpha={runsheet.get('alpha', ALPHA)}",
        result_path=json_path,
        stamped_plots=stamped_plots,
        observation=f'regime_counts={regime_counts}',
        conclusion='the Stage 12 pilot measures whether coarse-grained envelope basins persist on the frozen linear collective architecture',
    )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for rel in stamped_plots:
        print(rel)
    print(f'regime_counts={regime_counts}')


if __name__ == '__main__':
    main()
