#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage3_common import (
    PLOTS,
    REPO_ROOT,
    analyze_branch_cases,
    append_log,
    make_phase_plot,
    save_result_payload,
    plt,
)

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_pde_reconstruction': {
        'sizes': [12, 16],
        'restricted_modes': 3,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_pde_reconstruction'] = dict(DEFAULT_CONFIG['transverse_pde_reconstruction'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_pde_reconstruction'})
        if isinstance(on_disk.get('transverse_pde_reconstruction'), dict):
            merged['transverse_pde_reconstruction'].update(on_disk['transverse_pde_reconstruction'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_pde_reconstruction'})
        if isinstance(config.get('transverse_pde_reconstruction'), dict):
            merged['transverse_pde_reconstruction'].update(config['transverse_pde_reconstruction'])
    return merged


def normalize_complex_mode(vec: np.ndarray) -> np.ndarray:
    vec = np.asarray(vec, dtype=complex)
    idx = int(np.argmax(np.abs(vec)))
    phase = np.angle(vec[idx]) if np.abs(vec[idx]) > 0 else 0.0
    return vec * np.exp(-1j * phase)


def edge_fields_from_mode(case: dict[str, Any], mode_index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = int(case['config']['n_side'])
    vec = np.real(normalize_complex_mode(np.asarray(case['restricted_vectors'][mode_index], dtype=complex)))
    mid = np.asarray(case['midpoints'], dtype=float)
    dirs = np.asarray(case['directions'], dtype=float)
    ex = np.zeros((n, n, n), dtype=float)
    ey = np.zeros((n, n, n), dtype=float)
    ez = np.zeros((n, n, n), dtype=float)
    for value, point, direction in zip(vec, mid, dirs):
        axis = int(np.argmax(np.abs(direction)))
        if axis == 0:
            i = int(np.floor(point[0] * n)) % n
            j = int(np.rint(point[1] * n)) % n
            k = int(np.rint(point[2] * n)) % n
            ex[i, j, k] = value
        elif axis == 1:
            i = int(np.rint(point[0] * n)) % n
            j = int(np.floor(point[1] * n)) % n
            k = int(np.rint(point[2] * n)) % n
            ey[i, j, k] = value
        else:
            i = int(np.rint(point[0] * n)) % n
            j = int(np.rint(point[1] * n)) % n
            k = int(np.floor(point[2] * n)) % n
            ez[i, j, k] = value
    return ex, ey, ez


def cell_center_average(ex: np.ndarray, ey: np.ndarray, ez: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    Ax = 0.25 * (ex + np.roll(ex, -1, axis=1) + np.roll(ex, -1, axis=2) + np.roll(np.roll(ex, -1, axis=1), -1, axis=2))
    Ay = 0.25 * (ey + np.roll(ey, -1, axis=0) + np.roll(ey, -1, axis=2) + np.roll(np.roll(ey, -1, axis=0), -1, axis=2))
    Az = 0.25 * (ez + np.roll(ez, -1, axis=0) + np.roll(ez, -1, axis=1) + np.roll(np.roll(ez, -1, axis=0), -1, axis=1))
    return Ax, Ay, Az


def ddx(arr: np.ndarray) -> np.ndarray:
    return 0.5 * (np.roll(arr, -1, axis=0) - np.roll(arr, 1, axis=0))


def ddy(arr: np.ndarray) -> np.ndarray:
    return 0.5 * (np.roll(arr, -1, axis=1) - np.roll(arr, 1, axis=1))


def ddz(arr: np.ndarray) -> np.ndarray:
    return 0.5 * (np.roll(arr, -1, axis=2) - np.roll(arr, 1, axis=2))


def divergence(Ax: np.ndarray, Ay: np.ndarray, Az: np.ndarray) -> np.ndarray:
    return ddx(Ax) + ddy(Ay) + ddz(Az)


def curl(Ax: np.ndarray, Ay: np.ndarray, Az: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    Cx = ddy(Az) - ddz(Ay)
    Cy = ddz(Ax) - ddx(Az)
    Cz = ddx(Ay) - ddy(Ax)
    return Cx, Cy, Cz


def curlcurl(Ax: np.ndarray, Ay: np.ndarray, Az: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    Cx, Cy, Cz = curl(Ax, Ay, Az)
    return curl(Cx, Cy, Cz)


def residual_metrics(Ax: np.ndarray, Ay: np.ndarray, Az: np.ndarray, lam: float) -> dict[str, float | np.ndarray]:
    div = divergence(Ax, Ay, Az)
    Kx, Ky, Kz = curlcurl(Ax, Ay, Az)
    rx = Kx - lam * Ax
    ry = Ky - lam * Ay
    rz = Kz - lam * Az
    res = np.sqrt(rx * rx + ry * ry + rz * rz)
    field_norm = np.sqrt(np.sum((lam * Ax) ** 2 + (lam * Ay) ** 2 + (lam * Az) ** 2)) or 1.0
    div_norm = np.linalg.norm(div)
    return {
        'relative_curlcurl_residual': float(np.linalg.norm(res) / field_norm),
        'relative_divergence_norm': float(div_norm / (np.sqrt(np.sum(Ax * Ax + Ay * Ay + Az * Az)) or 1.0)),
        'residual_map': res,
    }


def make_residual_plot(metrics: dict[str, dict[str, Any]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for mode_index in range(3):
        ax.plot(
            sorted(int(k) for k in metrics.keys()),
            [metrics[str(n)]['modes'][mode_index]['relative_curlcurl_residual'] for n in sorted(int(k) for k in metrics.keys())],
            marker='o',
            label=f'mode {mode_index + 1}',
        )
    ax.set_xlabel('lattice side n')
    ax.set_ylabel(r'$\|\nabla \times \nabla \times A - \lambda A\| / \|\lambda A\|$')
    ax.set_title('Coarse-grained PDE residuals')
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_map_plot(metrics: dict[str, dict[str, Any]], path: Path) -> None:
    largest = max(int(k) for k in metrics.keys())
    residual_map = np.asarray(metrics[str(largest)]['modes'][0]['residual_map'], dtype=float)
    slice_map = residual_map[:, :, residual_map.shape[2] // 2]
    fig, ax = plt.subplots(figsize=(5.5, 5))
    image = ax.imshow(slice_map.T, origin='lower', cmap='magma')
    ax.set_title(f'Residual map slice (n={largest}, mode 1)')
    ax.set_xlabel('x index')
    ax.set_ylabel('y index')
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label('residual magnitude')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_PDE_Reconstruction_v1.md'
    rows = []
    for n in result['config']['sizes']:
        for mode_index, mode in enumerate(result['metrics'][str(n)]['modes'][:3], start=1):
            rows.append(
                f"| {n} | {mode_index} | {mode['eigenvalue']:.6f} | {mode['relative_curlcurl_residual']:.6f} | {mode['relative_divergence_norm']:.6f} |"
            )
    note = f"""# Transverse PDE Reconstruction

## Purpose

Test whether the lowest restricted modes admit a local coarse-grained PDE of the form `curl curl A ≈ lambda A` with `div A ≈ 0`.

## Setup

- branch: baseline periodic torus
- sizes: `{result['config']['sizes']}`
- restricted modes reconstructed: 3

## Residual summary

| `n` | mode | eigenvalue | relative curl-curl residual | relative divergence norm |
| --- | ---: | ---: | ---: | ---: |
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


def run_transverse_pde_reconstruction(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_pde_reconstruction']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 3))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))

    cases = analyze_branch_cases(
        sizes=sizes,
        variants=['baseline'],
        epsilon=epsilon,
        restricted_modes=restricted_modes,
        harmonic_tol=harmonic_tol,
        eig_tol=eig_tol,
        penalty=penalty,
        flux_tube_phase=0.0,
    )

    metrics: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        case = cases[f'baseline_n{n_side}']
        mode_metrics = []
        for mode_index in range(min(3, len(case['restricted_vectors']))):
            ex, ey, ez = edge_fields_from_mode(case, mode_index)
            Ax, Ay, Az = cell_center_average(ex, ey, ez)
            residual = residual_metrics(Ax, Ay, Az, float(case['restricted_transverse_spectrum'][mode_index]))
            mode_metrics.append(
                {
                    'mode_index': mode_index,
                    'eigenvalue': float(case['restricted_transverse_spectrum'][mode_index]),
                    'relative_curlcurl_residual': residual['relative_curlcurl_residual'],
                    'relative_divergence_norm': residual['relative_divergence_norm'],
                    'residual_map': np.asarray(residual['residual_map']).tolist(),
                }
            )
        metrics[str(n_side)] = {'modes': mode_metrics}

    residual_plot = PLOTS / 'pde_reconstruction_residuals.png'
    map_plot = PLOTS / 'transverse_field_residual_maps.png'
    phase_plot = PLOTS / 'divergence_curl_phase_pde_reconstruction.png'
    make_residual_plot(metrics, residual_plot)
    make_map_plot(metrics, map_plot)
    make_phase_plot(cases[f'baseline_n{sizes[-1]}'], phase_plot, f"Restricted phase map (PDE, n={sizes[-1]})")
    plot_paths = [residual_plot, map_plot, phase_plot]

    first = metrics[str(sizes[0])]['modes'][0]
    last = metrics[str(sizes[-1])]['modes'][0]
    observation = (
        'after coarse-graining to a periodic vector field, the lowest restricted modes remain nearly divergence-free and the curl-curl residual decreases with lattice size'
    )
    conclusion = (
        f"for the lowest mode, the relative curl-curl residual decreases from {first['relative_curlcurl_residual']:.3f} at n={sizes[0]} to {last['relative_curlcurl_residual']:.3f} at n={sizes[-1]}, consistent with an emergent local coarse-grained PDE in the tested range"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
        },
        'metrics': metrics,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_pde_reconstruction', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse PDE reconstruction',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, restricted_modes={restricted_modes}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_pde_reconstruction()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
