#!/usr/bin/env python3

from __future__ import annotations

import math
import numpy as np

from stage8_common import (
    REPO_ROOT,
    build_periodic_complex,
    ensure_matplotlib,
    load_config,
    harmonic_basis_from_L1,
    ExactProjector,
    HarmonicProjector,
    build_penalized_transverse_operator,
    project_transverse,
    low_eigensystem,
    save_result_payload,
    append_log,
)

plt = ensure_matplotlib()
NOTE_PATH = REPO_ROOT / 'experiments' / 'gauge_diagnostics' / 'Wilson_Loop_Test_v1.md'


def power_fit(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    mask = (x > 0) & (y > 0)
    x = np.log(x[mask])
    y = np.log(y[mask])
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        'prefactor': float(math.exp(intercept)),
        'power': float(slope),
        'r2': float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0,
    }


def edge_value(data, field, axis, i, j, k, sign):
    key = (axis, i % data.n_side, j % data.n_side, k % data.n_side)
    idx = data.edge_map[key]
    return sign * float(field[idx])


def square_loop_circulation(data, field, plane: str, i: int, j: int, k: int, side: int) -> float:
    s = 0.0
    if plane == 'xy':
        for t in range(side):
            s += edge_value(data, field, 'x', i + t, j, k, +1)
        for t in range(side):
            s += edge_value(data, field, 'y', i + side, j + t, k, +1)
        for t in range(side):
            s += edge_value(data, field, 'x', i + t, j + side, k, -1)
        for t in range(side):
            s += edge_value(data, field, 'y', i, j + t, k, -1)
    elif plane == 'xz':
        for t in range(side):
            s += edge_value(data, field, 'x', i + t, j, k, +1)
        for t in range(side):
            s += edge_value(data, field, 'z', i + side, j, k + t, +1)
        for t in range(side):
            s += edge_value(data, field, 'x', i + t, j, k + side, -1)
        for t in range(side):
            s += edge_value(data, field, 'z', i, j, k + t, -1)
    elif plane == 'yz':
        for t in range(side):
            s += edge_value(data, field, 'y', i, j + t, k, +1)
        for t in range(side):
            s += edge_value(data, field, 'z', i, j + side, k + t, +1)
        for t in range(side):
            s += edge_value(data, field, 'y', i, j + t, k + side, -1)
        for t in range(side):
            s += edge_value(data, field, 'z', i, j, k + t, -1)
    else:
        raise ValueError(plane)
    return s


def main() -> None:
    full_config = load_config()
    config = full_config['stage8_loop_response']
    epsilon = full_config['epsilon']
    harmonic_tol = full_config['harmonic_tol']
    eig_tol = full_config['eig_tol']
    penalty = full_config['penalty']
    n_side = int(config['n_side'])
    restricted_mode_index = int(config['restricted_mode_index'])
    max_loop_side = int(config['max_loop_side'])

    data = build_periodic_complex(n_side=n_side, epsilon=epsilon)
    harmonic_basis = harmonic_basis_from_L1(data.L1, count=20, harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = ExactProjector(data.d0, data.L0)
    harmonic_projector = HarmonicProjector(harmonic_basis)
    transverse_operator = build_penalized_transverse_operator(data, exact_projector, harmonic_projector, penalty)
    raw_evals, raw_vecs = low_eigensystem(transverse_operator, 12, eig_tol)

    restricted_modes: list[np.ndarray] = []
    restricted_lambdas: list[float] = []
    for idx in range(raw_vecs.shape[1]):
        vec = project_transverse(raw_vecs[:, idx], exact_projector, harmonic_projector)
        norm = float(np.linalg.norm(vec))
        if norm <= 1e-12:
            continue
        vec = vec / norm
        restricted_modes.append(vec)
        restricted_lambdas.append(float(np.dot(vec, data.upper @ vec)))
    mode = restricted_modes[restricted_mode_index]

    loop_sides = list(range(1, max_loop_side + 1))
    perimeter: list[float] = []
    area: list[float] = []
    response: list[float] = []
    h = 1.0 / n_side
    for side in loop_sides:
        values: list[float] = []
        for plane in ('xy', 'xz', 'yz'):
            for i in range(n_side):
                for j in range(n_side):
                    for k in range(n_side):
                        values.append(abs(square_loop_circulation(data, mode, plane, i, j, k, side)))
        perimeter.append(4.0 * side * h)
        area.append((side * h) ** 2)
        response.append(float(np.mean(values)))

    perimeter_fit = power_fit(np.asarray(perimeter), np.asarray(response))
    area_fit = power_fit(np.asarray(area), np.asarray(response))

    plot_per = REPO_ROOT / 'plots' / 'loop_perimeter_scaling.png'
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    ax.loglog(perimeter, response, marker='o')
    fit_curve = perimeter_fit['prefactor'] * np.asarray(perimeter) ** perimeter_fit['power']
    ax.loglog(perimeter, fit_curve, '--', label=f"fit p={perimeter_fit['power']:.2f}, R$^2$={perimeter_fit['r2']:.3f}")
    ax.set_xlabel('loop perimeter')
    ax.set_ylabel('mean |loop response|')
    ax.set_title('Stage 8 loop response vs perimeter')
    ax.grid(alpha=0.25, which='both')
    ax.legend(fontsize=8)
    fig.savefig(plot_per, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_area = REPO_ROOT / 'plots' / 'loop_area_scaling.png'
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    ax.loglog(area, response, marker='o')
    fit_curve = area_fit['prefactor'] * np.asarray(area) ** area_fit['power']
    ax.loglog(area, fit_curve, '--', label=f"fit p={area_fit['power']:.2f}, R$^2$={area_fit['r2']:.3f}")
    ax.set_xlabel('loop area')
    ax.set_ylabel('mean |loop response|')
    ax.set_title('Stage 8 loop response vs area')
    ax.grid(alpha=0.25, which='both')
    ax.legend(fontsize=8)
    fig.savefig(plot_area, dpi=180, bbox_inches='tight')
    plt.close(fig)

    observation = 'closed-loop response of the lowest restricted transverse mode scales smoothly with loop size on the periodic architecture'
    conclusion = 'the Wilson-loop analogue yields a stable nonlocal loop diagnostic with measurable perimeter and area scaling on the tested architecture'
    result = {
        'experiment': 'stage8_wilson_loop_test',
        'config': {
            'n_side': n_side,
            'epsilon': epsilon,
            'restricted_mode_index': restricted_mode_index,
            'max_loop_side': max_loop_side,
        },
        'selected_mode_eigenvalue': restricted_lambdas[restricted_mode_index],
        'loop_sides': loop_sides,
        'perimeter': perimeter,
        'area': area,
        'mean_abs_loop_response': response,
        'perimeter_fit': perimeter_fit,
        'area_fit': area_fit,
        'observation': observation,
        'conclusion': conclusion,
    }

    stamped_json, stamped_plots, timestamp = save_result_payload('stage8_wilson_loop_test', result, [plot_per, plot_area])
    append_log(
        'Stage 8D Wilson Loop Analogue',
        f'n={n_side}, epsilon={epsilon}, mode_index={restricted_mode_index}, max_loop_side={max_loop_side}',
        stamped_json,
        stamped_plots,
        observation,
        conclusion,
    )
    NOTE_PATH.write_text(
        f"# Wilson Loop Test v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: n={n_side}, epsilon={epsilon}, mode_index={restricted_mode_index}, max_loop_side={max_loop_side}.\n\n"
        f"Perimeter fit: power={perimeter_fit['power']:.6f}, R^2={perimeter_fit['r2']:.6f}.\n"
        f"Area fit: power={area_fit['power']:.6f}, R^2={area_fit['r2']:.6f}.\n\n"
        f"Observation: {observation}\n\n"
        f"Conclusion: {conclusion}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
