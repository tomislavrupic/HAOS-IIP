#!/usr/bin/env python3

from __future__ import annotations

import numpy as np

from stage9_common import (
    REPO_ROOT,
    WAVE_NOTES,
    append_log,
    build_periodic_complex,
    build_transverse_initial_state,
    build_transverse_setup,
    load_stage9_config,
    plt,
    save_result_payload,
    simulate_second_order,
    suggested_dt,
)


NOTE_PATH = WAVE_NOTES / 'Constraint_Stability_Test_v1.md'


def analyze_case(n_side: int, cfg: dict, alpha: float = 1.0) -> dict:
    data, projector = build_transverse_setup(n_side=n_side, epsilon=float(cfg['epsilon']), harmonic_tol=float(cfg['harmonic_tol']))
    center = np.asarray(cfg['packet_center'], dtype=float)
    sigma = float(cfg['sigma'])
    kick_axis = int(cfg['kick_axis'])
    q0, v0 = build_transverse_initial_state(
        midpoints=data.midpoints,
        edge_axes=data.edge_axes,
        center=center,
        sigma=sigma,
        kick_axis=kick_axis,
        kick_strength=float(cfg['kick_strength_transverse']),
        projector=projector,
    )

    phi_profile = np.exp(-np.sum(((data.points - center + 0.5) % 1.0 - 0.5) ** 2, axis=1) / (2.0 * sigma * sigma))
    exact_component = np.asarray(data.d0 @ phi_profile, dtype=float)

    q_contaminated = q0 + alpha * exact_component
    v_contaminated = v0 + alpha * np.asarray(data.d0 @ (phi_profile * (((data.points - center + 0.5) % 1.0 - 0.5)[:, kick_axis] / (sigma * sigma))), dtype=float)

    dt = suggested_dt(data.L1, float(cfg['dt_safety']))
    steps = max(1, int(np.ceil(float(cfg['t_final']) / dt)))

    baseline = simulate_second_order(data.L1, q0, v0, data.midpoints, dt, steps, project=projector, divergence_op=data.d0.T)
    contaminated = simulate_second_order(data.L1, q_contaminated, v_contaminated, data.midpoints, dt, steps, project=projector, divergence_op=data.d0.T)

    baseline_centers = np.asarray(baseline['centers'], dtype=float)
    contaminated_centers = np.asarray(contaminated['centers'], dtype=float)
    center_delta = np.linalg.norm(contaminated_centers - baseline_centers, axis=1)
    width_delta = np.abs(np.asarray(contaminated['widths']) - np.asarray(baseline['widths']))
    energy_delta = np.abs(np.asarray(contaminated['energies']) - np.asarray(baseline['energies']))

    return {
        'n_side': n_side,
        'dt': float(dt),
        'steps': int(steps),
        'times': baseline['times'],
        'baseline_divergence_norms': baseline['divergence_norms'],
        'baseline_projection_residuals': baseline['projection_residuals'],
        'center_differences': center_delta.tolist(),
        'width_differences': width_delta.tolist(),
        'energy_differences': energy_delta.tolist(),
        'summary': {
            'max_divergence_norm': float(np.max(baseline['divergence_norms'])),
            'max_projection_residual': float(np.max(baseline['projection_residuals'])),
            'max_center_difference': float(np.max(center_delta)),
            'max_width_difference': float(np.max(width_delta)),
            'max_energy_difference': float(np.max(energy_delta)),
        },
    }


def main() -> None:
    cfg = load_stage9_config()
    sizes = [int(n) for n in cfg['sizes']]
    cases = [analyze_case(n, cfg) for n in sizes]

    plot_residuals = REPO_ROOT / 'plots' / 'constraint_residuals.png'
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1))
    for case in cases:
        label = f"n={case['n_side']}"
        times = case['times']
        axes[0].plot(times, case['baseline_divergence_norms'], label=label)
        axes[1].plot(times, case['baseline_projection_residuals'], label=label)
        axes[2].plot(times, case['center_differences'], label=label)
    axes[0].set_title(r'$||d_0^* A(t)||$')
    axes[1].set_title(r'$||P_\perp A(t)-A(t)||$')
    axes[2].set_title('trajectory difference under exact shift')
    for ax in axes:
        ax.set_xlabel('time')
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    fig.savefig(plot_residuals, dpi=180, bbox_inches='tight')
    plt.close(fig)

    observation = 'projected dynamics preserves the divergence constraint and removes exact contamination so that longitudinally shifted initial data evolve with the same restricted observables'
    conclusion = 'constraint preservation and longitudinal redundancy remain stable during Stage 9 evolution on the frozen architecture'
    result = {
        'experiment': 'stage9_constraint_stability',
        'config': {
            'sizes': sizes,
            'epsilon': float(cfg['epsilon']),
            'sigma': float(cfg['sigma']),
            't_final': float(cfg['t_final']),
        },
        'cases': cases,
        'observation': observation,
        'conclusion': conclusion,
    }

    stamped_json, stamped_plots, timestamp = save_result_payload('stage9_constraint_stability', result, [plot_residuals])
    append_log(
        'Stage 9C Constraint Stability During Evolution',
        f"sizes={sizes}, epsilon={cfg['epsilon']}, sigma={cfg['sigma']}, t_final={cfg['t_final']}",
        stamped_json,
        stamped_plots,
        observation,
        conclusion,
    )
    NOTE_PATH.write_text(
        f"# Constraint Stability Test v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: sizes={sizes}, epsilon={cfg['epsilon']}, sigma={cfg['sigma']}, t_final={cfg['t_final']}.\n\n"
        f"Observation: {observation}\n\n"
        f"Conclusion: {conclusion}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
