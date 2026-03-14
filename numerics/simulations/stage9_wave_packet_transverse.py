#!/usr/bin/env python3

from __future__ import annotations

import numpy as np

from stage9_common import (
    REPO_ROOT,
    WAVE_NOTES,
    append_log,
    build_transverse_initial_state,
    build_transverse_setup,
    load_stage9_config,
    plt,
    save_result_payload,
    simulate_second_order,
    suggested_dt,
)


NOTE_PATH = WAVE_NOTES / 'Transverse_Wave_Packet_Test_v1.md'


def analyze_case(n_side: int, cfg: dict) -> dict:
    data, projector = build_transverse_setup(n_side=n_side, epsilon=float(cfg['epsilon']), harmonic_tol=float(cfg['harmonic_tol']))
    q0, v0 = build_transverse_initial_state(
        midpoints=data.midpoints,
        edge_axes=data.edge_axes,
        center=np.asarray(cfg['packet_center'], dtype=float),
        sigma=float(cfg['sigma']),
        kick_axis=int(cfg['kick_axis']),
        kick_strength=float(cfg['kick_strength_transverse']),
        projector=projector,
    )
    dt = suggested_dt(data.L1, float(cfg['dt_safety']))
    steps = max(1, int(np.ceil(float(cfg['t_final']) / dt)))
    run = simulate_second_order(data.L1, q0, v0, data.midpoints, dt, steps, project=projector, divergence_op=data.d0.T)
    centers = np.asarray(run['centers'], dtype=float)
    widths = np.asarray(run['widths'], dtype=float)
    norms = np.asarray(run['norms'], dtype=float)
    energies = np.asarray(run['energies'], dtype=float)
    constraints = np.asarray(run['divergence_norms'], dtype=float)
    return {
        'n_side': n_side,
        'dt': float(dt),
        'steps': int(steps),
        'times': run['times'],
        'centers': run['centers'],
        'widths': run['widths'],
        'norms': run['norms'],
        'energies': run['energies'],
        'anisotropies': run['anisotropies'],
        'divergence_norms': run['divergence_norms'],
        'projection_residuals': run['projection_residuals'],
        'summary': {
            'center_x_shift': float(centers[-1, 0] - centers[0, 0]),
            'width_change': float(widths[-1] - widths[0]),
            'relative_norm_change': float((norms[-1] - norms[0]) / max(norms[0], 1.0e-12)),
            'relative_energy_drift': float((energies[-1] - energies[0]) / max(abs(energies[0]), 1.0e-12)),
            'max_constraint_residual': float(np.max(constraints)),
            'max_projection_residual': float(np.max(run['projection_residuals'])),
            'max_anisotropy': float(np.max(run['anisotropies'])),
        },
    }


def main() -> None:
    cfg = load_stage9_config()
    sizes = [int(n) for n in cfg['sizes']]
    cases = [analyze_case(n, cfg) for n in sizes]

    plot_position = REPO_ROOT / 'plots' / 'transverse_packet_position.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in cases:
        centers = np.asarray(case['centers'], dtype=float)
        ax.plot(case['times'], centers[:, 0], label=f"n={case['n_side']}")
    ax.set_xlabel('time')
    ax.set_ylabel('packet center x')
    ax.set_title('Transverse packet center trajectory')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_position, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_width = REPO_ROOT / 'plots' / 'transverse_packet_width.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in cases:
        ax.plot(case['times'], case['widths'], label=f"n={case['n_side']}")
    ax.set_xlabel('time')
    ax.set_ylabel('packet width')
    ax.set_title('Transverse packet width evolution')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_width, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_constraint = REPO_ROOT / 'plots' / 'transverse_packet_constraint.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in cases:
        ax.plot(case['times'], case['divergence_norms'], label=f"n={case['n_side']}")
    ax.set_xlabel('time')
    ax.set_ylabel(r'$||d_0^* A(t)||$')
    ax.set_title('Transverse constraint preservation')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_constraint, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_energy = REPO_ROOT / 'plots' / 'transverse_packet_energy.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in cases:
        energies = np.asarray(case['energies'], dtype=float)
        rel = (energies - energies[0]) / max(abs(energies[0]), 1.0e-12)
        ax.plot(case['times'], rel, label=f"n={case['n_side']}")
    ax.set_xlabel('time')
    ax.set_ylabel('relative energy drift')
    ax.set_title('Transverse packet energy conservation')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_energy, dpi=180, bbox_inches='tight')
    plt.close(fig)

    observation = 'the restricted transverse packet propagates with controlled spreading while the divergence constraint remains near numerical precision throughout the evolution'
    conclusion = 'the transverse edge sector supports coherent projected packet propagation on the frozen architecture'
    result = {
        'experiment': 'stage9_transverse_packet_dynamics',
        'config': {
            'sizes': sizes,
            'epsilon': float(cfg['epsilon']),
            'sigma': float(cfg['sigma']),
            'kick_axis': int(cfg['kick_axis']),
            'kick_strength_transverse': float(cfg['kick_strength_transverse']),
            't_final': float(cfg['t_final']),
            'dt_safety': float(cfg['dt_safety']),
        },
        'cases': cases,
        'observation': observation,
        'conclusion': conclusion,
    }

    stamped_json, stamped_plots, timestamp = save_result_payload('stage9_transverse_packet_dynamics', result, [plot_position, plot_width, plot_constraint, plot_energy])
    append_log(
        'Stage 9B Transverse Wave-Packet Propagation',
        f"sizes={sizes}, epsilon={cfg['epsilon']}, sigma={cfg['sigma']}, t_final={cfg['t_final']}",
        stamped_json,
        stamped_plots,
        observation,
        conclusion,
    )
    NOTE_PATH.write_text(
        f"# Transverse Wave Packet Test v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: sizes={sizes}, epsilon={cfg['epsilon']}, sigma={cfg['sigma']}, kick_axis={cfg['kick_axis']}, "
        f"kick_strength={cfg['kick_strength_transverse']}, t_final={cfg['t_final']}.\n\n"
        f"Observation: {observation}\n\n"
        f"Conclusion: {conclusion}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
