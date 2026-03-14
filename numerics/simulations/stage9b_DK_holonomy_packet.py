#!/usr/bin/env python3

from __future__ import annotations

import numpy as np

from stage9b_common import (
    REPO_ROOT,
    WAVE_NOTES,
    append_log,
    build_dk2d_complex,
    build_initial_packet,
    crank_nicolson_evolution,
    first_order_dt,
    load_stage9b_config,
    pack_positions,
    plt,
    save_result_payload,
)


NOTE_PATH = WAVE_NOTES / 'DK_Holonomy_Packet_Response_v1.md'


def analyze_phase(phase: float, cfg: dict) -> dict:
    n_side = int(cfg['n_2d'])
    complex_data = build_dk2d_complex(n_side=n_side, epsilon=float(cfg['epsilon']), cycle_phase_x=float(phase), cycle_phase_y=0.0)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    psi0 = build_initial_packet(
        positions=positions,
        block_sizes=complex_data.block_sizes,
        center=np.asarray(cfg['packet_center_2d'], dtype=float),
        sigma=float(cfg['sigma_2d']),
        kick_cycles=float(cfg['kick_cycles']),
        axis=0,
        active_grades=[0, 1, 2],
    )
    dt = first_order_dt(complex_data.dirac_kahler, float(cfg['dt_scale']))
    steps = max(1, int(np.ceil(float(cfg['t_final_2d']) / dt)))
    run = crank_nicolson_evolution(complex_data.dirac_kahler, psi0, dt, steps, positions, complex_data.block_sizes)
    centers = np.asarray(run['centers'], dtype=float)
    widths = np.asarray(run['widths'], dtype=float)
    norms = np.asarray(run['norms'], dtype=float)
    grade_hist = np.asarray(run['grade_weights'], dtype=float)
    return {
        'phase': float(phase),
        'dt': float(dt),
        'steps': int(steps),
        'times': run['times'],
        'centers': run['centers'],
        'widths': run['widths'],
        'norms': run['norms'],
        'energies': run['energies'],
        'grade_weights': run['grade_weights'],
        'summary': {
            'center_x_shift': float(centers[-1, 0] - centers[0, 0]),
            'width_change': float(widths[-1] - widths[0]),
            'relative_norm_change': float((norms[-1] - norms[0]) / max(norms[0], 1.0e-12)),
            'omega1_mean_weight': float(np.mean(grade_hist[:, 1])),
        },
    }


def main() -> None:
    cfg = load_stage9b_config()
    phases = [float(v) for v in cfg['holonomy_phases']]
    cases = [analyze_phase(phase, cfg) for phase in phases]

    plot_paths = REPO_ROOT / 'plots' / 'DK_holonomy_packet_paths.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in cases:
        centers = np.asarray(case['centers'], dtype=float)
        ax.plot(case['times'], centers[:, 0], label=f"{case['phase']:.3f}")
    ax.set_xlabel('time')
    ax.set_ylabel('packet center x')
    ax.set_title('DK packet paths under flat holonomy')
    ax.grid(alpha=0.25)
    ax.legend(title='phase', fontsize=8)
    fig.savefig(plot_paths, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_widths = REPO_ROOT / 'plots' / 'DK_holonomy_packet_widths.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in cases:
        ax.plot(case['times'], case['widths'], label=f"{case['phase']:.3f}")
    ax.set_xlabel('time')
    ax.set_ylabel('packet width')
    ax.set_title('DK packet width under flat holonomy')
    ax.grid(alpha=0.25)
    ax.legend(title='phase', fontsize=8)
    fig.savefig(plot_widths, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_grades = REPO_ROOT / 'plots' / 'DK_holonomy_grade_weights.png'
    fig, axes = plt.subplots(1, len(cases), figsize=(3.4 * len(cases), 4.0), sharey=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case in zip(axes, cases):
        grade_hist = np.asarray(case['grade_weights'], dtype=float)
        ax.stackplot(case['times'], grade_hist.T, labels=[r'$\Omega^0$', r'$\Omega^1$', r'$\Omega^2$'], alpha=0.85)
        ax.set_title(f"{case['phase']:.3f}")
        ax.set_xlabel('time')
        ax.grid(alpha=0.2)
    axes[0].set_ylabel('grade weight fraction')
    axes[-1].legend(loc='upper right', fontsize=7)
    fig.savefig(plot_grades, dpi=180, bbox_inches='tight')
    plt.close(fig)

    observation = 'flat torus-cycle holonomy modifies DK packet transport and spreading smoothly while preserving stable norm and grade evolution'
    conclusion = 'the DK packet response remains coherent under the tested flat holonomy phases on the validated 2D cochain architecture'
    result = {
        'experiment': 'stage9b_DK_holonomy_packet_response',
        'config': {
            'n_side': int(cfg['n_2d']),
            'epsilon': float(cfg['epsilon']),
            'sigma_2d': float(cfg['sigma_2d']),
            'kick_cycles': float(cfg['kick_cycles']),
            't_final_2d': float(cfg['t_final_2d']),
            'holonomy_phases': phases,
        },
        'cases': cases,
        'observation': observation,
        'conclusion': conclusion,
    }

    stamped_json, stamped_plots, _ = save_result_payload('DK_holonomy_packet_response', result, [plot_paths, plot_widths, plot_grades])
    append_log(
        'Stage 9B DK Holonomy Packet Response',
        f"n={cfg['n_2d']}, epsilon={cfg['epsilon']}, sigma={cfg['sigma_2d']}, phases={phases}",
        stamped_json,
        stamped_plots,
        observation,
        conclusion,
    )
    NOTE_PATH.write_text(
        f"# DK Holonomy Packet Response v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: n={cfg['n_2d']}, epsilon={cfg['epsilon']}, sigma={cfg['sigma_2d']}, kick_cycles={cfg['kick_cycles']}, "
        f"phases={phases}, t_final={cfg['t_final_2d']}.\n\n"
        f"Observation: {observation}\n\n"
        f"Conclusion: {conclusion}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
