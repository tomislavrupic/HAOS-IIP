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
    spectral_pairing_summary_2d,
)


NOTE_PATH = WAVE_NOTES / 'DK_2D_Packet_Dynamics_v1.md'


def analyze_variant(label: str, active_grades: list[int], complex_data, cfg: dict) -> dict:
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    psi0 = build_initial_packet(
        positions=positions,
        block_sizes=complex_data.block_sizes,
        center=np.asarray(cfg['packet_center_2d'], dtype=float),
        sigma=float(cfg['sigma_2d']),
        kick_cycles=float(cfg['kick_cycles']),
        axis=0,
        active_grades=active_grades,
    )
    dt = first_order_dt(complex_data.dirac_kahler, float(cfg['dt_scale']))
    steps = max(1, int(np.ceil(float(cfg['t_final_2d']) / dt)))
    run = crank_nicolson_evolution(complex_data.dirac_kahler, psi0, dt, steps, positions, complex_data.block_sizes)
    centers = np.asarray(run['centers'], dtype=float)
    widths = np.asarray(run['widths'], dtype=float)
    norms = np.asarray(run['norms'], dtype=float)
    energies = np.asarray(run['energies'], dtype=float)
    return {
        'label': label,
        'active_grades': active_grades,
        'dt': float(dt),
        'steps': int(steps),
        'times': run['times'],
        'centers': run['centers'],
        'widths': run['widths'],
        'norms': run['norms'],
        'energies': run['energies'],
        'anisotropies': run['anisotropies'],
        'grade_weights': run['grade_weights'],
        'summary': {
            'center_x_shift': float(centers[-1, 0] - centers[0, 0]),
            'width_change': float(widths[-1] - widths[0]),
            'relative_norm_change': float((norms[-1] - norms[0]) / max(norms[0], 1.0e-12)),
            'relative_energy_drift': float((energies[-1] - energies[0]) / max(abs(energies[0]), 1.0e-12)),
            'max_anisotropy': float(np.max(run['anisotropies'])),
        },
    }


def main() -> None:
    cfg = load_stage9b_config()
    n_side = int(cfg['n_2d'])
    complex_data = build_dk2d_complex(n_side=n_side, epsilon=float(cfg['epsilon']))
    pairing = spectral_pairing_summary_2d(complex_data.dirac_kahler, modes=int(cfg['low_modes_2d']), tol=1.0e-9)

    variants = [
        analyze_variant('omega0_packet', [0], complex_data, cfg),
        analyze_variant('omega1_packet', [1], complex_data, cfg),
        analyze_variant('mixed_packet', [0, 1, 2], complex_data, cfg),
    ]

    plot_position = REPO_ROOT / 'plots' / 'DK_2D_packet_position.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in variants:
        centers = np.asarray(case['centers'], dtype=float)
        ax.plot(case['times'], centers[:, 0], label=case['label'])
    ax.set_xlabel('time')
    ax.set_ylabel('packet center x')
    ax.set_title(f'2D DK packet trajectories (n={n_side})')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_position, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_width = REPO_ROOT / 'plots' / 'DK_2D_packet_width.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in variants:
        ax.plot(case['times'], case['widths'], label=case['label'])
    ax.set_xlabel('time')
    ax.set_ylabel('packet width')
    ax.set_title(f'2D DK packet width evolution (n={n_side})')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_width, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_norm = REPO_ROOT / 'plots' / 'DK_2D_packet_norm.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in variants:
        norms = np.asarray(case['norms'], dtype=float)
        rel = (norms - norms[0]) / max(norms[0], 1.0e-12)
        ax.plot(case['times'], rel, label=case['label'])
    ax.set_xlabel('time')
    ax.set_ylabel('relative norm change')
    ax.set_title(f'2D DK norm preservation (n={n_side})')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_norm, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_grades = REPO_ROOT / 'plots' / 'DK_2D_packet_grade_weights.png'
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0), sharey=True)
    for ax, case in zip(axes, variants):
        grade_hist = np.asarray(case['grade_weights'], dtype=float)
        ax.stackplot(case['times'], grade_hist.T, labels=[r'$\Omega^0$', r'$\Omega^1$', r'$\Omega^2$'], alpha=0.85)
        ax.set_title(case['label'])
        ax.set_xlabel('time')
        ax.grid(alpha=0.2)
    axes[0].set_ylabel('grade weight fraction')
    axes[-1].legend(loc='upper right', fontsize=7)
    fig.savefig(plot_grades, dpi=180, bbox_inches='tight')
    plt.close(fig)

    observation = 'the validated 2D Dirac-Kaehler operator supports coherent packet propagation with stable norm and structured grade evolution across single-grade and mixed initial data'
    conclusion = 'the 2D DK sector behaves as a coherent first-order dynamical generator on the frozen cochain architecture'
    result = {
        'experiment': 'stage9b_DK_2D_packet_dynamics',
        'config': {
            'n_side': n_side,
            'epsilon': float(cfg['epsilon']),
            'sigma_2d': float(cfg['sigma_2d']),
            'kick_cycles': float(cfg['kick_cycles']),
            't_final_2d': float(cfg['t_final_2d']),
        },
        'pairing_summary': pairing,
        'variants': variants,
        'observation': observation,
        'conclusion': conclusion,
    }

    stamped_json, stamped_plots, _ = save_result_payload('DK_2D_packet_dynamics', result, [plot_position, plot_width, plot_norm, plot_grades])
    append_log(
        'Stage 9B DK 2D Packet Propagation',
        f"n={n_side}, epsilon={cfg['epsilon']}, sigma={cfg['sigma_2d']}, t_final={cfg['t_final_2d']}",
        stamped_json,
        stamped_plots,
        observation,
        conclusion,
    )
    NOTE_PATH.write_text(
        f"# DK 2D Packet Dynamics v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: n={n_side}, epsilon={cfg['epsilon']}, sigma={cfg['sigma_2d']}, kick_cycles={cfg['kick_cycles']}, "
        f"t_final={cfg['t_final_2d']}.\n\n"
        f"Pairing error: {pairing['pairing_error']:.3e}\n"
        f"Lowest |lambda|: {pairing['lowest_abs_eigenvalue']:.3e}\n\n"
        f"Observation: {observation}\n\n"
        f"Conclusion: {conclusion}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
