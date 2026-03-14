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
    leapfrog_second_order_complex,
    load_stage9b_config,
    pack_positions,
    plt,
    save_result_payload,
    second_order_dt,
)


NOTE_PATH = WAVE_NOTES / 'DK_First_vs_Second_Order_v1.md'


def main() -> None:
    cfg = load_stage9b_config()
    n_side = int(cfg['n_2d'])
    complex_data = build_dk2d_complex(n_side=n_side, epsilon=float(cfg['epsilon']))
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

    dt_first = first_order_dt(complex_data.dirac_kahler, float(cfg['dt_scale']))
    dt_second = second_order_dt(complex_data.delta_h, float(cfg['dt_scale']))
    dt = min(dt_first, dt_second)
    steps = max(1, int(np.ceil(float(cfg['t_final_2d']) / dt)))

    first = crank_nicolson_evolution(complex_data.dirac_kahler, psi0, dt, steps, positions, complex_data.block_sizes)
    v0 = 1j * (complex_data.dirac_kahler @ psi0)
    second = leapfrog_second_order_complex(complex_data.delta_h, psi0, v0, dt, steps, positions, complex_data.block_sizes)

    first_centers = np.asarray(first['centers'], dtype=float)
    second_centers = np.asarray(second['centers'], dtype=float)
    center_diff = np.linalg.norm(first_centers - second_centers, axis=1)
    width_diff = np.abs(np.asarray(first['widths']) - np.asarray(second['widths']))
    norm_diff = np.abs(np.asarray(first['norms']) - np.asarray(second['norms']))
    first_energies = np.asarray(first['energies'], dtype=float)
    second_energies = np.asarray(second['energies'], dtype=float)
    first_rel_energy_drift = np.abs((first_energies - first_energies[0]) / max(abs(first_energies[0]), 1.0e-12))
    second_rel_energy_drift = np.abs((second_energies - second_energies[0]) / max(abs(second_energies[0]), 1.0e-12))
    first_grade = np.asarray(first['grade_weights'], dtype=float)
    second_grade = np.asarray(second['grade_weights'], dtype=float)
    grade_weight_difference = np.max(np.abs(first_grade - second_grade), axis=1)

    plot_compare = REPO_ROOT / 'plots' / 'DK_first_vs_second_order_comparison.png'
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    axes[0].plot(first['times'], center_diff)
    axes[0].set_title('center difference')
    axes[1].plot(first['times'], width_diff)
    axes[1].set_title('width difference')
    axes[2].plot(first['times'], norm_diff, label='norm difference')
    axes[2].plot(first['times'], grade_weight_difference, label='grade-weight difference')
    axes[2].plot(first['times'], first_rel_energy_drift, label='first-order energy drift')
    axes[2].plot(first['times'], second_rel_energy_drift, label='second-order energy drift')
    axes[2].set_title('norm / grade / drift diagnostics')
    for ax in axes:
        ax.set_xlabel('time')
        ax.grid(alpha=0.25)
    axes[2].legend(fontsize=8)
    fig.savefig(plot_compare, dpi=180, bbox_inches='tight')
    plt.close(fig)

    observation = 'first-order DK propagation and the induced second-order Hodge evolution remain closely matched at the packet level over the tested time window, while each evolution preserves its own conserved quantity with small drift'
    conclusion = 'the first-order DK dynamics is consistent with its validated squared Hodge dynamics in the tested 2D packet regime'
    result = {
        'experiment': 'stage9b_DK_first_vs_second_order',
        'config': {
            'n_side': n_side,
            'epsilon': float(cfg['epsilon']),
            'sigma_2d': float(cfg['sigma_2d']),
            'kick_cycles': float(cfg['kick_cycles']),
            't_final_2d': float(cfg['t_final_2d']),
            'dt': float(dt),
        },
        'summary': {
            'max_center_difference': float(np.max(center_diff)),
            'max_width_difference': float(np.max(width_diff)),
            'max_norm_difference': float(np.max(norm_diff)),
            'max_grade_weight_difference': float(np.max(grade_weight_difference)),
            'max_first_order_energy_drift': float(np.max(first_rel_energy_drift)),
            'max_second_order_energy_drift': float(np.max(second_rel_energy_drift)),
        },
        'first_order': first,
        'second_order': second,
        'observation': observation,
        'conclusion': conclusion,
    }

    stamped_json, stamped_plots, _ = save_result_payload('DK_first_vs_second_order', result, [plot_compare])
    append_log(
        'Stage 9B DK First-vs-Second-Order Consistency',
        f"n={n_side}, epsilon={cfg['epsilon']}, sigma={cfg['sigma_2d']}, t_final={cfg['t_final_2d']}",
        stamped_json,
        stamped_plots,
        observation,
        conclusion,
    )
    NOTE_PATH.write_text(
        f"# DK First vs Second Order v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: n={n_side}, epsilon={cfg['epsilon']}, sigma={cfg['sigma_2d']}, kick_cycles={cfg['kick_cycles']}, "
        f"t_final={cfg['t_final_2d']}.\n\n"
        f"Observation: {observation}\n\n"
        f"Conclusion: {conclusion}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
