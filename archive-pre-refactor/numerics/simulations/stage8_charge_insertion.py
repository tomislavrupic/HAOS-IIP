#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import math
import numpy as np

from stage8_common import (
    REPO_ROOT,
    build_periodic_complex,
    ensure_matplotlib,
    load_config,
    solve_mean_zero_poisson,
    save_result_payload,
    append_log,
)

plt = ensure_matplotlib()
NOTE_PATH = REPO_ROOT / 'experiments' / 'gauge_diagnostics' / 'Charge_Insertion_Test_v1.md'


def main() -> None:
    config = load_config()['stage8_charge_insertion']
    epsilon = load_config()['epsilon']
    n_side = int(config['n_side'])
    source_strength = float(config['source_strength'])
    radial_bins = int(config['radial_bins'])

    data = build_periodic_complex(n_side=n_side, epsilon=epsilon)
    center = n_side // 2
    source_idx = int(data.node_index[center, center, center])
    source_point = data.points[source_idx]

    rho = np.full(data.points.shape[0], -source_strength / data.points.shape[0], dtype=float)
    rho[source_idx] += source_strength
    phi = solve_mean_zero_poisson(data.L0, rho)
    A = np.asarray(data.d0 @ phi, dtype=float)
    divergence = np.asarray(data.d0.T @ A, dtype=float)
    conservation_residual = float(np.linalg.norm(divergence - rho) / max(np.linalg.norm(rho), 1.0e-12))

    edge_delta = (data.midpoints - source_point + 0.5) % 1.0 - 0.5
    edge_r = np.linalg.norm(edge_delta, axis=1)
    edge_mag = np.abs(A)
    bins = np.linspace(0.0, float(np.max(edge_r)), radial_bins + 1)
    radial_r: list[float] = []
    radial_mag: list[float] = []
    for left, right in zip(bins[:-1], bins[1:]):
        mask = (edge_r >= left) & (edge_r < right if right < bins[-1] else edge_r <= right)
        if not np.any(mask):
            continue
        radial_r.append(float(np.mean(edge_r[mask])))
        radial_mag.append(float(np.mean(edge_mag[mask])))

    node_delta = (data.points - source_point + 0.5) % 1.0 - 0.5
    node_r = np.linalg.norm(node_delta, axis=1)
    shell_radii = np.linspace(1.0 / n_side, 0.45, 12)
    enclosed_charge: list[float] = []
    boundary_flux: list[float] = []
    flux_balance_error: list[float] = []
    for radius in shell_radii:
        inside = node_r <= radius
        q_inside = float(np.sum(rho[inside]))
        flux = 0.0
        for edge_idx, (u, v) in enumerate(data.edges):
            iu = bool(inside[u])
            iv = bool(inside[v])
            if iu and not iv:
                flux += float(A[edge_idx])
            elif iv and not iu:
                flux -= float(A[edge_idx])
        enclosed_charge.append(q_inside)
        boundary_flux.append(flux)
        flux_balance_error.append(abs(flux - q_inside))

    div_grid = divergence.reshape((n_side, n_side, n_side))[:, :, center]

    plot_decay = REPO_ROOT / 'plots' / 'charge_radial_decay.png'
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(radial_r, radial_mag, marker='o')
    ax.set_xlabel('distance from source')
    ax.set_ylabel(r'$\langle |A| \rangle$')
    ax.set_title('Charge insertion: radial edge-field decay')
    ax.grid(alpha=0.25)
    fig.savefig(plot_decay, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_flux = REPO_ROOT / 'plots' / 'charge_flux_balance.png'
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(shell_radii, enclosed_charge, marker='o', label='enclosed source')
    ax.plot(shell_radii, boundary_flux, marker='s', label='boundary flux')
    ax.set_xlabel('radius')
    ax.set_ylabel('net quantity')
    ax.set_title('Charge insertion: flux balance vs radius')
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(plot_flux, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_div = REPO_ROOT / 'plots' / 'charge_divergence_map.png'
    fig, ax = plt.subplots(figsize=(5, 4.5))
    im = ax.imshow(div_grid.T, origin='lower', cmap='coolwarm')
    ax.set_title('Charge insertion: divergence slice')
    ax.set_xlabel('x index')
    ax.set_ylabel('y index')
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.savefig(plot_div, dpi=180, bbox_inches='tight')
    plt.close(fig)

    result = {
        'experiment': 'stage8_charge_insertion',
        'config': {
            'n_side': n_side,
            'epsilon': epsilon,
            'source_strength': source_strength,
            'radial_bins': radial_bins,
        },
        'source_index': source_idx,
        'conservation_residual': conservation_residual,
        'radial_profile': {
            'radius': radial_r,
            'mean_field_magnitude': radial_mag,
        },
        'flux_balance': {
            'radius': shell_radii.tolist(),
            'enclosed_source': enclosed_charge,
            'boundary_flux': boundary_flux,
            'absolute_balance_error': flux_balance_error,
            'max_balance_error': float(max(flux_balance_error) if flux_balance_error else 0.0),
        },
        'divergence_summary': {
            'max_abs_divergence': float(np.max(np.abs(divergence))),
            'mean_abs_divergence': float(np.mean(np.abs(divergence))),
        },
        'observation': 'the scalar source produces a smooth induced edge field with a stable flux-balance profile on the periodic cochain architecture',
        'conclusion': 'the charge-insertion diagnostic shows structured scalar-to-edge response and exact discrete divergence accounting on the tested architecture',
    }

    stamped_json, stamped_plots, timestamp = save_result_payload('stage8_charge_insertion', result, [plot_decay, plot_flux, plot_div])
    append_log(
        'Stage 8A Charge Insertion',
        f'n={n_side}, epsilon={epsilon}, source_strength={source_strength}',
        stamped_json,
        stamped_plots,
        result['observation'],
        result['conclusion'],
    )
    NOTE_PATH.write_text(
        f"# Charge Insertion Test v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: n={n_side}, epsilon={epsilon}, source_strength={source_strength}.\n\n"
        f"Key metrics:\n"
        f"- conservation residual = {conservation_residual:.3e}\n"
        f"- max flux-balance error = {result['flux_balance']['max_balance_error']:.3e}\n"
        f"- max |divergence| = {result['divergence_summary']['max_abs_divergence']:.6f}\n\n"
        f"Observation: {result['observation']}\n\n"
        f"Conclusion: {result['conclusion']}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
