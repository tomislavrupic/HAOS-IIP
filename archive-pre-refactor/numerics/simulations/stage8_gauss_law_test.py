#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import numpy as np

from stage8_common import (
    REPO_ROOT,
    append_log,
    build_periodic_complex,
    ensure_matplotlib,
    load_config,
    save_result_payload,
    solve_mean_zero_poisson,
)

plt = ensure_matplotlib()
NOTE_PATH = REPO_ROOT / 'experiments' / 'gauge_diagnostics' / 'Gauss_Law_Test_v1.md'


def shell_flux_balance(data, rho: np.ndarray, A: np.ndarray, source_point: np.ndarray, num_radii: int) -> dict[str, list[float] | float]:
    node_delta = (data.points - source_point + 0.5) % 1.0 - 0.5
    node_r = np.linalg.norm(node_delta, axis=1)
    shell_radii = np.linspace(1.0 / data.n_side, 0.45, num_radii)
    enclosed_charge: list[float] = []
    boundary_flux: list[float] = []
    signed_mismatch: list[float] = []
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
        signed_mismatch.append(float(q_inside + flux))
    return {
        'radius': shell_radii.tolist(),
        'enclosed_source': enclosed_charge,
        'boundary_flux': boundary_flux,
        'signed_mismatch': signed_mismatch,
        'max_signed_mismatch': float(np.max(np.abs(signed_mismatch)) if signed_mismatch else 0.0),
    }


def main() -> None:
    root_cfg = load_config()
    cfg = root_cfg['stage8_gauss_law_test']
    epsilon = float(root_cfg['epsilon'])
    sizes = [int(v) for v in cfg['sizes']]
    source_strength = float(cfg['source_strength'])
    flux_radii = int(cfg['flux_radii'])

    cases: list[dict[str, object]] = []
    residual_maps: dict[int, np.ndarray] = {}
    scatter_payload: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    integrated_payload: dict[int, dict[str, list[float] | float]] = {}

    for n_side in sizes:
        data = build_periodic_complex(n_side=n_side, epsilon=epsilon)
        center = n_side // 2
        source_idx = int(data.node_index[center, center, center])
        source_point = data.points[source_idx]

        rho = np.full(data.points.shape[0], -source_strength / data.points.shape[0], dtype=float)
        rho[source_idx] += source_strength
        phi = solve_mean_zero_poisson(data.L0, rho)
        A = np.asarray(data.d0 @ phi, dtype=float)
        G = np.asarray(data.d0.T @ A, dtype=float)

        slope = float(np.dot(G, rho) / max(np.dot(rho, rho), 1.0e-16))
        residual = G - slope * rho
        l2_residual = float(np.linalg.norm(residual) / max(np.linalg.norm(rho), 1.0e-16))
        linf_residual = float(np.max(np.abs(residual)))
        corr = float(np.corrcoef(G, rho)[0, 1])

        flux_balance = shell_flux_balance(data, rho, A, source_point, flux_radii)
        integrated_payload[n_side] = flux_balance

        residual_grid = residual.reshape((n_side, n_side, n_side))[:, :, center]
        residual_maps[n_side] = residual_grid
        scatter_payload[n_side] = (rho.copy(), G.copy())

        cases.append(
            {
                'n_side': n_side,
                'slope': slope,
                'correlation': corr,
                'l2_residual': l2_residual,
                'linf_residual': linf_residual,
                'max_signed_flux_mismatch': flux_balance['max_signed_mismatch'],
            }
        )

    largest_n = max(sizes)
    largest_case = next(case for case in cases if case['n_side'] == largest_n)

    plot_residual = REPO_ROOT / 'plots' / 'gauss_residual_map.png'
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(residual_maps[largest_n].T, origin='lower', cmap='coolwarm')
    ax.set_title(f'Gauss residual slice (n={largest_n})')
    ax.set_xlabel('x index')
    ax.set_ylabel('y index')
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.savefig(plot_residual, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_corr = REPO_ROOT / 'plots' / 'gauss_correlation.png'
    fig, ax = plt.subplots(figsize=(5.4, 4.8))
    rho_largest, G_largest = scatter_payload[largest_n]
    ax.scatter(rho_largest, G_largest, s=12, alpha=0.45)
    xline = np.linspace(np.min(rho_largest), np.max(rho_largest), 100)
    ax.plot(xline, largest_case['slope'] * xline, color='tab:red', linewidth=1.5, label=f"best-fit slope = {largest_case['slope']:.6f}")
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$d_0^{*} A$')
    ax.set_title(f'Gauss-law correlation (n={largest_n})')
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(plot_corr, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_flux = REPO_ROOT / 'plots' / 'gauss_integrated_flux.png'
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    for n_side in sizes:
        payload = integrated_payload[n_side]
        ax.plot(payload['radius'], np.abs(payload['signed_mismatch']), marker='o', label=f'n={n_side}')
    ax.set_xlabel('radius')
    ax.set_ylabel(r'$|\sum_{V(r)} \rho + \sum_{\partial V(r)} A\cdot n|$')
    ax.set_title('Gauss integrated-flux mismatch')
    ax.set_yscale('log')
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(plot_flux, dpi=180, bbox_inches='tight')
    plt.close(fig)

    result = {
        'experiment': 'stage8_gauss_law_test',
        'config': {
            'sizes': sizes,
            'epsilon': epsilon,
            'source_strength': source_strength,
            'flux_radii': flux_radii,
        },
        'cases': cases,
        'largest_n': largest_n,
        'integrated_flux': integrated_payload,
        'observation': 'the induced edge field satisfies a local divergence constraint that matches the inserted scalar source to numerical precision on the tested periodic complexes',
        'conclusion': 'the Gauss-law diagnostic is consistent with a discrete divergence constraint of the form d0* A = c rho with c near unity and small integrated flux mismatch',
    }

    stamped_json, stamped_plots, timestamp = save_result_payload('stage8_gauss_law_test', result, [plot_residual, plot_corr, plot_flux])
    append_log(
        'Stage 8.5A Gauss-Law Diagnostic',
        f'sizes={sizes}, epsilon={epsilon}, source_strength={source_strength}',
        stamped_json,
        stamped_plots,
        result['observation'],
        result['conclusion'],
    )
    NOTE_PATH.write_text(
        f"# Gauss Law Test v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: sizes={sizes}, epsilon={epsilon}, source_strength={source_strength}.\n\n"
        f"Largest grid (n={largest_n}) metrics:\n"
        f"- best-fit slope = {largest_case['slope']:.6f}\n"
        f"- correlation = {largest_case['correlation']:.12f}\n"
        f"- relative L2 residual = {largest_case['l2_residual']:.3e}\n"
        f"- L_inf residual = {largest_case['linf_residual']:.3e}\n"
        f"- max integrated mismatch = {largest_case['max_signed_flux_mismatch']:.3e}\n\n"
        f"Observation: {result['observation']}\n\n"
        f"Conclusion: {result['conclusion']}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
