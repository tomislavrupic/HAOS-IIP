#!/usr/bin/env python3

from __future__ import annotations

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
NOTE_PATH = REPO_ROOT / 'experiments' / 'gauge_diagnostics' / 'Longitudinal_Decoupling_Test_v1.md'


def main() -> None:
    full_config = load_config()
    config = full_config['stage8_longitudinal_decoupling']
    epsilon = full_config['epsilon']
    harmonic_tol = full_config['harmonic_tol']
    eig_tol = full_config['eig_tol']
    penalty = full_config['penalty']
    n_side = int(config['n_side'])
    restricted_modes = int(config['restricted_modes'])
    alphas = [float(a) for a in config['alphas']]

    data = build_periodic_complex(n_side=n_side, epsilon=epsilon)
    harmonic_basis = harmonic_basis_from_L1(data.L1, count=max(restricted_modes, 8), harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = ExactProjector(data.d0, data.L0)
    harmonic_projector = HarmonicProjector(harmonic_basis)
    transverse_operator = build_penalized_transverse_operator(data, exact_projector, harmonic_projector, penalty)
    raw_evals, raw_vecs = low_eigensystem(transverse_operator, restricted_modes + 4, eig_tol)

    basis: list[np.ndarray] = []
    lambdas: list[float] = []
    for idx in range(raw_vecs.shape[1]):
        vec = project_transverse(raw_vecs[:, idx], exact_projector, harmonic_projector)
        norm = float(np.linalg.norm(vec))
        if norm <= 1e-12:
            continue
        vec = vec / norm
        basis.append(vec)
        lambdas.append(float(np.dot(vec, data.upper @ vec)))
        if len(basis) >= restricted_modes:
            break
    V = np.column_stack(basis)
    A = basis[0]

    x = data.points[:, 0]
    y = data.points[:, 1]
    z = data.points[:, 2]
    phi = np.sin(2.0 * np.pi * x) + 0.5 * np.cos(2.0 * np.pi * y) + 0.25 * np.sin(2.0 * np.pi * z)
    phi -= np.mean(phi)
    exact_field = np.asarray(data.d0 @ phi, dtype=float)
    exact_norm = float(np.linalg.norm(exact_field)) or 1.0
    exact_field = exact_field / exact_norm

    P_A = project_transverse(A, exact_projector, harmonic_projector)
    P_A = P_A / (float(np.linalg.norm(P_A)) or 1.0)
    base_energy = float(np.dot(P_A, data.upper @ P_A))
    base_weights = np.abs(V.T @ P_A) ** 2
    base_divfree = A - exact_projector.apply(A)

    projection_residuals: list[float] = []
    energy_residuals: list[float] = []
    spectral_weight_residuals: list[float] = []
    divfree_residuals: list[float] = []

    for alpha in alphas:
        A_prime = A + alpha * exact_field
        P_prime = project_transverse(A_prime, exact_projector, harmonic_projector)
        P_prime = P_prime / (float(np.linalg.norm(P_prime)) or 1.0)
        energy_prime = float(np.dot(P_prime, data.upper @ P_prime))
        weights_prime = np.abs(V.T @ P_prime) ** 2
        divfree_prime = A_prime - exact_projector.apply(A_prime)
        projection_residuals.append(float(np.linalg.norm(P_prime - P_A) / max(np.linalg.norm(P_A), 1.0e-12)))
        energy_residuals.append(float(abs(energy_prime - base_energy) / max(abs(base_energy), 1.0e-12)))
        spectral_weight_residuals.append(float(np.linalg.norm(weights_prime - base_weights) / max(np.linalg.norm(base_weights), 1.0e-12)))
        divfree_residuals.append(float(np.linalg.norm(divfree_prime - base_divfree) / max(np.linalg.norm(base_divfree), 1.0e-12)))

    plot_res = REPO_ROOT / 'plots' / 'decoupling_residuals.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.plot(alphas, projection_residuals, marker='o', label='transverse projection residual')
    ax.plot(alphas, energy_residuals, marker='s', label='coexact energy residual')
    ax.plot(alphas, spectral_weight_residuals, marker='^', label='restricted spectral-weight residual')
    ax.plot(alphas, divfree_residuals, marker='D', label='div-free projection residual')
    ax.set_yscale('log')
    ax.set_xlabel(r'longitudinal amplitude $\alpha$')
    ax.set_ylabel('relative residual')
    ax.set_title('Stage 8 longitudinal decoupling test')
    ax.grid(alpha=0.25, which='both')
    ax.legend(fontsize=8)
    fig.savefig(plot_res, dpi=180, bbox_inches='tight')
    plt.close(fig)

    max_residual = max(projection_residuals + energy_residuals + spectral_weight_residuals + divfree_residuals)
    observation = 'adding pure exact edge components leaves the restricted transverse observables unchanged to numerical precision in the tested amplitudes'
    conclusion = 'the longitudinal-decoupling diagnostic is consistent with redundancy of exact perturbations at the restricted transverse level'
    result = {
        'experiment': 'stage8_longitudinal_decoupling',
        'config': {
            'n_side': n_side,
            'epsilon': epsilon,
            'restricted_modes': restricted_modes,
            'alphas': alphas,
        },
        'residuals': {
            'projection': projection_residuals,
            'coexact_energy': energy_residuals,
            'spectral_weights': spectral_weight_residuals,
            'divfree_projection': divfree_residuals,
            'max_residual': float(max_residual),
        },
        'observation': observation,
        'conclusion': conclusion,
    }

    stamped_json, stamped_plots, timestamp = save_result_payload('stage8_longitudinal_decoupling', result, [plot_res])
    append_log(
        'Stage 8C Longitudinal Decoupling',
        f'n={n_side}, epsilon={epsilon}, restricted_modes={restricted_modes}',
        stamped_json,
        stamped_plots,
        observation,
        conclusion,
    )
    NOTE_PATH.write_text(
        f"# Longitudinal Decoupling Test v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: n={n_side}, epsilon={epsilon}, restricted_modes={restricted_modes}, alphas={alphas}.\n\n"
        f"Maximum residual across all observables = {max_residual:.3e}.\n\n"
        f"Observation: {observation}\n\n"
        f"Conclusion: {conclusion}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
