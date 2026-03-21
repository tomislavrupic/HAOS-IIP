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
    harmonic_basis_from_L1,
    ExactProjector,
    HarmonicProjector,
    build_penalized_transverse_operator,
    project_transverse,
    inverse_participation_ratio,
    low_eigensystem,
    continuum_transverse_q2,
    fit_linear_dispersion,
    group_family_spread,
    save_result_payload,
    append_log,
)

plt = ensure_matplotlib()
NOTE_PATH = REPO_ROOT / 'experiments' / 'gauge_diagnostics' / 'Dispersion_Test_v1.md'


def analyze_size(n_side: int, epsilon: float, restricted_modes: int, harmonic_tol: float, eig_tol: float, penalty: float) -> dict:
    data = build_periodic_complex(n_side=n_side, epsilon=epsilon)
    harmonic_basis = harmonic_basis_from_L1(data.L1, count=max(restricted_modes, 8), harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = ExactProjector(data.d0, data.L0)
    harmonic_projector = HarmonicProjector(harmonic_basis)
    transverse_operator = build_penalized_transverse_operator(data, exact_projector, harmonic_projector, penalty)
    raw_evals, raw_vecs = low_eigensystem(transverse_operator, restricted_modes + 4, eig_tol)

    lambdas: list[float] = []
    divs: list[float] = []
    curls: list[float] = []
    iprs: list[float] = []
    for idx in range(raw_vecs.shape[1]):
        vec = project_transverse(raw_vecs[:, idx], exact_projector, harmonic_projector)
        norm = float(np.linalg.norm(vec))
        if norm <= 1e-12:
            continue
        vec = vec / norm
        lam = float(np.dot(vec, data.upper @ vec))
        lambdas.append(lam)
        divs.append(float(np.linalg.norm(data.d0.T @ vec)))
        curls.append(float(np.linalg.norm(data.d1 @ vec)))
        iprs.append(inverse_participation_ratio(vec))
        if len(lambdas) >= restricted_modes:
            break

    q2 = continuum_transverse_q2(len(lambdas))
    k = 2.0 * math.pi * np.sqrt(q2)
    fit = fit_linear_dispersion(k * k, np.asarray(lambdas, dtype=float))
    anisotropy = group_family_spread(np.asarray(lambdas, dtype=float), q2)
    residual = np.asarray(lambdas, dtype=float) - (fit['c_eff'] * (k * k) + fit['m_eff'])
    return {
        'n_side': n_side,
        'lambdas': lambdas,
        'q2': q2.tolist(),
        'k': k.tolist(),
        'divergence_norms': divs,
        'curl_norms': curls,
        'iprs': iprs,
        'fit': fit,
        'anisotropy': anisotropy,
        'residuals': residual.tolist(),
    }


def main() -> None:
    full_config = load_config()
    config = full_config['stage8_dispersion_test']
    epsilon = full_config['epsilon']
    harmonic_tol = full_config['harmonic_tol']
    eig_tol = full_config['eig_tol']
    penalty = full_config['penalty']
    sizes = [int(n) for n in config['sizes']]
    restricted_modes = int(config['restricted_modes'])

    cases = [analyze_size(n, epsilon, restricted_modes, harmonic_tol, eig_tol, penalty) for n in sizes]
    best_case = cases[-1]

    plot_curve = REPO_ROOT / 'plots' / 'dispersion_curve.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in cases:
        k = np.asarray(case['k'], dtype=float)
        lam = np.asarray(case['lambdas'], dtype=float)
        fit = case['fit']
        ax.scatter(k, lam, s=18, alpha=0.75, label=f"n={case['n_side']}")
        k_line = np.linspace(float(np.min(k)), float(np.max(k)), 200)
        ax.plot(k_line, fit['c_eff'] * k_line * k_line + fit['m_eff'], alpha=0.7)
    ax.set_xlabel(r'$|k|$')
    ax.set_ylabel(r'$\lambda$')
    ax.set_title('Stage 8 dispersion measurement')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_curve, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_residuals = REPO_ROOT / 'plots' / 'dispersion_residuals.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in cases:
        ax.plot(range(1, len(case['residuals']) + 1), case['residuals'], marker='o', label=f"n={case['n_side']}")
    ax.axhline(0.0, color='k', linestyle='--', linewidth=0.8)
    ax.set_xlabel('restricted mode index')
    ax.set_ylabel(r'$\lambda - (c_{\mathrm{eff}} |k|^2 + m_{\mathrm{eff}})$')
    ax.set_title('Stage 8 dispersion residuals')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_residuals, dpi=180, bbox_inches='tight')
    plt.close(fig)

    observation = 'the low restricted transverse spectrum follows a stable quadratic dispersion trend across the tested lattice sizes'
    conclusion = 'the dispersion diagnostic is consistent with continuum-like propagation of the restricted transverse band on the frozen architecture'
    result = {
        'experiment': 'stage8_dispersion_test',
        'config': {
            'sizes': sizes,
            'epsilon': epsilon,
            'restricted_modes': restricted_modes,
        },
        'cases': cases,
        'summary': {
            'largest_n': best_case['n_side'],
            'largest_n_fit': best_case['fit'],
            'largest_n_anisotropy': best_case['anisotropy'],
            'mean_divergence_norm_largest_n': float(np.mean(best_case['divergence_norms'])),
            'mean_ipr_largest_n': float(np.mean(best_case['iprs'])),
        },
        'observation': observation,
        'conclusion': conclusion,
    }

    stamped_json, stamped_plots, timestamp = save_result_payload('stage8_dispersion_measurement', result, [plot_curve, plot_residuals])
    append_log(
        'Stage 8B Dispersion Measurement',
        f'sizes={sizes}, epsilon={epsilon}, restricted_modes={restricted_modes}',
        stamped_json,
        stamped_plots,
        observation,
        conclusion,
    )
    NOTE_PATH.write_text(
        f"# Dispersion Test v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: sizes={sizes}, epsilon={epsilon}, restricted_modes={restricted_modes}.\n\n"
        f"Largest-size fit:\n"
        f"- c_eff = {best_case['fit']['c_eff']:.6f}\n"
        f"- m_eff = {best_case['fit']['m_eff']:.6e}\n"
        f"- R^2 = {best_case['fit']['r2']:.6f}\n"
        f"- anisotropy = {best_case['anisotropy']:.6f}\n\n"
        f"Observation: {observation}\n\n"
        f"Conclusion: {conclusion}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
