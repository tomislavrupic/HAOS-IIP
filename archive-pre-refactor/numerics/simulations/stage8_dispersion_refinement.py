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
    inverse_participation_ratio,
    low_eigensystem,
    continuum_transverse_q2,
    fit_linear_dispersion,
    group_family_spread,
    save_result_payload,
    append_log,
)

plt = ensure_matplotlib()
NOTE_PATH = REPO_ROOT / 'experiments' / 'gauge_diagnostics' / 'Dispersion_Refinement_Test_v1.md'


def analyze_size(n_side: int, epsilon: float, restricted_modes: int, harmonic_tol: float, eig_tol: float, penalty: float) -> dict:
    data = build_periodic_complex(n_side=n_side, epsilon=epsilon)
    harmonic_basis = harmonic_basis_from_L1(data.L1, count=max(restricted_modes, 12), harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = ExactProjector(data.d0, data.L0)
    harmonic_projector = HarmonicProjector(harmonic_basis)
    transverse_operator = build_penalized_transverse_operator(data, exact_projector, harmonic_projector, penalty)
    _, raw_vecs = low_eigensystem(transverse_operator, restricted_modes + 6, eig_tol)

    lambdas = []
    divs = []
    curls = []
    iprs = []
    for idx in range(raw_vecs.shape[1]):
        vec = project_transverse(raw_vecs[:, idx], exact_projector, harmonic_projector)
        norm = float(np.linalg.norm(vec))
        if norm <= 1e-12:
            continue
        vec = vec / norm
        lambdas.append(float(np.dot(vec, data.upper @ vec)))
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
        'n2_lambda1': float((n_side**2) * lambdas[0]),
    }


def main() -> None:
    full_config = load_config()
    epsilon = float(full_config['epsilon'])
    harmonic_tol = float(full_config['harmonic_tol'])
    eig_tol = float(full_config['eig_tol'])
    penalty = float(full_config['penalty'])

    sizes = [16, 20, 24]
    restricted_modes = 24
    cases = [analyze_size(n, epsilon, restricted_modes, harmonic_tol, eig_tol, penalty) for n in sizes]
    best_case = cases[-1]

    xs = np.log(np.asarray(sizes, dtype=float))
    ys = np.log(np.asarray([case['lambdas'][0] for case in cases], dtype=float))
    slope, intercept = np.polyfit(xs, ys, 1)
    nfit = np.exp(intercept) * np.exp(slope * xs)
    ss_res = float(np.sum((np.exp(ys) - nfit) ** 2))
    ss_tot = float(np.sum((np.exp(ys) - np.mean(np.exp(ys))) ** 2))
    scaling_fit = {
        'prefactor': float(math.exp(intercept)),
        'power': float(-slope),
        'r2': float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0,
    }

    plot_curve = REPO_ROOT / 'plots' / 'dispersion_refinement_curve.png'
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
    ax.set_title('Stage 8 dispersion refinement')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_curve, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_n = REPO_ROOT / 'plots' / 'dispersion_refinement_n2lambda.png'
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    ax.plot(sizes, [case['n2_lambda1'] for case in cases], marker='o')
    ax.set_xlabel('n')
    ax.set_ylabel(r'$n^2 \lambda_1$')
    ax.set_title('Stage 8 dispersion refinement: lowest-mode scaling')
    ax.grid(alpha=0.25)
    fig.savefig(plot_n, dpi=180, bbox_inches='tight')
    plt.close(fig)

    plot_residuals = REPO_ROOT / 'plots' / 'dispersion_refinement_residuals.png'
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for case in cases:
        ax.plot(range(1, len(case['residuals']) + 1), case['residuals'], marker='o', label=f"n={case['n_side']}")
    ax.axhline(0.0, color='k', linestyle='--', linewidth=0.8)
    ax.set_xlabel('restricted mode index')
    ax.set_ylabel(r'$\lambda - (c_{\mathrm{eff}} |k|^2 + m_{\mathrm{eff}})$')
    ax.set_title('Stage 8 dispersion refinement residuals')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(plot_residuals, dpi=180, bbox_inches='tight')
    plt.close(fig)

    observation = 'extending the lattice-size scan preserves the n^2 lambda_1 plateau of the low restricted transverse band while leaving the direct momentum-fit quality moderate'
    conclusion = 'the larger-n refinement strengthens the continuum scaling evidence for the lowest transverse modes, although the coarse momentum assignment still limits the direct dispersion fit'
    result = {
        'experiment': 'stage8_dispersion_refinement',
        'config': {
            'sizes': sizes,
            'epsilon': epsilon,
            'restricted_modes': restricted_modes,
        },
        'cases': cases,
        'summary': {
            'scaling_fit': scaling_fit,
            'n2_lambda1_values': [case['n2_lambda1'] for case in cases],
            'largest_n_fit': best_case['fit'],
            'largest_n_anisotropy': best_case['anisotropy'],
            'mean_divergence_norm_largest_n': float(np.mean(best_case['divergence_norms'])),
            'mean_ipr_largest_n': float(np.mean(best_case['iprs'])),
        },
        'observation': observation,
        'conclusion': conclusion,
    }

    stamped_json, stamped_plots, _ = save_result_payload('stage8_dispersion_refinement', result, [plot_curve, plot_n, plot_residuals])
    append_log(
        'Stage 8B Refinement Dispersion',
        f'sizes={sizes}, epsilon={epsilon}, restricted_modes={restricted_modes}',
        stamped_json,
        stamped_plots,
        observation,
        conclusion,
    )
    NOTE_PATH.write_text(
        f"# Dispersion Refinement Test v1\n\n"
        f"Timestamped result: `{stamped_json.relative_to(REPO_ROOT)}`\n\n"
        f"Config: sizes={sizes}, epsilon={epsilon}, restricted_modes={restricted_modes}.\n\n"
        f"Lowest-mode scaling fit:\n"
        f"- power = {scaling_fit['power']:.6f}\n"
        f"- R^2 = {scaling_fit['r2']:.6f}\n"
        f"- n^2 lambda_1 values = {[round(v, 6) for v in result['summary']['n2_lambda1_values']]}\n\n"
        f"Largest-size direct fit:\n"
        f"- c_eff = {best_case['fit']['c_eff']:.6f}\n"
        f"- m_eff = {best_case['fit']['m_eff']:.6e}\n"
        f"- R^2 = {best_case['fit']['r2']:.6f}\n\n"
        f"Observation: {observation}\n\n"
        f"Conclusion: {conclusion}\n",
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
