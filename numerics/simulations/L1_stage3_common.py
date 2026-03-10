#!/usr/bin/env python3

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from L1_stage2_common import (
    PLOTS,
    REPO_ROOT,
    RESULTS,
    EXPERIMENT_LOG,
    save_result_payload,
    append_log,
    ensure_matplotlib,
)
from L1_transverse_band_test import analyze_case as analyze_band_case

plt = ensure_matplotlib()

VARIANT_LABELS = {
    'baseline': 'baseline periodic torus',
    'puncture': 'single cubic puncture',
    'line_defect': 'line defect',
    'flux_tube': 'flux-tube defect',
}

VARIANT_MARKERS = {
    'baseline': 'o',
    'puncture': 's',
    'line_defect': '^',
    'flux_tube': 'D',
}


def continuum_transverse_q2(k: int, max_m: int = 8) -> np.ndarray:
    values: list[float] = []
    for mx in range(-max_m, max_m + 1):
        for my in range(-max_m, max_m + 1):
            for mz in range(-max_m, max_m + 1):
                if mx == 0 and my == 0 and mz == 0:
                    continue
                q2 = float(mx * mx + my * my + mz * mz)
                values.extend([q2, q2])
    values.sort()
    if len(values) < k:
        raise ValueError('Increase max_m to generate enough continuum modes')
    return np.asarray(values[:k], dtype=float)


def linear_gap_fit(q2: np.ndarray, lambdas: np.ndarray) -> dict[str, float]:
    x = np.asarray(q2, dtype=float)
    y = np.asarray(lambdas, dtype=float)
    A = np.column_stack([x, np.ones_like(x)])
    coeffs, *_ = np.linalg.lstsq(A, y, rcond=None)
    c_eff, m_eff = coeffs
    pred = A @ coeffs
    resid = y - pred
    ss_res = float(np.sum(resid * resid))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        'c_eff': float(c_eff),
        'm_eff': float(m_eff),
        'r2': float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0,
    }


def mode_spacing(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return np.zeros(0, dtype=float)
    return np.diff(values)


def relative_spread(curves: list[np.ndarray], count: int) -> float:
    stack = np.asarray([curve[:count] for curve in curves], dtype=float)
    mean_curve = np.mean(stack, axis=0)
    return float(np.mean(np.std(stack, axis=0) / np.maximum(np.abs(mean_curve), 1e-12)))


def periodic_displacement(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    anchor = np.asarray(anchor, dtype=float)
    delta = points - anchor
    return (delta + 0.5) % 1.0 - 0.5


def analyze_branch_cases(
    sizes: list[int],
    variants: list[str],
    epsilon: float,
    restricted_modes: int,
    harmonic_tol: float,
    eig_tol: float,
    penalty: float,
    flux_tube_phase: float,
) -> dict[str, dict[str, Any]]:
    cases: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        for variant in variants:
            cases[f'{variant}_n{n_side}'] = analyze_band_case(
                n_side=n_side,
                epsilon=epsilon,
                variant=variant,
                restricted_modes=restricted_modes,
                harmonic_tol=harmonic_tol,
                eig_tol=eig_tol,
                penalty=penalty,
                flux_tube_phase=flux_tube_phase,
            )
    return cases


def make_phase_plot(case: dict[str, Any], path: Path, title: str) -> None:
    records = case['restricted_transverse_modes']
    x = [record['divergence_norm'] for record in records]
    y = [record['curl_norm'] for record in records]
    c = [record['eigenvalue'] for record in records]
    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(x, y, c=c, cmap='viridis', s=28)
    ax.set_xlabel(r'$||d_0^* a_k||$')
    ax.set_ylabel(r'$||d_1 a_k||$')
    ax.set_title(title)
    ax.grid(alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r'$\lambda_k$')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)
