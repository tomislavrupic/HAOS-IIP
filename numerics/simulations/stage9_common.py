#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
import scipy.sparse as sp

from stage8_common import (
    REPO_ROOT,
    RESULTS,
    PLOTS,
    append_log,
    build_periodic_complex,
    ensure_matplotlib,
    harmonic_basis_from_L1,
    ExactProjector,
    HarmonicProjector,
    load_config as load_stage8_config,
    project_transverse,
    save_result_payload,
)

plt = ensure_matplotlib()

WAVE_NOTES = REPO_ROOT / 'experiments' / 'wave_dynamics'
WAVE_NOTES.mkdir(exist_ok=True)

DEFAULT_STAGE9: dict[str, Any] = {
    'sizes': [16, 20, 24],
    'epsilon': 0.2,
    'sigma': 0.08,
    'packet_center': [0.25, 0.5, 0.5],
    'kick_axis': 0,
    'kick_strength_scalar': 1.0,
    'kick_strength_transverse': 1.0,
    't_final': 0.45,
    'dt_safety': 0.45,
    'harmonic_tol': 1.0e-8,
    'penalty': 10.0,
}


def load_stage9_config(config_path: Path | None = None) -> dict[str, Any]:
    base = load_stage8_config(config_path)
    merged = json.loads(json.dumps(DEFAULT_STAGE9))
    merged['epsilon'] = float(base.get('epsilon', merged['epsilon']))
    merged['harmonic_tol'] = float(base.get('harmonic_tol', merged['harmonic_tol']))
    merged['penalty'] = float(base.get('penalty', merged['penalty']))
    stage9_on_disk = base.get('stage9', {})
    if isinstance(stage9_on_disk, dict):
        merged.update(stage9_on_disk)
    return merged


def periodic_displacement(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    anchor = np.asarray(anchor, dtype=float)
    delta = points - anchor
    return (delta + 0.5) % 1.0 - 0.5


def gaussian_packet(points: np.ndarray, center: np.ndarray, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    delta = periodic_displacement(points, center)
    r2 = np.sum(delta * delta, axis=1)
    profile = np.exp(-r2 / (2.0 * sigma * sigma))
    return profile, delta


def circular_center(points: np.ndarray, weights: np.ndarray) -> np.ndarray:
    weights = np.asarray(weights, dtype=float)
    total = float(np.sum(weights))
    if total <= 0.0:
        return np.zeros(3, dtype=float)
    center = np.zeros(3, dtype=float)
    for axis in range(3):
        angles = 2.0 * math.pi * np.asarray(points[:, axis], dtype=float)
        s = float(np.sum(weights * np.sin(angles)))
        c = float(np.sum(weights * np.cos(angles)))
        center[axis] = ((math.atan2(s, c) / (2.0 * math.pi)) + 1.0) % 1.0
    return center


def packet_width(points: np.ndarray, weights: np.ndarray, center: np.ndarray) -> float:
    total = float(np.sum(weights))
    if total <= 0.0:
        return 0.0
    delta = periodic_displacement(points, center)
    r2 = np.sum(delta * delta, axis=1)
    return float(np.sqrt(np.sum(weights * r2) / total))


def anisotropy_ratio(points: np.ndarray, weights: np.ndarray, center: np.ndarray) -> float:
    total = float(np.sum(weights))
    if total <= 0.0:
        return 0.0
    delta = periodic_displacement(points, center)
    cov = (delta.T * weights) @ delta / total
    evals = np.linalg.eigvalsh(cov)
    smallest = float(max(np.min(evals), 1.0e-12))
    largest = float(np.max(evals))
    return largest / smallest


def energy(operator: sp.csr_matrix, q: np.ndarray, v: np.ndarray) -> float:
    kinetic = float(np.vdot(v, v).real)
    potential = float(np.vdot(q, operator @ q).real)
    return 0.5 * (kinetic + potential)


def estimate_lambda_max(operator: sp.csr_matrix, iterations: int = 40) -> float:
    dim = operator.shape[0]
    rng = np.random.default_rng(9)
    vec = rng.normal(size=dim)
    vec /= np.linalg.norm(vec)
    value = 0.0
    for _ in range(iterations):
        w = np.asarray(operator @ vec, dtype=float)
        norm = float(np.linalg.norm(w))
        if norm <= 0.0:
            return 0.0
        vec = w / norm
        value = float(np.dot(vec, np.asarray(operator @ vec, dtype=float)))
    return max(value, 1.0e-12)


def suggested_dt(operator: sp.csr_matrix, safety: float) -> float:
    lam_max = estimate_lambda_max(operator)
    return float(safety * 2.0 / math.sqrt(lam_max))


def simulate_second_order(
    operator: sp.csr_matrix,
    q0: np.ndarray,
    v0: np.ndarray,
    positions: np.ndarray,
    dt: float,
    steps: int,
    project: Callable[[np.ndarray], np.ndarray] | None = None,
    divergence_op: sp.csr_matrix | None = None,
    reference_q: np.ndarray | None = None,
) -> dict[str, Any]:
    q = np.asarray(q0, dtype=float).copy()
    v = np.asarray(v0, dtype=float).copy()
    if project is not None:
        q = np.asarray(project(q), dtype=float)
        v = np.asarray(project(v), dtype=float)

    times: list[float] = []
    centers: list[list[float]] = []
    widths: list[float] = []
    norms: list[float] = []
    energies: list[float] = []
    anisotropies: list[float] = []
    divergence_norms: list[float] = []
    projection_residuals: list[float] = []
    reference_differences: list[float] = []

    def record(t: float) -> None:
        weights = np.abs(q) ** 2
        center = circular_center(positions, weights)
        times.append(float(t))
        centers.append(center.tolist())
        widths.append(packet_width(positions, weights, center))
        norms.append(float(np.linalg.norm(q)))
        energies.append(energy(operator, q, v))
        anisotropies.append(anisotropy_ratio(positions, weights, center))
        if divergence_op is not None:
            divergence_norms.append(float(np.linalg.norm(divergence_op @ q)))
        if project is not None:
            projection_residuals.append(float(np.linalg.norm(np.asarray(project(q), dtype=float) - q)))
        if reference_q is not None:
            reference_differences.append(float(np.linalg.norm(q - reference_q)))

    record(0.0)
    for step in range(steps):
        acc = -np.asarray(operator @ q, dtype=float)
        if project is not None:
            acc = np.asarray(project(acc), dtype=float)
        v_half = v + 0.5 * dt * acc
        q = q + dt * v_half
        if project is not None:
            q = np.asarray(project(q), dtype=float)
        acc_new = -np.asarray(operator @ q, dtype=float)
        if project is not None:
            acc_new = np.asarray(project(acc_new), dtype=float)
        v = v_half + 0.5 * dt * acc_new
        if project is not None:
            v = np.asarray(project(v), dtype=float)
        record((step + 1) * dt)

    return {
        'times': times,
        'centers': centers,
        'widths': widths,
        'norms': norms,
        'energies': energies,
        'anisotropies': anisotropies,
        'divergence_norms': divergence_norms,
        'projection_residuals': projection_residuals,
        'reference_differences': reference_differences,
    }


def build_scalar_initial_state(points: np.ndarray, center: np.ndarray, sigma: float, kick_axis: int, kick_strength: float) -> tuple[np.ndarray, np.ndarray]:
    profile, delta = gaussian_packet(points, center, sigma)
    q0 = profile / max(np.linalg.norm(profile), 1.0e-12)
    v0 = kick_strength * (delta[:, kick_axis] / (sigma * sigma)) * profile
    v0 = v0 / max(np.linalg.norm(v0), 1.0e-12)
    return q0, v0


def build_transverse_setup(n_side: int, epsilon: float, harmonic_tol: float) -> tuple[Any, Callable[[np.ndarray], np.ndarray]]:
    data = build_periodic_complex(n_side=n_side, epsilon=epsilon)
    harmonic_basis = harmonic_basis_from_L1(data.L1, count=8, harmonic_tol=harmonic_tol, tol=1.0e-8)
    exact_projector = ExactProjector(data.d0, data.L0)
    harmonic_projector = HarmonicProjector(harmonic_basis)

    def projector(vec: np.ndarray) -> np.ndarray:
        return project_transverse(np.asarray(vec, dtype=float), exact_projector, harmonic_projector)

    return data, projector


def build_transverse_initial_state(
    midpoints: np.ndarray,
    edge_axes: list[str],
    center: np.ndarray,
    sigma: float,
    kick_axis: int,
    kick_strength: float,
    projector: Callable[[np.ndarray], np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    profile, delta = gaussian_packet(midpoints, center, sigma)
    axis_mask = np.asarray([1.0 if axis == 'y' else 0.0 for axis in edge_axes], dtype=float)
    raw = profile * axis_mask
    q0 = projector(raw)
    q0 = q0 / max(np.linalg.norm(q0), 1.0e-12)

    raw_v = kick_strength * (delta[:, kick_axis] / (sigma * sigma)) * raw
    v0 = projector(raw_v)
    v0 = v0 / max(np.linalg.norm(v0), 1.0e-12)
    return q0, v0


__all__ = [
    'REPO_ROOT',
    'RESULTS',
    'PLOTS',
    'WAVE_NOTES',
    'append_log',
    'build_periodic_complex',
    'build_scalar_initial_state',
    'build_transverse_initial_state',
    'build_transverse_setup',
    'ensure_matplotlib',
    'load_stage9_config',
    'plt',
    'save_result_payload',
    'simulate_second_order',
    'suggested_dt',
]
