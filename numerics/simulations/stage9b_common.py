#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from stage9_common import REPO_ROOT, WAVE_NOTES, append_log, plt, save_result_payload
from DK_stage6_common import build_dk2d_complex, low_eigs as low_eigs_2d, pairing_error as pairing_error_2d
from DK_stage7_common import build_dk3d_complex, low_eigs as low_eigs_3d, pairing_error as pairing_error_3d


DEFAULT_STAGE9B: dict[str, Any] = {
    'epsilon': 0.2,
    'n_2d': 16,
    'n_3d': 6,
    'sigma_2d': 0.08,
    'sigma_3d': 0.12,
    'packet_center_2d': [0.25, 0.5],
    'packet_center_3d': [0.25, 0.5, 0.5],
    'kick_cycles': 1.0,
    't_final_2d': 0.6,
    't_final_3d': 0.35,
    'dt_scale': 0.35,
    'holonomy_phases': [0.0, 0.2617993877991494, 0.5235987755982988, 0.7853981633974483],
    'low_modes_2d': 40,
    'low_modes_3d': 32,
}


def load_stage9b_config(config_path: Path | None = None) -> dict[str, Any]:
    path = config_path or (REPO_ROOT / 'config.json')
    merged = json.loads(json.dumps(DEFAULT_STAGE9B))
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged['epsilon'] = float(on_disk.get('epsilon', merged['epsilon']))
        if isinstance(on_disk.get('stage9b'), dict):
            merged.update(on_disk['stage9b'])
    return merged


def circular_center(points: np.ndarray, weights: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    weights = np.asarray(weights, dtype=float)
    total = float(np.sum(weights))
    if total <= 0.0:
        return np.zeros(points.shape[1], dtype=float)
    center = np.zeros(points.shape[1], dtype=float)
    for axis in range(points.shape[1]):
        angles = 2.0 * math.pi * points[:, axis]
        s = float(np.sum(weights * np.sin(angles)))
        c = float(np.sum(weights * np.cos(angles)))
        center[axis] = ((math.atan2(s, c) / (2.0 * math.pi)) + 1.0) % 1.0
    return center


def periodic_displacement(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    anchor = np.asarray(anchor, dtype=float)
    delta = points - anchor
    return (delta + 0.5) % 1.0 - 0.5


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


def complex_power_lambda_max(operator: sp.csr_matrix, iterations: int = 40) -> float:
    dim = operator.shape[0]
    rng = np.random.default_rng(17)
    vec = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    vec /= np.linalg.norm(vec)
    value = 0.0
    for _ in range(iterations):
        w = operator @ vec
        norm = float(np.linalg.norm(w))
        if norm <= 0.0:
            return 0.0
        vec = w / norm
        value = float(np.vdot(vec, operator @ vec).real)
    return max(value, 1.0e-12)


def first_order_dt(operator: sp.csr_matrix, scale: float) -> float:
    lam_max = complex_power_lambda_max(operator @ operator)
    return float(scale / math.sqrt(lam_max))


def second_order_dt(operator: sp.csr_matrix, scale: float) -> float:
    lam_max = complex_power_lambda_max(operator)
    return float(scale * 2.0 / math.sqrt(lam_max))


def pack_positions(*blocks: np.ndarray) -> np.ndarray:
    return np.vstack([np.asarray(block, dtype=float) for block in blocks])


def graded_weights(vec: np.ndarray, block_sizes: tuple[int, ...]) -> list[float]:
    weights: list[float] = []
    start = 0
    total = float(np.sum(np.abs(vec) ** 2)) or 1.0
    for size in block_sizes:
        part = vec[start:start + size]
        weights.append(float(np.sum(np.abs(part) ** 2) / total))
        start += size
    return weights


def packet_observables(vec: np.ndarray, positions: np.ndarray, block_sizes: tuple[int, ...]) -> dict[str, Any]:
    weights = np.abs(vec) ** 2
    center = circular_center(positions, weights)
    return {
        'center': center.tolist(),
        'width': packet_width(positions, weights, center),
        'anisotropy': anisotropy_ratio(positions, weights, center),
        'grade_weights': graded_weights(vec, block_sizes),
        'norm': float(np.linalg.norm(vec)),
    }


def crank_nicolson_evolution(
    operator: sp.csr_matrix,
    psi0: np.ndarray,
    dt: float,
    steps: int,
    positions: np.ndarray,
    block_sizes: tuple[int, ...],
) -> dict[str, Any]:
    psi = np.asarray(psi0, dtype=complex).copy()
    ident = sp.identity(operator.shape[0], dtype=complex, format='csr')
    lhs = (ident - 0.5j * dt * operator).tocsc()
    rhs_op = (ident + 0.5j * dt * operator).tocsr()
    solve = spla.factorized(lhs)
    energy0 = float(np.vdot(psi, operator @ psi).real)

    times: list[float] = []
    centers: list[list[float]] = []
    widths: list[float] = []
    norms: list[float] = []
    energies: list[float] = []
    anisotropies: list[float] = []
    grade_hist: list[list[float]] = []

    def record(t: float) -> None:
        obs = packet_observables(psi, positions, block_sizes)
        times.append(float(t))
        centers.append(obs['center'])
        widths.append(float(obs['width']))
        norms.append(float(obs['norm']))
        energies.append(float(np.vdot(psi, operator @ psi).real))
        anisotropies.append(float(obs['anisotropy']))
        grade_hist.append(obs['grade_weights'])

    record(0.0)
    for step in range(steps):
        psi = solve(rhs_op @ psi)
        record((step + 1) * dt)

    return {
        'times': times,
        'centers': centers,
        'widths': widths,
        'norms': norms,
        'energies': energies,
        'anisotropies': anisotropies,
        'grade_weights': grade_hist,
        'energy0': energy0,
    }


def leapfrog_second_order_complex(
    operator: sp.csr_matrix,
    q0: np.ndarray,
    v0: np.ndarray,
    dt: float,
    steps: int,
    positions: np.ndarray,
    block_sizes: tuple[int, ...],
) -> dict[str, Any]:
    q = np.asarray(q0, dtype=complex).copy()
    v = np.asarray(v0, dtype=complex).copy()

    times: list[float] = []
    centers: list[list[float]] = []
    widths: list[float] = []
    norms: list[float] = []
    energies: list[float] = []
    anisotropies: list[float] = []
    grade_hist: list[list[float]] = []

    def record(t: float) -> None:
        obs = packet_observables(q, positions, block_sizes)
        times.append(float(t))
        centers.append(obs['center'])
        widths.append(float(obs['width']))
        norms.append(float(obs['norm']))
        kinetic = float(np.vdot(v, v).real)
        potential = float(np.vdot(q, operator @ q).real)
        energies.append(0.5 * (kinetic + potential))
        anisotropies.append(float(obs['anisotropy']))
        grade_hist.append(obs['grade_weights'])

    record(0.0)
    for step in range(steps):
        acc = -(operator @ q)
        v_half = v + 0.5 * dt * acc
        q = q + dt * v_half
        acc_new = -(operator @ q)
        v = v_half + 0.5 * dt * acc_new
        record((step + 1) * dt)

    return {
        'times': times,
        'centers': centers,
        'widths': widths,
        'norms': norms,
        'energies': energies,
        'anisotropies': anisotropies,
        'grade_weights': grade_hist,
    }


def gaussian_profile(points: np.ndarray, center: np.ndarray, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    delta = periodic_displacement(points, center)
    r2 = np.sum(delta * delta, axis=1)
    profile = np.exp(-r2 / (2.0 * sigma * sigma))
    return profile, delta


def build_initial_packet(positions: np.ndarray, block_sizes: tuple[int, ...], center: np.ndarray, sigma: float, kick_cycles: float, axis: int, active_grades: list[int]) -> np.ndarray:
    parts: list[np.ndarray] = []
    start = 0
    for grade_idx, size in enumerate(block_sizes):
        coords = positions[start:start + size]
        profile, delta = gaussian_profile(coords, center, sigma)
        phase = np.exp(1j * 2.0 * math.pi * kick_cycles * delta[:, axis])
        part = profile * phase if grade_idx in active_grades else np.zeros(size, dtype=complex)
        parts.append(np.asarray(part, dtype=complex))
        start += size
    psi = np.concatenate(parts)
    return psi / max(np.linalg.norm(psi), 1.0e-12)


def spectral_pairing_summary_2d(operator: sp.csr_matrix, modes: int, tol: float) -> dict[str, Any]:
    evals, _ = low_eigs_2d(operator, k=modes, sigma=0.0, tol=tol)
    return {
        'pairing_error': pairing_error_2d(evals),
        'lowest_abs_eigenvalue': float(np.min(np.abs(evals))),
        'eigenvalues': evals.tolist(),
    }


def spectral_pairing_summary_3d(operator: sp.csr_matrix, modes: int, tol: float) -> dict[str, Any]:
    evals, _ = low_eigs_3d(operator, k=modes, sigma=0.0, tol=tol)
    return {
        'pairing_error': pairing_error_3d(evals),
        'lowest_abs_eigenvalue': float(np.min(np.abs(evals))),
        'eigenvalues': evals.tolist(),
    }


__all__ = [
    'REPO_ROOT',
    'WAVE_NOTES',
    'append_log',
    'build_dk2d_complex',
    'build_dk3d_complex',
    'build_initial_packet',
    'crank_nicolson_evolution',
    'first_order_dt',
    'gaussian_profile',
    'leapfrog_second_order_complex',
    'load_stage9b_config',
    'pack_positions',
    'packet_observables',
    'plt',
    'save_result_payload',
    'second_order_dt',
    'spectral_pairing_summary_2d',
    'spectral_pairing_summary_3d',
]
