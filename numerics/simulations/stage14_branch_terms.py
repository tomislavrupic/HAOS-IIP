#!/usr/bin/env python3

from __future__ import annotations

import numpy as np

from stage12_coarse_field_structure import edge_field_to_grid, gaussian_smooth_periodic


def smooth_edge_profile(midpoints: np.ndarray, edge_values: np.ndarray, n_side: int, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    grid = edge_field_to_grid(midpoints, np.asarray(edge_values, dtype=float), n_side)
    env = np.maximum(gaussian_smooth_periodic(grid, sigma), 0.0)
    max_env = float(np.max(env))
    if max_env > 0.0:
        env = env / max_env
    coords = np.floor(midpoints * n_side).astype(int) % n_side
    profile = env[coords[:, 0], coords[:, 1], coords[:, 2]]
    return profile.astype(float), env


def phase_lock_retention_profile(midpoints: np.ndarray, component_qs: list[np.ndarray], n_side: int, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    amps = np.stack([np.abs(q) for q in component_qs], axis=0)
    denom = np.maximum(np.sum(amps, axis=0), 1.0e-12)
    alignment = np.abs(np.sum(component_qs, axis=0)) / denom
    overlap = np.clip((np.sum(amps, axis=0) - np.max(amps, axis=0)) / denom, 0.0, 1.0)
    raw = alignment * overlap
    return smooth_edge_profile(midpoints, raw, n_side, sigma)


def envelope_coupling_profiles(midpoints: np.ndarray, component_qs: list[np.ndarray], n_side: int, sigma: float) -> tuple[list[np.ndarray], list[np.ndarray]]:
    amps = [np.abs(q) for q in component_qs]
    total_amp = np.sum(amps, axis=0)
    profiles: list[np.ndarray] = []
    envs: list[np.ndarray] = []
    for amp in amps:
        other_amp = np.clip(total_amp - amp, 0.0, None)
        profile, env = smooth_edge_profile(midpoints, other_amp, n_side, sigma)
        profiles.append(profile)
        envs.append(env)
    return profiles, envs
