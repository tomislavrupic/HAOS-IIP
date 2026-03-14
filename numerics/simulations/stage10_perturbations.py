#!/usr/bin/env python3

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy.sparse as sp

from stage10_common import Stage10Complex


@dataclass(frozen=True)
class PerturbedOperators:
    axis: str
    strength: float
    seed: int
    d0: sp.csr_matrix
    d1: sp.csr_matrix
    L0: sp.csr_matrix
    L1: sp.csr_matrix
    edge_scales: np.ndarray
    face_scales: np.ndarray
    metadata: dict[str, Any]


def _face_scales_from_edges(data: Stage10Complex, edge_scales: np.ndarray) -> np.ndarray:
    d1 = data.d1.tocsr()
    scales = np.ones(d1.shape[0], dtype=float)
    for face_idx in range(d1.shape[0]):
        start = d1.indptr[face_idx]
        end = d1.indptr[face_idx + 1]
        cols = d1.indices[start:end]
        if cols.size:
            scales[face_idx] = float(np.mean(edge_scales[cols]))
    return scales


def _operators_from_edge_scales(
    data: Stage10Complex,
    edge_scales: np.ndarray,
    axis: str,
    strength: float,
    seed: int,
    metadata: dict[str, Any],
) -> PerturbedOperators:
    edge_scales = np.asarray(edge_scales, dtype=float)
    edge_scales = np.clip(edge_scales, 1.0e-3, None)
    face_scales = _face_scales_from_edges(data, edge_scales)

    d0 = (sp.diags(edge_scales) @ data.d0).tocsr()
    d1 = (sp.diags(face_scales) @ data.d1).tocsr()
    L0 = (d0.T @ d0).tocsr()
    L1 = (d0 @ d0.T + d1.T @ d1).tocsr()

    meta = {
        **metadata,
        'edge_scale_min': float(np.min(edge_scales)),
        'edge_scale_max': float(np.max(edge_scales)),
        'face_scale_min': float(np.min(face_scales)),
        'face_scale_max': float(np.max(face_scales)),
    }
    return PerturbedOperators(
        axis=axis,
        strength=float(strength),
        seed=int(seed),
        d0=d0,
        d1=d1,
        L0=L0,
        L1=L1,
        edge_scales=edge_scales,
        face_scales=face_scales,
        metadata=meta,
    )


def perturb_edge_weight_noise(data: Stage10Complex, strength: float, seed: int) -> PerturbedOperators:
    rng = np.random.default_rng(seed)
    weight_factors = 1.0 + rng.uniform(-strength, strength, size=len(data.edge_axes))
    edge_scales = np.sqrt(np.clip(weight_factors, 1.0e-3, None))
    return _operators_from_edge_scales(
        data=data,
        edge_scales=edge_scales,
        axis='edge_weight_noise',
        strength=strength,
        seed=seed,
        metadata={'distribution': 'uniform', 'weight_factor_mean': float(np.mean(weight_factors))},
    )


def perturb_sparse_defects(
    data: Stage10Complex,
    density: float,
    seed: int,
    weakened_weight_factor: float = 0.1,
) -> PerturbedOperators:
    rng = np.random.default_rng(seed)
    count = max(1, int(round(density * len(data.edge_axes))))
    chosen = rng.choice(len(data.edge_axes), size=count, replace=False)
    weight_factors = np.ones(len(data.edge_axes), dtype=float)
    weight_factors[chosen] = weakened_weight_factor
    edge_scales = np.sqrt(weight_factors)
    return _operators_from_edge_scales(
        data=data,
        edge_scales=edge_scales,
        axis='sparse_graph_defects',
        strength=density,
        seed=seed,
        metadata={
            'defect_count': int(count),
            'defect_density': float(density),
            'weakened_weight_factor': float(weakened_weight_factor),
        },
    )


def perturb_anisotropic_bias(
    data: Stage10Complex,
    factor: float,
    seed: int,
    preferred_axis: str = 'x',
) -> PerturbedOperators:
    del seed
    edge_scales = np.ones(len(data.edge_axes), dtype=float)
    edge_axes = np.asarray(data.edge_axes)
    edge_scales[edge_axes == preferred_axis] = math.sqrt(factor)
    return _operators_from_edge_scales(
        data=data,
        edge_scales=edge_scales,
        axis='anisotropic_bias',
        strength=factor,
        seed=0,
        metadata={'preferred_axis': preferred_axis, 'bias_factor': float(factor)},
    )


def perturb_kernel_width_jitter(
    data: Stage10Complex,
    epsilon: float,
    strength: float,
    seed: int,
) -> PerturbedOperators:
    rng = np.random.default_rng(seed)
    if data.boundary_type == 'periodic':
        h = 1.0 / data.n_side
    else:
        h = 1.0 / max(data.n_side - 1, 1)
    local_factors = 1.0 + rng.uniform(-strength, strength, size=len(data.edge_axes))
    local_epsilon = np.clip(epsilon * local_factors, 1.0e-6, None)
    base_weight = math.exp(-(h * h) / (2.0 * epsilon))
    local_weights = np.exp(-(h * h) / (2.0 * local_epsilon))
    weight_factors = np.clip(local_weights / base_weight, 1.0e-3, None)
    edge_scales = np.sqrt(weight_factors)
    return _operators_from_edge_scales(
        data=data,
        edge_scales=edge_scales,
        axis='kernel_width_jitter',
        strength=strength,
        seed=seed,
        metadata={'epsilon': float(epsilon), 'jitter_strength': float(strength)},
    )


def build_perturbed_operators(
    data: Stage10Complex,
    perturbation_axis: str,
    perturbation_strength: float,
    seed: int,
    epsilon: float,
) -> PerturbedOperators:
    if perturbation_axis == 'edge_weight_noise':
        return perturb_edge_weight_noise(data, perturbation_strength, seed)
    if perturbation_axis == 'sparse_graph_defects':
        return perturb_sparse_defects(data, perturbation_strength, seed)
    if perturbation_axis == 'anisotropic_bias':
        return perturb_anisotropic_bias(data, perturbation_strength, seed)
    if perturbation_axis == 'kernel_width_jitter':
        return perturb_kernel_width_jitter(data, epsilon, perturbation_strength, seed)
    raise ValueError(f'unsupported perturbation axis: {perturbation_axis}')


def build_paired_perturbed_operators(
    data: Stage10Complex,
    perturbations: list[dict[str, Any]],
    epsilon: float,
    pair_name: str,
) -> PerturbedOperators:
    if len(perturbations) != 2:
        raise ValueError('Atlas-2 paired perturbations require exactly two component perturbations')

    components: list[PerturbedOperators] = []
    combined_scales = np.ones(len(data.edge_axes), dtype=float)
    for spec in perturbations:
        item = build_perturbed_operators(
            data=data,
            perturbation_axis=str(spec['axis']),
            perturbation_strength=float(spec['strength']),
            seed=int(spec.get('seed', 0)),
            epsilon=epsilon,
        )
        components.append(item)
        combined_scales *= item.edge_scales

    metadata = {
        'pair_name': pair_name,
        'components': [
            {
                'axis': item.axis,
                'strength': item.strength,
                'seed': item.seed,
                'metadata': item.metadata,
            }
            for item in components
        ],
    }
    return _operators_from_edge_scales(
        data=data,
        edge_scales=combined_scales,
        axis='paired_perturbation',
        strength=1.0,
        seed=0,
        metadata=metadata,
    )


__all__ = ['PerturbedOperators', 'build_perturbed_operators', 'build_paired_perturbed_operators']
