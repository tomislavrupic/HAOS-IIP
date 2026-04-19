#!/usr/bin/env python3

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from L1_transverse_band_test import (
    VARIANT_LABELS,
    VARIANT_MARKERS,
    build_periodic_complex,
    defect_distances,
    local_support_fraction,
    radial_profile,
)
from stage8_common import (
    ExactProjector,
    HarmonicProjector,
    build_penalized_transverse_operator,
    harmonic_basis_from_L1,
    inverse_participation_ratio,
    project_transverse,
)


def as_real_sparse(matrix: sp.csr_matrix) -> sp.csr_matrix:
    return matrix.real.tocsr()


def build_real_defect_complex(n_side: int, epsilon: float, variant: str) -> Any:
    data = build_periodic_complex(n_side=n_side, epsilon=epsilon, variant=variant, flux_tube_phase=0.0)
    return SimpleNamespace(
        points=np.asarray(data.points, dtype=float),
        edges=list(data.edges),
        directions=np.asarray(data.directions, dtype=float),
        midpoints=np.asarray(data.midpoints, dtype=float),
        edge_weights=np.asarray(data.edge_weights, dtype=float),
        face_weights=np.asarray(data.face_weights, dtype=float),
        d0=as_real_sparse(data.d0),
        d1=as_real_sparse(data.d1),
        lower=as_real_sparse(data.lower),
        upper=as_real_sparse(data.upper),
        L1=as_real_sparse(data.L1),
        node_laplacian=as_real_sparse(data.node_laplacian),
        n_side=int(data.n_side),
        variant=str(data.variant),
        removed_nodes=int(data.removed_nodes),
        removed_requested_edges=int(data.removed_requested_edges),
    )


def analyze_case(
    n_side: int,
    epsilon: float,
    variant: str,
    phase_modes: int,
    restricted_modes: int,
    harmonic_tol: float,
    eig_tol: float,
    penalty: float,
) -> dict[str, Any]:
    data = build_real_defect_complex(n_side=n_side, epsilon=epsilon, variant=variant)
    harmonic_basis, full_evals, full_evecs = harmonic_basis_from_L1(data.L1, count=max(phase_modes, restricted_modes, 8), harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = ExactProjector(data.d0, data.node_laplacian)
    harmonic_projector = HarmonicProjector(harmonic_basis)

    phase_count = min(phase_modes, len(full_evals))
    phase_records: list[dict[str, Any]] = []
    for idx in range(phase_count):
        vec = np.asarray(full_evecs[:, idx], dtype=float)
        div_norm = float(np.linalg.norm(data.d0.T @ vec))
        curl_norm = float(np.linalg.norm(data.d1 @ vec)) if data.d1.shape[0] else 0.0
        exact_part = exact_projector.apply(vec)
        harmonic_part = harmonic_projector.apply(vec)
        norm_sq = float(np.dot(vec, vec)) or 1.0
        exact_fraction = float(np.dot(exact_part, exact_part) / norm_sq)
        harmonic_fraction = float(np.dot(harmonic_part, harmonic_part) / norm_sq)
        phase_records.append(
            {
                'mode_index': idx,
                'eigenvalue': float(full_evals[idx]),
                'divergence_norm': div_norm,
                'curl_norm': curl_norm,
                'exact_fraction': exact_fraction,
                'harmonic_fraction': harmonic_fraction,
                'coexact_fraction': float(max(0.0, 1.0 - exact_fraction - harmonic_fraction)),
            }
        )

    transverse_operator = build_penalized_transverse_operator(data, exact_projector, harmonic_projector, penalty=penalty)
    k_eff = min(restricted_modes, data.upper.shape[0] - 2) if data.upper.shape[0] > 2 else 1
    restricted_evals, restricted_vecs = spla.eigsh(transverse_operator, k=k_eff, which='SA', tol=eig_tol)
    order = np.argsort(restricted_evals.real)
    restricted_evals = np.asarray(restricted_evals[order], dtype=float)
    restricted_vecs = np.asarray(restricted_vecs[:, order], dtype=float)

    restricted_spectrum: list[float] = []
    restricted_records: list[dict[str, Any]] = []
    distances = defect_distances(data.midpoints, variant=variant, n_side=n_side)
    for idx in range(len(restricted_evals)):
        vec = project_transverse(restricted_vecs[:, idx], exact_projector, harmonic_projector)
        norm = float(np.linalg.norm(vec))
        if norm <= 1.0e-12:
            continue
        vec = vec / norm
        lam = float(np.dot(vec, data.upper @ vec))
        exact_part = exact_projector.apply(vec)
        harmonic_part = harmonic_projector.apply(vec)
        exact_fraction = float(np.dot(exact_part, exact_part))
        harmonic_fraction = float(np.dot(harmonic_part, harmonic_part))
        restricted_spectrum.append(lam)
        restricted_records.append(
            {
                'mode_index': len(restricted_spectrum) - 1,
                'eigenvalue': lam,
                'divergence_norm': float(np.linalg.norm(data.d0.T @ vec)),
                'curl_norm': float(np.linalg.norm(data.d1 @ vec)) if data.d1.shape[0] else 0.0,
                'coexact_fraction': float(max(0.0, 1.0 - exact_fraction - harmonic_fraction)),
                'ipr': inverse_participation_ratio(vec),
                'near_defect_fraction': local_support_fraction(vec, distances, variant=variant, n_side=n_side),
                'radial_profile': radial_profile(vec, distances) if variant in {'puncture', 'line_defect'} else {'radius': [], 'density': []},
            }
        )

    return {
        'label': f'{variant}_n{n_side}',
        'config': {
            'variant': variant,
            'variant_label': VARIANT_LABELS[variant],
            'n_side': int(n_side),
            'nodes': int(len(data.points)),
            'edges': int(len(data.edges)),
            'faces': int(len(data.face_weights)),
            'removed_nodes': int(data.removed_nodes),
            'removed_requested_edges': int(data.removed_requested_edges),
            'epsilon': float(epsilon),
        },
        'dimensions': {
            'harmonic': int(harmonic_basis.shape[1]),
        },
        'phase_modes': phase_records,
        'restricted_transverse_spectrum': restricted_spectrum,
        'restricted_transverse_modes': restricted_records,
    }
