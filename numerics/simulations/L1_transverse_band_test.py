#!/usr/bin/env python3

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
for path in (RESULTS, PLOTS):
    path.mkdir(exist_ok=True)

VARIANT_LABELS = {
    "baseline": "baseline periodic torus",
    "puncture": "single cubic puncture",
    "line_defect": "line defect",
    "flux_tube": "flux-tube defect",
    "mild_disorder": "smooth mild-disorder torus",
}

VARIANT_MARKERS = {
    "baseline": "o",
    "puncture": "s",
    "line_defect": "^",
    "flux_tube": "D",
    "mild_disorder": "v",
}


@dataclass(frozen=True)
class ComplexData:
    points: np.ndarray
    edges: list[tuple[int, int]]
    edge_axes: list[str]
    directions: np.ndarray
    midpoints: np.ndarray
    edge_weights: np.ndarray
    face_weights: np.ndarray
    d0: sp.csr_matrix
    d1: sp.csr_matrix
    lower: sp.csr_matrix
    upper: sp.csr_matrix
    L1: sp.csr_matrix
    node_laplacian: sp.csr_matrix
    n_side: int
    variant: str
    removed_nodes: int
    removed_requested_edges: int


def central_pair(n_side: int) -> tuple[int, int]:
    return (n_side // 2 - 1, n_side // 2)


def periodic_delta(values: np.ndarray, center: float) -> np.ndarray:
    delta = np.asarray(values, dtype=float) - float(center)
    return np.minimum(np.abs(delta), 1.0 - np.abs(delta))


def disorder_profile(midpoint: np.ndarray, axis: str) -> float:
    x, y, z = np.asarray(midpoint, dtype=float)
    axis_phase = {"x": 0.0, "y": 2.0 * math.pi / 3.0, "z": 4.0 * math.pi / 3.0}[axis]
    return float(
        0.45 * math.sin(2.0 * math.pi * (x + 0.37 * y) + axis_phase)
        + 0.35 * math.cos(2.0 * math.pi * (y + 0.23 * z) - 0.5 * axis_phase)
        + 0.20 * math.sin(2.0 * math.pi * (z + 0.19 * x + 0.11 * y) + 0.3 * axis_phase)
    )


def build_periodic_complex(
    n_side: int,
    epsilon: float,
    variant: str,
    flux_tube_phase: float = 0.0,
    disorder_strength: float = 0.12,
    cycle_holonomy_phase: float = 0.0,
) -> ComplexData:
    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))
    raw_points = np.array(
        [[i / n_side, j / n_side, k / n_side] for i in range(n_side) for j in range(n_side) for k in range(n_side)],
        dtype=float,
    )
    active_mask = np.ones((n_side, n_side, n_side), dtype=bool)
    removed_edge_keys: set[tuple[str, int, int, int]] = set()

    if variant == "puncture":
        ci0, ci1 = central_pair(n_side)
        for i in (ci0, ci1):
            for j in (ci0, ci1):
                for k in (ci0, ci1):
                    active_mask[i, j, k] = False
    elif variant == "line_defect":
        ci = n_side // 2
        cj = n_side // 2
        for k in range(n_side):
            removed_edge_keys.add(("z", ci, cj, k))
    elif variant not in {"baseline", "flux_tube", "mild_disorder"}:
        raise ValueError(f"Unsupported variant: {variant}")

    active_old_indices = [int(node_index[i, j, k]) for i in range(n_side) for j in range(n_side) for k in range(n_side) if active_mask[i, j, k]]
    old_to_new = {old: new for new, old in enumerate(active_old_indices)}
    points = raw_points[active_old_indices]

    h = 1.0 / n_side
    edge_weight = math.exp(-(h * h) / (2.0 * epsilon))
    basis = {
        "x": np.array([1.0, 0.0, 0.0]),
        "y": np.array([0.0, 1.0, 0.0]),
        "z": np.array([0.0, 0.0, 1.0]),
    }

    def ux(i: int, j: int, k: int) -> complex:
        return np.exp(1j * cycle_holonomy_phase) if i == n_side - 1 else 1.0 + 0.0j

    def uy(i: int, j: int, k: int) -> complex:
        return 1.0 + 0.0j

    def uz(i: int, j: int, k: int) -> complex:
        return 1.0 + 0.0j

    phase_fn = {"x": ux, "y": uy, "z": uz}

    edges: list[tuple[int, int]] = []
    edge_axes: list[str] = []
    directions: list[np.ndarray] = []
    midpoints: list[np.ndarray] = []
    edge_weights: list[float] = []
    edge_map: dict[tuple[str, int, int, int], int] = {}
    b0_rows: list[int] = []
    b0_cols: list[int] = []
    b0_data: list[complex] = []

    def node_active(i: int, j: int, k: int) -> bool:
        return bool(active_mask[i % n_side, j % n_side, k % n_side])

    def add_edge(axis: str, i: int, j: int, k: int, ni: int, nj: int, nk: int) -> None:
        if (axis, i, j, k) in removed_edge_keys:
            return
        if not node_active(i, j, k) or not node_active(ni, nj, nk):
            return
        u_old = int(node_index[i, j, k])
        v_old = int(node_index[ni % n_side, nj % n_side, nk % n_side])
        u = old_to_new[u_old]
        v = old_to_new[v_old]
        edge_idx = len(edges)
        edge_map[(axis, i, j, k)] = edge_idx
        edges.append((u, v))
        edge_axes.append(axis)
        directions.append(basis[axis])
        midpoint = np.array(
            [
                ((i + 0.5) / n_side) % 1.0 if axis == "x" else i / n_side,
                ((j + 0.5) / n_side) % 1.0 if axis == "y" else j / n_side,
                ((k + 0.5) / n_side) % 1.0 if axis == "z" else k / n_side,
            ],
            dtype=float,
        )
        midpoints.append(midpoint)
        local_edge_weight = edge_weight
        if variant == "mild_disorder":
            local_edge_weight *= 1.0 + float(disorder_strength) * disorder_profile(midpoint, axis)
        edge_weights.append(local_edge_weight)
        link_phase = phase_fn[axis](i, j, k)
        b0_rows.extend([edge_idx, edge_idx])
        b0_cols.extend([u, v])
        b0_data.extend([-1.0 + 0.0j, link_phase])

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_edge("x", i, j, k, (i + 1) % n_side, j, k)
                add_edge("y", i, j, k, i, (j + 1) % n_side, k)
                add_edge("z", i, j, k, i, j, (k + 1) % n_side)

    n_nodes = len(points)
    n_edges = len(edges)
    B0 = sp.coo_matrix((b0_data, (b0_rows, b0_cols)), shape=(n_edges, n_nodes), dtype=complex).tocsr()

    c_rows: list[int] = []
    c_cols: list[int] = []
    c_data: list[complex] = []
    face_weights: list[float] = []
    face_idx = 0
    ci = n_side // 2
    cj = n_side // 2

    def add_face(face_type: str, boundary_keys: list[tuple[str, int, int, int, int, complex]]) -> None:
        nonlocal face_idx
        oriented_edges: list[tuple[int, complex]] = []
        local_weights: list[float] = []
        phase_multiplier = 1.0 + 0.0j
        if variant == "flux_tube" and face_type == "xy":
            i0, j0, _, _, _ = boundary_keys[0][1:6]
            if i0 == ci and j0 == cj:
                phase_multiplier = np.exp(1j * flux_tube_phase)
        for boundary_pos, (axis, i, j, k, sign, phase) in enumerate(boundary_keys):
            key = (axis, i % n_side, j % n_side, k % n_side)
            if key not in edge_map:
                return
            edge_id = edge_map[key]
            coeff = complex(sign) * phase
            if boundary_pos == 0:
                coeff *= phase_multiplier
            oriented_edges.append((edge_id, coeff))
            local_weights.append(edge_weights[edge_id])
        for edge_id, coeff in oriented_edges:
            c_rows.append(face_idx)
            c_cols.append(edge_id)
            c_data.append(coeff)
        face_weights.append(float(np.mean(local_weights)))
        face_idx += 1

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_face(
                    "xy",
                    [
                        ("x", i, j, k, +1, 1.0 + 0.0j),
                        ("y", (i + 1) % n_side, j, k, +1, np.conj(ux(i, j, k))),
                        ("x", i, (j + 1) % n_side, k, -1, np.conj(uy(i, j, k))),
                        ("y", i, j, k, -1, 1.0 + 0.0j),
                    ],
                )
                add_face(
                    "xz",
                    [
                        ("x", i, j, k, +1, 1.0 + 0.0j),
                        ("z", (i + 1) % n_side, j, k, +1, np.conj(ux(i, j, k))),
                        ("x", i, j, (k + 1) % n_side, -1, np.conj(uz(i, j, k))),
                        ("z", i, j, k, -1, 1.0 + 0.0j),
                    ],
                )
                add_face(
                    "yz",
                    [
                        ("y", i, j, k, +1, 1.0 + 0.0j),
                        ("z", i, (j + 1) % n_side, k, +1, np.conj(uy(i, j, k))),
                        ("y", i, j, (k + 1) % n_side, -1, np.conj(uz(i, j, k))),
                        ("z", i, j, k, -1, 1.0 + 0.0j),
                    ],
                )

    n_faces = face_idx
    C = sp.coo_matrix((c_data, (c_rows, c_cols)), shape=(n_faces, n_edges), dtype=complex).tocsr()

    edge_weights_arr = np.asarray(edge_weights, dtype=float)
    face_weights_arr = np.asarray(face_weights, dtype=float)
    d0 = (sp.diags(np.sqrt(edge_weights_arr)) @ B0).tocsr()
    d1 = (sp.diags(np.sqrt(face_weights_arr)) @ C).tocsr() if n_faces else sp.csr_matrix((0, n_edges), dtype=complex)
    lower = (d0 @ d0.conj().T).tocsr()
    upper = (d1.conj().T @ d1).tocsr() if n_faces else sp.csr_matrix((n_edges, n_edges), dtype=complex)
    L1 = (lower + upper).tocsr()
    node_laplacian = (d0.conj().T @ d0).tocsr()

    return ComplexData(
        points=points,
        edges=edges,
        edge_axes=edge_axes,
        directions=np.asarray(directions, dtype=float),
        midpoints=np.asarray(midpoints, dtype=float),
        edge_weights=edge_weights_arr,
        face_weights=face_weights_arr,
        d0=d0,
        d1=d1,
        lower=lower,
        upper=upper,
        L1=L1,
        node_laplacian=node_laplacian,
        n_side=n_side,
        variant=variant,
        removed_nodes=int(np.size(active_mask) - np.count_nonzero(active_mask)),
        removed_requested_edges=int(len(removed_edge_keys)),
    )


class ExactProjector:
    def __init__(self, d0: sp.csr_matrix, node_laplacian: sp.csr_matrix):
        self.d0 = d0.tocsr()
        reduced = node_laplacian.tocsc()[1:, 1:]
        self.solve = spla.factorized(reduced)
        self.node_dim = node_laplacian.shape[0]

    def apply(self, vec: np.ndarray) -> np.ndarray:
        rhs = np.asarray(self.d0.conj().T @ vec, dtype=complex)
        coeffs = np.zeros(self.node_dim, dtype=complex)
        coeffs[1:] = self.solve(rhs[1:])
        return np.asarray(self.d0 @ coeffs, dtype=complex)


class HarmonicProjector:
    def __init__(self, basis: np.ndarray):
        self.basis = np.asarray(basis, dtype=complex)

    def apply(self, vec: np.ndarray) -> np.ndarray:
        if self.basis.size == 0:
            return np.zeros_like(vec)
        return self.basis @ (self.basis.conj().T @ vec)


def inverse_participation_ratio(vec: np.ndarray) -> float:
    abs_sq = np.abs(vec) ** 2
    norm_sq = float(np.sum(abs_sq))
    if norm_sq == 0.0:
        return 0.0
    return float(np.sum(abs_sq * abs_sq) / (norm_sq * norm_sq))


def low_eigensystem(operator: sp.csr_matrix | spla.LinearOperator, k: int, tol: float) -> tuple[np.ndarray, np.ndarray]:
    dim = operator.shape[0]
    k_eff = min(k, dim - 2) if dim > 2 else 1
    evals, evecs = spla.eigsh(operator, k=k_eff, which="SA", tol=tol)
    order = np.argsort(evals.real)
    return np.asarray(evals[order], dtype=float), np.asarray(evecs[:, order], dtype=complex)


def harmonic_basis_from_reference(reference: ComplexData, count: int, harmonic_tol: float, tol: float) -> np.ndarray:
    evals, evecs = low_eigensystem(reference.L1, k=count, tol=tol)
    mask = evals < harmonic_tol
    basis = np.asarray(evecs[:, mask], dtype=complex)
    if basis.size:
        basis, _ = np.linalg.qr(basis, mode="reduced")
    return basis


def build_penalized_transverse_operator(
    analysis: ComplexData,
    exact_projector: ExactProjector,
    harmonic_projector: HarmonicProjector,
    penalty: float,
) -> spla.LinearOperator:
    edge_dim = analysis.upper.shape[0]

    def matvec(vec: np.ndarray) -> np.ndarray:
        vec = np.asarray(vec, dtype=complex)
        return np.asarray(analysis.upper @ vec, dtype=complex) + penalty * exact_projector.apply(vec) + penalty * harmonic_projector.apply(vec)

    return spla.LinearOperator((edge_dim, edge_dim), matvec=matvec, dtype=complex)


def project_transverse(vec: np.ndarray, exact_projector: ExactProjector, harmonic_projector: HarmonicProjector) -> np.ndarray:
    projected = vec - exact_projector.apply(vec)
    projected = projected - harmonic_projector.apply(projected)
    return projected


def defect_distances(midpoints: np.ndarray, edge_axes: list[str], variant: str, n_side: int) -> np.ndarray:
    center = 0.5 - 0.5 / n_side
    if variant in {"line_defect", "flux_tube"}:
        dx = periodic_delta(midpoints[:, 0], center)
        dy = periodic_delta(midpoints[:, 1], center)
        return np.sqrt(dx * dx + dy * dy)
    if variant == "mild_disorder":
        strength = np.asarray([abs(disorder_profile(midpoint, axis)) for midpoint, axis in zip(midpoints, edge_axes)], dtype=float)
        max_strength = float(np.max(strength)) or 1.0
        return 1.0 - strength / max_strength
    dx = periodic_delta(midpoints[:, 0], center)
    dy = periodic_delta(midpoints[:, 1], center)
    dz = periodic_delta(midpoints[:, 2], center)
    return np.sqrt(dx * dx + dy * dy + dz * dz)


def local_support_fraction(mode: np.ndarray, distances: np.ndarray, variant: str, n_side: int) -> float:
    h = 1.0 / n_side
    if variant == "puncture":
        mask = distances <= 2.5 * h
    elif variant in {"line_defect", "flux_tube"}:
        mask = distances <= 1.5 * h
    elif variant == "mild_disorder":
        threshold = float(np.quantile(distances, 0.2)) if len(distances) else 0.0
        mask = distances <= threshold
    else:
        mask = distances <= 2.0 * h
    norm = float(np.sum(np.abs(mode) ** 2)) or 1.0
    return float(np.sum(np.abs(mode[mask]) ** 2) / norm)


def radial_profile(mode: np.ndarray, distances: np.ndarray, bins: int = 20) -> dict[str, list[float]]:
    if len(distances) == 0:
        return {"radius": [], "density": []}
    order = np.argsort(distances)
    d_sorted = distances[order]
    weights = np.abs(mode[order]) ** 2
    edges = np.linspace(float(np.min(d_sorted)), float(np.max(d_sorted)), bins + 1)
    radii: list[float] = []
    density: list[float] = []
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (d_sorted >= left) & (d_sorted < right if right < edges[-1] else d_sorted <= right)
        if not np.any(mask):
            continue
        radii.append(float(np.mean(d_sorted[mask])))
        density.append(float(np.sum(weights[mask]) / np.sum(mask)))
    return {"radius": radii, "density": density}


def dominant_axis(mode: np.ndarray, edge_axes: list[str]) -> str:
    axis_energy = {
        axis: float(np.sum(np.abs(mode[[idx for idx, value in enumerate(edge_axes) if value == axis]]) ** 2))
        for axis in ("x", "y", "z")
    }
    return max(axis_energy, key=axis_energy.get)


def classify_support(mode: np.ndarray, edge_axes: list[str], exact_fraction: float, harmonic_fraction: float, coexact_fraction: float, divergence_norm: float, curl_norm: float, ipr: float) -> str:
    axis = dominant_axis(mode, edge_axes)
    if harmonic_fraction > 0.8 and divergence_norm < 1.0e-8 and curl_norm < 1.0e-8:
        return f"harmonic {axis}-cycle"
    if coexact_fraction > 0.8 and harmonic_fraction < 0.2 and curl_norm >= divergence_norm:
        return f"coexact {axis}-biased circulation"
    if coexact_fraction > 0.5 and harmonic_fraction > 0.2:
        return f"harmonic/coexact {axis}-mix"
    if exact_fraction > 0.7:
        return f"exact {axis}-biased gradient"
    if ipr > 0.06:
        return f"localized {axis}-biased edge concentration"
    return f"mixed {axis}-biased support"


def analyze_full_mode(
    vec: np.ndarray,
    eigenvalue: float,
    reference: ComplexData,
    analysis: ComplexData,
    exact_projector: ExactProjector,
    harmonic_projector: HarmonicProjector,
    mode_index: int,
    distances: np.ndarray,
    variant: str,
    n_side: int,
) -> dict[str, Any]:
    norm_sq = float(np.real(np.vdot(vec, vec))) or 1.0
    exact_part = exact_projector.apply(vec)
    harmonic_part = harmonic_projector.apply(vec)
    exact_fraction = float(np.real(np.vdot(vec, exact_part)) / norm_sq)
    harmonic_fraction = float(np.real(np.vdot(vec, harmonic_part)) / norm_sq)
    coexact_fraction = float(max(0.0, 1.0 - exact_fraction - harmonic_fraction))
    divergence_norm = float(np.linalg.norm(reference.d0.conj().T @ vec))
    curl_norm = float(np.linalg.norm(analysis.d1 @ vec)) if analysis.d1.shape[0] else 0.0
    ipr = inverse_participation_ratio(vec)
    return {
        "mode_index": int(mode_index),
        "eigenvalue": float(eigenvalue),
        "exact_fraction": exact_fraction,
        "harmonic_fraction": harmonic_fraction,
        "coexact_fraction": coexact_fraction,
        "divergence_norm": divergence_norm,
        "curl_norm": curl_norm,
        "ipr": ipr,
        "near_defect_fraction": local_support_fraction(vec, distances, variant=variant, n_side=n_side),
        "support_pattern": classify_support(vec, analysis.edge_axes, exact_fraction, harmonic_fraction, coexact_fraction, divergence_norm, curl_norm, ipr),
    }


def analyze_case(
    n_side: int,
    epsilon: float,
    variant: str,
    restricted_modes: int,
    harmonic_tol: float,
    eig_tol: float,
    penalty: float,
    flux_tube_phase: float,
    disorder_strength: float = 0.12,
    full_mode_scan_count: int | None = None,
    cycle_holonomy_phase: float = 0.0,
) -> dict[str, Any]:
    reference_variant = "baseline" if variant == "flux_tube" else variant
    reference = build_periodic_complex(
        n_side=n_side,
        epsilon=epsilon,
        variant=reference_variant,
        flux_tube_phase=0.0,
        disorder_strength=disorder_strength,
        cycle_holonomy_phase=0.0,
    )
    analysis = build_periodic_complex(
        n_side=n_side,
        epsilon=epsilon,
        variant=variant,
        flux_tube_phase=flux_tube_phase,
        disorder_strength=disorder_strength,
        cycle_holonomy_phase=cycle_holonomy_phase,
    )
    distances = defect_distances(analysis.midpoints, analysis.edge_axes, variant=variant, n_side=n_side)

    harmonic_basis = harmonic_basis_from_reference(reference, count=max(restricted_modes, 12), harmonic_tol=harmonic_tol, tol=eig_tol)
    exact_projector = ExactProjector(reference.d0, reference.node_laplacian)
    harmonic_projector = HarmonicProjector(harmonic_basis)

    full_count = max(int(full_mode_scan_count or 0), restricted_modes, 12)
    full_evals, full_vecs = low_eigensystem(analysis.L1, k=full_count, tol=eig_tol)
    full_modes = [
        analyze_full_mode(
            full_vecs[:, idx],
            float(full_evals[idx]),
            reference,
            analysis,
            exact_projector,
            harmonic_projector,
            idx,
            distances,
            variant,
            n_side,
        )
        for idx in range(len(full_evals))
    ]

    transverse_operator = build_penalized_transverse_operator(analysis, exact_projector, harmonic_projector, penalty=penalty)
    k_eff = min(restricted_modes + 4, analysis.upper.shape[0] - 2) if analysis.upper.shape[0] > 2 else 1
    restricted_spectrum: list[float] = []
    restricted_records: list[dict[str, Any]] = []
    restricted_vectors: list[np.ndarray] = []

    projected_rays: list[tuple[float, np.ndarray]] = []
    raw_evals, raw_vecs_full = spla.eigsh(transverse_operator, k=k_eff, which="SA", tol=eig_tol)
    order = np.argsort(raw_evals.real)
    for idx in order:
        projected_rays.append((float(np.real(raw_evals[idx])), np.asarray(raw_vecs_full[:, idx], dtype=complex)))

    for _, raw_vec in projected_rays:
        vec = project_transverse(raw_vec, exact_projector, harmonic_projector)
        norm = float(np.linalg.norm(vec))
        if norm <= 1.0e-12:
            continue
        vec = vec / norm
        lam = float(np.real(np.vdot(vec, analysis.upper @ vec)))
        if lam < 1.0e-10 and len(restricted_spectrum) > 0:
            continue
        norm_sq = float(np.real(np.vdot(vec, vec))) or 1.0
        exact_part = exact_projector.apply(vec)
        harmonic_part = harmonic_projector.apply(vec)
        exact_fraction = float(np.real(np.vdot(vec, exact_part)) / norm_sq)
        harmonic_fraction = float(np.real(np.vdot(vec, harmonic_part)) / norm_sq)
        coexact_fraction = float(max(0.0, 1.0 - exact_fraction - harmonic_fraction))
        restricted_spectrum.append(lam)
        restricted_vectors.append(vec)
        restricted_records.append(
            {
                "mode_index": len(restricted_spectrum) - 1,
                "eigenvalue": lam,
                "exact_fraction": exact_fraction,
                "harmonic_fraction": harmonic_fraction,
                "coexact_fraction": coexact_fraction,
                "divergence_norm": float(np.linalg.norm(reference.d0.conj().T @ vec)),
                "curl_norm": float(np.linalg.norm(analysis.d1 @ vec)) if analysis.d1.shape[0] else 0.0,
                "ipr": inverse_participation_ratio(vec),
                "near_defect_fraction": local_support_fraction(vec, distances, variant=variant, n_side=n_side),
                "radial_profile": radial_profile(vec, distances) if variant in {"puncture", "line_defect", "flux_tube"} else {"radius": [], "density": []},
            }
        )
        if len(restricted_spectrum) >= restricted_modes:
            break

    return {
        "label": f"{variant}_n{n_side}",
        "config": {
            "variant": variant,
            "variant_label": VARIANT_LABELS[variant],
            "n_side": int(n_side),
            "nodes": int(len(analysis.points)),
            "edges": int(len(analysis.edges)),
            "faces": int(len(analysis.face_weights)),
                "removed_nodes": int(analysis.removed_nodes),
                "removed_requested_edges": int(analysis.removed_requested_edges),
                "epsilon": float(epsilon),
                "flux_tube_phase": float(flux_tube_phase if variant == "flux_tube" else 0.0),
                "disorder_strength": float(disorder_strength if variant == "mild_disorder" else 0.0),
                "cycle_holonomy_phase": float(cycle_holonomy_phase),
            },
        "dimensions": {
            "harmonic": int(harmonic_basis.shape[1]),
        },
        "full_spectrum": [float(value) for value in full_evals[: max(restricted_modes, 12)]],
        "full_modes": full_modes,
        "restricted_transverse_spectrum": restricted_spectrum,
        "restricted_transverse_modes": restricted_records,
        "restricted_vectors": restricted_vectors,
        "midpoints": analysis.midpoints,
        "directions": analysis.directions,
    }


def plot_edge_mode(ax: Any, midpoints: np.ndarray, directions: np.ndarray, mode: np.ndarray, title: str) -> None:
    from matplotlib import cm

    magnitudes = np.abs(mode)
    threshold = np.quantile(magnitudes, 0.85) if len(magnitudes) > 12 else 0.0
    select = magnitudes >= threshold
    if not np.any(select):
        select = np.ones_like(magnitudes, dtype=bool)
    max_mag = float(np.max(magnitudes[select])) or 1.0
    vectors = directions[select] * (np.real(mode[select]) / max_mag)[:, None] * 0.18
    phases = np.angle(mode[select])
    colors = cm.twilight((phases + math.pi) / (2.0 * math.pi))
    ax.scatter(
        midpoints[select, 0],
        midpoints[select, 1],
        midpoints[select, 2],
        c=phases,
        cmap="twilight",
        s=22 + 150 * magnitudes[select] / max_mag,
        vmin=-math.pi,
        vmax=math.pi,
    )
    ax.quiver(
        midpoints[select, 0],
        midpoints[select, 1],
        midpoints[select, 2],
        vectors[:, 0],
        vectors[:, 1],
        vectors[:, 2],
        colors=colors,
        linewidth=0.8,
        arrow_length_ratio=0.25,
    )
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
