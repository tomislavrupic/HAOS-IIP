from __future__ import annotations

import math
from collections import Counter, deque
from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy.sparse as sp


@dataclass(frozen=True)
class DK2DComplex:
    n_side: int
    epsilon: float
    cycle_phase_x: float
    cycle_phase_y: float
    points: np.ndarray
    edge_midpoints: np.ndarray
    face_centers: np.ndarray
    edge_axes: np.ndarray
    edge_weights: np.ndarray
    face_weights: np.ndarray
    d0: sp.csr_matrix
    d1: sp.csr_matrix
    delta1: sp.csr_matrix
    delta2: sp.csr_matrix
    delta_h: sp.csr_matrix
    dirac_kahler: sp.csr_matrix
    block_sizes: tuple[int, int, int]


def periodic_displacement(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    anchor = np.asarray(anchor, dtype=float)
    delta = points - anchor
    return (delta + 0.5) % 1.0 - 0.5


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


def pack_positions(*blocks: np.ndarray) -> np.ndarray:
    return np.vstack([np.asarray(block, dtype=float) for block in blocks])


def edge_phase_x(i: int, n_side: int, cycle_phase_x: float) -> complex:
    if cycle_phase_x == 0.0:
        return 1.0 + 0.0j
    return np.exp(1j * cycle_phase_x) if i == n_side - 1 else 1.0 + 0.0j


def edge_phase_y(j: int, n_side: int, cycle_phase_y: float) -> complex:
    if cycle_phase_y == 0.0:
        return 1.0 + 0.0j
    return np.exp(1j * cycle_phase_y) if j == n_side - 1 else 1.0 + 0.0j


def build_dk2d_complex(
    n_side: int,
    epsilon: float,
    cycle_phase_x: float = 0.0,
    cycle_phase_y: float = 0.0,
) -> DK2DComplex:
    h = 1.0 / n_side
    edge_weight = math.exp(-(h * h) / (2.0 * epsilon))
    face_weight = edge_weight
    root_edge = math.sqrt(edge_weight)
    root_face = math.sqrt(face_weight)

    node_index = np.arange(n_side * n_side).reshape((n_side, n_side))
    points = np.array([[i / n_side, j / n_side] for i in range(n_side) for j in range(n_side)], dtype=float)

    edge_rows: list[int] = []
    edge_cols: list[int] = []
    edge_data: list[complex] = []
    edge_midpoints: list[np.ndarray] = []
    edge_axes: list[int] = []
    edge_map: dict[tuple[str, int, int], int] = {}

    edge_id = 0
    for i in range(n_side):
        for j in range(n_side):
            u = int(node_index[i, j])
            v = int(node_index[(i + 1) % n_side, j])
            phase = edge_phase_x(i, n_side, cycle_phase_x)
            edge_map[("x", i, j)] = edge_id
            edge_rows.extend([edge_id, edge_id])
            edge_cols.extend([u, v])
            edge_data.extend([-root_edge + 0.0j, root_edge * phase])
            edge_midpoints.append(np.array([((i + 0.5) / n_side) % 1.0, j / n_side], dtype=float))
            edge_axes.append(0)
            edge_id += 1
    for i in range(n_side):
        for j in range(n_side):
            u = int(node_index[i, j])
            v = int(node_index[i, (j + 1) % n_side])
            phase = edge_phase_y(j, n_side, cycle_phase_y)
            edge_map[("y", i, j)] = edge_id
            edge_rows.extend([edge_id, edge_id])
            edge_cols.extend([u, v])
            edge_data.extend([-root_edge + 0.0j, root_edge * phase])
            edge_midpoints.append(np.array([i / n_side, ((j + 0.5) / n_side) % 1.0], dtype=float))
            edge_axes.append(1)
            edge_id += 1

    n_nodes = n_side * n_side
    n_edges = edge_id
    d0 = sp.coo_matrix((edge_data, (edge_rows, edge_cols)), shape=(n_edges, n_nodes), dtype=complex).tocsr()

    face_rows: list[int] = []
    face_cols: list[int] = []
    face_data: list[complex] = []
    face_centers: list[np.ndarray] = []
    face_id = 0
    for i in range(n_side):
        for j in range(n_side):
            ex0 = edge_map[("x", i, j)]
            ey1 = edge_map[("y", (i + 1) % n_side, j)]
            ex1 = edge_map[("x", i, (j + 1) % n_side)]
            ey0 = edge_map[("y", i, j)]
            ux = edge_phase_x(i, n_side, cycle_phase_x)
            uy = edge_phase_y(j, n_side, cycle_phase_y)
            entries = [
                (ex0, +root_face),
                (ey1, +root_face * ux),
                (ex1, -root_face * uy),
                (ey0, -root_face),
            ]
            for edge_idx, value in entries:
                face_rows.append(face_id)
                face_cols.append(edge_idx)
                face_data.append(value)
            face_centers.append(np.array([((i + 0.5) / n_side) % 1.0, ((j + 0.5) / n_side) % 1.0], dtype=float))
            face_id += 1

    n_faces = face_id
    d1 = sp.coo_matrix((face_data, (face_rows, face_cols)), shape=(n_faces, n_edges), dtype=complex).tocsr()
    delta1 = d0.getH().tocsr()
    delta2 = d1.getH().tocsr()
    delta0_op = (delta1 @ d0).tocsr()
    delta1_op = (d0 @ delta1 + delta2 @ d1).tocsr()
    delta2_op = (d1 @ delta2).tocsr()
    delta_h = sp.block_diag((delta0_op, delta1_op, delta2_op), format="csr")

    zero00 = sp.csr_matrix((n_nodes, n_nodes), dtype=complex)
    zero01 = sp.csr_matrix((n_nodes, n_faces), dtype=complex)
    zero10 = sp.csr_matrix((n_edges, n_edges), dtype=complex)
    zero20 = sp.csr_matrix((n_faces, n_nodes), dtype=complex)
    zero22 = sp.csr_matrix((n_faces, n_faces), dtype=complex)
    dirac_kahler = sp.bmat(
        [
            [zero00, delta1, zero01],
            [d0, zero10, delta2],
            [zero20, d1, zero22],
        ],
        format="csr",
    )

    return DK2DComplex(
        n_side=n_side,
        epsilon=epsilon,
        cycle_phase_x=cycle_phase_x,
        cycle_phase_y=cycle_phase_y,
        points=points,
        edge_midpoints=np.asarray(edge_midpoints, dtype=float),
        face_centers=np.asarray(face_centers, dtype=float),
        edge_axes=np.asarray(edge_axes, dtype=int),
        edge_weights=np.full(n_edges, edge_weight, dtype=float),
        face_weights=np.full(n_faces, face_weight, dtype=float),
        d0=d0,
        d1=d1,
        delta1=delta1,
        delta2=delta2,
        delta_h=delta_h,
        dirac_kahler=dirac_kahler,
        block_sizes=(n_nodes, n_edges, n_faces),
    )


def build_graph(config: dict[str, Any]) -> Any:
    kind = str(config.get("kind", "dk2d_periodic"))
    if kind == "dk2d_periodic":
        return build_dk2d_complex(
            n_side=int(config["n_side"]),
            epsilon=float(config["epsilon"]),
            cycle_phase_x=float(config.get("cycle_phase_x", 0.0)),
            cycle_phase_y=float(config.get("cycle_phase_y", 0.0)),
        )
    raise ValueError(f"unsupported graph kind: {kind}")


def packet_state_with_profile(
    positions: np.ndarray,
    block_sizes: tuple[int, ...],
    center: np.ndarray,
    sigma: float,
    amplitude: float,
    phase_offset: float,
    kick_vector: np.ndarray,
    kick_cycles: float,
    grade0_scale: float,
    grade1_scale: float,
    anisotropy: tuple[float, float],
) -> np.ndarray:
    axis_x = max(sigma * anisotropy[0], 1.0e-12)
    axis_y = max(sigma * anisotropy[1], 1.0e-12)
    parts: list[np.ndarray] = []
    start = 0
    for grade_idx, size in enumerate(block_sizes):
        coords = positions[start:start + size]
        delta = periodic_displacement(coords, center)
        scaled_r2 = (delta[:, 0] / axis_x) ** 2 + (delta[:, 1] / axis_y) ** 2
        profile = np.exp(-0.5 * scaled_r2)
        base_phase = 2.0 * math.pi * kick_cycles * (delta @ kick_vector)
        phase = np.exp(1j * (base_phase + phase_offset))
        if grade_idx == 0:
            part = float(grade0_scale) * profile * phase
        elif grade_idx == 1:
            part = float(grade1_scale) * profile * phase
        else:
            part = np.zeros(size, dtype=complex)
        parts.append(np.asarray(part, dtype=complex))
        start += size
    psi = np.concatenate(parts)
    psi /= max(np.linalg.norm(psi), 1.0e-12)
    return float(amplitude) * psi


def local_grade_mix(vec: np.ndarray, block_sizes: tuple[int, ...]) -> np.ndarray:
    n0, n1, _ = block_sizes
    q0 = vec[:n0]
    q1 = vec[n0:n0 + n1]
    out = np.zeros_like(vec)
    m = min(n0, n1)
    out[:m] += q1[:m]
    out[n0:n0 + m] += q0[:m]
    return out


def evolve_signed(
    base_operator: sp.csr_matrix,
    block_sizes: tuple[int, ...],
    state0: np.ndarray,
    dt: float,
    steps: int,
    beta: float,
    direction: float,
) -> list[np.ndarray]:
    state = state0.copy()
    states = [state.copy()]
    base_norm = float(np.linalg.norm(state0))
    for _ in range(steps):
        acc = -direction * (base_operator @ state)
        if beta > 0.0:
            acc += direction * beta * local_grade_mix(state, block_sizes)
        state = state + dt * acc
        norm = float(np.linalg.norm(state))
        if not np.isfinite(norm) or norm <= 1.0e-12:
            raise RuntimeError("numerical instability indicator triggered")
        state = state / norm * base_norm
        states.append(state.copy())
    return states


def evolve_unsigned(
    base_operator: sp.csr_matrix,
    block_sizes: tuple[int, ...],
    state0: np.ndarray,
    dt: float,
    steps: int,
    beta: float,
) -> list[np.ndarray]:
    return evolve_signed(base_operator, block_sizes, state0, dt, steps, beta, direction=1.0)


def evolve_grade_locked(
    base_operator: sp.csr_matrix,
    block_sizes: tuple[int, ...],
    state0: np.ndarray,
    dt: float,
    steps: int,
    beta: float,
    lock_alpha: float,
) -> list[np.ndarray]:
    state = state0.copy()
    states = [state.copy()]
    base_norm = float(np.linalg.norm(state0))
    n0, n1, _ = block_sizes
    init_energy = np.asarray(
        [
            float(np.sum(np.abs(state0[:n0]) ** 2)),
            float(np.sum(np.abs(state0[n0:n0 + n1]) ** 2)),
            float(np.sum(np.abs(state0[n0 + n1:]) ** 2)),
        ],
        dtype=float,
    )
    init_fractions = init_energy / max(float(np.sum(init_energy)), 1.0e-12)
    for _ in range(steps):
        acc = -(base_operator @ state)
        if beta > 0.0:
            acc += float(beta) * local_grade_mix(state, block_sizes)
        state = state + dt * acc
        if lock_alpha > 0.0:
            parts = [state[:n0].copy(), state[n0:n0 + n1].copy(), state[n0 + n1:].copy()]
            total = float(np.sum(np.abs(state) ** 2))
            for idx, part in enumerate(parts):
                current = float(np.sum(np.abs(part) ** 2))
                target = (1.0 - float(lock_alpha)) * current + float(lock_alpha) * init_fractions[idx] * total
                if current > 1.0e-12 and target > 0.0:
                    part *= math.sqrt(target / current)
            state = np.concatenate(parts)
        norm = float(np.linalg.norm(state))
        if not np.isfinite(norm) or norm <= 1.0e-12:
            raise RuntimeError("numerical instability indicator triggered")
        state = state / norm * base_norm
        states.append(state.copy())
    return states


def run_transport(graph: Any, kernel_config: dict[str, Any]) -> list[np.ndarray]:
    mode = str(kernel_config.get("mode", "signed"))
    state0 = np.asarray(kernel_config["state0"], dtype=complex)
    dt = float(kernel_config["dt"])
    steps = int(kernel_config["steps"])
    beta = float(kernel_config.get("beta", 0.0))
    if mode == "signed":
        return evolve_signed(graph.dirac_kahler, graph.block_sizes, state0, dt, steps, beta, float(kernel_config.get("direction", 1.0)))
    if mode == "unsigned":
        return evolve_unsigned(graph.dirac_kahler, graph.block_sizes, state0, dt, steps, beta)
    if mode == "grade_locked":
        return evolve_grade_locked(
            graph.dirac_kahler,
            graph.block_sizes,
            state0,
            dt,
            steps,
            beta,
            float(kernel_config.get("lock_alpha", 0.0)),
        )
    raise ValueError(f"unsupported transport mode: {mode}")


def field_grid_2d(positions: np.ndarray, vec: np.ndarray, n_side: int) -> np.ndarray:
    grid = np.zeros((n_side, n_side), dtype=float)
    coords = np.floor(np.asarray(positions) * n_side).astype(int) % n_side
    weights = np.abs(vec) ** 2
    for idx, coord in enumerate(coords):
        grid[coord[0], coord[1]] += float(weights[idx])
    return grid


def local_peak_candidates(grid: np.ndarray) -> list[tuple[float, tuple[int, int]]]:
    n_side = grid.shape[0]
    peaks: list[tuple[float, tuple[int, int]]] = []
    for i in range(n_side):
        for j in range(n_side):
            value = float(grid[i, j])
            if value <= 0.0:
                continue
            is_peak = True
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    if di == 0 and dj == 0:
                        continue
                    ni = (i + di) % n_side
                    nj = (j + dj) % n_side
                    if grid[ni, nj] > value:
                        is_peak = False
                        break
                if not is_peak:
                    break
            if is_peak:
                peaks.append((value, (i, j)))
    if not peaks:
        flat_idx = int(np.argmax(grid))
        peaks.append((float(grid.flat[flat_idx]), np.unravel_index(flat_idx, grid.shape)))
    peaks.sort(key=lambda item: item[0], reverse=True)
    return peaks


def top_peak_positions(grid: np.ndarray, count: int) -> tuple[list[np.ndarray], int]:
    peaks = local_peak_candidates(grid)
    raw_count = len(peaks)
    taken = set()
    positions: list[np.ndarray] = []
    n_side = grid.shape[0]
    for _, (i, j) in peaks:
        if len(positions) >= count:
            break
        taken.add((i, j))
        positions.append(np.array([(i + 0.5) / n_side, (j + 0.5) / n_side], dtype=float))
    if len(positions) < count:
        flat_order = np.argsort(grid.ravel())[::-1]
        for flat_idx in flat_order:
            ij = np.unravel_index(int(flat_idx), grid.shape)
            if ij in taken:
                continue
            positions.append(np.array([(ij[0] + 0.5) / n_side, (ij[1] + 0.5) / n_side], dtype=float))
            taken.add(ij)
            if len(positions) >= count:
                break
    while len(positions) < count:
        positions.append(positions[-1].copy())
    return positions[:count], raw_count


def pairing_cost(points_a: list[np.ndarray], points_b: list[np.ndarray]) -> float:
    total = 0.0
    for a, b in zip(points_a, points_b):
        total += float(np.linalg.norm(periodic_displacement(np.asarray([b]), np.asarray(a))[0]))
    return total


def ordered_peaks(
    grid: np.ndarray,
    count: int,
    anchors: list[np.ndarray] | None,
    previous: list[np.ndarray] | None,
) -> tuple[list[np.ndarray], int]:
    from itertools import permutations

    peaks, raw_count = top_peak_positions(grid, count)
    if anchors is not None:
        best = min(permutations(peaks, len(peaks)), key=lambda perm: pairing_cost(list(anchors), list(perm)))
        return [np.asarray(p, dtype=float) for p in best], raw_count
    if previous is not None:
        best = min(permutations(peaks, len(peaks)), key=lambda perm: pairing_cost(list(previous), list(perm)))
        return [np.asarray(p, dtype=float) for p in best], raw_count
    peaks.sort(key=lambda p: (p[0], p[1]))
    return peaks, raw_count


def mean_pair_distance(centers: list[np.ndarray]) -> float:
    from itertools import combinations

    values = []
    for i, j in combinations(range(len(centers)), 2):
        delta = periodic_displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float))[0]
        values.append(float(np.linalg.norm(delta)))
    return float(np.mean(values)) if values else 0.0


def min_pair_distance(centers: list[np.ndarray]) -> float:
    from itertools import combinations

    values = []
    for i, j in combinations(range(len(centers)), 2):
        delta = periodic_displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float))[0]
        values.append(float(np.linalg.norm(delta)))
    return float(np.min(values)) if values else 0.0


def mean_pair_separation_series(center_histories: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    mean_vals: list[float] = []
    min_vals: list[float] = []
    steps = len(center_histories[0]) if center_histories else 0
    for idx in range(steps):
        centers = [history[idx] for history in center_histories]
        mean_vals.append(mean_pair_distance(centers))
        min_vals.append(min_pair_distance(centers))
    return np.asarray(mean_vals, dtype=float), np.asarray(min_vals, dtype=float)


def pair_phase_difference(state: np.ndarray, packet_states: list[np.ndarray]) -> float:
    refs = [float(np.angle(np.vdot(packet, state))) for packet in packet_states]
    if len(refs) < 2:
        return 0.0
    diff = refs[1] - refs[0]
    return float((diff + math.pi) % (2.0 * math.pi) - math.pi)


def overlap_phase_alignment(state: np.ndarray, packet_states: list[np.ndarray]) -> float:
    phases = [float(np.angle(np.vdot(packet, state))) for packet in packet_states]
    if len(phases) < 2:
        return 0.0
    anchor = phases[0]
    diffs = []
    for phase in phases[1:]:
        diff = (phase - anchor + math.pi) % (2.0 * math.pi) - math.pi
        diffs.append(abs(diff))
    mean_abs = float(np.mean(diffs)) if diffs else 0.0
    return float(np.clip(1.0 - mean_abs / math.pi, 0.0, 1.0))


def edge_grid(edge_midpoints: np.ndarray, edge_values: np.ndarray, n_side: int) -> np.ndarray:
    grid = np.zeros((n_side, n_side), dtype=float)
    coords = np.floor(np.asarray(edge_midpoints) * n_side).astype(int) % n_side
    for idx, coord in enumerate(coords):
        grid[coord[0], coord[1]] += float(edge_values[idx])
    return grid


def connected_components(mask: np.ndarray) -> int:
    seen = np.zeros_like(mask, dtype=bool)
    n0, n1 = mask.shape
    count = 0
    for i in range(n0):
        for j in range(n1):
            if not mask[i, j] or seen[i, j]:
                continue
            count += 1
            stack = [(i, j)]
            seen[i, j] = True
            while stack:
                x, y = stack.pop()
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    nx = x + dx
                    ny = y + dy
                    if 0 <= nx < n0 and 0 <= ny < n1 and mask[nx, ny] and not seen[nx, ny]:
                        seen[nx, ny] = True
                        stack.append((nx, ny))
    return count


def grade_exchange_signal(grade_hist: np.ndarray) -> tuple[np.ndarray, float]:
    signal = np.linalg.norm(grade_hist - grade_hist[0][None, :], axis=1)
    peak = float(np.max(signal)) if signal.size else 0.0
    if peak <= 1.0e-12:
        return signal, 0.0
    coherence = float(np.mean(signal) / peak)
    return signal, max(0.0, min(1.0, coherence))


def flow_concentration(grid: np.ndarray) -> float:
    total = float(np.sum(grid))
    if total <= 1.0e-12:
        return 0.0
    flat = np.sort(grid.ravel())[::-1]
    top_k = max(1, int(math.ceil(0.1 * flat.size)))
    return float(np.sum(flat[:top_k]) / total)


def count_channels(grid: np.ndarray) -> int:
    threshold = max(0.15 * float(np.max(grid)), float(np.mean(grid) + np.std(grid)))
    return connected_components(grid >= threshold)


def sigma_pair(base_sigma: float, ratio_a_to_b: float) -> tuple[float, float]:
    root = math.sqrt(max(ratio_a_to_b, 1.0e-12))
    return base_sigma * root, base_sigma / root


def topology_survival_time(metrics: dict[str, Any], dt: float) -> float:
    topology = str(metrics["topology_class"])
    trace = [str(label) for label in metrics["topology_trace"]]
    return float(sum(label == topology for label in trace) * dt)


def separation_oscillation_indicator(center_histories: list[list[list[float]]]) -> float:
    histories = [np.asarray(track, dtype=float) for track in center_histories]
    mean_pair_distances, _ = mean_pair_separation_series(histories)
    if mean_pair_distances.size < 6:
        return 0.0
    tail = mean_pair_distances[len(mean_pair_distances) // 3 :]
    centered = tail - float(np.mean(tail))
    amplitude = float(np.max(np.abs(centered))) if centered.size else 0.0
    if amplitude <= 1.0e-12:
        return 0.0
    zero_crossings = int(np.sum(centered[1:] * centered[:-1] < 0.0))
    return float(zero_crossings / max(1, tail.size - 1))


def dominant_label(labels: list[str], anchor_label: str) -> str:
    counts = Counter(labels)
    max_count = max(counts.values())
    winners = [label for label, count in counts.items() if count == max_count]
    if anchor_label in winners:
        return anchor_label
    priority = ["braid_like_exchange", "transfer_smeared", "localized_capture", "dispersive_pass"]
    return min(winners, key=lambda label: priority.index(label) if label in priority else len(priority))


def derive_phase_label(
    topology_class: str,
    refinement_flag: bool,
    weak_coupling_flag: bool,
    bidirectional_flag: bool,
) -> str:
    stable = bool(refinement_flag and weak_coupling_flag and bidirectional_flag)
    if topology_class == "braid_like_exchange" and stable:
        return "stable_braid_phase"
    if topology_class == "transfer_smeared" and stable:
        return "smeared_transfer_phase"
    if topology_class == "localized_capture" and stable:
        return "localized_encounter_phase"
    return "transient_mixed_phase"


def build_operator_scaled_01(base_operator: sp.csr_matrix, block_sizes: tuple[int, ...], scale: float) -> sp.csr_matrix:
    if abs(float(scale) - 1.0) <= 1.0e-12:
        return base_operator
    n0, n1, _ = block_sizes
    lil = base_operator.tocsr().astype(complex).tolil(copy=True)
    for row_start, row_stop, col_start, col_stop in (
        (0, n0, n0, n0 + n1),
        (n0, n0 + n1, 0, n0),
    ):
        for row in range(row_start, row_stop):
            cols = lil.rows[row]
            data = lil.data[row]
            for idx, col in enumerate(cols):
                if col_start <= col < col_stop:
                    data[idx] *= float(scale)
    return lil.tocsr()


def apply_selector(state: dict[str, Any], selector_config: dict[str, Any]) -> Any:
    kind = str(selector_config.get("kind", "threshold_rule"))
    if kind == "threshold_rule":
        value = float(state[selector_config["feature"]])
        threshold = float(selector_config["threshold"])
        direction = str(selector_config.get("direction", "<="))
        return value <= threshold if direction == "<=" else value >= threshold
    if kind == "conjunction":
        return all(bool(state.get(feature)) == bool(expected) for feature, expected in selector_config["features"].items())
    if kind == "phase_label":
        return derive_phase_label(
            topology_class=str(state["topology_class"]),
            refinement_flag=bool(state["refinement_stability_flag"]),
            weak_coupling_flag=bool(state["weak_coupling_stability_flag"]),
            bidirectional_flag=bool(state["bidirectional_stability_flag"]),
        )
    if kind == "topology_classification":
        packet_count = int(selector_config["packet_count"])
        braid_flag = bool(selector_config.get("braid_flag", False))
        concentration = float(state["flow_concentration_index"])
        coherence = float(state["grade_exchange_coherence"])
        amplitude = float(state["grade_exchange_amplitude"])
        recurrence = float(state.get("recurrence_indicator", 0.0))
        close_fraction = float(state.get("close_fraction", 0.0))
        thresholds = selector_config.get("thresholds", {})
        if packet_count == 2 and braid_flag and concentration >= float(thresholds.get("braid_concentration", 0.885)) and coherence >= float(thresholds.get("braid_coherence", 0.40)):
            return "braid_like_exchange"
        if close_fraction >= float(thresholds.get("capture_close_fraction", 0.58)) and recurrence >= float(thresholds.get("capture_recurrence", 0.08)) and concentration >= float(thresholds.get("capture_concentration", 0.72)):
            return "localized_capture"
        if coherence >= float(thresholds.get("smeared_coherence", 0.33)) and amplitude >= float(thresholds.get("smeared_amplitude", 0.08)):
            return "transfer_smeared"
        return str(selector_config.get("fallback", "dispersive_pass"))
    raise ValueError(f"unsupported selector kind: {kind}")


def compute_invariants(state: dict[str, Any]) -> dict[str, Any]:
    invariants: dict[str, Any] = {}
    if "edge_grid" in state:
        grid = np.asarray(state["edge_grid"], dtype=float)
        invariants["flow_concentration_index"] = flow_concentration(grid)
        invariants["channel_count"] = count_channels(grid)
    if "grade_hist" in state:
        signal, coherence = grade_exchange_signal(np.asarray(state["grade_hist"], dtype=float))
        invariants["grade_exchange_trace"] = signal.tolist()
        invariants["grade_exchange_coherence"] = coherence
        invariants["grade_exchange_amplitude"] = float(np.max(signal)) if signal.size else 0.0
    if "metrics" in state and "dt" in state:
        invariants["topology_survival_time"] = topology_survival_time(state["metrics"], float(state["dt"]))
    if "center_histories" in state:
        invariants["separation_oscillation_indicator"] = separation_oscillation_indicator(state["center_histories"])
    return invariants
