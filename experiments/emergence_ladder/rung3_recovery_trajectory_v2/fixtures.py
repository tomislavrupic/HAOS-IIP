from __future__ import annotations

from dataclasses import dataclass

import numpy as np


PERTURBATION_FAMILIES = (
    "localized_block_shock",
    "sparse_node_corruption",
    "relational_twist",
)


@dataclass(frozen=True)
class GridFixture:
    n_side: int
    incidence: np.ndarray
    edges: tuple[tuple[int, int, str], ...]
    probes: np.ndarray


def node_id(i: int, j: int, n_side: int) -> int:
    return (i % n_side) * n_side + (j % n_side)


def make_fixture(n_side: int = 8) -> GridFixture:
    edges: list[tuple[int, int, str]] = []
    for i in range(n_side):
        for j in range(n_side):
            source = node_id(i, j, n_side)
            edges.extend(
                [
                    (source, node_id(i + 1, j, n_side), "horizontal"),
                    (source, node_id(i, j + 1, n_side), "vertical"),
                    (source, node_id(i + 1, j + 1, n_side), "diagonal_positive"),
                    (source, node_id(i + 1, j - 1, n_side), "diagonal_negative"),
                ]
            )
    incidence = np.zeros((len(edges), n_side * n_side), dtype=float)
    for index, (source, target, _) in enumerate(edges):
        incidence[index, source] = -1.0
        incidence[index, target] = 1.0
    coords = np.array([(i / n_side, j / n_side) for i in range(n_side) for j in range(n_side)], dtype=float)
    probes = np.vstack(
        [
            np.cos(2.0 * np.pi * coords[:, 0]),
            np.sin(2.0 * np.pi * coords[:, 0]),
            np.cos(2.0 * np.pi * coords[:, 1]),
            np.sin(2.0 * np.pi * coords[:, 1]),
        ]
    )
    probes /= np.linalg.norm(probes, axis=1, keepdims=True)
    return GridFixture(n_side, incidence, tuple(edges), probes)


def initial_state(fixture: GridFixture, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = fixture.n_side
    coords = np.array([(i / n, j / n) for i in range(n) for j in range(n)], dtype=float)
    coefficients = rng.normal(size=6)
    state = (
        coefficients[0] * np.cos(2.0 * np.pi * coords[:, 0])
        + coefficients[1] * np.sin(2.0 * np.pi * coords[:, 0])
        + coefficients[2] * np.cos(2.0 * np.pi * coords[:, 1])
        + coefficients[3] * np.sin(2.0 * np.pi * coords[:, 1])
        + 0.45 * coefficients[4] * np.cos(2.0 * np.pi * (coords[:, 0] + coords[:, 1]))
        + 0.45 * coefficients[5] * np.sin(2.0 * np.pi * (coords[:, 0] - coords[:, 1]))
    )
    state -= np.mean(state)
    state /= max(float(np.std(state)), 1.0e-12)
    return state


def perturbation(fixture: GridFixture, family: str, magnitude: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed + {name: index * 1009 for index, name in enumerate(PERTURBATION_FAMILIES)}[family])
    n = fixture.n_side
    values = np.zeros(n * n, dtype=float)
    if family == "localized_block_shock":
        i0 = int(rng.integers(0, n))
        j0 = int(rng.integers(0, n))
        for di in range(2):
            for dj in range(2):
                values[node_id(i0 + di, j0 + dj, n)] = rng.choice((-1.0, 1.0))
    elif family == "sparse_node_corruption":
        selected = rng.choice(n * n, size=max(4, (n * n) // 6), replace=False)
        values[selected] = rng.choice((-1.0, 1.0), size=len(selected))
    elif family == "relational_twist":
        phase = float(rng.uniform(0.0, 2.0 * np.pi))
        coords = np.array([(i / n, j / n) for i in range(n) for j in range(n)], dtype=float)
        values = np.sin(4.0 * np.pi * (coords[:, 0] - coords[:, 1]) + phase)
    else:
        raise ValueError(f"unknown perturbation family: {family}")
    values -= np.mean(values)
    values /= max(float(np.sqrt(np.mean(values**2))), 1.0e-12)
    return float(magnitude) * values


def spanning_tree_incidence(fixture: GridFixture) -> np.ndarray:
    n = fixture.n_side
    selected: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(n - 1):
            selected.append((node_id(i, j, n), node_id(i, j + 1, n)))
    for i in range(n - 1):
        selected.append((node_id(i, 0, n), node_id(i + 1, 0, n)))
    incidence = np.zeros((len(selected), n * n), dtype=float)
    for index, (source, target) in enumerate(selected):
        incidence[index, source] = -1.0
        incidence[index, target] = 1.0
    return incidence


def rewired_incidence(fixture: GridFixture, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    node_count = fixture.n_side * fixture.n_side
    incidence = np.zeros_like(fixture.incidence)
    for index in range(len(incidence)):
        source, target = rng.choice(node_count, size=2, replace=False)
        incidence[index, source] = -1.0
        incidence[index, target] = 1.0
    return incidence
