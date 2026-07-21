from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from experiments.emergence_ladder.rung3_recovery_trajectory_v2.fixtures import node_id


@dataclass(frozen=True)
class RelationalAlphabet:
    """Forty-eight local order relations grouped into sixteen 2x2 cells."""

    incidence: np.ndarray
    blocks: tuple[tuple[int, int, int], ...]
    records: tuple[tuple[int, int, str, int], ...]
    n_side: int

    @classmethod
    def build(cls, n_side: int = 8) -> "RelationalAlphabet":
        if n_side % 2:
            raise ValueError("the relational tiling requires an even side length")
        rows: list[np.ndarray] = []
        blocks: list[tuple[int, int, int]] = []
        records: list[tuple[int, int, str, int]] = []
        for i in range(0, n_side, 2):
            for j in range(0, n_side, 2):
                anchor = node_id(i, j, n_side)
                pairs = (
                    (anchor, node_id(i + 1, j, n_side), "horizontal"),
                    (anchor, node_id(i, j + 1, n_side), "vertical"),
                    (anchor, node_id(i + 1, j + 1, n_side), "diagonal"),
                )
                indices = []
                block_id = len(blocks)
                for source, target, name in pairs:
                    row = np.zeros(n_side * n_side, dtype=float)
                    row[source] = -1.0
                    row[target] = 1.0
                    indices.append(len(rows))
                    rows.append(row)
                    records.append((source, target, name, block_id))
                blocks.append(tuple(indices))
        return cls(np.vstack(rows), tuple(blocks), tuple(records), n_side)

    def encode(self, state: np.ndarray) -> np.ndarray:
        values = self.incidence @ np.asarray(state, dtype=float)
        return (values >= 0.0).astype(np.uint8)

    def block_error_counts(self, left: np.ndarray, right: np.ndarray) -> np.ndarray:
        errors = np.asarray(left, dtype=np.uint8) ^ np.asarray(right, dtype=np.uint8)
        return np.asarray([int(np.sum(errors[list(block)])) for block in self.blocks], dtype=int)

    def local_nodes(self, block_id: int) -> tuple[int, ...]:
        nodes: set[int] = set()
        for index in self.blocks[block_id]:
            source, target, _, _ = self.records[index]
            nodes.update((source, target))
        return tuple(sorted(nodes))


def symbol_agreement(expected: np.ndarray, observed: np.ndarray) -> float:
    return float(np.mean(np.asarray(expected, dtype=np.uint8) == np.asarray(observed, dtype=np.uint8)))

