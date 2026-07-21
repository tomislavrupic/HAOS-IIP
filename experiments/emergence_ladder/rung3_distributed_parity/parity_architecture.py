from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np

from .alphabet import RelationalAlphabet


LOCAL_H = np.asarray(((1, 1, 0), (0, 1, 1)), dtype=np.uint8)
PROHIBITED_FIELD_FRAGMENTS = ("continuous", "reference", "checkpoint", "future", "function", "score", "label")


@dataclass(frozen=True)
class DistributedParityMemory:
    """Two parity bits per block; no continuous values or primary symbols."""

    parity_by_block: tuple[tuple[int, int], ...]
    owner_by_block: tuple[int, ...]
    architecture_id: str = "local_3_1_3_coset_v1"

    @classmethod
    def encode(cls, symbols: np.ndarray, alphabet: RelationalAlphabet) -> "DistributedParityMemory":
        z = np.asarray(symbols, dtype=np.uint8)
        parity = []
        for block in alphabet.blocks:
            value = (LOCAL_H @ z[list(block)]) % 2
            parity.append((int(value[0]), int(value[1])))
        return cls(tuple(parity), tuple(alphabet.local_nodes(i)[0] for i in range(len(alphabet.blocks))))

    @property
    def bit_count(self) -> int:
        return 2 * len(self.parity_by_block)

    def parity_array(self) -> np.ndarray:
        return np.asarray(self.parity_by_block, dtype=np.uint8)


def assert_memory_schema() -> None:
    names = {field.name.lower() for field in fields(DistributedParityMemory)}
    for fragment in PROHIBITED_FIELD_FRAGMENTS:
        if any(fragment in name for name in names):
            raise ValueError(f"prohibited parity-memory field: {fragment}")


def global_check_matrix(blocks: int = 16) -> np.ndarray:
    matrix = np.zeros((2 * blocks, 3 * blocks), dtype=np.uint8)
    for block in range(blocks):
        matrix[2 * block : 2 * block + 2, 3 * block : 3 * block + 3] = LOCAL_H
    return matrix


def syndrome(symbols: np.ndarray, memory: DistributedParityMemory, alphabet: RelationalAlphabet) -> np.ndarray:
    z = np.asarray(symbols, dtype=np.uint8)
    observed = []
    for block_id, block in enumerate(alphabet.blocks):
        current = (LOCAL_H @ z[list(block)]) % 2
        observed.append(current ^ np.asarray(memory.parity_by_block[block_id], dtype=np.uint8))
    return np.asarray(observed, dtype=np.uint8)


def information_accounting() -> dict[str, object]:
    matrix = global_check_matrix()
    rank = int(np.linalg.matrix_rank(matrix.astype(float)))
    nullity = int(matrix.shape[1] - rank)
    return {
        "continuous_state_variables": 64,
        "state_precision_bits": 64,
        "full_checkpoint_bits": 4096,
        "primary_relational_symbols": 48,
        "stored_parity_bits": 32,
        "rt02_memory_bits_upper_bound": 256,
        "maximum_local_decoder_visible_bits": 5,
        "global_binary_check_rank": rank,
        "global_binary_nullity": nullity,
        "symbol_vectors_per_parity_coset": 2**nullity,
        "parity_to_checkpoint_ratio": 32 / 4096,
        "parity_to_rt02_ratio": 32 / 256,
        "single_component_has_complete_parity": False,
        "continuous_state_reconstruction_possible": False,
    }

