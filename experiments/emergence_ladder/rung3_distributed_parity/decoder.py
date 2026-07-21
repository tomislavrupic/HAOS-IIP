from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .alphabet import RelationalAlphabet
from .parity_architecture import DistributedParityMemory, syndrome


SYNDROME_TO_LOCAL_ERROR = {(1, 0): 0, (1, 1): 1, (0, 1): 2}


@dataclass(frozen=True)
class DecodeResult:
    target_symbols: np.ndarray
    estimated_error: np.ndarray
    syndrome_before: np.ndarray
    syndrome_after: np.ndarray
    block_status: tuple[str, ...]
    visible_bits_per_block: tuple[int, ...]

    @property
    def syndrome_eliminated(self) -> bool:
        return bool(not np.any(self.syndrome_after))


def decode_local(
    current_symbols: np.ndarray,
    memory: DistributedParityMemory,
    alphabet: RelationalAlphabet,
    available_checks: np.ndarray | None = None,
) -> DecodeResult:
    current = np.asarray(current_symbols, dtype=np.uint8)
    before = syndrome(current, memory, alphabet)
    checks = np.ones_like(before, dtype=bool) if available_checks is None else np.asarray(available_checks, dtype=bool)
    if checks.shape != before.shape:
        raise ValueError("available-check mask shape mismatch")
    estimated = np.zeros_like(current)
    statuses: list[str] = []
    visibility: list[int] = []
    for block_id, block in enumerate(alphabet.blocks):
        visible = checks[block_id]
        visibility.append(3 + int(np.sum(visible)))
        if not np.all(visible):
            statuses.append("unresolved_missing_check")
            continue
        pair = tuple(int(value) for value in before[block_id])
        if pair == (0, 0):
            statuses.append("no_error_or_undetected")
            continue
        local_index = SYNDROME_TO_LOCAL_ERROR[pair]
        estimated[block[local_index]] = 1
        statuses.append("unique_radius_one_correction")
    target = current ^ estimated
    after = syndrome(target, memory, alphabet)
    return DecodeResult(target, estimated, before, after, tuple(statuses), tuple(visibility))

