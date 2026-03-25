"""Deterministic diagnostics for geometry-like structure."""

from .diagnostics import (
    compute_effective_dimension,
    compute_flow_bending_index,
    compute_neighborhood_persistence,
    compute_path_distortion,
    compute_recoverability,
    compute_transport_arrival_times,
)

__all__ = [
    "compute_effective_dimension",
    "compute_flow_bending_index",
    "compute_neighborhood_persistence",
    "compute_path_distortion",
    "compute_recoverability",
    "compute_transport_arrival_times",
]

