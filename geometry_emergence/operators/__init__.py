"""Operator primitives for the geometry emergence sandbox."""

from .sandbox import (
    InteractionGraph,
    KernelField,
    TransportOperator,
    apply_transport_step,
    build_clustered_graph,
    build_random_interaction_graph,
    build_transport_operator,
)

__all__ = [
    "InteractionGraph",
    "KernelField",
    "TransportOperator",
    "apply_transport_step",
    "build_clustered_graph",
    "build_random_interaction_graph",
    "build_transport_operator",
]

