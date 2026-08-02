"""Presentation semantics that never convert negative evidence into run failure."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ResultSemantics:
    execution_health: str
    evidence_outcome: str
    research_state: str
    display_tone: str
    interpretation: str


_VERIFIED_NEGATIVES = {
    "TERMINAL_NEGATIVE_MODE_SAMPLING_CURRENT_REGIME",
    "VERIFIED_NEGATIVE_SYNTHETIC_GENERATIVE_LINE_PAUSED",
}
_OPEN_NEGATIVES = {
    "OPEN_NO_DIRECTED_YIELD_ADVANTAGE",
    "OPEN_NO_UNIQUE_DIRECTED_ADVANTAGE",
    "OPEN_NO_STRUCTURAL_INTERVENTION_ADVANTAGE",
}
_POSITIVE_BOUNDED = {
    "PASS_BOUNDED_DIRECTED_YIELD_ADVANTAGE",
    "PASS_BOUNDED_UNIQUE_DIRECTED_ADVANTAGE",
    "STRUCTURAL_EXPANSION_CONFIRMED",
}
_EXECUTION_FAILURES = {
    "EXECUTION_FAILED",
    "INVALID_ARTIFACT",
    "CONTRACT_VIOLATION",
    "REPRODUCIBILITY_FAILED",
}


def interpret_result_status(status: str) -> ResultSemantics:
    """Separate run validity from the scientific direction of its result."""

    if status in _VERIFIED_NEGATIVES:
        return ResultSemantics(
            execution_health="VERIFIED",
            evidence_outcome="NEGATIVE",
            research_state="TERMINAL_BOUNDARY",
            display_tone="BLUE",
            interpretation="The experiment completed correctly and established a bounded negative result.",
        )
    if status in _OPEN_NEGATIVES:
        return ResultSemantics(
            execution_health="VERIFIED",
            evidence_outcome="NEGATIVE",
            research_state="OPEN",
            display_tone="AMBER",
            interpretation="The experiment completed correctly; the tested advantage was not supported.",
        )
    if status in _POSITIVE_BOUNDED:
        return ResultSemantics(
            execution_health="VERIFIED",
            evidence_outcome="POSITIVE_BOUNDED",
            research_state="OPEN_REPLICATION",
            display_tone="GREEN",
            interpretation="A bounded positive result passed its declared gates and still requires replication.",
        )
    if status in _EXECUTION_FAILURES:
        return ResultSemantics(
            execution_health="FAILED",
            evidence_outcome="UNKNOWN",
            research_state="INVALID_RUN",
            display_tone="RED",
            interpretation="The run or evidence contract failed; no scientific conclusion is licensed.",
        )
    return ResultSemantics(
        execution_health="UNKNOWN",
        evidence_outcome="UNKNOWN",
        research_state="UNCLASSIFIED",
        display_tone="GRAY",
        interpretation="The status has not been assigned an evidence-semantic classification.",
    )
