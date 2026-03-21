#!/usr/bin/env python3

from __future__ import annotations


REGIME_ORDER = {
    'Fragmenting': 0,
    'Chaotic or irregular': 1,
    'Diffusive': 2,
    'Localized': 3,
    'Oscillatory trapped': 3,
    'Ballistic dispersive': 4,
    'Metastable structured': 4,
    'Ballistic coherent': 5,
}


def classify_transition_type(
    baseline_regime: str,
    perturbed_regime: str,
    baseline_width_ratio: float,
    perturbed_width_ratio: float,
    coherence_delta: float,
) -> str:
    if perturbed_regime == baseline_regime:
        return 'stable'
    if perturbed_regime == 'Fragmenting':
        return 'fragmentation_induced'
    if perturbed_regime == 'Chaotic or irregular':
        return 'chaos_induced'

    baseline_rank = REGIME_ORDER.get(baseline_regime, 0)
    perturbed_rank = REGIME_ORDER.get(perturbed_regime, 0)
    width_growth = float(perturbed_width_ratio) - float(baseline_width_ratio)

    if perturbed_rank < baseline_rank or coherence_delta < -0.02 or width_growth > 0.08:
        return 'regime_softening'
    return 'regime_shift'


__all__ = ['classify_transition_type']
