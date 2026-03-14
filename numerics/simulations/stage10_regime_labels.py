#!/usr/bin/env python3

from __future__ import annotations

from typing import Any


def classify_secondary_sector(operator_sector: str, max_leakage: float, max_constraint: float) -> str:
    if max_constraint > 1.0e-6:
        return 'Constraint-failing'
    if max_leakage > 1.0e-6:
        return 'Leakage-dominant'
    if operator_sector == 'scalar':
        return 'Scalar-dominant'
    if operator_sector == 'transverse':
        return 'Transverse-dominant'
    return 'Mixed-sector'


def classify_primary_regime(summary: dict[str, float], max_leakage: float, max_constraint: float) -> str:
    shift = float(summary['center_shift'])
    width_ratio = float(summary['width_ratio'])
    coherence_ratio = float(summary['coherence_ratio'])
    path_eff = float(summary['path_efficiency'])
    rel_energy_drift = abs(float(summary['relative_energy_drift']))

    if max_constraint > 1.0e-4 or max_leakage > 1.0e-4:
        return 'Chaotic or irregular'
    if coherence_ratio < 0.45 and width_ratio > 1.35:
        return 'Fragmenting'
    if shift < 0.01 and width_ratio < 1.06:
        return 'Localized'
    if shift < 0.02 and width_ratio >= 1.06:
        return 'Oscillatory trapped'
    if shift < 0.05 and width_ratio > 1.18:
        return 'Diffusive'
    if path_eff >= 0.9 and width_ratio <= 1.12 and coherence_ratio >= 0.65 and rel_energy_drift < 0.02:
        return 'Ballistic coherent'
    if shift >= 0.02 and path_eff >= 0.75 and width_ratio <= 1.35:
        return 'Ballistic dispersive'
    if coherence_ratio >= 0.65 and rel_energy_drift < 0.03:
        return 'Metastable structured'
    return 'Chaotic or irregular'


def proto_spacetime_score(summary: dict[str, float], max_leakage: float, max_constraint: float, reproducible_under_refinement: bool) -> tuple[int, dict[str, int]]:
    stable_propagation = int(abs(float(summary['relative_energy_drift'])) < 0.02 and float(summary['path_efficiency']) > 0.7)
    low_leakage = int(max_leakage < 1.0e-6)
    bounded_constraint = int(max_constraint < 1.0e-6)
    low_directional_distortion = int(float(summary['max_anisotropy']) < 1.75)
    refinement = int(bool(reproducible_under_refinement))
    components = {
        'stable_propagation': stable_propagation,
        'low_leakage': low_leakage,
        'bounded_constraint_growth': bounded_constraint,
        'low_directional_distortion': low_directional_distortion,
        'reproducible_under_refinement': refinement,
    }
    return sum(components.values()), components


def classify_run(operator_sector: str, summary: dict[str, float], max_leakage: float, max_constraint: float, reproducible_under_refinement: bool) -> dict[str, Any]:
    primary = classify_primary_regime(summary, max_leakage, max_constraint)
    secondary = classify_secondary_sector(operator_sector, max_leakage, max_constraint)
    score, components = proto_spacetime_score(summary, max_leakage, max_constraint, reproducible_under_refinement)
    return {
        'regime_label': primary,
        'sector_label': secondary,
        'proto_spacetime_score': score,
        'proto_spacetime_components': components,
    }


__all__ = ['classify_run']
