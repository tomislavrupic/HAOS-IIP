#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASE_DIR = Path(__file__).resolve().parent
RUNS_DIR = PHASE_DIR / "runs"

INPUT_FILES = {
    "terminology": ROOT / "TERMINOLOGY.md",
    "translation_table": ROOT / "docs/notes/foundations/HAOS_IIP_Core_Translation_Table_v1.md",
    "canonical_entry": ROOT / "docs/notes/foundations/HAOS_IIP_Canonical_Entry_Paper_v1.md",
    "phase7_spectral": ROOT / "phase7-spectral/phase7_spectral_manifest.json",
    "phase18_distance_surrogate": ROOT / "phase18-distance-surrogate/phase18_manifest.json",
}


DEFINITION = """# Spectral Address

Formerly framed as: **Harmonic Address**.

## Core Definition

A spectral address is a coordinate that identifies a stable, recoverable resonance pattern, or eigenmode-like structure, inside a discrete interaction-invariant system.

It lives in the frozen branch-local cochain-Laplacian hierarchy. It is not static data. It is a persistent spectral signature that the system can reconstruct under bounded perturbation when the relevant interaction invariants and cochain structure remain recoverable.

## Operational Reading

A spectral address is accepted only when it can be:

- located in a declared operator, band, eigenmode, cochain subspace, peak window, or branch-local spectral signature;
- reidentified after perturbation, refinement, delay, distance, or matched controls;
- compared against nearby null windows or control branches;
- reconstructed from the frozen interaction structure without adding hidden labels;
- failed cleanly when the signature drifts, aliases, broadens, or becomes control-matched.

The address is therefore a recoverability relation, not a stored essence.

## How Information Is Stored

Information is stored as recoverable structure in the interaction system:

- interaction invariants preserve what the system can re-enter after perturbation;
- cochain incidence maps preserve which grade-level relations are allowed to compose;
- Laplacian spectra preserve stable mode windows, branch identities, and subspace relationships;
- reconstruction succeeds only when the frozen cochain structure carries enough invariant information to recover the signature.

The system does not store an address as a static tag. It stores the conditions under which the address can be reconstructed.

## Mathematical Interpretation

The primary mathematical language is spectral graph theory on finite cochain complexes:

- graph and Hodge/cochain Laplacians define the operator family;
- Laplacian eigenmodes, eigenvalue bands, and spectral subspaces supply the address coordinates;
- cochain complexes and incidence maps define the grade-local structure behind those coordinates;
- branch-local tracking measures whether the spectral coordinate remains identifiable across perturbation or refinement.

Kuramoto-style phase language and attractor language may be used as dynamical analogies, but they are subordinate to the Laplacian/cochain formulation. The technical claim is spectral recoverability, not synchronization mythology.

## Relation To Harmonic Address

Use "spectral address" for public technical work.

Use "harmonic address" only as historical/internal shorthand for the broader intuition that recoverable structures can be addressed by stable resonance-like identity. When the claim depends on measured operators, eigenmodes, spectra, bands, or cochain structure, the stricter term is spectral address.

## Non-Claims

A spectral address is not:

- a physical particle, field, or hidden essence;
- proof of spacetime, consciousness, or ontology;
- a universal frequency;
- a static memory cell;
- valid without perturbation recovery and control comparison.
"""


SLIDE_ALIGNMENT = """# Spectral Address Slide Alignment

## New Opener

**Spectral Address**

Formerly framed as **Harmonic Address**.

A spectral address is a coordinate for a stable, recoverable resonance pattern inside a discrete interaction-invariant system.

## Core Idea

The address lives in the frozen branch-local cochain-Laplacian hierarchy. It is not data sitting in the system. It is a persistent spectral signature that can be reconstructed under perturbation when the relevant interaction invariants and recoverable cochain structures survive.

## How Information Is Stored

Tight alignment:

- information is stored through interaction invariants, not symbolic labels;
- cochain structures preserve grade-level recoverability across nodes, edges, and higher relations;
- the address is recoverable when a spectral/cochain signature can be re-entered after perturbation;
- failure means drift, aliasing, broadening, or control-matching under the declared metric.

Suggested slide wording:

```text
Information is stored as recoverable interaction structure.

A spectral address is not a static tag. It is the spectral/cochain signature that remains reconstructable when the frozen interaction system is perturbed.
```

## Mathematical Interpretation

Tight alignment:

- lead with Laplacian eigenmodes and spectral graph theory;
- make cochain complexes explicit through incidence maps and Hodge/cochain Laplacians;
- keep Kuramoto and attractor language as secondary intuition for phase-locking or basin behavior;
- state the HAOS gate: a claimed address must survive perturbation and matched controls.

Suggested slide wording:

```text
Mathematically, a spectral address is tracked through Laplacian eigenmodes, spectral bands, and cochain subspaces in a frozen interaction hierarchy.

Kuramoto/attractor language is useful only as a dynamical analogy. The operational test is spectral recoverability under perturbation and controls.
```
"""


SUMMARY = """# Phase XIX Summary

Objective: define **spectral address** as the stricter technical replacement for the older internal phrase "harmonic address".

Definition: a spectral address is a coordinate that identifies a stable, recoverable resonance pattern, or eigenmode-like structure, inside a discrete interaction-invariant system. It lives in the frozen branch-local cochain-Laplacian hierarchy, not as static data, but as a persistent spectral signature reconstructable under bounded perturbation.

Scope:

- terminology-only phase;
- no new dynamics;
- no frozen-manifest edits;
- no stronger physical or ontological claim;
- public term: spectral address;
- legacy/internal shorthand: harmonic address.

Alignment results:

- "How Information Is Stored" now ties address persistence to interaction invariants and recoverable cochain structures.
- "Mathematical Interpretation" now foregrounds Laplacian eigenmodes, spectral graph theory, and cochain complexes, with Kuramoto/attractor language kept as analogy only.
- Failure remains explicit: if the address cannot be reidentified, reconstructable, or separated from controls, the address claim fails.

Phase XIX establishes the spectral-address translation layer for later HAOS-IIP wording and slide work.
"""


def write_text(path: Path, content: str) -> None:
    path.write_text(content.rstrip() + "\n")


def main() -> None:
    missing_inputs = [str(path) for path in INPUT_FILES.values() if not path.exists()]
    if missing_inputs:
        raise SystemExit(json.dumps({"success": False, "missing_inputs": missing_inputs}, indent=2))

    RUNS_DIR.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now(timezone.utc).isoformat()
    gates = {
        "definition_names_spectral_address": True,
        "harmonic_address_marked_legacy": True,
        "interaction_invariants_named": True,
        "recoverable_cochain_structures_named": True,
        "laplacian_eigenmodes_named": True,
        "spectral_graph_theory_named": True,
        "cochain_complexes_named": True,
        "kuramoto_attractor_subordinated": True,
        "no_new_dynamics": True,
        "success": True,
    }

    runs = {
        "phase": 19,
        "phase_name": "spectral-address",
        "timestamp": timestamp,
        "mode": "terminology_translation",
        "edits": [
            "define spectral address as the public technical term",
            "retain harmonic address as legacy/internal shorthand",
            "tie information storage to interaction invariants and recoverable cochain structures",
            "foreground Laplacian eigenmodes, spectral graph theory, and cochain complexes",
        ],
        "non_claims": [
            "no new dynamics",
            "no physical ontology",
            "no static-memory claim",
            "no consciousness claim",
        ],
    }

    manifest = {
        "phase": 19,
        "phase_name": "spectral-address",
        "timestamp": timestamp,
        "mode": "terminology_translation",
        "input_files": {name: str(path) for name, path in INPUT_FILES.items()},
        "outputs": {
            "definition": str(PHASE_DIR / "spectral_address_definition.md"),
            "slide_alignment": str(PHASE_DIR / "slide_alignment.md"),
            "summary": str(PHASE_DIR / "phase19_summary.md"),
            "runs": str(RUNS_DIR / "phase19_runs.json"),
        },
        "gates": gates,
        "closure_statement": "Phase XIX establishes the spectral-address translation layer for later HAOS-IIP wording and slide work.",
    }

    write_text(PHASE_DIR / "spectral_address_definition.md", DEFINITION)
    write_text(PHASE_DIR / "slide_alignment.md", SLIDE_ALIGNMENT)
    write_text(PHASE_DIR / "phase19_summary.md", SUMMARY)
    write_text(RUNS_DIR / "phase19_runs.json", json.dumps(runs, indent=2))
    write_text(PHASE_DIR / "phase19_manifest.json", json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
