# API Contract

Status: contract only  
Scope: artifact-driven operator, initialization, control, and telemetry API for the frozen HAOS-IIP stack  
Authority: this document does not change any earlier frozen phase; it defines the target callable surface for later implementation work

## Purpose

This contract freezes one thin implementation surface for later runners and probes.

The design rule is:

- artifact-driven
- minimal
- deterministic
- no phase-source imports
- no operator-class drift

The contract is anchored to the already frozen scientific layers:

- operator definition: [phase6-operator/phase6_operator_manifest.json](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/phase6-operator/phase6_operator_manifest.json)
- localized branch candidate: [phaseX-proto-particle/phaseX_integrated_manifest.json](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/phaseX-proto-particle/phaseX_integrated_manifest.json)
- protection and persistence logic: [phase11-protection/phase11_manifest.json](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/phase11-protection/phase11_manifest.json)
- interaction and collective sector continuity: [phase12-interactions/phase12_manifest.json](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/phase12-interactions/phase12_manifest.json), [phase13-sector-formation/phase13_manifest.json](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/phase13-sector-formation/phase13_manifest.json), [phase14-collective-dynamics/phase14_manifest.json](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/phase14-collective-dynamics/phase14_manifest.json)
- ordering and causal telemetry: [phase16-temporal-ordering/phase16_manifest.json](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/phase16-temporal-ordering/phase16_manifest.json), [phase17-causal-closure/phase17_manifest.json](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/phase17-causal-closure/phase17_manifest.json)

## Module Layout

Target files:

- `operators/cochain_laplacian.py`
- `initialization/frozen_branch.py`
- `controls/altered_connectivity.py`
- `telemetry/metrics.py`

These paths are part of the contract.

## Shared Type Assumptions

The signatures below intentionally avoid implementation-specific sparse backends.

- `"SparseMatrix"` means the project’s chosen sparse matrix type.
- `np.ndarray` means a dense state vector or dense observable history.
- Metadata dictionaries must stay JSON-serializable.

Suggested imports:

```python
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
```

## 1. Operator Layer

File path:

- `operators/cochain_laplacian.py`

Frozen object:

- full block cochain Laplacian `delta_h`
- periodic DK2D complex
- `h = 1 / n_side`
- fixed `epsilon = 0.2`
- deterministic level family across refinement

This is the single authoritative operator-definition layer.

Exact signatures:

```python
from dataclasses import dataclass


@dataclass(frozen=True)
class OperatorSpec:
    n_side: int
    epsilon: float = 0.2
    boundary: str = "periodic"
    operator_name: str = "delta_h"
    hierarchy_h: float = 0.0


def build_delta_h(spec: OperatorSpec) -> tuple["SparseMatrix", dict]:
    """
    Returns:
        delta_h: frozen branch-local periodic DK2D block cochain Laplacian
        meta: {
            "n_side": ...,
            "h": ...,
            "dimension": ...,
            "nnz": ...,
            "nullspace_estimate": ...,
            "spectral_radius_sample": ...,
        }
    """
```

Contract rules:

- exactly one constructor
- exactly one operator class
- no alternate normalization
- no boundary-family expansion inside this file
- `meta["h"]` must equal `1.0 / spec.n_side`

## 2. Frozen Branch Initialization

File path:

- `initialization/frozen_branch.py`

Frozen default branch family:

- `low_mode_localized_wavepacket`

This is the validated default because the integrated proto-particle authority path kept this as the only aggregate-persistent branch candidate across the frozen repeat levels.

Exact signatures:

```python
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class BranchInitSpec:
    n_side: int
    family: str = "low_mode_localized_wavepacket"
    seed: int = 1303
    amplitude_norm: float = 1.0
    placement: str = "single_centered"  # or "dilute_layout"
    min_physical_separation: float | None = None
    population_size: int | None = None


def build_frozen_branch_state(spec: BranchInitSpec, delta_h) -> tuple[np.ndarray, dict]:
    """
    Returns:
        psi0: deterministic frozen branch state
        meta: {
            "family": ...,
            "seed": ...,
            "population_size": ...,
            "layout_rule": ...,
            "candidate_centers": ...,
        }
    """
```

Contract rules:

- default family must remain `low_mode_localized_wavepacket`
- deterministic seed handling only
- `placement="single_centered"` must produce one validated centered branch state
- `placement="dilute_layout"` must reuse already frozen deterministic layout families rather than inventing new placement logic
- this file must not redefine persistence criteria

## 3. Control Construction

File path:

- `controls/altered_connectivity.py`

Control exposure rule:

- one public interface
- named variants allowed underneath that interface
- same node counts as the matched branch level

Named control families preserved by contract:

- `periodic_diagonal_augmented_control`
- `open_grid_scalar_block`

Exact signatures:

```python
from dataclasses import dataclass


@dataclass(frozen=True)
class ControlSpec:
    n_side: int
    control_name: str = "altered_connectivity"
    variant: str = "periodic_diagonal_augmented_control"  # or "open_grid_scalar_block"
    seed: int = 1303
    epsilon: float = 0.2


def build_control_operator(spec: ControlSpec) -> tuple["SparseMatrix", dict]:
    """
    Returns:
        control_delta_h: deterministic matched control operator
        meta: {
            "control_name": ...,
            "variant": ...,
            "node_count_match": True,
            "connectivity_texture_changed": True,
        }
    """
```

Contract rules:

- control construction must stay parallel to the branch operator layer
- same `n_side` means same matched total node count
- no random topology generation outside deterministic seed policy
- this file must expose controls as contrasts, not as new primary operators

## 4. Telemetry Layer

File path:

- `telemetry/metrics.py`

The telemetry contract has three frozen levels:

- local mode metrics
- survival and persistence metrics
- temporal-ordering and causal-graph metrics

### 4.1 Local Mode Metrics

Exact signatures:

```python
import numpy as np


def overlap(reference: np.ndarray, state: np.ndarray) -> float:
    """Normalized absolute inner-product overlap."""


def localization_width(state: np.ndarray, coords: np.ndarray) -> float:
    """Second-moment width of |state|^2 around its center."""


def concentration_retention(state: np.ndarray, mask: np.ndarray) -> float:
    """Mass retained inside frozen local support region."""


def participation_ratio(state: np.ndarray) -> float:
    """Inverse participation ratio style compactness diagnostic."""


def recovery_score(
    reference: np.ndarray,
    state: np.ndarray,
    coords: np.ndarray,
    mask: np.ndarray,
) -> float:
    """
    One frozen structural score.
    Keep it deterministic and composed only from already-frozen ingredients.
    """
```

Contract rules:

- no competing recovery scores
- no plotting logic
- no threshold logic inside these primitive metric functions

### 4.2 Survival and Persistence

Exact signatures:

```python
from dataclasses import dataclass


@dataclass(frozen=True)
class SurvivalThresholds:
    max_width_growth: float
    min_concentration: float
    max_participation_growth: float
    min_overlap: float
    min_recovery_score: float


def classify_single_mode(
    reference,
    state,
    coords,
    mask,
    thresholds: SurvivalThresholds,
) -> str:
    """Returns: persistent | diffusive | unstable"""


def persistence_time(
    history: list[dict],
    tau_grid: list[float],
    thresholds: SurvivalThresholds,
) -> float:
    """Largest tau for which the state remains inside the frozen gate."""
```

Contract rules:

- classification vocabulary is fixed to `persistent | diffusive | unstable` at this layer
- higher-level regime classes such as merger or decoherence remain phase-runner logic, not primitive telemetry
- `tau_grid` must be interpreted against the frozen grid already established in the later phase stack

### 4.3 Temporal Ordering and Causal Graph

Exact signatures:

```python
import numpy as np


def first_threshold_crossing(signal: np.ndarray, threshold: float) -> int | None:
    """First index where observable crosses frozen threshold."""


def front_arrival_order(
    probe_histories: dict[str, np.ndarray],
    threshold: float,
) -> dict[str, int | None]:
    """Arrival ordering map for a disturbance family."""


def reconstruct_influence_edges(
    event_times: dict[str, int | None],
) -> list[tuple[str, str]]:
    """Directed edges from source to responding probes under frozen policy."""


def acyclicity_score(
    edges: list[tuple[str, str]],
    nodes: list[str],
) -> float:
    """1.0 for DAG-like, lower when cycles appear."""


def causal_depths(
    edges: list[tuple[str, str]],
    source: str,
) -> dict[str, int]:
    """Shortest directed path depth from source."""


def order_compatibility(edge_set, arrival_ordering) -> float:
    """Mismatch rate between inferred edges and front-arrival order."""
```

Contract rules:

- these functions must remain sufficient for the frozen Phase XVI and Phase XVII telemetry stack
- they must not trigger new simulations
- they must operate on observable histories or reconstructed event maps only

## Implementation Order

If implemented later, the order is part of the intended discipline:

1. freeze `build_delta_h`
2. freeze `build_frozen_branch_state`
3. freeze one matched `build_control_operator`
4. freeze the telemetry layer
5. only then write runners

## Explicit Non-Goals

This contract does not:

- introduce new operators
- redefine frozen persistence thresholds
- redefine frozen trace or coefficient contracts
- add new causal or geometric ontology
- authorize broad repo refactors

If coherence cannot be recovered, shrink the system.
