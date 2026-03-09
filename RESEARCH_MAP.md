# RESEARCH_MAP

## Core pipeline

```text
one substrate
    ->
one kernel
    ->
weighted discrete geometry
    ->
three operator sectors
```

## Operator hierarchy

```text
node lift   -> L0  -> scalar / geometry sector
edge lift   -> L1  -> gauge / vector sector
Dirac lift  -> D_H -> fermion / spinor sector
```

## Repository navigation by stage

### 1. Interaction Kernel

Core files:

- [docs/canon/KERNEL_DEFINITION.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/KERNEL_DEFINITION.md)
- [theory/kernel/HAOS_IIP_Kernel_v1.docx](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/theory/kernel/HAOS_IIP_Kernel_v1.docx)
- [theory/kernel/HAOS_Spectral_Recoverability_Principle_v1.docx](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/theory/kernel/HAOS_Spectral_Recoverability_Principle_v1.docx)

Question:

- what interaction rule is primary?

### 2. Operator hierarchy

Core files:

- [docs/canon/OPERATOR_STRUCTURE.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/OPERATOR_STRUCTURE.md)
- [docs/canon/HAOS_IIP_CORE_THEORY.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/HAOS_IIP_CORE_THEORY.md)
- [theory/operators/Three_Operator_Architecture_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/theory/operators/Three_Operator_Architecture_v1.md)

Question:

- what lifts of the same weighted geometry are actually being used?

### 3. Emergent geometry

Core files:

- [docs/canon/EMERGENT_GEOMETRY.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/EMERGENT_GEOMETRY.md)
- [theory/geometry/HAOS_IIP_Spectral_Geometry_Note_v1.docx](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/theory/geometry/HAOS_IIP_Spectral_Geometry_Note_v1.docx)
- [theory/geometry/HAOS_IIP_Spectral_Gravity_Note_v1.docx](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/theory/geometry/HAOS_IIP_Spectral_Gravity_Note_v1.docx)

Question:

- when does `L0` behave like a scalar geometric operator?

### 4. Gauge program

Core files:

- [docs/canon/GAUGE_PROGRAM.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/GAUGE_PROGRAM.md)
- [theory/gauge/HAOS_IIP_Emergent_Gauge_Sector_Note_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/theory/gauge/HAOS_IIP_Emergent_Gauge_Sector_Note_v1.md)
- [experiments/gauge_tests/README.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/experiments/gauge_tests/README.md)

Question:

- can an edge/Hodge branch `L1` produce a genuine gauge/vector sector?

### 5. Particle program

Core files:

- [theory/particles/Rebuilding_Dudas_Picture_in_HAOS_IIP_v1.docx](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/theory/particles/Rebuilding_Dudas_Picture_in_HAOS_IIP_v1.docx)
- [docs/canon/OPEN_PROBLEMS.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/canon/OPEN_PROBLEMS.md)

Question:

- what Dirac-type or chiral operator structure is required before particle-like claims are legitimate?

## Experiment map

### Current explicit experiment

- [experiments/eigenmodes/haos_iip_3d_low_mode_study/HAOS_IIP_3D_Low_Mode_Study_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/experiments/eigenmodes/haos_iip_3d_low_mode_study/HAOS_IIP_3D_Low_Mode_Study_v1.md)

Direct result:

- scalar geometry-like low modes of `L0` are visible
- connection sensitivity appears under phase dressing
- no credible vector or fermion sector is present yet

### Next experiment slots

- `experiments/gauge_tests/`: weighted edge/Hodge, holonomy, and transverse-mode experiments
- `experiments/eigenmodes/`: defect and Dirac-type low-mode studies
- `papers/drafts/`: paper-ready synthesis after experiment stabilization
