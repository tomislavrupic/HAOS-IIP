# PB-02 Observation Map

Status: frozen draft

Target

- real-world power-grid recovery and cascading failure behavior
- public dataset candidate: PowerGraph

Prediction target

- post-perturbation recovery ranking
- secondary targets: recovery time, residual error, return probability,
  post-recovery similarity, and recovered/failed label where the dataset
  supports them

Feature families to report

- baseline only
- HAOS only
- baseline + HAOS
- ablated HAOS
- matched null

Frozen split policy

- development: tune conventional baselines and lock the HAOS feature bundle
- calibration: choose exactly one HAOS feature bundle and stop
- holdout: untouched scoring split; no feature changes after inspection
- replication: second external dataset only after PB-02 holdout is frozen

Holdout rule

- no new features after holdout inspection
- no threshold movement after holdout inspection
- no feature selection using holdout outcomes

Controls

- shuffled target labels
- topology-destroyed graph
- degree-preserving rewiring
- weight shuffling
- parameter-matched null
- perturbation-free baseline
- seed repeat
- intentional leakage positive control

Success criterion

- best baseline + frozen HAOS features improves holdout ranking over the best
  conventional baseline on the same frozen split, and the gain survives matched
  controls and replication

Failure criterion

- no stable holdout improvement
- control invalidation
- target leakage
- underpowered comparison

Boundary

This is an execution map, not a claim of empirical success.
