# Phase XV Summary

Phase XV tests whether the frozen dilute collective sector supports bounded propagation-level structure under deterministic disturbances.

## Frozen Inputs

- Candidate label: `low_mode_localized_wavepacket`
- Refinement levels: `[48, 60, 72, 84]`
- Ensemble sizes: `[3, 5, 7]`
- Layout seeds: `[1301, 1302, 1303]`
- Tau grid: `[0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 3.6, 4.2]`

## Aggregate Metrics

- Branch effective-speed drift: `0.016755292342`
- Control effective-speed drift: `0.041888230853`
- Branch mean effective speed: `0.041363141386`
- Control mean effective speed: `0.046222397796`
- Branch influence-order mean: `0.848677248677`
- Control influence-order mean: `0.848677248677`
- Branch low-k drift: `0.013983435184`
- Control low-k drift: `0.021355474473`
- Branch shell-latency mean span: `2.333333333333`
- Control shell-latency mean span: `2.266666666667`

## Probe Means

- Density pulse branch/control speed means: `0.041363141386` / `0.046222397796`
- Candidate removal branch/control speed means: `0.046192032634` / `0.053100775884`
- Bias onset branch/control speed means: `0.199043712902` / `0.199073838022`

## Success Flags

- Propagation drift bounded: `True`
- Effective speed band compact: `True`
- Influence ordering stable: `True`
- Low-k transport coherent: `True`
- Control hierarchy different: `True`

## Control Flags

- Control propagation drift bounded: `False`
- Control effective speed band compact: `False`
- Control influence ordering stable: `True`
- Control low-k transport coherent: `False`

Phase XV establishes effective propagation-structure feasibility for the frozen operator hierarchy.
