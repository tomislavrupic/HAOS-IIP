# RP-01 Functional Target

Status: frozen before calibration and validation.

## Function

The function is preservation of a four-component low-frequency input/output
response:

```text
F(x) = [<x, cos(2 pi x)>, <x, sin(2 pi x)>,
        <x, cos(2 pi y)>, <x, sin(2 pi y)>]
```

This is the same independent probe family used to expose missing functional
restoration in RT-02. It is not parity satisfaction, syndrome magnitude,
symbol agreement, state distance, or the experiment label.

Functional recovery is:

```text
1 - ||F(x_final)-F(x_pre)|| / ||F(x_post)-F(x_pre)||
```

with zero assigned when the perturbation does not measurably alter the function.

## Frozen Gates

- calibration-derived functional-recovery threshold, lower bounded by `0.75`;
- target functional recovery rate must exceed passive and frozen RT-02;
- paired run-level confidence interval versus both baselines must exclude zero;
- function must recover in more than one perturbation family;
- recovery must be stronger inside the declared parity radius than above it.

The function matters because Rung 3 requires return of behavior, not only
repair of discrete consistency.
