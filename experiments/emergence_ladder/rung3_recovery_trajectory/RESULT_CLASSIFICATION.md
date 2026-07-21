# EL-R3-RT-01 Result Classification

Classification: `NEGATIVE_RESULT`

## Decisive Values

- Mean branch recovery gain: `-0.046654`.
- Branch recovery-gain 95% CI: `[-0.049307, -0.044018]`.
- Mean passive recovery gain: `0.000000`.
- Mean operator-only recovery gain: `-0.030471`.
- Mean topology-altered recovery gain: `-0.046528`.
- Mean grade-signature error: `0.360626`.
- Recovered-run fraction: `0.000000`.

The branch fails every official recovery gate except rejection of the trivial
attractor. The trivial control creates apparent convergence but has only
`0.187132` identity with the branch reference, so it is correctly rejected as
recovery to the wrong state.

## Frozen Interpretation

The tested grade-locked cochain transport supports persistence diagnostics in
earlier phases but does not recover its unperturbed trajectory after the
predeclared orthogonal disruptions. No parameter, threshold, or control will be
tuned in EL-R3-RT-01.

Rung 3 requires a new recovery mechanism, not another score over this one.
