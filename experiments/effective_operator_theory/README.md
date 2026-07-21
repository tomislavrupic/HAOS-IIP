# Effective Operator Theory Scaffold

Status: synthetic diagnostic scaffold.

Lifecycle note: EOT-01 is retained as a supporting calibration.
`EOT-02-ARTIFACT-DERIVED` is deferred until lower recovery and closure rungs are
supported and enough independent refinements exist for a real untouched
holdout. Any successor must use unseen refinements and predeclared conventional
baselines; it may not select terms after holdout inspection. See the
[branch lifecycle registry](../../docs/branch_governance/branch_lifecycle_summary.md).

This bundle imports a restrained lesson from effective field theory: declare a
scale, allowed terms, cutoffs, controls, and claim ceilings before interpreting
coarse dynamics. It does not derive EFT quantum gravity, physical field theory,
spacetime, Lorentz invariance, matter content, constants, empirical physics, or
ontology.

The first benchmark is deliberately boring. It uses a deterministic one
dimensional ring graph Laplacian with a diffusion-like leading term and a
suppressed higher-order correction:

```text
rate(lambda) ~= c1 * lambda + c2 * lambda^2
```

The benchmark asks whether the generated synthetic mode decay rates recover
that hierarchy under frozen rules and whether controls detect missing coupling
or unstable sign structure.

## Run

```bash
uv run python experiments/effective_operator_theory/run_effective_operator_expansion.py
```

## Outputs

- `results/precommitment_contract.json`
- `results/allowed_terms.json`
- `results/coefficient_sweep.csv`
- `results/control_results.csv`
- `results/effective_operator_report.md`
- `results/effective_operator_result.json`

Current frozen result hash:

```text
effective_operator_a740ce933dd1dbd479931847
```

## Interpretation Boundary

A `PASS` means only:

- the synthetic operator hierarchy recovered its predeclared effective terms;
- the leading Laplacian-like term was detected;
- the higher-order correction remained suppressed;
- the long-wavelength residual stayed below the declared tolerance;
- the identity and unstable-sign controls behaved as expected.

It does not mean HAOS-IIP has derived a physical field theory or quantum
gravity. This is a methodological bridge: a way to make future coarse-grained
operator claims more disciplined, smaller, and easier to falsify.
