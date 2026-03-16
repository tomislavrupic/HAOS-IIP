# Stage 23.11 -- DK Minimal Effective Model Extraction

## Role

This stage closes the descriptive braid-atlas pass by extracting a compact effective surrogate for the clustered Dirac-Kahler braid / smear sector.

No new dynamical ingredient is allowed.
This is a reduction step, not a discovery scan.

## Context

Phase III established on the frozen HAOS-IIP architecture that:

- clustered DK encounters support a reproducible braid-like exchange window
- a finite smeared midpoint band exists in phase space
- no persistence or metastable transport regime opens
- non-clustered families do not recover the braid sector as a general mechanism

The correct next move is therefore:

compress the clustered braid texture into the smallest workable effective state model.

## Objective

Construct a two-state discrete surrogate that:

1. uses only reduced observables extracted from the frozen clustered DK branch
2. reproduces the braid / smear / unresolved class structure on a small anchor lattice
3. preserves the coarse detuning-threshold structure
4. can be drawn as a low-parameter phase-plane system

## Allowed ingredients

Use only:

- frozen Dirac-Kahler propagator structure
- clustered packet seed template
- bounded clustered-seed grade detuning inside the existing seed definition
- scalar observables already logged or derivable from the same runs:
  - signed separation envelope
  - grade-exchange amplitude
  - flow concentration index
  - topology class label

Do not add:

- memory
- kernel adaptation
- harmonic-address forcing
- new fields
- long-range coupling
- persistence logic

## Effective state

Use:

- `X_n`: signed separation envelope
- `Y_n`: net grade-exchange amplitude

Fit a discrete map of the form

`X_{n+1} = X_n + f(Y_n, phase_offset)`

`Y_{n+1} = Y_n + g(X_n, detuning_delta)`

with:

- only constant, linear, and quadratic terms
- no dependence of `f` on detuning
- no dependence of `g` on phase

## Calibration lattice

The source brief listed more axis values than a 9-run matrix can contain.
To keep the stage finite and aligned with the already resolved clustered braid window, use the 3 x 3 anchor lattice:

- `phase_offset_fraction_of_pi in {0.375, 0.500, 0.625}`
- `detuning_delta in {0.0, 0.5, 0.75}`

Interpret `detuning_delta` here as bounded initial grade-balance detuning of the clustered DK seed, not as a new operator term.

## Per-run outputs

For each of the 9 runs produce:

- reduced phase-plane trajectory plot
- observed vs predicted topology label
- scalar error diagnostics

## Summary outputs

Produce:

- effective phase diagram
- fitted parameter table
- failure map of surrogate breakdown
- short interpretation note

Target note:

`Stage_III_C2_Minimal_Effective_Model_Extraction_v1.md`

## Validation rule

Set `effective_model_valid = TRUE` only if:

- topology agreement fraction >= 0.75
- detuning-threshold error <= 0.1
- the predicted phase structure contains a single smeared band on the low-detuning slice

## Interpretation rule

If valid:

- state that the clustered DK braid texture admits a compact effective surrogate on the tested anchor lattice
- do not upgrade this into a persistence law or a family-general mechanism

If invalid:

- state that the braid / smear sector remains intrinsically higher-dimensional at this resolution
- motivate a later sector change rather than more surrogate forcing
