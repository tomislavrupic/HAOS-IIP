# HAOS-IIP Prediction Protocol

## Purpose

The prediction workbench converts HAOS-IIP ideas into bounded, falsifiable records.

The central format is:

> In system X, under perturbation P, recoverability observable R should degrade before visible or conventional failure marker F.

This is a diagnostic prediction format, not an ontology claim.

## Admission Gate

A prediction may enter the ledger only if it has:

- a concrete system scope
- a controlled perturbation
- measurable observables
- explicit pass criteria
- explicit fail criteria
- at least one falsifier
- scope limitations
- non-claims

If a statement cannot be perturbed, measured, failed, or revised, it belongs outside this ledger.

## Prediction Levels

- `P0_toy_internal`: Supported only by toy or in-repo demonstrations.
- `P1_cross_system`: Supported across more than one toy or internal system.
- `P2_external_dataset`: Ready for external data or third-party benchmark tests.
- `P3_field_facing`: A field-facing diagnostic prediction with clear controls.
- `P4_physics_candidate`: A candidate physics-facing prediction. This requires external replication before any strong language is allowed.

## Status Values

- `candidate_unvalidated`: Stated clearly, not yet supported.
- `toy_supported_not_external`: Supported by toy or internal artifacts only.
- `external_test_ready`: Has a specified external test plan.
- `externally_supported`: Supported by external data or independent replication.
- `failed`: Failed under stated criteria.
- `retired`: Removed from active use without being treated as proven false.

## Novelty Boundary

The workbench may say:

- "not currently tracked by standard observables"
- "field-facing diagnostic candidate"
- "recoverability precursor prediction"
- "external validation required"

The workbench must not say:

- "HAOS-IIP proves new physics"
- "HAOS-IIP explains biology"
- "HAOS-IIP replaces existing field theory"
- "failure to refute is confirmation"

## Required Falsification Form

Every record must include at least one falsifier in this form:

> If X is measured under control Y and recoverability does not move before marker Z, the prediction fails for that system class.

## Evidence Rule

Internal demos can support only `P0_toy_internal` or `P1_cross_system` status.

External claims require external datasets, independent controls, or third-party replication. Until then, the correct language is "candidate" or "toy-supported diagnostic."

## Shrink Rule

If a prediction fails, shrink scope before discarding the entire line:

1. Reduce the system class.
2. Tighten perturbation conditions.
3. Replace vague observables with direct measurements.
4. Move the record to `failed` if recoverability cannot be recovered under the revised scope.

