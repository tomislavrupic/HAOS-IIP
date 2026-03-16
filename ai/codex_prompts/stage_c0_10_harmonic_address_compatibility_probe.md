# Stage C0.10 -- Harmonic Address Compatibility Probe

HAOS / IIP Pre-Metric Combinatorial Program

## Objective

Test whether a purely relational harmonic-address labeling rule on the combinatorial kernel graph produces:

- topology-selective persistence effects
- exchange-texture bias
- or protection / destabilization of previously observed braid-like flow families

without introducing:

- metric distance
- geometric embedding
- explicit new field dynamics
- external stabilizer terms

This stage asks whether address compatibility alone can bias encounter phenomenology.

## Conceptual rule

Each packet component carries a discrete harmonic address label:

`a in Z_N`

Address compatibility modulates the existing static combinatorial kernel contribution:

- exact match: `W_addr = 1`
- near match: `W_addr = eta_match`
- mismatch: `W_addr = eta_mismatch`

so that:

`K_eff(i, j) = W_addr(i, j) * K_base(i, j)`

No coordinate distance, Gaussian falloff, or Euclidean embedding distance is added.

## Experimental questions

1. Does harmonic-address matching increase local persistence?
2. Does mismatch destroy braid-like exchange topology?
3. Does near-match create intermediate exchange classes?
4. Is there evidence of address-selective encounter sorting?
5. Does address symmetry breaking open new flow channels?

## Baseline configuration

Use the clustered Dirac-Kahler mixed-encounter family from Phase III as seed:

- tight clustered separation baseline
- previous clustered baseline width
- projected-transverse state sector
- static kernel
- no memory
- no curvature trigger

## Minimal run design

Run a 9-case matrix on one fixed clustered seed:

- Group A: fully compatible address assignments
- Group B: near-compatible address assignments
- Group C: incompatible assignments plus randomized-address control

Include:

- one symmetry-broken address distribution
- one randomized address-field control

## Required observables

Per run compute:

- persistence_gain
- topology_class
- braid_exchange_indicator
- transfer_smear_index
- flow_concentration_index
- encounter_sorting_score
- address_selectivity_index

Derived classifications:

- address_protected_braid
- address_fragile_braid
- address_induced_smear
- address_sorted_encounter
- no_address_effect

## Numerical discipline

- keep graph resolution fixed
- keep timestep fixed across the 9-run matrix
- keep packet initialization envelope fixed
- keep the Dirac-Kahler update rule fixed
- vary only the harmonic-address configuration

## Required outputs

Produce:

- stamped JSON
- stamped CSV
- classification table
- topology vs address heatmap
- persistence distribution panel
- short interpretive note

Target note:

`Stage_C0_10_Harmonic_Address_Compatibility_Probe_v1.md`
