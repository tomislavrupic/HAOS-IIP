# Stage 23.4 - Dirac-Kahler Metastable Transport Window Scan

HAOS (Harmonic Address Operating System) · IIP (Interaction Invariance Physics)
Pre-Geometry Atlas · Phase III continuation

## Purpose

Stage 23.3 established:
- clustered DK encounters form a long-lived mixed regime
- coherence can increase without capture
- minimal stabilizers do not produce bounded composites

Stage 23.4 tests a different hypothesis:

The clustered DK family may define a metastable transport window
where grade-exchange enables structured energy redistribution
without persistence promotion.

The goal is classification, not promotion.

## Frozen architecture constraints

Do not modify:
- kernel form
- graph topology
- operator discretization
- damping / drive baseline
- persistence classifier thresholds

Allowed variation axes only:
1. packet separation scale
2. packet width ratio
3. initial phase offset
4. mild grade-coupling modulation
5. weak anisotropic kernel stretch

No new operators.
No new fields.

## Observable targets

For each run compute:
- persistence_gain
- collision_class
- grade_exchange_amplitude
- grade_exchange_coherence
- transport_span
- post-collision energy asymmetry
- recurrence_indicator

Derived labels:
- localized encounter
- metastable transport
- dispersive washout

## Minimal run family (concept)

Three structured triplets:

### Triplet A - separation scaling
- tight clustered baseline
- medium clustered
- near-overlapping seed

### Triplet B - width contrast
- symmetric widths
- probe-A wider
- probe-B wider

### Triplet C - phase-transport tuning
- in-phase
- quarter-phase offset
- half-phase offset

Total target: 9 runs

## Classification logic

After simulation:

If
- transport_span increases
- grade_exchange_coherence stays high
- persistence_gain remains weak

-> classify as metastable transport

If
- recurrence_indicator rises
- oscillatory separation emerges

-> candidate proto-bound composite

If
- transport_span collapses
- coherence drops

-> dispersive decay

## Output protocol

Generate:
1. JSON runsheet log
2. CSV metric table
3. transport_span vs time panel
4. grade_coherence heatmap
5. short snapshot note

Snapshot note must answer:
- Does DK define a reproducible transport regime?
- Is transport sensitive to phase geometry or width contrast?
- Is there any signal of delayed persistence promotion?

## Freeze criteria

Stage 23.4 freezes as:
- positive exploratory layer
  if a stable metastable-transport classification family appears
- negative boundary
  if all clustered encounters remain unresolved mixed events

No persistence claims allowed.

## Commit tag suggestion

Freeze Stage 23.4 DK metastable transport scan

## Continuation hint

If transport turns out real, a natural continuation is:
- Stage 23.5 as energy-flow topology mapping
