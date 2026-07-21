# EL-R3-RP-01 Theory Interpretation

## What RP-01 Added Beyond RT-02

RT-02 stores independent one-bit edge orientations. RP-01 adds sparse parity
relations among three local order symbols, explicit syndrome detection, and a
finite one-error-per-block correction prediction. The memory is distributed:
32 stored bits, two at each of sixteen owners, with five bits visible to any
decoder. A float64 checkpoint would require 4096 bits.

## What Worked

- Every constructed one-symbol error is corrected by the local decoder.
- Two errors in a block can be shown to miscorrect, matching distance three.
- Distinct continuous states share the same relational symbols and parity.
- Within-radius validation cases recover the correct relational symbols.
- Within-radius median identity and operator restoration are `1.0`.
- The continuous bridge improves the relation-defined state region.
- Signal-blind correction with the same energy does not reproduce that repair.

Restorative information is stored in distributed parity bits. Syndrome
comparison detects relational inconsistency; local bit flipping estimates the
error; bounded inequality feedback performs the continuous correction.

## What Failed

The discrete relations do not determine the independent low-frequency
function. Inside the correction radius, median functional restoration is
`0.0` even when decoding, identity, state-region repair, and operator-region
repair succeed. The continuous bridge restores inequalities, not the
amplitudes and phase relations needed by the four probe responses.

Outside the radius, syndrome elimination is not reliable evidence of identity:
the decoder often reaches the wrong member of the parity coset. This is the
expected ambiguity of a local `[3,1,3]` block, not a hidden recovery basin.

Parity redundancy therefore supplies information missing from RT-02 only at
the discrete-consistency level. It does not supply the mesoscopic amplitude,
phase, or transition information required for function.

## Theory-Space Reduction

Three tested ideas are now separated:

1. grade locking supports persistence but not restoration;
2. independent orientation memory supports identity-preserving adaptation but
   not function;
3. local parity over order relations supports finite-radius relation repair but
   not functional organizational recovery.

The next honest action is theory revision, not more decoder tuning. A future
Rung 3 proposal must explain how function-relevant mesoscopic information is
stored without becoming a continuous checkpoint. Candidate directions must be
distinct at the information representation level, such as dynamical phase
memory, learned local generative constraints frozen on development systems, or
multi-timescale latent reconstruction. None is authorized here.

No operational closure, cross-scale emergence, universality, continuum, or
physical claim follows.

