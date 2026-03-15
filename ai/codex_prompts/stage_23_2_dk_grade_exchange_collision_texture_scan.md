# Stage 23.2 - Dirac-Kaehler Grade-Exchange and Collision-Texture Scan

## Role
Stage 23.1 showed that the Dirac-Kaehler sector changes collision texture through widespread grade exchange, but it did not establish a persistence-level advantage over the matched scalar control.

Stage 23.2 therefore narrows the question:

Can the DK sector support a reproducible family of grade-exchange / collision-texture behaviors that remains structured under refinement and phase variation, even if persistence is still weak?

## Discipline
Keep fixed:
- same periodic regular lattice family
- same kernel baseline
- same packet bandwidth and amplitude class
- no memory, reinforcement, or gating
- no extra fields

Change only:
- collision-family focus
- inclusion of refinement in the strongest geometry families

## Geometries
Use only the Phase III families that were most informative in Stage 23.1:
1. counter-propagating corridor pair
2. offset glancing collision
3. tight clustered pair

## Conditions
For each geometry:
- in phase
- out of phase

At each resolution:
- n = 12
- n = 24

Total minimal scan:
- 3 geometries x 2 phase relations x 2 resolutions = 12 DK runs

## Primary observables
- collision label
- persistence label
- composite lifetime
- binding persistence
- encounter dwell time
- grade-transfer amplitude
- grade-transfer onset time
- omega1 peak weight
- omega2 peak weight
- resolution stability of collision texture

## What counts as a positive Stage 23.2 signal
Not binding.

A positive 23.2 outcome means:
- grade-exchange strength is geometry-selective
- phase relation changes collision texture reproducibly
- at least one texture family remains stable under n = 12 -> 24 refinement
- DK requires new collision-language beyond scalar pass-through/dispersion

## Negative freeze rule
Freeze 23.2 negative if:
- grade exchange is ubiquitous but structureless
- phase sensitivity is weak or inconsistent
- refinement washes out the apparent collision families
- no geometry-specific texture family survives comparison

## Interpretation boundary
Even a positive 23.2 result implies only:
- the DK sector changes collision phenomenology in a structured way

It does not imply:
- fermionic matter
- exclusion law
- binding
- emergent geometry

## Deliverables
- `stage23_2_dk_grade_exchange_runs.json`
- `stage23_2_dk_grade_exchange_texture_scan.py`
- `Stage_23_2_DK_Grade_Exchange_Texture_Scan_v1.md`
- stamped JSON / CSV outputs
- grade-transfer summary plots
- texture stability summary plots
