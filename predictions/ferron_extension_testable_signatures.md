# HAOS-IIP Ferron Extension - Testable Signatures

This note converts Materials Bridge Line A into forward-facing, falsifiable
prediction records. It is not a new theory, not a proof of HAOS-IIP, and not a
claim that ferrons embody HAOS-IIP.

## Calibration Source

The calibration point is release `64.1`, the bounded ferron / NbOI2 materials
bridge over Choe et al. source data:

- real Figshare XLSX data loaded and hashed
- `9` time traces and `6` spectra or peak records parsed
- `15/15` target-window spectral records found near `3.13 THz`
- mean peak frequency: `3.13893 THz`
- mean spectral recoverability: `0.928455`
- published Fig. 2d target-band traces: `8` traces, `216` time points
- published group-velocity proxies: `58.8 km/s` at `96 nm`, `100 km/s` at `224 nm`
- computed-from-time-trace STFT envelope correlation vs. published: `0.995318`
- unified persistence proxy: `0.957192`
- raw STFT/time-frequency grid status: `NO_STFT_DATA_FOUND`

The raw-grid boundary remains active. Published target-band traces and computed
STFT diagnostics are not raw time-frequency maps.

## Prediction 1 - Polar-Mode Address Selection

In layered polar ferroelectrics or closely related polar oxyhalides with an
independently measured polar transverse optical mode, pulsed excitation should
produce a narrowband coherent readout concentrated in a pre-registered target
window around that mode. Non-polar or mode-mismatched controls should not show
the same target-window recoverability.

Suggested pass criteria:

- target-window peak found in at least `80%` of controlled records
- linewidth proxy below `0.1 THz` when spectra support linewidth estimation
- mean spectral recoverability at or above `0.85`
- non-polar or mode-mismatched controls fail at least one of those criteria

Falsifier:

If matched polar materials repeatedly show only broadband emission, or if
non-polar controls show the same target-window recoverability under matched
excitation, this prediction fails for that material class.

## Prediction 2 - Thickness / Boundary Propagation Dependence

In NbOI2-like or WO2Br2-like flakes, the target-band packet arrival time should
increase monotonically with pump-probe separation, and the group-velocity proxy
should change with thickness or confinement under matched excitation.

Suggested pass criteria:

- peak time versus distance fit has positive slope
- monotonic peak-delay fraction is at least `0.75`
- group velocity differs between at least two controlled thickness groups by
  more than the experimental uncertainty
- target-band recoverability remains above `0.70` before visible propagation
  failure

Falsifier:

If controlled thickness series show no distance-ordered peak delay and no
thickness dependence beyond uncertainty, while the target-band signal remains
measurable, this prediction fails for that material class.

## Prediction 3 - Re-Addressing After Secondary Perturbation

After a secondary pulse, controlled defect scattering, or mild boundary
perturbation below damage threshold, the coherent readout should preferentially
return to the original target frequency window before the full amplitude
envelope necessarily recovers.

Suggested pass criteria:

- post-perturbation peak frequency returns to the pre-registered target window
- frequency stability is at least `0.95`
- mean target-window spectral recoverability is at least `0.85`
- broadband background does not dominate the recovered signal

Falsifier:

If secondary perturbation repeatedly produces broadband thermalization or a
different dominant frequency with no preferential recovery of the original
target window, this prediction fails for that perturbation regime.

## Prediction 4 - Cross-Dataset Narrowband Control Discrimination

Across external datasets with pre-registered narrowband collective modes,
HAOS-style recoverability metrics should separate coherent directional
propagation records from broadband or non-propagating controls better than
amplitude-only ranking.

Suggested pass criteria:

- target-window recoverability separates coherent propagation records from
  controls in a held-out dataset
- the separation survives fixed thresholds chosen before inspection
- amplitude-only ranking performs worse or is less stable under matched controls

Falsifier:

If pre-registered target-window recoverability does not distinguish coherent
propagation records from broadband or non-propagating controls better than
amplitude-only baselines, this prediction fails for the selected dataset class.

## Boundary

These predictions are field-facing diagnostic candidates. They do not claim new
physics, a materials ontology, a universal harmonic grid, or a proof of
HAOS-IIP. The numbers above are provisional thresholds seeded by the ferron
audit, not universal constants.
