# Ferron Coherence Validation

## Status
PASS

## Data Provenance
- Source DOI: 10.6084/m9.figshare.31895293
- Source URL: https://doi.org/10.6084/m9.figshare.31895293
- Source data Fig.1.xlsx | sha256=05a38cc576816479f851c0f4dfbe6370bed3ca163043e2c386c546980b9e82f3 | source=https://ndownloader.figshare.com/files/63310870
- Source data Fig.2.xlsx | sha256=8cba21ae3cba830a63a71a997532d892669139e3ecc7ba73833b53227a81dc25 | source=https://ndownloader.figshare.com/files/63310873
- Source data Fig.4.xlsx | sha256=42ea381f87d9046fbacce64246324d7b33599e5289fb25430bbad2c74f7fa93e | source=https://ndownloader.figshare.com/files/63310876

## Parsed Signals
- Time traces: 9
- Frequency spectra or peak records: 6
- STFT amplitudes: 0

## HAOS Mapping
- carrier: ferroelectric polarization order
- perturbation: laser pulse, fluence, distance, temperature, sample thickness, boundary/disorder effects
- recoverable structure: coherent ferron / polarization-wave mode
- observable readout: narrow-band THz emission, 3.13 THz amplitude, transient reflectance oscillation, Fourier/STFT peak
- visible failure: loss of narrow-band signal, spectral broadening, amplitude collapse, incoherent background dominance, or propagation failure

## Results
- baseline_recoverability: 1
- k_star_level: None
- first_visible_failure_level: None
- early_detection: False
- confidence summary: {"max": 0.9999926836070802, "mean": 0.9798696276083381, "min": 0.8654037886340976}

## Spectral Feature Audit
- status: PARTIAL
- spectral records audited: 15
- target peak search window: 3.13 THz +/- 0.2 THz
- target peak records found: 15
- strongest / most stable peak frequency summary: mean 3.13893 THz
- mean spectral recoverability: 0.928455
- k_star remains absent: True
- visible failure detected: False
- Target-window peak detected in 15/15 audited spectral records near 3.13 THz.
- Detected peak-frequency range is 3.13236-3.192 THz (mean 3.13893 THz).
- Absolute drift from the 3.13 THz target ranges from 0.00236 to 0.062 THz.
- FWHM-like linewidth proxy is available for 14 records, with mean 0.0253178 THz.
- Peak-to-background ratio is available for 14 records, with mean 8194.98.
- Distance-series target peaks remain detectable; amplitude retention at maximum parsed distance is 1-1, while monotonic decay consistency is 0-0.333333.
- Spectral recoverability ranges from 0.742975 to 0.99705 across audited records.
- missing data notes:
  - 9 audited spectra are FFT-derived from parsed time traces.
  - Linewidth proxy unavailable for Source data Fig.4_Fig. 4c Im_r_pp_: invalid_linewidth_width
  - No STFT maps parsed in the current data load.
  - Peak-to-background unavailable for Source data Fig.4_Fig. 4c Im_r_pp_: no_background_points_outside_target_window
- bounded interpretation: Standout spectral feature: a target-window peak near 3.13 THz is detected in the parsed spectral records. Peak frequency remains stable under the available conditions, and no sustained k-star collapse is detected under the current proxy thresholds.

## Limitations
- external-data probe only
- no proof of HAOS-IIP
- no claim that ferrons embody HAOS-IIP
- no new materials theory
- no consciousness claim
- metric proxies are provisional
- results depend on available data quality
