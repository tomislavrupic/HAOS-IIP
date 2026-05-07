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

## STFT / Time-Frequency Audit
- status: NO_STFT_DATA_FOUND
- raw files inspected: 3
- candidate sheets inspected: 1
- usable STFT maps found: 0
- target frequency search window: 3.13 THz +/- 0.2 THz
- usable STFT data found: False
- No usable raw STFT/time-frequency map was found in the downloaded XLSX files. No STFT recoverability claim is made.
- rejected candidate sheets:
  - Source data Fig.2.xlsx#Fig.2d: long_form_frequency_column_missing; matrix_axes_not_numeric_enough
- missing data notes:
  - Source data Fig.2.xlsx#Fig.2d rejected: long_form_frequency_column_missing; matrix_axes_not_numeric_enough
- bounded interpretation: No raw STFT/time-frequency map was parsed. Existing spectral audit result remains unchanged.

## Published STFT Target-Band Intensity Trace Audit
- status: PASS
- mode: published_stft_target_band_intensity_trace
- raw STFT/time-frequency grid status remains: NO_STFT_DATA_FOUND
- traces found: 8
- time points audited: 216
- peak records: 8
- thickness groups: 2
- mean group velocity m/s: 79411.8
- mean velocity consistency with 1e5 m/s: 0.794118
- mean envelope recoverability: 0.625
- group velocity rows:
  - thickness=96 nm | slope=17 ps/um | velocity=58823.5 m/s | r_squared=0.955372
  - thickness=224 nm | slope=10 ps/um | velocity=100000 m/s | r_squared=0.987654
- bounded interpretation: Published target-band STFT intensity traces were parsed as post-processed 3.13 THz-band envelope data. Peak arrival time increases with distance in the parsed thickness groups, supporting a bounded propagation-readout audit. This does not alter the raw STFT-grid result: no raw time-frequency grid is claimed.

## Computed From Time-Trace STFT Diagnostic
- status: PASS
- mode: computed_from_time_trace_stft_diagnostic
- raw STFT/time-frequency grid status remains: NO_STFT_DATA_FOUND
- target frequency: 3.13 THz
- frequency tolerance: +/- 0.05 THz
- window / step: 10 ps / 5 ps
- traces computed: 8
- time points computed: 216
- mean computed group velocity m/s: 101754
- mean envelope correlation vs published: 0.995318
- mean absolute peak-time delta vs published ps: 5.84163
- sheet thickness assignment note: Fig.2a/Fig.2b headers do not include thickness labels; assignment follows the published Fig.2d propagation peak-time ordering.
- computed group velocity rows:
  - thickness=224 nm | slope=7.5 ps/um | velocity=133333 m/s | r_squared=1
  - thickness=96 nm | slope=14.25 ps/um | velocity=70175.4 m/s | r_squared=0.998157
- bounded interpretation: Computed STFT target-band envelopes were derived from the raw Fig. 2 time traces using a 10 ps window and 5 ps step. This is a derived diagnostic, not raw published STFT-grid parsing. It can be compared to the published Fig. 2d target-band STFT intensity traces as a bounded reproducibility check.

## Limitations
- external-data probe only
- no proof of HAOS-IIP
- no claim that ferrons embody HAOS-IIP
- no new materials theory
- no consciousness claim
- metric proxies are provisional
- results depend on available data quality
