# Ferron Harmonic-Address Persistence Summary

## Status
- status: PASS
- mode: ferron_harmonic_address_persistence_summary
- target_frequency_THz: 3.13
- combined_persistence_proxy: 0.957192
- raw_stft_time_frequency_grid_status: NO_STFT_DATA_FOUND

## Components
- target_peak_presence: 1 (included; source=spectral_feature_summary.json)
- spectral_recoverability: 0.928455 (included; source=spectral_feature_summary.json)
- no_sustained_k_star: 1 (included; source=spectral_feature_summary.json)
- published_propagation_monotonicity: 1 (included; source=published_stft_trace_summary.json)
- published_velocity_consistency: 0.794118 (included; source=published_stft_trace_summary.json)
- computed_envelope_match: 0.995318 (included; source=computed_stft_summary.json)
- computed_velocity_consistency: 0.982456 (included; source=computed_stft_summary.json)
- raw_stft_grid_boundary: None (not included; source=stft_feature_summary.json)

## Bounded Interpretation
The ferron sidecar shows a stable target-window spectral feature near 3.13 THz, published target-band STFT intensity traces with monotonic propagation delay, and computed-from-time-trace STFT envelopes that closely match the published traces. The raw STFT time-frequency grid remains unavailable in the downloaded XLSX files; that boundary is preserved. This is a bounded external-data persistence summary, not a proof or materials-theory claim.

## Claim Boundary
- external-data sidecar summary only
- not a proof of HAOS-IIP
- not a claim that ferrons embody HAOS-IIP
- not a new materials theory
- not a raw STFT-grid claim
- computed STFT outputs are derived diagnostics, not published raw STFT maps
