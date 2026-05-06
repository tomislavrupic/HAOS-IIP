# Materials Bridge Line A - Ferron Coherence Recoverability Probe

This is an external-data sidecar experiment for HAOS-IIP.

It targets published ferron / NbOI2 coherent polarization-wave data and asks a
bounded software question: can HAOS-style recoverability metrics be applied to
real materials coherence readouts when suitable public data files are available?

This does not claim ferrons prove HAOS-IIP. It does not claim ferrons embody
HAOS-IIP. It does not modify the frozen HAOS-IIP core, frozen phases,
telemetry definitions, theory files, existing experiments, or paper releases.
It is not a new materials theory and it is not a quantum consciousness claim.

## Sources

Primary target:

- Choe et al., "Observation of coherent ferron emission and propagation,"
  Nature Materials, 2026, DOI `10.1038/s41563-026-02597-4`
- Source data DOI: `10.6084/m9.figshare.31895293`

Secondary reference if primary public data is unavailable:

- Wang et al., "Ultrafast dynamics of ferroelectric polarization of NbOI2
  captured with femtosecond electron diffraction," Nature Communications,
  2025, DOI `10.1038/s41467-025-63533-9`

## Run

```bash
uv run --with numpy --with matplotlib python experiments/materials_bridge/ferron_coherence_demo/run_ferron_coherence_demo.py
```

The runner first attempts to fetch the Figshare source-data record. If public
data cannot be downloaded automatically, it creates a small synthetic smoke-test
dataset solely to verify the pipeline. Smoke-test output is not scientific
validation and must not be cited as external evidence.

To run only the software smoke test:

```bash
uv run --with numpy --with matplotlib python experiments/materials_bridge/ferron_coherence_demo/run_ferron_coherence_demo.py --force-smoke-test
```

## Outputs

- `data_manifest.json`
- `outputs/data_inspection.json`
- `outputs/results.csv`
- `outputs/summary.json`
- `outputs/real_data_metrics.csv` when real data parses
- `outputs/smoke_test_metrics.csv` when only smoke-test data is available
- plot PNGs under `outputs/`
- `ferron_coherence_validation.md`

Raw downloaded files, if available, are stored untouched under:

```text
outputs/raw/
```

Archives are extracted only into derived storage for parsing:

```text
outputs/derived_extracted/
```

## HAOS Mapping

- carrier: ferroelectric polarization order
- perturbation: laser pulse, fluence, distance, temperature, sample thickness,
  boundary/disorder effects
- recoverable structure: coherent ferron / polarization-wave mode
- observable readout: narrow-band THz emission, 3.13 THz amplitude, transient
  reflectance oscillation, Fourier/STFT peak
- visible failure: loss of narrow-band signal, spectral broadening, amplitude
  collapse, incoherent background dominance, or propagation failure

## Validation Boundary

`PASS` is allowed only when real external data is loaded, provenance is recorded,
at least one ferron-related signal is parsed, at least one recoverability curve
is computed, and k-star / safety-margin logic runs deterministically.

`DATA_UNAVAILABLE_SMOKE_TEST_ONLY` means the code path ran on a synthetic
pipeline-check dataset. It is software verification only.

`FAIL` means raw data existed but could not be parsed into supported signal
forms, or no recoverability curve could be computed.
