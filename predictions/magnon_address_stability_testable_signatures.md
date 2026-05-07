# HAOS-IIP Magnon Address Stability - Testable Signatures

**Version**: `0.1`

**Seeded by**: Materials Bridge Line A ferron audit (`0.957192` persistence
proxy) plus user-reported alpha-Fe2O3 magnon audit (`0.902` persistence proxy).

This note creates the parallel spin-sector prediction layer. It is not a new
spin-wave theory, not an external magnon proof of HAOS-IIP, and not a claim
that ferrons and magnons share a literal substrate.

## Provenance Boundary

Local repo provenance:

- ferron sidecar: committed and reproducible
- ferron persistence proxy: `0.957192`
- ferron spectral recoverability: `0.928455`
- ferron target-window peak audit: `15/15` near `3.13 THz`
- ferron raw STFT/time-frequency grid status: `NO_STFT_DATA_FOUND`

User-provided magnon anchor:

- alpha-Fe2O3 audit persistence proxy: `0.902`
- raw public STFT-grid boundary: `NO_PUBLIC_STFT_GRID_FOUND`
- audit artifact location in this repo: not found at note creation time

The `0.902` value is therefore recorded as a user-provided audit anchor until
the corresponding magnon audit files are committed.

## Candidate External Anchors

These references motivate audit targets. They do not by themselves validate
HAOS-IIP recoverability metrics.

- Li et al., "Anisotropic Coherent Propagation of Sub-Terahertz Magnons in
  Altermagnetic alpha-Fe2O3," Advanced Optical Materials, 2026,
  DOI `10.1002/adom.202503604`.
- Yang et al., "Spin-orbit torque manipulation of sub-terahertz magnons in
  antiferromagnetic alpha-Fe2O3," Nature Communications, 2024,
  DOI `10.1038/s41467-024-48431-w`.
- Sun et al., "Observation of Chiral Magnon Band Splitting in Altermagnetic
  Hematite," Physical Review Letters, 2025.
- Choe et al., "Long-lived zone-boundary magnons in an antiferromagnet,"
  Nature Communications, 2025, DOI `10.1038/s41467-025-60287-2`.
- Yuan et al., "Inelastic Neutron Scattering Study of Field Dependence of
  Magnetic Excitations in CoTiO3," Physical Review B, 2024,
  DOI `10.1103/PhysRevB.109.174440`.
- University of Vienna press release, May 4, 2026, reporting magnon lifetimes
  up to `18 us`.

## Core Prediction

Discrete harmonic addresses in altermagnetic and antiferromagnetic substrates
should appear as long-lived, narrowband coherent magnons with crystal-axis
anisotropy and high spectral recoverability when excitation is matched to the
spin-sector address.

The spin-sector HAOS mapping is:

- carrier: magnetic order / spin precession
- perturbation: optical pulse, strain pulse, spin-orbit torque, magnetic field,
  temperature, thickness, boundary condition, or controlled defect scattering
- recoverable structure: coherent magnon / spin-wave mode
- observable readout: narrowband Raman, THz, neutron, magneto-optic, or
  spin-torque magnon feature
- visible failure: broadband thermalization, rapid damping, amplitude collapse,
  loss of directionality, or loss of target-window recovery

## HAOS-PRED-010 - Crystal-Orientation-Dependent Narrowband Magnon Address

Prediction: sub-THz exchange magnons should show highest coherence and
narrowest linewidth only for specific wavevector / crystal-axis alignments,
especially in Dzyaloshinskii-Moriya-active or altermagnetic directions.

Expected signature:

- linewidth proxy below `0.1 THz` where instrumental resolution permits
- target-window frequency stability at or above `0.95`
- target-window recoverability at or above `0.85`
- propagation beyond a few microns when imaging data are available

Falsification:

If matched orientation scans show isotropic broadband response, or if
pre-registered crystal-axis alignment does not improve recoverability over
misaligned controls, this prediction fails for that material class.

## HAOS-PRED-011 - Thickness / Boundary-Dependent Group Velocity

Prediction: magnon packet velocity should scale with thickness, confinement, or
boundary condition, analogous in form to the ferron thickness dependence but not
assumed to share the same mechanism.

Expected signature:

- group-velocity proxy changes systematically across controlled thickness or
  boundary groups
- peak arrival time or phase delay varies monotonically with propagation
  distance along at least one crystal direction
- target-window recoverability remains above `0.70` before visible propagation
  failure

Falsification:

If controlled-thickness or boundary series show no systematic velocity scaling,
no distance-ordered delay, and no recoverability change beyond uncertainty while
the target-window signal remains measurable, this prediction fails for that
material class.

## HAOS-PRED-012 - Target-Address Recovery After Perturbation

Prediction: after a secondary pulse, field quench, strain perturbation, or
spin-orbit torque below irreversible damage threshold, the system should
preferentially recover the original narrowband magnon address.

Expected signature:

- post-perturbation peak frequency returns to the pre-registered target window
- frequency stability is at least `0.95`
- target-window spectral recoverability is at least `0.85`
- magnetic-field tuning may shift the packet within the registered family
  without broadband thermalization

Falsification:

If secondary perturbation repeatedly causes complete thermalization, a different
dominant mode, or no preferential return to the original target window under
non-damaging controls, this prediction fails for that perturbation regime.

## HAOS-PRED-013 - Hybrid Magnon-Phonon / Cross-Sector Locking

Prediction: in van der Waals magnetic hybrids or strongly coupled magnetic
heterostructures, magnon addresses should lock to compatible phonon or polar
modes, producing enhanced lifetime, cleaner linewidth, or more directional
propagation than either uncoupled control.

Candidate systems include CrSBr, MnBi2Te4, FePSe2-like vdW magnets, and
engineered magnon-phonon cavities.

Expected signature:

- hybrid mode has higher recoverability than magnon-only and phonon-only
  controls
- hybrid linewidth or lifetime improves under matched coupling conditions
- target-window peak remains stable under a perturbation where uncoupled
  controls broaden or collapse

Falsification:

If controlled hybridization produces no lifetime, linewidth, directionality, or
recoverability benefit beyond uncoupled controls, this prediction fails for that
hybrid class.

## Boundary

These predictions are field-facing diagnostic candidates. They do not claim new
physics, a universal harmonic grid, or a replacement for spin-wave theory. The
next required step is to commit the alpha-Fe2O3 audit artifact or build a new
magnon data bridge that loads public spectra / traces and runs the same bounded
recoverability audit used for ferrons.
