# Periodic / Twisted L1 Flux Scan

## Purpose

Map the low periodic `L1` branch under a wider range of background flux values and test whether the mostly coexact family found in the earlier binary-flux experiment develops into a stable transverse band or remains dominated by harmonic/topological structure.

The direct questions are:

- does the low coexact branch move continuously under flux?
- does it remain mostly coexact?
- is the behavior robust across lattice sizes?
- does it still look harmonic/topological, or is there any sign of a cleaner propagating vector band?

## Setup

- substrate: periodic cubic lattice (`3`-torus)
- comparison baseline: open cubic lattice
- kernel weight:

$$
w_e = \exp\left(-\frac{h^2}{2\varepsilon}\right), \qquad \varepsilon = 0.2
$$

- periodic sizes: `n = 4, 5`
- flux scan: `m = 0, 1, 2, 3, 4`
- plaquette angle:

$$
\phi = \frac{2\pi m}{n^2}
$$

- operator: covariant edge/Hodge branch

$$
L_{1,A} = d_{0,A} d_{0,A}^* + d_{1,A}^* d_{1,A}
$$

- diagnostics per low mode:
  - eigenvalue
  - exact fraction
  - coexact fraction
  - `||d_0^* a||`
  - `||d_1 a||`
  - inverse participation ratio
  - qualitative support pattern

The selected low branch is defined as the lowest mode among the first four full-spectrum modes whose coexact fraction is at least `0.8`. If no such mode exists, the scan falls back to the best available low mode.

## Open vs Periodic Baseline

The open-box comparison remains exact-gradient dominated at the bottom of the spectrum. The periodic branch does not.

This reproduces the earlier result and remains the baseline reference for the flux scan:

- open low `L1` modes: exact-gradient leakage
- periodic `m = 0` low `L1` modes: harmonic torus-cycle structure

So the scan starts from an already non-scalar periodic branch.

## Low-Branch Flow

### Size `n = 5`

The selected low branch stays strongly coexact through the whole scan:

| `m` | `\phi` | `\lambda` | coexact fraction |
| --- | ---: | ---: | ---: |
| 0 | 0.000000 | 0.000000 | 1.0000 |
| 1 | 0.251327 | 0.013057 | 0.9755 |
| 2 | 0.502655 | 0.018837 | 0.9837 |
| 3 | 0.753982 | 0.024179 | 0.9876 |
| 4 | 1.005310 | 0.032351 | 0.9870 |

This branch is nearly smooth from `m = 1` onward. The support pattern changes from harmonic torus-cycle structure at `m = 0` to divergence-free circulation at nonzero flux.

### Size `n = 4`

The selected low branch also stays mostly coexact, but the ordering is less smooth:

| `m` | `\phi` | `\lambda` | coexact fraction |
| --- | ---: | ---: | ---: |
| 0 | 0.000000 | 0.000000 | 1.0000 |
| 1 | 0.392699 | 0.014110 | 0.9798 |
| 2 | 0.785398 | 0.171805 | 0.9228 |
| 3 | 1.178097 | 0.049192 | 0.9821 |
| 4 | 1.570796 | 0.275223 | 0.9077 |

This indicates that a mostly coexact family survives, but its spectral ordering re-shuffles enough that it does not look like one clean band.

## Projected Spectrum

The coexact/divergence-free projected branch confirms that the low family is not exact leakage.

For `n = 5`, the first nonzero projected eigenvalue moves as:

$$
0.000000 \to 0.027359 \to 0.035496 \to 0.041094 \to 0.058335
$$

For `n = 4`, the first nonzero projected eigenvalue moves as:

$$
0.000000 \to 0.032130 \to 0.197545 \to 0.083671 \to 0.327827
$$

The `n = 5` branch is smoother. The `n = 4` branch still shows strong reordering.

## Interpretation

### What is stable

- the low periodic branch remains mostly coexact across the full scan
- the branch persists for both tested lattice sizes
- the low family responds to flux instead of collapsing back into exact gradients
- at nonzero flux, the support pattern shifts away from pure torus-cycle zero modes toward divergence-free circulation structure

### What is not yet stable

- the low branch does not trace one clean smooth band across the full scan
- the smaller `n = 4` lattice still shows strong reordering under flux
- no isolated propagating Maxwell-like transverse band appears

So the scan strengthens the earlier result but does not yet upgrade it to a clean propagating vector sector.

## Direct Answers

### Does the low coexact branch move continuously under flux?

Not cleanly across the full scan.

For `n = 5` the branch is close to smooth after the zero-flux point, but for `n = 4` the branch reorders strongly enough that it does not behave like a single clean band.

### Does it stay mostly coexact across the scan?

Yes.

The selected low branch stays above `0.90` coexact fraction for both tested sizes and all scanned flux values.

### Is the behavior robust across lattice sizes?

Yes in the qualitative sense.

Both `n = 4` and `n = 5` retain a low mostly coexact branch with nontrivial flux response. The difference is in spectral smoothness, not in whether the branch exists.

### Does the branch still look harmonic/topological, or is there any sign of a more propagating vector band?

There is some departure from pure harmonic torus-cycle behavior once flux is turned on, because the low support patterns become divergence-free circulation modes rather than exact zero cycles.

But there is still no clean evidence for a propagating Maxwell-like band. The branch remains best interpreted as flux-deformed topological/coexact edge structure.

## Final verdict

The wider periodic flux scan preserves the main positive result:

- a low mostly coexact `L1` branch survives across the scanned flux values and lattice sizes
- the branch is genuinely flux-responsive
- the branch is not reducible to scalar-gradient leakage

But the scan also sharpens the current limitation:

- the low branch does not yet organize into one clean smooth transverse band
- the present evidence still points to harmonic/topological edge structure with flux-lifted circulation rather than a clean emergent photon-like sector

## Current verdict

- what worked: the low periodic `L1` branch stayed mostly coexact across `m = 0, 1, 2, 3, 4` and remained present for `n = 4, 5`
- what did not appear: a single smooth propagating transverse band separated cleanly from branch reordering and residual topological structure
- what must be tried next: puncture or defect backgrounds, plus larger periodic sizes, to test whether the low coexact family separates from torus-cycle bookkeeping into a more local vector sector
