# Gauss Law Test v1

Timestamped result: `data/20260310_192921_stage8_gauss_law_test.json`

Config: sizes=[12, 16, 20], epsilon=0.2, source_strength=1.0.

Largest grid (n=20) metrics:
- best-fit slope = 1.000000
- correlation = 1.000000000000
- relative L2 residual = 4.141e-15
- L_inf residual = 3.854e-15
- max integrated mismatch = 3.128e-03

Observation: the induced edge field satisfies a local divergence constraint that matches the inserted scalar source to numerical precision on the tested periodic complexes

Conclusion: the Gauss-law diagnostic is consistent with a discrete divergence constraint of the form d0* A = c rho with c near unity and small integrated flux mismatch
