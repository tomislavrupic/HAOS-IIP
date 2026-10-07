# PIX-7 candidate 01: nonlocal quasitopological gravity

Date: 2026-09-08. Status: LEAD FOR AUDIT; ALL-GATES PASS NOT_DEMONSTRATED.
This is a targeted literature screen and an elementary metric calculation, not a comprehensive review, independent replication, or HAOS-IIP physics result.

## Selection

Nominate **nonlocal quasitopological gravity (NLQT)** for the first adversarial audit. This is a research choice, not a finding that the theory can meet the whole specification.

[Bueno, Cano, Hennigar and Murcia, arXiv:2607.07790v1](https://arxiv.org/html/2607.07790v1) give a covariant nonlocal completion preserving regular spherical vacuum solutions. Equations (10)–(12) define the construction; select their entire function Ω(z)=γ²z², γ>0, and geometric-series curvature couplings α_n=α^(n−1), α>0. The demonstrated ghost argument covers maximally symmetric backgrounds and the spherical perturbation argument covers spherical vacuum backgrounds. The polynomial construction uses D≥5; four-dimensional nonpolynomial alternatives are mentioned. Exact nonlinear matter solutions remain an outstanding task. Do not import local-theory matter-collapse results into this completion.

**Candidate identity is only partially frozen:** the selected couplings and form factor still require an explicit choice of covariant curvature densities, operator domain, boundary prescription and matter action. In particular, a four-dimensional realization is not established merely by substituting D=4 into a metric formula. Completing this identity is milestone M0.

The related four-dimensional action construction of [Borissova and Carballo-Rubio, arXiv:2602.16773v2](https://arxiv.org/html/2602.16773v2), section IV.A, yields the Hayward geometry. Its second-order property is restricted to spherical symmetry. Treat this as a separate local benchmark; neither its existence nor a common metric transfers NLQT's stability results to it.

## Ten-gate audit

Statuses below are our assessment against PIX-7, not the authors' classifications. OPEN means evidence is missing, not that a counterexample has been proved.

| Gate | Current audit status | Required next evidence |
| --- | --- | --- |
| C1 All observables bounded | PARTIAL: metric check below only | Uniform estimates for matter, tidal probes and needed derivatives on a preregistered evolution domain |
| C2 Complete evolution | OPEN | Global continuation and predictability through horizons for the selected full law |
| C3 Recover tested GR | OPEN at full gate | Four-dimensional field equations, solution limits and quantitative precision-test comparison |
| C4 Recover Newton | OPEN for physical four-dimensional candidate | Calibrated Poisson limit after fixing the four-dimensional action and matter sector |
| C5 Conservation closure | STRUCTURAL ROUTE, not full pass | Vary the full action; verify the matter identity and constraint propagation under the chosen nonlocal prescription |
| C6 Stable, causal, well-posed | PARTIAL | Nonspherical perturbations, physical kinetic signs, characteristic analysis, and nonlinear initial-value control |
| C7 Dynamical regularization | PARTIAL | Smooth matter collapse in the same completed law, not an imposed static solution |
| C8 Independent admissibility | OPEN | Freeze data constraints and function spaces without excluding eventual divergences by definition |
| C9 Distinguishing observation | OPEN | Derive a measurable response with noise and nuisance parameters; an altered metric alone is insufficient |
| C10 Explain GR singularity | PARTIAL: matched metric limit below | Repeat the truncation comparison for actual coupled evolution |

No numeric score is assigned. Six partial successes cannot compensate for one fatal pathology.

## Independent analytic check: the four-dimensional Hayward benchmark

This calculation tests a metric shared by related constructions. It is not a derivation of the NLQT equations. Use G=c=1, M>0, ℓ>0, r≥0:

\[
ds^2=-f(r)dt^2+f(r)^{-1}dr^2+r^2d\Omega^2,
\qquad f(r)=1-\frac{2Mr^2}{r^3+2M\ell^2}.
\]

Writing ψ=(1−f)/r² gives, by algebra,

\[
\frac{\psi}{1-\ell^2\psi}=\frac{2M}{r^3},\qquad
\psi=\frac{2M}{r^3+2M\ell^2}\le\ell^{-2}.
\]

Here the inverse relation limits the curvature proxy. The left-hand response itself grows without bound at ψ=ℓ⁻². Thus this mechanism does not implement the earlier optional demand for a bounded response function; it targets bounded physical curvature. Insisting on both would be an additional gate that this benchmark does not meet.

For this single-function spherical metric the Kretschmann scalar is

\[
K=(f'')^2+4(f'/r)^2+4[(1-f)/r^2]^2.
\]

Set x=r³/(2Mℓ²). Direct differentiation yields

\[
\frac{1-f}{r^2}=\frac{1}{\ell^2(1+x)},\quad
\frac{f'}r=\frac{x-2}{\ell^2(1+x)^2},\quad
f''=\frac{-2x^2+14x-2}{\ell^2(1+x)^3}.
\]

For x≥0 their absolute values are at most ℓ⁻², 2ℓ⁻² and 7ℓ⁻². The last loose bound follows from
2x²+14x+2 ≤ 7(1+x)³. Consequently,

\[
0\le K\le69\ell^{-4}.
\]

This is a conservative bound uniform in r and positive M for this metric family. It does not establish bounds for other solutions or freely falling boosted measurements. In the usual curvature-sign convention,

\[
\lim_{r\to0}R=12/\ell^2,\qquad
\lim_{r\to0}K=24/\ell^4.
\]

At large radius,

\[
f=1-\frac{2M}{r}+\frac{4M^2\ell^2}{r^4}+O(r^{-7}).
\]

At fixed r>0, sending ℓ→0 recovers Schwarzschild; its K=48M²/r⁶ diverges as r→0. The expansion is uncontrolled when 2Mℓ²/r³ is order one. This demonstrates a nonuniform GR limit in a static family, not recovery of all Einstein phenomenology. M=0 is a separate Minkowski limit. Negative mass is outside this analytic check, not thereby excluded from the full theory's admissible domain.

## Immediate adversarial checks

1. **Four-dimensional action and nonlocal definition.** Fix one action, including branches of nonpolynomial invariants, a matter Lagrangian, and the functional definition of the entire operator. Verify that it is differentiable around the intended low-curvature backgrounds. A shared spherical equation is insufficient.
2. **Break spherical symmetry.** Compute physical odd/even nonspherical perturbations of its black-hole solution. A new ghost, loss of well-posedness, or uncontrolled growth rejects that realization before expensive observable fitting.
3. **Perturb the inner horizon with smooth matter.** The benchmark's horizon equation is r³−2Mr²+2Mℓ²=0. Its double root occurs at r=√3ℓ, M=3√3ℓ/4. For larger masses two positive horizons occur, with |κ_h|=|r_h²−3ℓ²|/(2r_h³), generically nonzero. This identifies a place to test amplification of infalling disturbances; it is not itself a proof of instability in the completed theory.
4. **Track matter independently.** A regular gravitational response to a distributional point source does not make that source's local density finite. Use smooth finite-energy data and evolve the source, rather than infer matter regularity from the metric.

The inner-horizon literature is not unanimous about growth rates: [Bonanno, Panassiti and Saueressig](https://arxiv.org/abs/2507.03581) report softened growth in some effective models. That does not give this selected action a boundedness proof. Likewise, [Bueno et al., matter-coupled QT study](https://arxiv.org/abs/2603.10110) explicitly separates geometric regularity from matter regularity. Neither source supplies a counterexample for an independently fixed NLQT smooth-data domain.

## Why the other screened routes were not selected first

- [Asymptotically free mimetic gravity](https://arxiv.org/abs/1905.01343) closely matches the limiting-curvature motivation. Its [stability extension](https://arxiv.org/html/2103.12442v3) addresses scalar gradient trouble using higher spatial derivatives but reports problems with the primordial spectra in the simplest model. It is a useful comparator, not an automatic C6/C9 pass.
- [Born–Infeld-inspired bouncing cosmologies](https://arxiv.org/abs/1707.08953) have tensor-instability obstacles in studied models. This does not reject every Born–Infeld theory.
- [Modified loop quantum cosmologies](https://arxiv.org/abs/1812.08937) distinguish removal of strong singularities from bounded curvature; pressure-related divergences remain relevant. Symmetry-reduced cosmology alone cannot pass the full gravity specification.

Decision: **NLQT is the lead research direction from this screen; the four-dimensional Hayward action is the concrete local benchmark. Neither is certified to pass C1–C10.** Do not combine separate papers' favorable results into a nonexistent single validated theory. The next executable scientific milestone is M0 plus the nonspherical stability audit, not an assertion that a new theory has been found.
