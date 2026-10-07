# PIX-7 constructed candidate: finite-capacity relational exchange v1

Created 2026-09-08 in response to the user's request to build a candidate from HAOS-IIP without searching the web.

Classification: **CONSTRUCTED_FINITE_GRAPH_MODEL / GR_RECOVERY_NOT_DEMONSTRATED**.
Originality: a construction made in this session using familiar mathematical ingredients; historical novelty is unknown. Neither an exhaustive theory search nor a probability measure over theories has been supplied.

## 1. The proposed mechanism

Treat a site's load as one component of a finite internal state. Neighboring states exchange load through conservative rotations. A second internal state mediates an attractive response. Concentration changes the orientations and therefore the exchange currents; no component is clipped after evolution.

The intended distinction is between a finite physical state and an unbounded tangent approximation to it. This model tests whether that distinction can coexist with attraction, conservation and reversible dynamics. It does not yet explain spacetime singularities.

HAOS starting point: the quadratic and relational Dirichlet forms in `HAOS_Quadratic_Defect_Theorem_T1_v1.md` and `HAOS_Relational_Dirichlet_Theorem_T2_v1.md`. We **assume positive edge weights explicitly**, rather than infer them from global positivity alone. The latter does not in general force every pair coefficient to be positive. No existing theorem or frozen result is edited by this sidecar.

## 2. A bounded ansatz search, not a search of all theories

Restrict attention to two three-component unit vectors per site, pairwise isotropic bilinear exchange within each species, and an onsite bilinear cross-coupling. These restrictions are design assumptions, not consequences of HAOS-IIP.

Write the cross-coupling as sᵀCt. Requiring invariance under rotations of s around its z-axis for every s,t eliminates the x and y rows of C. The remaining vector can be aligned with the mediator's x-axis by a change of its internal basis, since its exchange is isotropic. Thus the permitted onsite bilinear term has the form −g s_z t_x. In this restricted ansatz, charge conservation selects the coupling's shape rather than an aesthetic choice of correction function.

This is a small algebraic reduction of candidate space. It says nothing about other state spaces, nonlinear couplings, gauge theories or tensor sectors.

## 3. State space and interaction law

Let G be a finite undirected connected graph with positive symmetric edge weights w_ij and weighted Laplacian L. Let

\[
s_i,t_i\in S^2,\qquad n_i=(1+s_i^z)/2\in[0,1].
\]

The vectors are internal states, not directions in physical space. The compact manifold is an assumption about the substrate, not a derived law of nature. Define fixed parameters J,K>0, g≥0 and a background m∈[−1,1]. The neutral sector has ∑s_i^z=Nm. The Hamiltonian is

\[
H=\frac J2\sum_{\{i,j\}}w_{ij}|s_i-s_j|^2
 +\frac K2\sum_{\{i,j\}}w_{ij}|t_i-t_j|^2
 -g\sum_i(s_i^z-m)t_i^x.
\]

Choose the spin Poisson structure with evolution

\[
\dot s_i=\frac{\partial H}{\partial s_i}\times s_i,
\qquad
\dot t_i=\frac{\partial H}{\partial t_i}\times t_i.
\]

Equivalently,

\[
\dot s_i=[J(Ls)_i-g t_i^x e_z]\times s_i,
\qquad
\dot t_i=[K(Lt)_i-g(s_i^z-m)e_x]\times t_i.
\]

This is a local first-order law on the supplied graph. In each run m is fixed once from the neutral initial charge; it is not recomputed as a global feedback signal. The completeness argument also applies off the neutral sector with a fixed m, but the massless static Poisson equation below then lacks its zero-mean compatibility condition.

## 4. What can actually be proved

### A. Invariant capacity and complete finite-graph evolution

By orthogonality of a cross product,

\[
\frac d{dt}|s_i|^2=\frac d{dt}|t_i|^2=0.
\]

Therefore n_i stays in [0,1] without clipping, including initial states at n_i=0 or 1. The vector field is polynomial in ambient coordinates and smooth on the compact manifold (S²)^(2N). Local uniqueness plus compactness gives a unique solution for every finite initial state and every t∈R. No divergent vector field is hidden at the poles.

This is a global ordinary-differential-equation result on each fixed finite graph. All fixed smooth observables on this compact state space are bounded. It does not prove that every real physical observable has such a representation: a singular readout such as 1/(1−n) would defeat that inference, and the physical readout map is still missing.

### B. Conservation and finite transfer

Let Q=∑n_i. The coupling produces no z-torque on s; the exchange terms cancel pairwise. Hence dQ/dt=0. The current into i from j is

\[
I_{j\to i}=\frac{Jw_{ij}}2(s_i\times s_j)_z,
\quad I_{i\to j}=-I_{j\to i},\quad
|I_{j\to i}|\le Jw_{ij}/2.
\]

At n_i=0 or 1 its instantaneous load current vanishes because s_i has no transverse component. Its other components can still rotate, allowing subsequent transfer; the state is not artificially frozen.

Furthermore,

\[
\dot H=\sum_i\nabla_{s_i}H\cdot(\nabla_{s_i}H\times s_i)
 +\nabla_{t_i}H\cdot(\nabla_{t_i}H\times t_i)=0.
\]

With d_i=∑_j w_ij,

\[
|\dot s_i|\le2Jd_i+g,\qquad |\dot t_i|\le2Kd_i+2g.
\]

The energy is bounded below by −2gN and above by 2(J+K)∑edges w+2gN. These are finite-graph statements; quantum ghost freedom and relativistic energy positivity have not been established.

### C. An attractive graph-Poisson sector

Near s=e_x, t=e_z, in the m=0 sector, use tangent coordinates

\[
s\simeq(1,y,q),\quad t\simeq(x,v,1).
\]

To quadratic order,

\[
H_2=\frac J2(q^TLq+y^TLy)+\frac K2(x^TLx+v^TLv)-gq^Tx.
\]

In the static, mean-zero mediator sector,

\[
KLx=gq,\qquad
H_{\mathrm{induced}}=-\frac{g^2}{2K}q^TL^+q.
\]

Thus like-signed contrasts have an attractive Green-function interaction in this approximation. This is a static elimination, not a claim that the undamped mediator relaxes automatically or follows arbitrary moving sources adiabatically.

If an independently justified three-dimensional continuum limit gave L→−∇², the carrier potential Φ=−gx would satisfy ∇²Φ=(g²/K)q. That is a conditional route to a Newton-shaped scalar equation. The actual numerical substrate here is two-dimensional, and no physical mass calibration, universal free fall, inverse-square measurement, or Einstein limit follows from this algebra.

### D. A controlled failure of the approximation

The tangent evolution is

\[
\dot q=-JLy,\quad \dot y=JLq-gx,\quad
\dot x=KLv,\quad \dot v=-KLx+gq.
\]

For J=K=1 and an L-eigenmode λ>0, the q=x branch satisfies

\[
\ddot q=\lambda(g-\lambda)q.
\]

If g>λ it grows exponentially. The full model remains compact. This demonstrates how a legitimate weak-amplitude approximation can acquire an unbounded continuation when used outside its regime.

It is **not a finite-time singularity**, and it is not the GR singularity mechanism. The uniform reference is linearly unstable in this regime by design; global boundedness must not be described as linear stability. What is shown is bounded nonlinear evolution after a growing perturbation.

## 5. Execution and evidence

Code: `experiments/physics_bridge/pix7_finite_capacity_exchange/run_probe.py`.
Graph obtained from the existing `haos_core.build_graph` with zero cycle phases. No public APIs changed. DOP853 integrates the ambient equations; unit norms are never projected during evolution.

Protocol was written before the first run: J=K=1, g=2, contrast amplitude 0.001, t=12, fixed edge weight exp(−1/2), sizes 4², 6², 8². The size sweep is not a continuum-convergence study. Three random initial states and a near-capacity state test outside the symmetric seed family. No failed seed or parameter trial was removed. The initial tolerance-only comparison used the same maximum step; verification was strengthened to halve that step and check the inverse flow. No model parameters changed.

| Sites | Maximum q magnitude, tangent | Maximum q magnitude, full |
| --- | --- | --- |
| 16 | 61.7790 | 0.97703 |
| 36 | 30.9122 | 0.90791 |
| 64 | 4.81772 | 0.80476 |

For the 36-site case, occupancy stays approximately between 0.04604 and 0.95396. The full trajectory turns around after concentration while the tangent model keeps growing. This is load redistribution in a graph, not a spacetime bounce.

The zero-coupling control stays within occupancy 0.4995–0.5005. Halving the weak initial amplitude gives roughly a fourfold reduction in relative tangent error: 4.7131×10⁻⁴, 1.1782×10⁻⁴, 2.9455×10⁻⁵ for amplitudes 0.04, 0.02, 0.01. This supports the derived weak-amplitude expansion over the tested short interval.

The tested random-state energy drift per site is below 4.4×10⁻¹³, and the instantaneous graph-Poisson residual is about 2.1×10⁻¹⁵. Exact proofs above and numerical checks are distinct evidence. Machine-readable results, source hash, controls and tolerance diagnostics are retained in `outputs/results.json`; arrays are in `outputs/trajectory.npz`; the figure is `outputs/comparison.png`.

## 6. Hostile PIX-7 verdict

| Gate | Verdict for this construction |
| --- | --- |
| C1 Bounded physical observables | Finite internal components/currents proved. Physical density, curvature and their readout map missing. |
| C2 Complete evolution | Complete finite-graph ODE proved. Geodesic completeness has no defined metric to test. |
| C3 Einstein recovery | NOT_DEMONSTRATED. No tensor gravitational sector; cannot call this a replacement for GR. |
| C4 Newton recovery | Conditional scalar graph-Poisson shape only. Full Newtonian phenomenology not recovered. |
| C5 Conservation closure | Exact toy Hamiltonian and total occupancy conservation. Covariant stress-energy closure undefined. |
| C6 Pathology exclusion | Smooth global ODE with bounded variables; homogeneous reference has a known growth instability. Relativistic causality and quantum modes undefined. |
| C7 Regularization mechanism | Finite load currents and state rotation; observed redistribution. No physical collapse calculation. |
| C8 Independent admissibility | Every point of the chosen compact state space is admitted, including poles. Compactness itself is an unverified substrate postulate. |
| C9 Observation | Full/tangent contrast is a simulation discriminator, not an astronomical prediction. |
| C10 Apparent GR singularity | Tangent extrapolation failure demonstrated; GR truncation not derived. |

**As a completed gravity theory this candidate does not qualify. As a finite conservative interaction mechanism it is concrete, executable and mathematically bounded.**

## 7. Conditions that stop promotion

- Finite occupancy is not finite density unless a physical cell volume is defined and controlled. Sending that volume to zero can restore divergent density. No minimum volume is derived here.
- A supplied graph, preferred time and internal axes are assumed. The experiment derives neither space, spatial dimension, Lorentz symmetry nor diffeomorphism invariance.
- Finite local time derivatives do not imply a strict light cone. Continuous-time lattice evolution can produce small tails at distant sites. Relativistic causality is not claimed.
- Compactness permits recurrence and strong nonlinear behavior; it does not guarantee equilibrium, relaxation or generic bounce outcomes.
- A graph Laplacian acting on scalar/internal vectors is not a massless spin-2 field with universal stress-energy coupling. The two transverse directions of a unit vector must not be renamed gravitational-wave polarizations.
- The analytic current bounds depend on degree and couplings. Refinement with diverging weights requires a separate uniform-limit theorem.

The next research decision is whether this finite exchange mechanism can coexist with a constrained tensor sector that reproduces linearized Einstein dynamics. Until such a construction and its constraint algebra exist, use it as a regularization laboratory only. Adding an Einstein solver beside it would impose recovery rather than derive it from HAOS-IIP.
