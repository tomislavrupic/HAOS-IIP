"""Internal finite-graph hypothesis probe; no web, clipping, or state projection."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from haos_core import build_graph


def operator(side):
    graph = build_graph({"kind": "dk2d_periodic", "n_side": side,
                         "epsilon": 1.0 / side**2})
    lap = (graph.d0.conj().T @ graph.d0).real.tocsr()
    return graph, lap


def rhs(_, state, lap, coupling, mean_z):
    s, t = state.reshape(2, -1, 3)
    hs = lap @ s
    ht = lap @ t
    hs[:, 2] -= coupling * t[:, 0]
    ht[:, 0] -= coupling * (s[:, 2] - mean_z)
    return np.stack((np.cross(hs, s), np.cross(ht, t))).ravel()


def energy(state, lap, coupling, mean_z):
    s, t = state.reshape(2, -1, 3)
    return float(0.5 * np.sum(s * (lap @ s))
                 + 0.5 * np.sum(t * (lap @ t))
                 - coupling * np.dot(s[:, 2] - mean_z, t[:, 0]))


def initial_mode(graph, amplitude):
    q = amplitude * np.cos(2 * np.pi * graph.points[:, 0])
    s = np.column_stack((np.sqrt(1 - q*q), np.zeros_like(q), q))
    t = np.column_stack((q, np.zeros_like(q), np.sqrt(1-q*q)))
    return np.stack((s, t)).ravel()


def evolve(initial, lap, coupling, end=12., tolerance=1e-10, max_step=.04):
    mean_z = float(initial.reshape(2, -1, 3)[0, :, 2].mean())
    times = np.linspace(0, end, 601)
    sol = solve_ivp(rhs, (0, end), initial, args=(lap, coupling, mean_z),
                    method="DOP853", t_eval=times, rtol=tolerance,
                    atol=tolerance*0.01, max_step=max_step)
    if not sol.success or sol.t.size != times.size:
        raise RuntimeError(sol.message)
    states = sol.y.T.reshape(-1, 2, lap.shape[0], 3)
    energies = np.array([energy(y, lap, coupling, mean_z) for y in sol.y.T])
    masses = ((1 + states[:, 0, :, 2])/2).sum(axis=1)
    measures = {
        "max_spin_norm_error": float(np.max(np.abs(np.linalg.norm(states, axis=3)-1))),
        "max_mass_drift": float(np.max(np.abs(masses-masses[0]))),
        "max_absolute_energy_drift": float(np.max(np.abs(energies-energies[0]))),
        "energy_drift_per_site": float(np.max(np.abs(energies-energies[0])) / lap.shape[0]),
        "min_occupancy": float(np.min((1+states[:, 0, :, 2])/2)),
        "max_occupancy": float(np.max((1+states[:, 0, :, 2])/2)),
        "max_mediator_abs_x": float(np.max(np.abs(states[:, 1, :, 0]))),
        "function_evaluations": sol.nfev,
    }
    return times, states, measures


def tangent_rhs(_, state, lap, coupling):
    q, y, x, v = state.reshape(4, -1)
    return np.stack((-lap@y, lap@q-coupling*x,
                     lap@v, -lap@x+coupling*q)).ravel()


def tangent(initial, lap, coupling, times):
    s, t = initial.reshape(2, -1, 3)
    init = np.stack((s[:, 2], s[:, 1], t[:, 0], t[:, 1])).ravel()
    sol = solve_ivp(tangent_rhs, (times[0], times[-1]), init,
                    args=(lap, coupling), t_eval=times, method="DOP853",
                    rtol=1e-11, atol=1e-13, max_step=.04)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T.reshape(-1, 4, lap.shape[0])


def run():
    out = Path(__file__).resolve().parent / "outputs"
    out.mkdir(exist_ok=True)
    results = {
        "classification": "INTERNAL_FINITE_GRAPH_PROBE",
        "physical_gravity_demonstrated": False,
        "novelty_established": False,
        "state_projection_or_clipping": False,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "fixed_protocol": {"sides": [4, 6, 8], "J": 1., "K": 1., "g": 2.,
                           "amplitude": .001, "end_time": 12., "random_seeds": [7, 19, 41]},
        "size_runs": [], "random_runs": [], "weak_amplitude_runs": [],
    }
    main = None
    for side in [4, 6, 8]:
        graph, lap = operator(side)
        init = initial_mode(graph, .001)
        times, states, metrics = evolve(init, lap, 2.)
        linear = tangent(init, lap, 2., times)
        lam = float(np.linalg.eigvalsh(lap.toarray())[1])
        results["size_runs"].append({"side": side, "nodes": side**2, **metrics,
            "lambda_first": lam,
            "predicted_tangent_growth_rate": float(np.sqrt(lam*(2-lam))),
            "tangent_max_abs_q": float(np.max(np.abs(linear[:, 0]))),
            "full_max_abs_q": float(np.max(np.abs(states[:, 0, :, 2])))})
        if side == 6:
            main = (times, states, linear, lap, init)
    times, states, linear, lap, init = main
    _, fine, fine_metrics = evolve(init, lap, 2., tolerance=1e-12, max_step=.02)
    results["tolerance_check"] = {**fine_metrics,
        "reference_max_step": .04, "refined_max_step": .02,
        "maximum_state_difference": float(np.max(np.abs(fine-states)))}
    mean = float(init.reshape(2, -1, 3)[0, :, 2].mean())
    reverse = solve_ivp(lambda time, state: -rhs(time, state, lap, 2., mean),
                        (0., 12.), states[-1].ravel(), method="DOP853",
                        rtol=1e-12, atol=1e-14, max_step=.02)
    if not reverse.success:
        raise RuntimeError(reverse.message)
    results["inverse_flow_error"] = float(np.max(np.abs(reverse.y[:, -1]-init)))
    _, _, no_coupling = evolve(init, lap, 0.)
    results["zero_coupling_control"] = no_coupling

    for seed in [7, 19, 41]:
        rng = np.random.default_rng(seed)
        z = rng.normal(size=(2, lap.shape[0], 3))
        z /= np.linalg.norm(z, axis=2, keepdims=True)  # initial data only
        _, _, metrics = evolve(z.ravel(), lap, 2., end=3.)
        results["random_runs"].append({"seed": seed, **metrics})

    # Explicit near-capacity initial data: not rejected for being near a pole.
    graph, _ = operator(6)
    near = initial_mode(graph, 1.-1e-8)
    _, _, results["near_capacity_run"] = evolve(near, lap, 2., end=3.)

    for amplitude in [.04, .02, .01]:
        weak = initial_mode(graph, amplitude)
        tm, st, _ = evolve(weak, lap, 2., end=.5)
        lin = tangent(weak, lap, 2., tm)
        full_tangent = np.stack((st[:, 0, :, 2], st[:, 0, :, 1],
                                 st[:, 1, :, 0], st[:, 1, :, 1]), axis=1)
        results["weak_amplitude_runs"].append({"amplitude": amplitude,
            "relative_tangent_error": float(np.linalg.norm(full_tangent-lin)/np.linalg.norm(lin))})

    # Structural checks independent of trajectory fit.
    rng = np.random.default_rng(101)
    z = rng.normal(size=(2, lap.shape[0], 3))
    z /= np.linalg.norm(z, axis=2, keepdims=True)
    mean = float(z[0, :, 2].mean())
    zdot = rhs(0, z.ravel(), lap, 2., mean).reshape(z.shape)
    grad_s, grad_t = lap@z[0], lap@z[1]
    grad_s[:, 2] -= 2*z[1, :, 0]
    grad_t[:, 0] -= 2*(z[0, :, 2]-mean)
    permutation = rng.permutation(lap.shape[0])
    lperm = lap[permutation][:, permutation]
    perm_dot = rhs(0, z[:, permutation].ravel(), lperm, 2., mean).reshape(z.shape)
    source = np.cos(2*np.pi*graph.points[:, 0])
    phi = 2*np.linalg.pinv(lap.toarray())@source
    results["structural_checks"] = {
        "instantaneous_energy_derivative": float(np.sum(grad_s*zdot[0])+np.sum(grad_t*zdot[1])),
        "instantaneous_total_occupancy_derivative": float(.5*np.sum(zdot[0, :, 2])),
        "instantaneous_norm_derivative_max": float(np.max(np.abs(2*np.sum(z*zdot, axis=2)))),
        "relabeling_equivariance_error": float(np.max(np.abs(perm_dot-zdot[:, permutation]))),
        "tangent_poisson_residual": float(np.linalg.norm(lap@phi-2*source)/np.linalg.norm(2*source)),
    }
    all_metrics = results["size_runs"]+results["random_runs"]+[results["near_capacity_run"]]
    checks = {
        "numerical_invariants": all(m["max_spin_norm_error"]<1e-7 and m["max_mass_drift"]<1e-7
                                    and m["energy_drift_per_site"]<1e-7 for m in all_metrics),
        "occupancy_range": all(m["min_occupancy"]>=-1e-7 and m["max_occupancy"]<=1+1e-7 for m in all_metrics),
        "tangent_outgrows_capacity": all(m["tangent_max_abs_q"]>1 for m in results["size_runs"]),
        "tolerance_agreement": results["tolerance_check"]["maximum_state_difference"]<1e-6,
        "inverse_flow": results["inverse_flow_error"]<1e-6,
        "structural_identities": all(abs(v)<1e-10 for v in results["structural_checks"].values()),
        "weak_limit_error_decreases": all(a["relative_tangent_error"]>b["relative_tangent_error"]
              for a,b in zip(results["weak_amplitude_runs"],results["weak_amplitude_runs"][1:])),
    }
    results["diagnostic_checks"] = checks
    results["diagnostic_status"] = "PASS_SCOPED" if all(checks.values()) else "INVESTIGATE"
    (out/"results.json").write_text(json.dumps(results, indent=2)+"\n")
    np.savez_compressed(out/"trajectory.npz", times=times, compact=states, tangent=linear)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(times, np.max(np.abs(linear[:, 0]), axis=1), label="Tangent approximation")
    axes[0].plot(times, np.max(np.abs(states[:, 0, :, 2]), axis=1), label="Full compact dynamics")
    axes[0].axhline(1, color="black", linestyle=":", label="State capacity")
    axes[0].set(yscale="log", xlabel="Dimensionless time", ylabel="Maximum |load contrast|",
                title="Same initial data, different high-amplitude behavior")
    axes[0].legend()
    axes[1].plot(times, np.max((1+states[:, 0, :, 2])/2, axis=1), label="Maximum occupancy")
    axes[1].plot(times, np.min((1+states[:, 0, :, 2])/2, axis=1), label="Minimum occupancy")
    axes[1].set(ylim=(-.03,1.03), xlabel="Dimensionless time", ylabel="Occupancy fraction",
                title="No clipping or normalization during evolution")
    axes[1].legend()
    fig.suptitle("PIX-7 finite graph experiment — not a spacetime or GR simulation")
    fig.tight_layout()
    fig.savefig(out/"comparison.png", dpi=170)
    plt.close(fig)
    print(json.dumps({"status": results["diagnostic_status"], "checks": checks,
                      "size_runs": results["size_runs"],
                      "tolerance": results["tolerance_check"],
                      "weak_limit": results["weak_amplitude_runs"]}, indent=2))


if __name__ == "__main__":
    run()
