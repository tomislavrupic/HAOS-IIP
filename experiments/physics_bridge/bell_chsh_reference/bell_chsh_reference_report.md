# Bell/CHSH Computational Reference Probe

Implemented fact: this sidecar runs a deterministic CHSH reference calculation and seeded finite-sample controls.
Design choice: CHSH uses `|E(a,b)+E(a,b')+E(a',b)-E(a',b')|` with the precommitted angle table.
Heuristic: sampled rows are finite-sample telemetry around analytic reference or classical-control generators.
Unverified hypothesis: no CST or HAOS Bell-recovery hypothesis is tested here.

## Verdict Labels
- BELL_REFERENCE_SANITY_PASS
- CST_BELL_STATUS_OPEN
- HAOS_DERIVATION_NOT_TESTED
- NO_PHYSICAL_EXPERIMENT

## CHSH Results
- deterministic local response max S: `2.000000000000`
- quantum analytic reference S: `2.828427124746`
- quantum sampled reference S: `2.820320000000` (95% CI `2.807888346077` to `2.832751653923`)
- local hidden-angle control sampled S: `2.010560000000` (95% CI `1.995404734708` to `2.025715265292`)
- uncorrelated control sampled S: `0.007880000000` (95% CI `0.000000000000` to `0.025410880805`)

## Boundary
- This is not a laboratory Bell experiment.
- This does not test loophole closure, detector effects, locality constraints, or real apparatus data.
- This does not derive Bell correlations from CST or HAOS-IIP.
- This does not change the frozen CST v0.2.2 checkpoint.
- The sampled local hidden-angle row is a finite-sample control; the exhaustive local enumeration is the actual local-bound check.
- `BELL_REFERENCE_SANITY_PASS` means only that the computational reference harness behaves as expected.
