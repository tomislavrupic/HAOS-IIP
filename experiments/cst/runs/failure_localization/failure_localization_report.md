# CST Failure Localization Report

Implemented fact: this hostile diagnostic uses the existing CST toy benchmark and frozen Phase 15-18 artifacts.
Design choice: all component comparisons use distance values where lower means closer.
Heuristic: effect size is Cohen d with positive direction meaning target-target is closer than target-control.
Unverified hypothesis: none added; this report tests where the current CST toy benchmark loses discrimination.

## Verdict Labels
- TARGET NOT DISTINCT
- METRIC INSUFFICIENT
- CONTROL INVALID
- UNDERPOWERED
- MIXED / OPEN

## H0
The target branch is no more recoverable than matched controls under the current signature and perturbation suite.

Decision: fails_to_reject_H0. Underpowered: True.
This is not proof of equivalence.

## Primary Failure Points
- periodic_diagonal_augmented_control / d_inv: control mean 0.16693897649 <= target mean 0.285768828366
- randomized_edge_control / d_inv: control mean 0 <= target mean 0.285768828366
- shuffled_structure_control / d_inv: control mean 0 <= target mean 0.285768828366
- periodic_diagonal_augmented_control / d_struct: control mean 0.234595959596 <= target mean 0.373737373737
- periodic_diagonal_augmented_control / d_spec: control mean 0.351237535169 <= target mean 0.434517473357
- randomized_edge_control / d_spec: control mean 0 <= target mean 0.434517473357
- shuffled_structure_control / d_spec: control mean 0 <= target mean 0.434517473357
- randomized_edge_control / d_causal: control mean 0 <= target mean 0
- shuffled_structure_control / d_causal: control mean 0 <= target mean 0
- periodic_diagonal_augmented_control / d_time: control mean 0.046768707483 <= target mean 0.164399104308

## Ablation Summary
- No leave-one-component-out variant rescues the current FAIL verdict for the focus controls.

## Timing Ablation
- interpolated_non_front_times: verdict=FAIL, target_mean=0.20973712996157232, focus_control_mean=0.17299244397993363
- temporal_component_disabled: verdict=FAIL, target_mean=0.2188047350922088, focus_control_mean=0.17748889195959386
- direct_observed_timing_only: verdict=FAIL, target_mean=0.21959477924350726, focus_control_mean=0.17263551071588729

## Classification
- Genuine absence of target-specific closure: Open, but current evidence fails to reject H0.
- Signature aliasing: Present; branch_signature_id uses compact address fields and omits several discriminating fields.
- Metric insensitivity: Present; shuffled/random controls preserve many tested marginals and remain close.
- Normalization artifacts: Present for address-floor and relative small-denominator risks; see normalization_audit.csv.
- Control leakage: Present for randomized/shuffled controls through retained labels and descriptor marginals.
- Seed insufficiency: Present; configured seed count is 2 and target-target pairs are 3.
- Heuristic timing contamination: Partial risk; timing mode changes distributions and is isolated in timing_ablation.csv.
