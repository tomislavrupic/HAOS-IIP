# Stage 24.2 / IV-A2 -- Smeared-Sector Effective Law Extraction

Context:
HAOS = Harmonic Address Operating System.
IIP = Interaction Invariance Physics.

Purpose

Run the first true effective-law extraction pass for Phase IV on the frozen clustered DK branch.

This stage attempts to extract the smallest valid effective law for the surviving smeared-transfer sector using only the accepted Phase IV observable core from Stage 24.1.

This is not a braid-recovery stage.
This is not a new mechanism-discovery stage.
This is not a family-expansion stage.

Its purpose is narrower:
1. determine whether the smeared sector admits a compact reduced law,
2. express that law only in terms of the accepted core observables,
3. keep corridor variables as external context coordinates rather than promoted state variables,
4. reject any model that only works by reintroducing previously excluded quantities,
5. freeze the minimal law if a genuine compression is achieved.

Authoritative upstream inputs

Treat the following stages as frozen and authoritative:
- `23.10 / III-C1`
- `23.11 / III-C2`
- `23.12 / III-C3`
- `23.13 / III-C4`
- `23.14 / III-C5`
- `23.15 / III-P1`
- `24.1 / IV-A1`

Accepted primary observables from `24.1`:
- `topology_survival_time`
- `refinement_stability_flag`
- `weak_coupling_stability_flag`
- `reverse_stability_flag`
- `flow_concentration_index`

Accepted context coordinates:
- `phase_corridor_position`
- `phase_corridor_width`

Accepted as secondary or derived only:
- `grade_exchange_coherence_level`
- `smeared_regime_membership_flag`

Explicitly rejected from the law core:
- `transfer_asymmetry`
- `coherence_decay_rate`
- `topology_return_error`
- `braid_like_exchange_indicator`
- `mechanism_sensitivity_score`

Core question

Can the smeared sector's persistence and stability be captured by a compact reduced law using only the accepted primary observables, while treating phase-corridor coordinates as external context rather than internal state variables?

Scientific constraints

Allowed core state variables:
- `topology_survival_time`
- `refinement_stability_flag`
- `weak_coupling_stability_flag`
- `reverse_stability_flag`
- `flow_concentration_index`

Allowed context coordinates:
- `phase_corridor_position`
- `phase_corridor_width`

Allowed secondary descriptors for diagnostics only:
- `grade_exchange_coherence_level`
- `smeared_regime_membership_flag`

Forbidden in the final law and selection logic:
- `transfer_asymmetry`
- `coherence_decay_rate`
- `topology_return_error`
- `braid_like_exchange_indicator`
- `mechanism_sensitivity_score`

Do not allow braid-labeled quantities back into the explanatory core through feature engineering, relabeling, hidden proxy variables, or rescue language.

Required outputs

Produce:
1. `Stage_IV_A2_Smeared_Sector_Effective_Law_Extraction_v1.md`
2. `stage24_2_smeared_sector_effective_law_extraction.py`
3. `stage24_2_smeared_sector_effective_law_runs.json`
4. `*_stage24_2_smeared_sector_effective_law_extraction.json`
5. `*_stage24_2_smeared_sector_effective_law_model_ledger.csv`
6. figures:
   - `*_stage24_2_law_family_comparison_table.png`
   - `*_stage24_2_minimal_law_decision_panel.png`
   - `*_stage24_2_prediction_vs_truth_panel.png`
   - optional `*_stage24_2_context_conditioning_panel.png`

Stage logic

Compare a small family of candidate effective-law forms and keep only the least complex form that actually works under the frozen branch constraints.

At minimum compare:
- Model family A -- threshold rule
- Model family B -- simple reduced surrogate
- Model family C -- context-conditioned threshold rule

You may add one additional minimal family if it materially clarifies selection.

Target outputs of the law

Use a branch-valid smeared-sector target such as:
- `stable_closed_smeared`
- `transient_mixed`

Do not let the target drift back into braid classification.

Success criteria

The stage succeeds only if it produces a compact law that:
- uses only the accepted primary observables as core inputs,
- treats corridor variables as external context only,
- predicts the selected smeared-sector target with good fidelity,
- compresses behavior beyond a trivial lookup table,
- remains interpretable and branch-local,
- does not require any rejected variable,
- yields a clear winning minimal model family.

Suggested success booleans:
- `effective_law_found`
- `minimal_model_selected`
- `uses_only_ledger_core`
- `context_variables_not_promoted`
- `rejected_variables_unused`
- `predictive_compression_valid`

Failure criteria

Mark the stage as failed if:
- the best-performing model requires a rejected variable,
- the best-performing model silently depends on braid-labeled structure,
- the law is really just a lookup table,
- the model cannot remain stable under the frozen branch context,
- corridor variables have to be promoted into state variables,
- or model complexity rises beyond what the branch justifies.

Consistency checks

Before accepting any law candidate, verify:
1. only the accepted primary observables are used as core inputs,
2. corridor variables appear only as context,
3. no forbidden variable appears directly or indirectly,
4. the law is simpler than a full empirical lookup table,
5. the law beats a naive baseline,
6. the law does not contradict the frozen smeared-dominant Phase III closure.

Preferred model-selection principle

Select the winner by:
1. branch validity
2. compliance with the observable ledger
3. simplicity
4. interpretability
5. predictive adequacy

Do not select a more complex model merely because it improves fit slightly.

Required final interpretation

If successful, end with language close to:

On the frozen clustered DK branch, the smeared sector admits a compact effective law built from the accepted Phase IV observable core. The selected law predicts branch-valid smeared-sector persistence/stability without reintroducing rejected braid- or asymmetry-based quantities. Phase-corridor variables remain external conditioning coordinates rather than promoted state variables.

If unsuccessful, say clearly that no compact law was found under the current observable ledger.
