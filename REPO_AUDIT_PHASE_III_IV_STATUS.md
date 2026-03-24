# Repo Audit - Phase III / IV Status

## Purpose
This audit fixes the Phase III / IV authority chain so the repo answers, without ambiguity:
- what is authoritative
- what is superseded
- what is frozen
- what is still active
- which run supports the current claim

## Current top-line status
- `23.10` through `23.15`: `frozen`
- `24.1` through `24.5`: `frozen`
- `24.4` earlier narrow read stamped `20260316_121048`: `superseded`
- `24.5` first freeze attempt stamped `20260316_122224`: `superseded`
- Phase III current truth is frozen by `23.14` and paperized by `23.15`
- Phase IV current truth is frozen by `24.5`

## Authority chain

| Stage | Status | Authoritative timestamp | Authoritative note | Notes |
| --- | --- | --- | --- | --- |
| `23.10 / III-C1` | `frozen` | `20260316_102709` | `Stage_23_10_DK_Braid_Mechanism_Isolation_Probe_v1.md` | Canonical non-generality read |
| `23.11 / III-C2` | `frozen` | `20260316_104146` | `Stage_III_C2_Minimal_Effective_Model_Extraction_v1.md` | Canonical Phase III reduction read |
| `23.12 / III-C3` | `frozen` | `20260316_105448` | `Stage_III_C3_Effective_Phase_Diagram_Closure_v1.md` | Canonical one-sided smeared closure |
| `23.13 / III-C4` | `frozen` | `20260316_110320` | `Stage_III_C4_Mechanism_Isolation_Probe_v1.md` | Canonical phase-corridor-primary mechanism read |
| `23.14 / III-C5` | `frozen` | `20260316_111201` | `Stage_III_C5_Final_Consolidation_and_Freeze_v1.md` | Canonical Phase III freeze |
| `23.15 / III-P1` | `frozen` | editorial package | `papers/drafts/Phase_III_Clustered_DK_Effective_Sector_Draft_v1.md` | Paper-draft stage, not a simulation stage |
| `24.1 / IV-A1` | `frozen` | `20260316_114317` | `Stage_IV_A1_Smeared_Sector_Observable_Ledger_v1.md` | Canonical observable ledger |
| `24.2 / IV-A2` | `frozen` | `20260316_115346` | `Stage_IV_A2_Smeared_Sector_Effective_Law_Extraction_v1.md` | Canonical selector law |
| `24.3 / IV-A3` | `frozen` | `20260316_120334` | `Stage_IV_A3_Stable_Smeared_Quasi_Invariant_Probe_v1.md` | Canonical local ordering read |
| `24.4 / IV-A4` | `frozen` | `20260316_121603` | `Stage_IV_A4_Stable_Smeared_Law_and_Ordering_Extension_Boundary_Test_v1.md` | Canonical bounded extension-boundary map |
| `24.5 / IV-A5` | `frozen` | `20260316_122247` | `Stage_IV_A5_Phase_IV_Consolidation_and_Freeze_v1.md` | Canonical Phase IV freeze |

## Supersession records

### `24.2`
- authoritative run: `20260316_115346`
- superseded run: `20260316_115312`
- reason: tie-break selection logic was sharpened so equal-fit threshold candidates prefer the transport-side continuous separator over a binary closure-support proxy
- interpretation change: yes
- authoritative claim impact: the frozen selector is `flow_concentration_index <= 0.884308`, not the earlier refinement-flag proxy

### `24.4`
- authoritative run: `20260316_121603`
- superseded run: `20260316_121048`
- reason: the earlier read only covered a narrower width-retreat manifold, while the authoritative rerun executes the full four-class extension lattice
- interpretation change: yes
- authoritative claim impact: `24.4` now supports the full extension-boundary statement used by `24.5`

### `24.5`
- authoritative run: `20260316_122247`
- superseded run: `20260316_122224`
- reason: the initial freeze attempt used an overly blunt portable-core contradiction check
- interpretation change: no
- authoritative claim impact: freeze validity is corrected to `true`; scientific meaning is unchanged

## Critical claim lineage
- `braid is not family-wide`
  - first frozen at `23.10`
  - carried by `23.14`
  - preserved by `24.5`
- `smeared sector closes one-sidedly`
  - first frozen at `23.12`
  - carried by `23.14`
  - preserved as background in `24.5`
- `stable_closed_smeared = 1 iff flow_concentration_index <= 0.884308`
  - frozen at `24.2`
  - carried through `24.4`
  - frozen into Phase IV closure at `24.5`
- `inverse_flow_survival_ordering`
  - frozen at `24.3`
  - tested under extension in `24.4`
  - frozen as local structure in `24.5`
- `motif perturbation is earliest ordering boundary`
  - frozen at `24.4` authoritative rerun `20260316_121603`
  - carried into `24.5`
- `moderate corridor geometry is first selector-failure boundary`
  - frozen at `24.4` authoritative rerun `20260316_121603`
  - carried into `24.5`

## Stage integrity findings
- simulation stages `23.10` to `24.5` have complete chains of prompt -> script -> runsheet -> note -> machine outputs -> plots
- `23.15` is a paper-draft stage and intentionally does not follow the simulation chain; its authoritative package is prompt -> draft -> abstract -> figure map -> outline
- no downstream Phase III / IV note currently cites the superseded `24.4` timestamp after the audit edits

## Plot hygiene findings
- `24.2` contains both `20260316_115312` and `20260316_115346` plot families; only `115346` is authoritative
- `24.4` contains both `20260316_121048` and `20260316_121603` plot families; only `121603` is authoritative
- `24.5` contains both `20260316_122224` and `20260316_122247` plot families; only `122247` is authoritative
- superseded plot families are preserved for history and must be read through the manifest, not by filename freshness alone

## Prompt drift findings
- stored authoritative prompts for `24.4` and `24.5` now align with the frozen interpretation
- the earlier narrower `24.4` read survives only as superseded outputs, not as live prompt logic

## Paper alignment findings
- the Phase III draft remains aligned to frozen Phase III claims only
- no Phase IV paper exists yet, so there is currently no downstream paper claim outrunning the Phase IV freeze

## Commit-unit recommendation
The clean current commit unit is:
- authoritative `24.5` freeze package
- repo audit note
- `phaseIII_phaseIV_authoritative_manifest.json`
- supersession annotations in the live `24.2`, `24.4`, and `24.5` notes

## Fast answers
- authoritative `24.4` run: `20260316_121603`
- superseded `24.4` run: `20260316_121048`
- exact `24.2` threshold: `0.884308`
- `24.3` selected ordering relation: `IVA3_Q3 = inverse_flow_survival_ordering`
- Phase III freeze files: `23.14` note/json/csv plus `23.15` paper package
- commit-ready uncommitted files now: `24.5` package + audit note + manifest + note supersession edits

## Audit conclusion
The Phase III / IV repo state is now reconstructable from written authority records rather than memory. A new reader can determine current truth, superseded runs, and frozen boundaries without replaying the conversation history.
