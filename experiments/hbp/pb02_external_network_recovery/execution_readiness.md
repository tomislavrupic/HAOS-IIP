# PB-02 Execution Readiness

Status: frozen draft readiness note

Purpose

Make the external benchmark runnable without changing the claim ceiling.

Required external input

- PowerGraph dataset as described in `dataset_manifest.json`
- a local, repository-relative extraction path for the dataset
- no changes to the frozen split manifests after inspection

Required frozen artifacts

- `precommitment_contract.json`
- `dataset_manifest.json`
- `execution_contract.json`
- `development_split.json`
- `calibration_split.json`
- `holdout_split.json`
- `observation_map.md`

Execution sequence

1. Acquire the PowerGraph dataset from the public source named in the manifest.
2. Verify the local extraction path and keep it outside the holdout decision loop.
3. Run the benchmark entrypoint only after the split manifests are frozen.
4. Score baseline-only, HAOS-only, baseline+HAOS, ablated HAOS, and matched null on the same split path.
5. Record holdout results without adding features or moving thresholds after inspection.

What this note is not

- not a claim that the benchmark has been executed
- not a claim of empirical predictive value
- not a claim of external validation

Boundary

This note exists to make the next real-world holdout run explicit and auditable.
