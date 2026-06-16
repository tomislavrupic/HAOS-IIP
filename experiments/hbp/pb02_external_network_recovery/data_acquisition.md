# PB-02 Data Acquisition

Status: frozen acquisition note

Purpose

Record the public source and local extraction convention needed before the
PB-02 holdout run.

Public source

- dataset landing page: <https://figshare.com/articles/dataset/PowerGraph/22820534>
- supplement: <https://github.com/PowerGraph-Datasets>
- archive size advertised by figshare: 5.23 GB

Local data path

- `/Volumes/Samsung T5/2026/HAOS/HAOS DOCS/DATA/Powergraph`

Expected dataset layout

- `dataset_pf_opf/ieee24`
- `dataset_pf_opf/ieee39`
- `dataset_pf_opf/uk`
- `dataset_pf_opf/ieee118`
- `dataset_pf_opf/texas`
- `dataset_cascades/ieee24`
- `dataset_cascades/ieee39`
- `dataset_cascades/uk`
- `dataset_cascades/ieee118`
- `raw.7z` for the Texas cascading-failure raw package

Local extraction rule

- unpack the downloaded archive into a repository-external data directory
- keep the extraction path out of the holdout decision loop
- do not modify any frozen split manifest after inspection
- do not add features after holdout inspection

Pre-run checks

- confirm the extracted tree matches the expected layout above
- confirm the split manifests remain unchanged
- confirm the execution contract is still the active one

What this note is not

- not a claim that the dataset has been downloaded here
- not a claim of external predictive validity
- not a result report

Boundary

This note exists to make the next real-world holdout run auditable once the
dataset is present locally.
