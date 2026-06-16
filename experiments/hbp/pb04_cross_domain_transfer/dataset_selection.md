# PB-04 Dataset Selection

PB-04 is frozen as a contract-first branch. This note freezes the current
cross-domain candidate direction so the next implementation step can be written
against a concrete source/target pairing.

## Selected Candidate Pair

- Source domain candidate: `experiments/biology/gene_network_demo`
- Target domain candidate: `/DATA/Powergraph/dataset_pf_opf`

## Why this candidate

- It is the clearest in-repo biology-line source artifact currently available.
- It gives a genuinely different target domain with real infrastructure data.
- It keeps the transfer question separate from the PB-02 holdout path.

## Limits

- The biology-side artifact is a synthetic gene-network demo, not biological
  tissue.
- This is therefore a proxy source domain, not a claim of real tissue transfer.
- The branch remains negative-until-proven.

## Freeze Statement

PB-04 will use the biology-line demo as the current source candidate and
PowerGraph PF/OPF data as the current target candidate unless a later explicit
precommitment supersedes this note.

