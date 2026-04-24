# 3D Tissue Demo

This is Biology Line C.

It is an experimental application layer on top of frozen HAOS-IIP. It models a simple 3D tissue-like interaction structure as a deterministic cell graph and tests early recoverability loss under local lesion perturbation.

The model is a toy 8 x 8 x 8 multicellular lattice. It does not modify HAOS core and makes no claims about real tissue biology.

Run from the repository root:

```bash
python experiments/biology/tissue_3d_demo/run_tissue_3d_demo.py
```

Outputs are written to `experiments/biology/tissue_3d_demo/outputs/`.
