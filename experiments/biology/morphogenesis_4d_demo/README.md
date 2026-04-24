# 4D Morphogenesis Demo

This is Biology Line D.

It is an experimental application layer on top of frozen HAOS-IIP. It models 4D morphogenesis as a 3D tissue state evolving through developmental time and tests early recoverability loss under time-dependent developmental perturbation.

The model is a toy deterministic 8 x 8 x 8 cell lattice with 24 developmental steps. It does not modify HAOS core and makes no claims about real developmental biology.

Run from the repository root:

```bash
python experiments/biology/morphogenesis_4d_demo/run_morphogenesis_4d_demo.py
```

Outputs are written to `experiments/biology/morphogenesis_4d_demo/outputs/`.
