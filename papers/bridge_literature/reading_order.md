# Bridge Reading Order

If you want the shortest path into the bridge corpus, read these first.

## 1. Spectral structure

1. Ulrike von Luxburg, *A Tutorial on Spectral Clustering*
2. Boaz Nadler, Stephane Lafon, Ronald R. Coifman, Ioannis G. Kevrekidis,
   *Diffusion Maps, Spectral Clustering and Eigenfunctions of Fokker-Planck operators*
3. Ronald R. Coifman and Matthew J. Hirn, *Diffusion maps for changing data*

Why first:

- this is the spectral-address backbone
- it gives the cleanest entry point into normalized Laplacians, Fiedler vectors,
  and diffusion geometry
- it is the easiest place to talk about recoverable structure without drifting
  into physics claims

## 2. Hodge structure

1. Mitchell Black and Amir Nayyeri, *Hodge Decomposition and General Laplacian
   Solvers for Embedded Simplicial Complexes*
2. T. Mitchell Roddenberry et al., *Hodgelets*
3. Xiaoqi Wei and Guo-Wei Wei, *Persistent Topological Laplacians -- a Survey*

Why second:

- this is the discrete sector-split backbone
- it maps cleanly onto exact / coexact / harmonic language
- it is the best route to the “constraints / degrees of freedom / residue” story

## 3. Curvature and flow

1. Jürgen Jost and Shiping Liu, *Ollivier's Ricci curvature, local clustering
   and curvature dimension inequalities on graphs*
2. R. P. Sreejith et al., *Forman curvature for complex networks*
3. Melanie Weber et al., *Characterizing Complex Networks with Forman-Ricci
   Curvature and Associated Geometric Flows*

Why third:

- this is the geometry-engine backbone
- it connects local overlap, bottlenecks, and flow-like smoothing
- it is where the bridge starts to feel physically suggestive without claiming
  anything physical

## 4. Transport and kernel distances

1. Marco Cuturi, *Sinkhorn Distances*
2. Gabriel Peyré and Marco Cuturi, *Computational Optimal Transport*
3. Arthur Gretton et al., *A Kernel Two-Sample Test*

Why fourth:

- these are the numerical comparison tools
- they are useful for turning bridge intuition into measurable diagnostics
- they help prevent the bridge from becoming pure metaphor

## 5. Physics-adjacent discrete geometry

1. Simon Catterall, *Dirac-Kähler fermions and exact lattice supersymmetry*
2. Vivien de Beauce, Samik Sen, James C. Sexton, *Chiral Dirac fermions on the
   lattice using Geometric Discretisation*
3. Roy R. Lederman and Bogdan Toader, *On Manifold Learning in Plato's Cave:
   Remarks on Manifold Learning and Physical Phenomena*

Why last:

- this is the most seductive layer
- it is also the easiest place to overread a metaphor as a derivation
- it is best approached after the structural pieces are already clear

## Practical rule

Read the spectral, Hodge, and curvature papers before you let yourself read
anything that sounds like a physics story.

That order keeps the bridge honest.
