# External Bibliography

This is a curated corpus of the most relevant external papers for HAOS-IIP's
spectral, Hodge, curvature, transport, and bridge-adjacent machinery.

The list is intentionally opinionated: it prioritizes papers that match the
current frozen operator stack and the physics-bridge discussion.

## 1. Spectral clustering and graph Laplacians

- Ulrike von Luxburg, [A Tutorial on Spectral Clustering](https://arxiv.org/abs/0711.0189) (2007). Core reference for normalized Laplacians, Fiedler vectors, and spectral partitioning.
- Boaz Nadler, Stephane Lafon, Ronald R. Coifman, Ioannis G. Kevrekidis, [Diffusion Maps, Spectral Clustering and Eigenfunctions of Fokker-Planck operators](https://arxiv.org/abs/math/0506090) (2005). Probabilistic interpretation of Laplacian eigenvectors and diffusion geometry.
- Boaz Nadler, Stephane Lafon, Ronald R. Coifman, Ioannis G. Kevrekidis, [Diffusion maps, spectral clustering and reaction coordinates of dynamical systems](https://arxiv.org/abs/math/0503445) (2005). Reaction-coordinate and dynamical-systems view of diffusion geometry.
- Ronald R. Coifman, Matthew J. Hirn, [Diffusion maps for changing data](https://arxiv.org/abs/1209.0245) (2012). Embedding evolution under parameter changes.
- Jean Gallier, [Notes on Elementary Spectral Graph Theory. Applications to Graph Clustering Using Normalized Cuts](https://arxiv.org/abs/1311.2492) (2013). Detailed normalized-cut and Laplacian notes.
- David K Hammond, Pierre Vandergheynst, Rémi Gribonval, [Wavelets on Graphs via Spectral Graph Theory](https://arxiv.org/abs/0912.3848) (2009). Graph wavelets from Laplacian spectra.
- Zhana Kuncheva, Giovanni Montana, [Multi-scale Community Detection in Temporal Networks Using Spectral Graph Wavelets](https://arxiv.org/abs/1708.04060) (2017). Multi-scale clustering with spectral wavelets.

## 2. Hodge decomposition, simplicial complexes, and higher-order spectra

- Mitchell Black, Amir Nayyeri, [Hodge Decomposition and General Laplacian Solvers for Embedded Simplicial Complexes](https://arxiv.org/abs/2205.02134) (2022). Practical Hodge decomposition and higher-order Laplacian solving.
- T. Mitchell Roddenberry, Florian Frantzen, Michael T. Schaub, Santiago Segarra, [Hodgelets: Localized Spectral Representations of Flows on Simplicial Complexes](https://arxiv.org/abs/2109.08728) (2021). Localized higher-order spectral representations that mirror the Hodge split.
- Shiquan Ren, Chengyuan Wu, Jie Wu, [Hodge Decompositions for Weighted Hypergraphs](https://arxiv.org/abs/1805.11331) (2018). Weighted generalization of the Hodge sector split.
- Carla Farsi, Elizabeth Gillaspy, Sooran Kang, Judith Packer, [Wavelets and graph C*-algebras](https://arxiv.org/abs/1601.00061) (2016). Useful for higher-rank graph wavelet language adjacent to Hodge structures.
- Simon Catterall, [Dirac-Kähler fermions and exact lattice supersymmetry](https://arxiv.org/abs/hep-lat/0509136) (2005). Strong lattice/Dirac-Kähler bridge tied to cochain-like discretizations.
- S. I. Kruglov, [Dirac-KÄhler Equation](https://arxiv.org/abs/hep-th/0110251) (2001). Dirac-Kähler formulation and symmetry discussion.
- S. I. Kruglov, [Dirac-Kähler Equation (Review)](https://arxiv.org/abs/hep-th/0110060) (2001). Review of the Dirac-Kähler formalism.
- Vivien de Beauce, Samik Sen, James C. Sexton, [Chiral Dirac fermions on the lattice using Geometric Discretisation](https://arxiv.org/abs/hep-lat/0309167) (2003). Discrete forms and lattice fermions.
- Xiaoqi Wei, Guo-Wei Wei, [Persistent Topological Laplacians -- a Survey](https://arxiv.org/abs/2312.07563) (2023). Survey of persistent Laplacians across simplicial complexes, hypergraphs, sheaves, and chains.
- Jian Liu, Jingyan Li, Jie Wu, [The algebraic stability for persistent Laplacians](https://arxiv.org/abs/2302.03902) (2023). Stability theory for persistent Laplacians.

## 3. Discrete curvature and curvature flow on graphs

- Jürgen Jost, Shiping Liu, [Ollivier's Ricci curvature, local clustering and curvature dimension inequalities on graphs](https://arxiv.org/abs/1103.4037) (2011). Connects graph curvature to neighborhood overlap and clustering.
- R. P. Sreejith, Karthikeyan Mohanraj, Jürgen Jost, Emil Saucan, Areejit Samal, [Forman curvature for complex networks](https://arxiv.org/abs/1603.00386) (2016). Fast combinatorial curvature for networks.
- Melanie Weber, Emil Saucan, Jürgen Jost, [Characterizing Complex Networks with Forman-Ricci Curvature and Associated Geometric Flows](https://arxiv.org/abs/1607.08654) (2016). Curvature flow as a graph-analysis tool.
- Pim van der Hoorn, William J. Cunningham, Gabor Lippner, Carlo Trugenberger, Dmitri Krioukov, [Ollivier curvature convergence in random geometric graphs](https://arxiv.org/abs/2008.01209) (2020). Continuum convergence result for graph curvature.
- Pim van der Hoorn, Gabor Lippner, Carlo Trugenberger, Dmitri Krioukov, [Ollivier curvature of random geometric graphs converges to Ricci curvature of their Riemannian manifolds](https://arxiv.org/abs/2009.04306) (2020). Stronger convergence result linking random geometric graphs to smooth Ricci curvature.

## 4. Optimal transport, Sinkhorn, and sliced Wasserstein

- Marco Cuturi, [Sinkhorn Distances: Lightspeed Computation of Optimal Transportation Distances](https://arxiv.org/abs/1306.0895) (2013). Canonical entropic OT / Sinkhorn reference.
- Jason Altschuler, Jonathan Weed, Philippe Rigollet, [Near-linear time approximation algorithms for optimal transport via Sinkhorn iteration](https://arxiv.org/abs/1705.09634) (2017). Complexity and approximation guarantees for Sinkhorn-style OT.
- Simone Di Marino, Augusto Gerolin, [Optimal Transport losses and Sinkhorn algorithm with general convex regularization](https://arxiv.org/abs/2007.00976) (2020). Generalized Sinkhorn framework.
- Yihe Dong, Yu Gao, Richard Peng, Ilya Razenshteyn, Saurabh Sawlani, [A Study of Performance of Optimal Transport](https://arxiv.org/abs/2005.01182) (2020). Practical OT algorithm comparisons.
- Jiqing Wu et al., [Sliced Wasserstein Generative Models](https://arxiv.org/abs/1706.02631) (2017). Sliced Wasserstein approximation route.
- Gabriel Peyré, Marco Cuturi, [Computational Optimal Transport](https://arxiv.org/abs/1803.00567) (2018). Canonical numerical OT reference.
- Luis Caicedo Torres, Luiz Manella Pereira, M. Hadi Amini, [A Survey on Optimal Transport for Machine Learning: Theory and Applications](https://arxiv.org/abs/2106.01963) (2021). High-level OT survey for ML use cases.
- Eduardo Fernandes Montesuma, Fred Ngolè Mboula, Antoine Souloumiac, [Recent Advances in Optimal Transport for Machine Learning](https://arxiv.org/abs/2306.16156) (2023). Recent OT survey.
- Young-Heon Kim, Brendan Pass, [Wasserstein Barycenters over Riemannian manifolds](https://arxiv.org/abs/1412.7726) (2014). Multi-measure Wasserstein geometry on manifolds.
- Sebastian Claici, Edward Chien, Justin Solomon, [Stochastic Wasserstein Barycenters](https://arxiv.org/abs/1802.05757) (2018). Practical barycenter computation.
- Ethan Anderes, Steffen Borgwardt, Jacob Miller, [Discrete Wasserstein Barycenters: Optimal Transport for Discrete Data](https://arxiv.org/abs/1507.07218) (2015). Discrete barycenter theory.

## 5. Kernel distances and two-sample testing

- Arthur Gretton et al., [A Kernel Two-Sample Test](https://jmlr.org/papers/v13/gretton12a.html) (2012). Canonical MMD / kernel two-sample test reference.

## 6. Physics-adjacent discrete geometry and lattice analogies

- On the lattice-geometry side, the Dirac-Kähler / staggered-fermion bridge is the main reason this folder keeps a physics-adjacent section at all.
  If you want the shortest path, follow Catterall, Kruglov, and the geometric
  discretization papers first.

## 7. Cautionary geometry / manifold-learning context

- Roy R. Lederman, Bogdan Toader, [On Manifold Learning in Plato's Cave: Remarks on Manifold Learning and Physical Phenomena](https://arxiv.org/abs/2304.14248) (2023). A useful warning that measurement geometry and source geometry can differ.

## 8. Optional broader context

- The classic Dirac–Kähler / staggered-fermion relationship is also summarized in standard references such as Kähler, Becher & Joos, and lattice QCD texts.

## Notes

- This is a curated bridge corpus, not a claim layer.
- Several of these papers are mathematically adjacent rather than directly
  causal for HAOS-IIP.
- The point is to keep the bridge conversation anchored to real literature.
