# The superbase of a uniform spanning tree on a torus

Paper and programs for a numerical study of the question of Lyons (Publ. Math. IHES 98 (2003), Section 12) and
Schramm (ICM 2006, Problem 6.2): the law of the homology classes carried by the dual of a uniform spanning tree on a
fine grid approximation of a torus, and its limit as the mesh tends to zero.

- `paper/` LaTeX source and PDF.
- `code/ust_torus_shear.py` Wilson's algorithm on the torus Z^2 / <(n1, 0), (s, n2)> (Numba); the classes of all
  fundamental cycles of each tree. `code/ust_torus_theta.py` the axis-parallel case.
- `code/ust_metal.swift` the same sampler on an Apple GPU, one tree per thread.
  Usage: `swiftc -O ust_metal.swift -o ust_metal && ./ust_metal n1 n2 s trees seed out.tsv`.
- `code/exact_small.py` the exact law on a small torus by listing every spanning tree; `code/vs_exact.py` the
  chi-squared comparison of a sample with it; `code/matrix_tree.py` the spanning-tree counts by the matrix-tree theorem.
- `code/classify.py`, `code/tally.py` tallies of superbases; `code/mesh_table.py` the law by mesh;
  `code/compare.py` two tori against each other, with a PSL2(Z) matrix carrying classes; `code/shape_fit.py`,
  `code/holdout.py`, `code/cand_fit.py` the fits of the weight of a superbase.
- `results/` the exact laws of the three small tori.
- The release `v1.0` carries `samples.tar.gz`: one file per torus with the class set of each of 10^6 trees (80 tori).

Python 3, NumPy, Numba; Swift with Metal on macOS.
