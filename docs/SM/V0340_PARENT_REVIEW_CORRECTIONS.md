# Candidate review corrections to v0.34.0

These changes correct documented search scope and generic API edge cases. They do not add a viable parent or expand the B1 search.

- `b1_scan` imposes colour product `a_c*b_c=0` and two individual mixed antitriplet blocks of index `-3`. These conditions define its search domain. Its 60/156 accepted assignments and finite trilemma concern that domain. Exact SM net content does not force the colour product to vanish: the regression counterexample has colour class `(1,1)`, index-level SM content, and degrees up to 12.
- An unsupported representation in `ParentCandidate.content` now makes the 8D anomaly verdict OPEN with its name recorded. It is not silently omitted.
- Identically zero colour/Higgs slope conditions and equal-stack blow-up slopes represent a whole positive-r or ample-cone locus.
- A nonzero scaling of the collinear blow-up flux leaves its D-flat locus at `eps/A_E=1/3` and ampleness threshold `r>1/3`. The default family has nine doublet pairs in the index ledger.
- Weak-triplet labels survive SM net-content normalisation.

The Standard-Model entries are net chiral indices and index-forced extra counts. Actual cohomology, Yukawa couplings, the global Green–Schwarz construction, and physical masses remain open.
