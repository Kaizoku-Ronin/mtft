# mtft v0.30.2 — SM-03: the up-type Yukawa tensor of M1 (2026-09-13)

`surface.rrspace` gains: reduction with composite matrices and exact derivatives (chain rule, checked against
finite differences to 4e-8), eta-quotient q-series with pole offsets, f_K = g1^2 w (weight-1 CM form g1 frozen
to q^400 in `_data/x0143_weight1_g1.json`), the 18-dimensional Higgs target space H^0(K(2 SigmaP)) from cubic
products with double vanishing at the 20 W143 points and P4 (rank 42 certified), and `up_yukawa_M1`: the
3 x 3 x 18 tensor with expansion residual 3e-12 — the products of the three-family sections lie in f_K times the
target space as the divisor bookkeeping demands.  Every Higgs direction gives a rank-3 mass matrix.
Open: W13 grading (flavour texture), compact-metric normalisation, down/lepton tensors, Higgs direction.
Full fast suite before build; clean tree.  Base: PyPI 0.30.0 + unpushed 0.30.1.  Pin four-way.
