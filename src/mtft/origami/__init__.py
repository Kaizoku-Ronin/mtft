"""mtft.origami — dimers, origami/t-embeddings, and the observable insertion calculus.

New in v0.20.0.  Three layers, all reusable:

  dimer      DimerGraph: APMs, boundary measurement, Plucker vector,
             Grassmann-Plucker certification, coarse-graining,
             ENSEMBLE CONSERVATION gate
  insertion  cumulants, Fisher metric, cubic Amari tensor, Brioschi and
             cumulant (E2) curvature, POTENTIAL PATH INDEPENDENCE gate
  perfect    Theta involution (9.2-9.3), perfect-system solver, cyclic
             orbit structure, equivariant Kasteleyn factorization
  instances  the certified (2,4) and (3,6) instances
  gates      the full battery, runnable end to end

Reference: P. Galashin, arXiv:2410.09574v2.
"""
from mtft.origami.dimer import DimerGraph, ensemble_conservation
from mtft.origami.insertion import (
    D_log, cumulants, fisher_metric, cubic_tensor,
    brioschi_curvature, cumulant_curvature, fisher_pack, path_independence,
)
from mtft.origami.perfect import (
    bracket, t_coefficients, Theta, winding, is_valid_pair,
    lambda_from_annihilator, perfect_residual, solve_perfect_branches,
    cyclic_matrix, orbit_structure, equivariant_kasteleyn_factor,
)
from mtft.origami.instances import (
    galashin_24, prism_36, PRISM_C, PRISM_LAMBDA0,
    t_embedding_24, mandelstams_24, closed_curvature_B,
)

__all__ = [
    "DimerGraph", "ensemble_conservation",
    "D_log", "cumulants", "fisher_metric", "cubic_tensor",
    "brioschi_curvature", "cumulant_curvature", "fisher_pack",
    "path_independence",
    "bracket", "t_coefficients", "Theta", "winding", "is_valid_pair",
    "lambda_from_annihilator", "perfect_residual", "solve_perfect_branches",
    "cyclic_matrix", "orbit_structure", "equivariant_kasteleyn_factor",
    "galashin_24", "prism_36", "PRISM_C", "PRISM_LAMBDA0",
    "t_embedding_24", "mandelstams_24", "closed_curvature_B",
]
