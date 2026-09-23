"""R2C-04 (2026-09-19): where an MTFT-native parent can live — exact representation-theoretic facts.

Irreducible 2-form x 6-form anomaly terms need a cubic Casimir (d_abc != 0): SU(N >= 3) (and U(N)) have one; SU(2), SO(N != 6), Sp(N),
G2, F4, E6, E7, E8 have none.  6D gauge anomalies factorize automatically when there is no independent quartic Casimir (tr F^4 = c (tr F^2)^2):
SU(2), SU(3), G2, F4, E6, E7, E8.  A simple parent from the second list with the three families as ONE complex block of a U(1) inside the
group, with flux degree 3 on S0 (h^0 = 3, h^1 = 0 by purity), realises "one parent, three internal modes":
  E7 -> E6 x U(1): 133 = 78_0 + 1_0 + 27_{+1} + 27bar_{-1} — one complex 27 block, degree 3 -> three 27's of E6 (each 27 = 16 + 10 + 1 of SO(10):
  a family, a Higgs 10 and a singlet).  E6 -> SO(10) x U(1): 78 = 45_0 + 1_0 + 16_{-3} + 16bar_{+3}: charge 3 needs a degree-1 flux (O(P)) for
  three 16's (S0(3P) is pure).  E8 -> E6 x SU(3) or SO(10) x SU(4): the family blocks carry the weights of a traceless factor, so the net
  family number sum_i d_i vanishes — no net families from E8's adjoint on a curve.
The net families of a complex block of charge q under a flux of degree d is q d (index of S0 (x) L^q).  Gravitational anomaly of a pure
gaugino theory: n_grav = dim(adj) must be 0 mod 4 for tensor cancellation (133, 78 fail); the (1,0) supersymmetric condition
H - V + 29 T = 273 gives integer hypermultiplet counts (E7, T = 1: H = 377).  Recorded as facts and gates; the E7 candidate M3 is a
research record, its Higgs sector, vacuum and scale untouched."""
from fractions import Fraction

CUBIC_CASIMIR = {"SU(2)": False, "SU(3)": True, "SU(N>=3)": True, "U(N>=3)": True, "SO(6)": True, "SO(N!=6)": False, "Sp(N)": False, "G2": False, "F4": False, "E6": False, "E7": False, "E8": False}
INDEPENDENT_QUARTIC_CASIMIR = {"SU(2)": False, "SU(3)": False, "SU(N>=4)": True, "SO(N>=7)": True, "SO(10)": True, "Sp(N>=2)": True, "G2": False, "F4": False, "E6": False, "E7": False, "E8": False}
ADJOINT_DIM = {"E6": 78, "E7": 133, "E8": 248, "SO(10)": 45, "SU(5)": 24}

def parent_admissible(group):
    """No cubic Casimir (no irreducible 2x6 terms) and no independent quartic Casimir (gauge anomaly factorizes)."""
    return {"group": group, "no_cubic_casimir": not CUBIC_CASIMIR.get(group, True), "gauge_anomaly_factorizes": not INDEPENDENT_QUARTIC_CASIMIR.get(group, True)}

def complex_block_families(charge, flux_degree): return int(charge) * int(flux_degree)

def e7_three_27_families(flux_degree=3):
    decomposition = {"78_0": 78, "1_0": 1, "27_+1": 27, "27bar_-1": 27}; assert sum(decomposition.values()) == 133
    return {"parent": "E7", "commutant": "E6 x U(1)", "decomposition": decomposition, "net_27_families": complex_block_families(1, flux_degree), "pure": flux_degree == 3,
            "so10_content_per_27": {"16": 1, "10": 1, "1": 1}, "admissible": parent_admissible("E7")}

def e6_three_16_families(flux_degree=1):
    decomposition = {"45_0": 45, "1_0": 1, "16_-3": 16, "16bar_+3": 16}; assert sum(decomposition.values()) == 78
    return {"parent": "E6", "commutant": "SO(10) x U(1)", "decomposition": decomposition, "net_16_families": complex_block_families(3, flux_degree), "flux_bundle": "O(P), degree 1 (S0(3P) pure)", "admissible": parent_admissible("E6")}

def e8_adjoint_net_families(weights_degrees):
    """Family blocks carrying the weights of a traceless factor (SU(3): 3 weights, SU(4): 4) have sum of degrees 0: net families vanish."""
    return {"net_families": sum(weights_degrees), "traceless": sum(weights_degrees) == 0}

def gaugino_p2_integrality(group):
    """Pure-gaugino p2 cancellation by chiral tensors needs dim(adj)/28 net tensors (CC-33: the tensor p2 coefficient is 28 Weyl units, not 4);
    E6 (78), E7 (133) and E8 (248) all fail.  R2C-04's E8 survivor is withdrawn."""
    n = ADJOINT_DIM[group]; need = Fraction(n, 28); return {"group": group, "n_grav": n, "net_tensors_needed": need, "integral": need.denominator == 1, "rule": "dim(adj) = 0 mod 28 (CC-33)"}

def susy_gravitational_condition(V, T=1): return {"V": V, "T": T, "H": 273 - 29 * T + V}

M3_RECORD = {"model_id": "M3", "kind": "6D E7 gauge parent on X0(143) with the E6-commutant U(1) flux L = O(P1+P2+P3)", "families": "three 27's of E6 (one complex block, degree 3, pure)",
             "higgs": "the 10 of SO(10) inside each 27; the 27 block is flux-charged (d = 3), so its vector modes are tachyonic at the compactification scale: 15 modes at m^2 = -1/8 (R2C-05, mesh-verified; CC-30 corrects an earlier -3/8) — the split E6 x U(1) background recombines; no flux-neutral Higgs in E7 -> E6 x U(1) — OPEN/negative",
             "anomalies": "no cubic Casimir; quartic factorizes; gravitational needs a SUSY completion (H = 377 for T = 1) — OPEN", "vacuum": "none", "scale": "none", "status": "research record (R2C-04)"}
