"""INT-07 (AXG-02 §6, AXG-03 §4, AXG-04 §5): flux/scalar gates, tensor transgression, factorization and integral lattices."""
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

def native_scalar_flux_gate(K, m):
    """Native periodic scalars require their charge row to annihilate the flux vector m: K m = 0 rowwise.  M1: K m = (12, 12);
    the difference row k_Delta = k_2 - k_1 = (3, -2, -1, -1, 1) has k_Delta . m = 0 (necessary, not sufficient)."""
    Km = K * sp.Matrix(m); kd = (K.row(1) - K.row(0)) if K.rows >= 2 else None
    return {"K_m": list(Km), "passing_rows": [i for i in range(K.rows) if Km[i] == 0], "k_delta": list(kd) if kd is not None else None, "k_delta_m": (kd * sp.Matrix(m))[0] if kd is not None else None}

def flux_transgression(K):
    """Integral tensor flux reduction induces 3 k_1 (not k_1): K_hat = [k_Delta; 3 k_1], SNF (1, 3) -> candidate Z3."""
    kd = K.row(1) - K.row(0); Kh = sp.Matrix.vstack(kd, 3 * K.row(0)); S = smith_normal_form(Kh, domain=sp.ZZ)
    return {"K_hat": Kh, "invariant_factors": [S[i, i] for i in range(min(S.shape))], "discrete_remnant": [S[i, i] for i in range(min(S.shape)) if S[i, i] not in (0, 1, -1)]}

def factorization_check(I, X, Y):
    """Symbolic identity I = X * Y for anomaly polynomials in commuting characteristic classes."""
    return sp.expand(I - sp.expand(X * Y)) == 0

def c3x_anomaly_polynomial():
    C, W, eta, h, x = sp.symbols("C W eta h x"); X4 = W + 9 * h ** 2; Y4 = 3 * C + W + eta - 27 * h ** 2 - 6 * x ** 2
    return {"symbols": (C, W, eta, h, x), "X4": X4, "Y4": Y4, "I8": sp.expand(X4 * Y4), "I6_flux_reduction": sp.expand(-12 * x * (W + 9 * h ** 2)), "BF_level": -12}

def integral_lattice_gate(Omega, vectors=None):
    """Even unimodular check of the tensor lattice, signature, and Omega-norms of source vectors (integrality witness)."""
    Om = sp.Matrix(Omega); det = Om.det(); even = all(Om[i, i] % 2 == 0 for i in range(Om.rows)); ev = [sp.re(v) for v in Om.eigenvals(multiple=True)]
    sig = (sum(1 for v in ev if v > 0), sum(1 for v in ev if v < 0)); norms = {k: (sp.Matrix(v).T * Om * sp.Matrix(v))[0] for k, v in (vectors or {}).items()}
    return {"det": det, "unimodular": abs(det) == 1, "even": even, "signature": sig, "norms": norms, "all_norms_integral": all(n.is_integer for n in norms.values()) if norms else None}

C3X_LATTICE = {"Omega": [[0, 1], [1, 0]], "vectors": {"a": (0, 2), "b_C": (0, -3), "b_W": (-1, -1), "b_hh": (18, -54), "b_XX": (0, -12), "b_hX": (0, 0)}}
