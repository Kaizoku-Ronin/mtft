"""INT-06 (AXG-01 §§4, 6; AXG-03 §4): shift kernel, Smith remnants, integer dressings and kinetic-normalised vector masses.

For a shift matrix K (axion charges Da = da + K A): the kernel is the massless U(1) subspace; SNF(K) records finite
remnants; an operator with stack-charge vector q_O is dressable iff q_O = K^T n for INTEGER n (rational solvability is
not enough: axion periodicity); the vector mass matrix with gauge kinetic matrix H and axion kinetic matrix G is
M^2 = H^{-1/2} K^T G K H^{-1/2}, of rank = rank K."""
import sympy as sp

def shift_kernel(K): return {"kernel": K.nullspace(), "rank": K.rank(), "massless_count": K.cols - K.rank()}

def smith_remnant(K):
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(K, domain=sp.ZZ); d = [S[i, i] for i in range(min(S.shape)) if S[i, i] != 0]
    return {"invariant_factors": d, "finite_remnant": [x for x in d if x not in (1, -1)]}

def integer_dressing(K, q_operator):
    """Solve K^T n = q_O over the integers (via Smith form).  Returns (dressable, n or None, rational_solvable)."""
    from sympy.matrices.normalforms import smith_normal_form
    A = K.T; b = sp.Matrix(q_operator)
    try: sol = A.gauss_jordan_solve(b)[0]; rational = True
    except ValueError: return {"dressable": False, "n": None, "rational_solvable": False}
    free = sol.free_symbols; sol0 = sol.subs({s_: 0 for s_ in free})
    # integer search on a bounded box of the free parameters (kernel of K^T is zero for primitive full-row-rank K)
    if not free: return {"dressable": all(x.is_integer for x in sol0), "n": sol0 if all(x.is_integer for x in sol0) else None, "rational_solvable": rational}
    import itertools
    for vals in itertools.product(range(-6, 7), repeat=len(free)):
        cand = sol.subs(dict(zip(free, vals)))
        if all(x.is_integer for x in cand): return {"dressable": True, "n": cand, "rational_solvable": True}
    return {"dressable": False, "n": None, "rational_solvable": True}

def canonical_vector_masses(K, H, G):
    """Eigenvalues of M^2 = H^{-1/2} K^T G K H^{-1/2} (exact), for positive diagonal H (gauge kinetic) and G (axion kinetic)."""
    Hm = sp.Matrix(H); Gm = sp.Matrix(G); Hs = sp.diag(*[1 / sp.sqrt(Hm[i, i]) for i in range(Hm.rows)])
    M2 = (Hs * K.T * Gm * K * Hs).applyfunc(sp.simplify); ev = M2.eigenvals()
    return {"M2": M2, "eigenvalues": {sp.simplify(k): v for k, v in ev.items()}, "rank": M2.rank()}
