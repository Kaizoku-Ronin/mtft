"""INT-05 (AXG-01 §3): 4D anomaly polynomial of a stack model with line fluxes, and its factorization through the shift matrix.

Conventions (AXG-01): stacks with ranks N_i, line degrees m_i; an oriented bifundamental (i, j-bar) has chiral index
m_i - m_j (charge conjugation accounted for by the orientation).  For U(1) charges q_i on the stacks:
  cubic     P(q)   = sum_{i<j} (m_i - m_j) N_i N_j (q_i - q_j)^3
  gravity   Agr(q) = sum_{i<j} (m_i - m_j) N_i N_j (q_i - q_j)
  mixed SU(N_k)^2:  A_k(q) = (1/2) sum_{j != k} N_j (m_k - m_j) (q_k - q_j)        (T(fund) = 1/2)
A(u,v,w) = tr_chiral(Q_u Q_v Q_w) is obtained by polarisation (directional derivatives / 3, / 2)."""
import itertools, sympy as sp

def stack_model(N, m, order=None):
    order = list(order or N); return {"order": order, "N": [sp.Integer(N[k]) for k in order], "m": [sp.Integer(m[k]) for k in order]}

def anomaly_polynomials(model, q=None):
    N, m, order = model["N"], model["m"], model["order"]; n = len(order); q = list(q) if q is not None else list(sp.symbols("q_0:%d" % n))
    P = sp.expand(sum((m[i] - m[j]) * N[i] * N[j] * (q[i] - q[j]) ** 3 for i in range(n) for j in range(i + 1, n)))
    Agr = sp.expand(sum((m[i] - m[j]) * N[i] * N[j] * (q[i] - q[j]) for i in range(n) for j in range(i + 1, n)))
    mixed = {order[k]: sp.expand(sum(sp.Rational(1, 2) * N[j] * (m[k] - m[j]) * (q[k] - q[j]) for j in range(n) if j != k)) for k in range(n) if N[k] >= 2}
    return {"q": q, "cubic": P, "gravity": Agr, "mixed": mixed}

def directional(p, q, v): return sp.expand(sum(sp.nsimplify(v[i]) * sp.diff(p, q[i]) for i in range(len(q))))

def polarised_trace(P, q, u, v, w):
    """tr_chiral(Q_u Q_v Q_w) from the cubic P: successive directional derivatives divided by 3, 2, 1."""
    return directional(directional(directional(P, q, u), q, v), q, w) / 6

def mixed_anomaly_matrix(model):
    """Rows: gradients of the mixed non-Abelian and gravitational anomalies (linear in q); its row space is spanned by the
    Stueckelberg-charged directions, its kernel by the anomaly-free U(1)s."""
    A = anomaly_polynomials(model); q = A["q"]; rows = [[sp.diff(p, z) for z in q] for p in list(A["mixed"].values()) + [A["gravity"]]]
    return sp.Matrix(rows)

def primitive_shift_matrix(model):
    """Primitive integer basis of the mixed-anomaly row lattice (Hermite normal form, zero rows dropped): the matrix K with
    Da = da + K A for the periodic axions.  M1: K = [[0, -2, 1, 1, 0], [3, -4, 0, 0, 1]] up to row operations."""
    from sympy.matrices.normalforms import hermite_normal_form
    M = mixed_anomaly_matrix(model); den = sp.lcm([sp.fraction(x)[1] for x in M]); Mi = (M * den).applyfunc(sp.Integer)
    Hn = hermite_normal_form(Mi.T).T if Mi.rank() else Mi
    rows = [list(Hn.row(i)) for i in range(Hn.rows) if any(Hn.row(i))]
    rows = [[x // sp.gcd(row) for x in row] for row in (list(map(sp.Integer, r)) for r in rows)]
    return sp.Matrix(rows)

def kernel_directions(K):
    """Rational kernel of K (anomaly-free U(1) directions), as column vectors."""
    return K.nullspace()
