"""mtft.surface.blowup — B3 of the parent-action ledger (PV-02, 2026-09-29): the blow-up of X0(143) x 143a1 at an arithmetic point, its
Neron–Severi lattice, the supersymmetric (topological) twist, chirality and slopes on it, and the exact D-flat Standard-Model search.

Lattice.  NS(S), S = X0(143) x E, contains U (+) <-2 delta>: f1 = [{x} x E], f2 = [X x {e}] (f1.f2 = 1, f1^2 = f2^2 = 0), and
gamma = Gamma - delta f1 - f2, the primitive part of the graph Gamma of the modular parametrisation X0(143) -> 143a1 of degree delta = 4
(Gamma^2 = -delta deg K_E = 0, Gamma.f1 = 1, Gamma.f2 = delta), gamma^2 = -2 delta = -8, gamma orthogonal to U.  Blowing up a point p adds the
exceptional curve E with E^2 = -1:  NS(S~) contains  U (+) <-8> (+) <-1>,  an ODD lattice.  A class is written (a, b, c, e) = a f1 + b f2 +
c gamma + e E, so that (a, b) is the package's bidegree (degree a on X, degree b on E), and
    (a,b,c,e).(a',b',c',e') = a b' + a' b - 8 c c' - e e',        K~ = pi* K_S + E = 24 f1 + E = (24, 0, 0, 1),        K~.L = 24 b - e.
K~ is not divisible by 2 (the lattice is odd): the blow-up carries no spin structure, so the spin twist of the parity theorem does not exist
on it and every twist is spin^c.

Twists and chirality.  With spinors (x) R, the net chirality of an adjoint block L is c1(L).(2R - K~) (Riemann–Roch, as in
product_surface.adjoint_net_chirality).  The supersymmetric topological twist (Beasley–Heckman–Vafa) is R = O: fermions are (0,q)-forms
valued in the bundle, and  net(L) = -c1(L).K~ = e - 24 b:  the Euler characteristic of X0(143) times the torus flux, shifted by the
exceptional flux.  On the product (e = 0) it is a multiple of 24; three families need e = 3 mod 24.  net is LINEAR in the class, so for
adjoint blocks L_i L_j^-1 it is a coboundary, net(i, j-bar) = phi_i - phi_j with phi(d) = e - 24 b: every chirality count of an adjoint
parent is a potential difference.  With hypercharge from the stack U(1)s (y_c = y_L + 1/6, u-type singlets y = y_L - 1/2, d-type singlets
y = y_L + 1/2) the exact Standard Model needs  phi_c - phi_L = 3,  sum_A (phi_A - phi_c) = 3,  sum_B (phi_B - phi_c) = 3,
sum_A (phi_A - phi_L) - sum_B (phi_B - phi_L) = 3,  sum_{A,B} (phi_B - phi_A) = 3,  which is consistent iff n_A = n_B + 1; the minimal
content is U(8) = U(3) x U(2) x U(1)_A1 x U(1)_A2 x U(1)_B — the same rank as PV-01's non-supersymmetric B1 class — with
phi_c - phi_L = 3, phi_B - phi_c = 3, phi_A1 + phi_A2 - 2 phi_c = 3.  Every such assignment carries at least six vector-like doublet pairs
(Higgs-type), because the doublet blocks total |phi_A1 + 3| + |6 - phi_A1| + 6 >= 15 for net 3.

Kähler classes and D-flatness.  omega~ = A_X f1 + A_E f2 - eps E (product polarisation pulled back, exceptional curve of area eps); ample iff
0 < eps < min(A_X, A_E): the Seshadri constant of a product polarisation at a point is the smaller fibre area, since a curve of multiplicity k
at p meets each fibre through p at least k times.  Slope: mu(L) = a A_E + b A_X + e eps (gamma is invisible to product classes).  D-flatness
of a split background = all stack slopes equal; the two ratios (r, eps/A_E) can now satisfy two independent conditions where the product's
single ratio zeroed one slope ray (the root of the trilemma).  `d_flat_sm_search` enumerates Kähler normals n = (A_E, A_X, eps) and the
classes orthogonal to them; `collinear_family` is the closed-form solution in which every stack is a power of one line bundle.

Status: EXACT (lattice arithmetic, Riemann–Roch, ampleness); the physics identifications (twist, D-term = slope) are standard and
conditional on a supersymmetric 8D gauge parent.  Cohomology dimensions on the blow-up (h^1, h^0(K (x) L): the Yukawa question) are NOT
computed here — only Euler characteristics are exact."""
import itertools
from fractions import Fraction
import sympy as sp

DELTA = 4                                                    # modular degree of X0(143) -> 143a1
K_TILDE = (24, 0, 0, 1)

def dot(u, v):
    """Intersection form on (a, b, c, e) = a f1 + b f2 + c gamma + e E."""
    u = _pad(u); v = _pad(v)
    return u[0] * v[1] + v[0] * u[1] - 2 * DELTA * u[2] * v[2] - u[3] * v[3]

def _pad(u):
    u = tuple(u)
    if len(u) == 2: return (u[0], u[1], 0, 0)
    if len(u) == 3: return (u[0], u[1], 0, u[2])               # (a, b, e): no gamma component
    return u

def self_intersection(L): return dot(L, L)

def is_spin():
    """K~ = (24, 0, 0, 1) is not 2-divisible: the blow-up is not spin."""
    return all(k % 2 == 0 for k in K_TILDE)

def net_chirality(L, twist=(0, 0, 0, 0)):
    """c1(L).(2R - K~) for the adjoint block L with spinors twisted by R; the SUSY topological twist R = O gives e - 24 b."""
    R = _pad(twist); L = _pad(L)
    return dot(L, tuple(2 * r - k for r, k in zip(R, K_TILDE)))

def phi(d):
    """The chirality potential phi(d) = e - 24 b; net(i, j-bar) = phi(d_i) - phi(d_j) for adjoint blocks."""
    d = _pad(d); return d[3] - 24 * d[1]

def chi(L):
    """Riemann–Roch on the blow-up: chi(O_S~) = chi(O_S) = chi(O_X) chi(O_E) = (1 - 13)(1 - 1) = 0 (Noether: (K^2 + e)/12 = (0 + 0)/12);
    chi(L) = (L^2 - L.K~)/2.  On the product this is (a - 12) b, product_surface.product_block's untwisted index."""
    L = _pad(L); return Fraction(dot(L, L) - dot(L, K_TILDE), 2)

def slope(L, A_X, A_E, eps):
    """mu(L) = a A_E + b A_X + e eps for the Kähler class A_X f1 + A_E f2 - eps E."""
    L = _pad(L); return L[0] * A_E + L[1] * A_X + L[3] * eps

def ample(A_X, A_E, eps):
    """0 < eps < min(A_X, A_E) (Seshadri constant of the product polarisation at a point = the smaller fibre area)."""
    return A_X > 0 and A_E > 0 and 0 < eps < min(A_X, A_E)

def polystable_locus(stacks):
    """All Kähler classes (A_X, A_E, eps) up to scale at which every stack slope is equal, as exact sympy solutions in the ratios
    r = A_X/A_E, t = eps/A_E, with the ampleness verdict.  stacks: ((name, rank, class), ...)."""
    r, t = sp.symbols("r t", positive=True)
    mus = [slope(d, r, 1, t) for _, _, d in stacks]
    eqs = [sp.Eq(m, mus[0]) for m in mus[1:]]
    sols = sp.solve(eqs, [r, t], dict=True) if eqs else [{}]
    out = []
    for s in sols:
        rr, tt = s.get(r, r), s.get(t, t)
        cond = None
        if rr.free_symbols or tt.free_symbols: cond = sp.And(tt > 0, tt < 1, tt < rr)
        else: cond = bool(ample(rr, 1, tt))
        out.append({"r": rr, "eps_over_A_E": tt, "ample": cond})
    return out


# ------------------------------------------------------------------ the supersymmetric Standard-Model conditions (adjoint parent)
def sm_phi_conditions(phis):
    """phis: dict with keys c, L, A1, A2, B (chirality potentials).  Returns the net counts of every SM-type block and the exact-SM verdict."""
    c, L, A1, A2, B = (phis[k] for k in ("c", "L", "A1", "A2", "B"))
    Q = c - L; u = [A1 - c, A2 - c]; d = [B - c]; lep = [A1 - L, A2 - L]; anti_lep = [B - L]; ec = [B - A1, B - A2]; sterile = [A1 - A2]
    net = {"Q": Q, "u^c": sum(u), "d^c": sum(d), "L": sum(lep) - sum(anti_lep), "e^c": sum(ec)}
    pairs = {"u^c": (sum(abs(x) for x in u) - abs(sum(u))) // 2, "e^c": (sum(abs(x) for x in ec) - abs(sum(ec))) // 2,
             "doublets": (sum(abs(x) for x in lep) + sum(abs(x) for x in anti_lep) - abs(net["L"])) // 2, "sterile_singlets": abs(sterile[0])}
    return {"net": net, "exact_SM": all(v == 3 for v in net.values()), "vector_like_pairs": pairs, "blocks": {"u": u, "d": d, "lep": lep, "anti_lep": anti_lep, "ec": ec}}

def sm_solution_space():
    """The general solution of the phi conditions: phi_L = phi_c - 3, phi_B = phi_c + 3, phi_A1 + phi_A2 = 2 phi_c + 3 (one free integer)."""
    pc, pa = sp.symbols("phi_c phi_A1", integer=True)
    return {"phi_L": pc - 3, "phi_B": pc + 3, "phi_A2": 2 * pc + 3 - pa, "free": (pc, pa), "min_doublet_pairs": 6, "note": "n_A = n_B + 1 is forced"}


# ------------------------------------------------------------------ the D-flat search
def classes_orthogonal(n, box_a=4, box_b=1, box_e=30):
    """Integer classes (a, b, e) with a n_a + b n_b + e n_e = 0 for the Kähler normal n = (A_E, A_X, eps) (slopes zero), in the box."""
    na, nb, ne = n; out = []
    for a in range(-box_a, box_a + 1):
        for b in range(-box_b, box_b + 1):
            num = -(a * na + b * nb)
            if num % ne == 0 and abs(num // ne) <= box_e: out.append((a, b, num // ne))
    return out

def d_flat_sm_search(n_box=8, box_a=4, box_b=1, box_e=30, max_pairs=None):
    """Exact finite search for D-flat (polystable) split backgrounds of a supersymmetric U(8) adjoint parent on the blow-up whose net
    chiral content is exactly the Standard Model.  Enumerates primitive integer Kähler normals n = (A_E, A_X, eps) with entries <= n_box
    and 0 < eps < min(A_X, A_E) (ampleness), the classes orthogonal to n (all stack slopes zero, colour fixed at the origin), and the
    assignments d_L (phi = -3), d_B (phi = 3), d_A1, d_A2 (phi sum 3), all distinct and nonzero.  Returns the solutions with their vector-like
    pair counts, sorted by total extras.  EXACT within the boxes."""
    from math import gcd
    sols = []; seen = set()
    for AE in range(1, n_box + 1):
        for AX in range(1, n_box + 1):
            for eps in range(1, min(AX, AE)):
                if gcd(gcd(AE, AX), eps) != 1: continue
                n = (AE, AX, eps); cls = classes_orthogonal(n, box_a, box_b, box_e)
                by_phi = {}
                for d in cls:
                    if d == (0, 0, 0): continue
                    by_phi.setdefault(phi(d), []).append(d)
                for dL in by_phi.get(-3, []):
                    for dB in by_phi.get(3, []):
                        for pA1, lst in by_phi.items():
                            for dA1 in lst:
                                for dA2 in by_phi.get(3 - pA1, []):
                                    if len({dL, dB, dA1, dA2}) < 4 or dA1 > dA2: continue
                                    rec = sm_phi_conditions({"c": 0, "L": -3, "A1": pA1, "A2": 3 - pA1, "B": 3})
                                    assert rec["exact_SM"]
                                    extras = sum(rec["vector_like_pairs"].values())
                                    if max_pairs is not None and extras > max_pairs: continue
                                    key = (n, dL, dB, dA1, dA2)
                                    if key in seen: continue
                                    seen.add(key)
                                    sols.append({"kahler_normal_(A_E,A_X,eps)": n, "r": Fraction(AX, AE), "eps_over_A_E": Fraction(eps, AE),
                                                 "stacks": {"c": (0, 0, 0), "L": dL, "B": dB, "A1": dA1, "A2": dA2}, "phi": {"c": 0, "L": -3, "B": 3, "A1": pA1, "A2": 3 - pA1},
                                                 "vector_like_pairs": rec["vector_like_pairs"], "extras": extras, "collinear": _collinear([dL, dB, dA1, dA2])})
    sols.sort(key=lambda s: (s["extras"], s["kahler_normal_(A_E,A_X,eps)"]))
    return {"n_solutions": len(sols), "solutions": sols, "boxes": {"n_box": n_box, "box_a": box_a, "box_b": box_b, "box_e": box_e}}

def _collinear(vs):
    vs = [tuple(v) for v in vs if any(v)]
    if len(vs) < 2: return True
    v0 = vs[0]
    for v in vs[1:]:
        # cross product zero?
        cx = (v0[1] * v[2] - v0[2] * v[1], v0[2] * v[0] - v0[0] * v[2], v0[0] * v[1] - v0[1] * v[0])
        if any(cx): return False
    return True

def collinear_family(a_L=1):
    """Closed form: with every stack a power of the single line bundle v = (a_L, 0, -3 a_L) — degree a_L on X0(143), exceptional flux -3 a_L —
    and charges k = (c: 0, L: 1, B: -1, A1: 2, A2: -3), the background is D-flat exactly on eps = a_L A_E / 3 for EVERY shape r > a_L/3
    (ampleness), the flux is traceless (sum n_i k_i = 0: an SU(8) background), and the net chirality of a block of relative charge q is
    -3 a_L q: three families for a_L = 1.  Vector-like extras: six doublet pairs, six u^c pairs, six e^c pairs, fifteen sterile singlets."""
    v = (a_L, 0, -3 * a_L); k = {"c": 0, "L": 1, "B": -1, "A1": 2, "A2": -3}; ranks = {"c": 3, "L": 2, "B": 1, "A1": 1, "A2": 1}
    stacks = {name: tuple(kk * x for x in v) for name, kk in k.items()}
    phis = {name: phi(d) for name, d in stacks.items()}
    rec = sm_phi_conditions(phis)
    r, t = sp.symbols("r t", positive=True)
    mus = {name: slope(d, r, 1, t) for name, d in stacks.items()}
    return {"v": v, "charges": k, "stacks": stacks, "traceless": sum(ranks[nm] * kk for nm, kk in k.items()) == 0, "phi": phis, "exact_SM": rec["exact_SM"],
            "vector_like_pairs": rec["vector_like_pairs"], "slopes": mus, "d_flat_locus": sp.solve(sp.Eq(mus["L"], mus["c"]), t)[0], "ample_iff": f"r > {a_L}/3",
            "net_per_unit_charge": phi(v)}

def report(n_box=8):
    fam = collinear_family(); srch = d_flat_sm_search(n_box=n_box)
    best = srch["solutions"][:5]
    return {"blow_up_is_spin": is_spin(), "net_chirality_susy_twist": "e - 24 b", "three_family_block_example": {"class": (-1, 0, 3), "net": net_chirality((-1, 0, 3)), "slope": slope((-1, 0, 3), sp.Symbol("A_X"), sp.Symbol("A_E"), sp.Symbol("eps"))},
            "collinear_family": {k: v for k, v in fam.items() if k in ("v", "charges", "traceless", "exact_SM", "vector_like_pairs", "d_flat_locus", "ample_iff", "net_per_unit_charge")},
            "search": {"n_solutions": srch["n_solutions"], "min_extras": best[0]["extras"] if best else None, "best": best, "boxes": srch["boxes"]}}
