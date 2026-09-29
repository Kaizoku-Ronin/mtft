"""PV-01 (2026-09-29): the ten-form anomaly polynomial of chiral eight-dimensional field content, exactly — the gate the
parent-action ledger (docs/SM/PARENT_ACTION_REQUIREMENTS.md, class B1) asks for — and the block/index census of complex
representations reduced on the split backgrounds of X0(143) x 143a1.

Conventions.  Per positive-chirality Weyl fermion  I_10 = [A-roof(T) ch(F)]_10,  A-roof = 1 - p1/24 + (7 p1^2 - 4 p2)/5760;
multiplicity and chirality enter as a signed count, and a negative-chirality fermion in R equals a positive one in R* (ch_odd flips).
Gauge data enter through the Chern character: for a U(N) factor A in the power sums  s_{A,k} = tr_fund (F/2 pi)^k  (ch_k(fund) = s_k/k!),
for a U(1) of charge q through q^k c1^k/k!.  Representations are built by the lambda-ring rules on the truncated Chern character:
  ch(V (+) W) = ch V + ch W,   ch(V (x) W) = ch V ch W,   ch(V*): odd degrees negated,
  ch(Lambda^2 V) = (ch(V)^2 - psi^2 ch V)/2,   ch(Sym^2 V) = (ch(V)^2 + psi^2 ch V)/2,   psi^2 ch_k = 2^k ch_k (Adams operation).
Hence  tr_{Lambda^2} F^5 = (N - 16) tr F^5 + 5 tr F tr F^4 + 10 tr F^2 tr F^3  (the 8D analogue of the 10D rule (N - 32)), and
tr_{Sym^2} F^5 = (N + 16) tr F^5 + 5 tr F tr F^4 + 10 tr F^2 tr F^3.  Real representations (adjoint, V (+) V*) have ch_odd = 0 and
contribute nothing; an 8D theory with fermions only in real representations is anomaly-free but its reduction on the spin surface
has even net chirality (product_surface.chirality_parity_theorem) — chirality on X0(143) x 143a1 costs complex representations,
and this module computes the cost.

Classification.  The primitive quintic s_5 of a U(N)/SU(N) factor with N >= 5 is IRREDUCIBLE: no Green–Schwarz counterterm removes
it.  Every other ten-form monomial is a product of lower invariants; an 8D Green–Schwarz mechanism cancels a factorised piece
X_2 X_8 (axion) or X_4 X_6 (2-form or 4-form) per form field.  For SU(N) (s_1 = 0) the whole reducible part is  s_3 (alpha s_2 + beta p1):
one 2-form with X_6 = tr F^3 suffices whenever the quintic coefficient vanishes — the cubic Casimir that obstructs 6D (R2C-04) is the
Green–Schwarz factor in 8D.  For N <= 4 the quintic is not independent (SU(2): s_3 = s_5 = 0; SU(3), SU(4): s_5 = 5 s_2 s_3 / 6);
`specialise` substitutes explicit Chern roots before deciding.  There is no pure gravitational anomaly in 8D (A-roof has no 10-form),
so gravitini and neutral tensors do not enter.

Status: EXACT (symbolic identities in commuting characteristic classes).  Global anomalies (Omega_9^Spin) and the choice of parent are
separate gates."""
import itertools
from fractions import Fraction
import sympy as sp

P1, P2 = sp.symbols("p1 p2")
DEG = 5                                                   # keep ch_0 .. ch_5 (ten-form = degree 5 in 2-forms)
A_ROOF = {0: sp.Integer(1), 2: -P1 / 24, 4: (7 * P1 ** 2 - 4 * P2) / 5760}      # by degree in 2-forms


class Ch:
    """Truncated Chern character: dict degree -> class (degree 0 = the rank)."""
    def __init__(self, parts): self.c = {k: sp.expand(parts.get(k, 0)) for k in range(DEG + 1)}
    def __add__(self, o): return Ch({k: self.c[k] + o.c[k] for k in self.c})
    def __mul__(self, o): return Ch({k: sum(self.c[i] * o.c[k - i] for i in range(k + 1)) for k in self.c})
    def scale(self, n): return Ch({k: n * v for k, v in self.c.items()})
    def dual(self): return Ch({k: (-1) ** k * v for k, v in self.c.items()})
    def adams2(self): return Ch({k: 2 ** k * v for k, v in self.c.items()})
    def antisym2(self): return (self * self + self.adams2().scale(-1)).scale(sp.Rational(1, 2))
    def sym2(self): return (self * self + self.adams2()).scale(sp.Rational(1, 2))
    def __repr__(self): return "Ch(" + ", ".join(f"{k}: {v}" for k, v in self.c.items() if v != 0) + ")"


def s(label, k): return sp.Symbol(f"s_{label}{k}")
def rank(label): return sp.Symbol(f"N_{label}")

def fund(label, N=None):
    """Fundamental of the U(N) factor `label`; N symbolic unless given."""
    return Ch({0: rank(label) if N is None else sp.Integer(N), **{k: s(label, k) / sp.factorial(k) for k in range(1, DEG + 1)}})

def antifund(label, N=None): return fund(label, N).dual()
def adjoint(label, N=None): return fund(label, N) * antifund(label, N)
def antisym2(label, N=None): return fund(label, N).antisym2()
def sym2(label, N=None): return fund(label, N).sym2()
def u1(q, label="c"):
    c1 = sp.Symbol(f"c1_{label}"); return Ch({k: sp.Integer(q) ** k * c1 ** k / sp.factorial(k) for k in range(DEG + 1)})
def tensor(*reps):
    out = Ch({0: 1})
    for r in reps: out = out * r
    return out


def weyl_I10(ch):
    """Ten-form anomaly of one positive-chirality Weyl fermion with Chern character `ch`."""
    return sp.expand(A_ROOF[0] * ch.c[5] + A_ROOF[2] * ch.c[3] + A_ROOF[4] * ch.c[1])

def I10(spectrum):
    """spectrum: iterable of (Ch, signed multiplicity).  Returns the total ten-form (expanded)."""
    return sp.expand(sum(sp.sympify(n) * weyl_I10(ch) for ch, n in spectrum))


def form_degree(monomial):
    d = 0
    for sym, e in monomial.as_powers_dict().items():
        name = str(sym)
        if name.startswith("s_"): d += 2 * int(name[-1]) * e
        elif name == "p1": d += 4 * e
        elif name == "p2": d += 8 * e
        elif name.startswith("c1_"): d += 2 * e
    return d

def classify(I):
    """Split a ten-form into its irreducible quintic parts (s_{A,5} for every factor A) and the reducible remainder, and try to write the
    remainder as X_a X_b (a single Green–Schwarz field)."""
    I = sp.expand(I); terms = I.as_coefficients_dict() if I != 0 else {}
    irreducible = {}; reducible = sp.Integer(0)
    for mon, coeff in terms.items():
        names = [str(x) for x in mon.free_symbols]
        if len(names) == 1 and names[0].startswith("s_") and names[0].endswith("5") and mon.as_powers_dict()[list(mon.free_symbols)[0]] == 1:
            irreducible[names[0]] = coeff
        else: reducible += coeff * mon
    fac = sp.factor_list(reducible) if reducible != 0 else (0, [])
    factors = [(sp.expand(f), e, form_degree(sp.expand(f).as_ordered_terms()[0])) for f, e in fac[1]] if reducible != 0 else []
    gs = None
    if reducible != 0:
        degs = sorted(d for _, e, d in factors for _ in range(e))
        if len(degs) >= 2 and sum(degs) == 10 and all(d in (2, 4, 6, 8) for d in degs):
            gs = degs
    return {"I10": I, "anomaly_free": I == 0, "irreducible_quintic": irreducible, "reducible": sp.expand(reducible),
            "reducible_factors": [(f, e, d) for f, e, d in factors], "green_schwarz_single_field": gs}


def specialise(I, ranks, su=True):
    """Substitute explicit Chern roots for every U(N) factor in `ranks` (label -> N); with su=True impose tracelessness x_N = -sum x_i.
    Returns the polynomial in the roots (zero iff the ten-form vanishes for that group)."""
    subs = {}
    for label, N in ranks.items():
        xs = list(sp.symbols(f"x_{label}1:{N + 1}"))
        if su: xs[-1] = -sum(xs[:-1])
        for k in range(1, DEG + 1): subs[s(label, k)] = sum(x ** k for x in xs)
        subs[rank(label)] = N
    return sp.expand(I.subs(subs))


# ------------------------------------------------------------------ the census of SU(N) / U(N) content in fund, Lambda^2, Sym^2
def sun_content_conditions(N=None):
    """For n_F fundamentals, n_A antisymmetric and n_S symmetric Weyl fermions (net, signed) of one U(N) factor: the ten-form in the
    s-basis and the exact linear conditions.  For SU(N >= 5): quintic  n_F + (N - 16) n_A + (N + 16) n_S = 0 (irreducible; must vanish);
    the reducible part is s_3 [ (n_A + n_S) s_2 / 12 - (n_F + (N - 4) n_A + (N + 4) n_S) p1 / 144 ].  Exact cancellation of everything
    forces n_F = n_A = n_S = 0: there is no anomaly-free chiral SU(N) content in these representations without Green–Schwarz."""
    nF, nA, nS = sp.symbols("n_F n_A n_S"); Nn = rank("A") if N is None else sp.Integer(N)
    I = I10([(fund("A", N), nF), (antisym2("A", N), nA), (sym2("A", N), nS)])
    Isu = sp.expand(I.subs(s("A", 1), 0))
    quintic = sp.expand(Isu.coeff(s("A", 5)))
    red = sp.expand(Isu - quintic * s("A", 5))
    alpha = sp.expand(red.coeff(s("A", 3)).coeff(s("A", 2))); beta = sp.expand(red.coeff(s("A", 3)).coeff(P1))
    full = sp.solve([quintic, alpha, beta], [nF, nA, nS], dict=True)
    return {"N": Nn, "I10_U(N)": I, "I10_SU(N)": Isu, "quintic_coefficient": quintic, "reducible": red,
            "reducible_is_s3_times_X4": sp.expand(red - s("A", 3) * (alpha * s("A", 2) + beta * P1)) == 0,
            "X4": alpha * s("A", 2) + beta * P1, "exact_cancellation_solutions": full,
            "green_schwarz_condition": sp.Eq(quintic, 0)}

def sun_gs_census(N, box=3):
    """For SU(N), N >= 5: every (n_A, n_S) in [-box, box]^2 not both zero, with the fundamental count n_F = -(N - 16) n_A - (N + 16) n_S that
    cancels the irreducible quintic; the remainder s_3 X_4 is Green–Schwarz-completable with one 2-form.  EXACT."""
    nF, nA, nS = sp.symbols("n_F n_A n_S"); cond = sun_content_conditions(N); out = []
    for b, c in itertools.product(range(-box, box + 1), repeat=2):
        if (b, c) == (0, 0): continue
        a = -(N - 16) * b - (N + 16) * c; sub = {nF: a, nA: b, nS: c}
        assert cond["quintic_coefficient"].subs(sub) == 0
        out.append({"n_F": a, "n_A": b, "n_S": c, "X4": sp.expand(cond["X4"].subs(sub)), "X6": s("A", 3)})
    return out

def gs_decomposition(I, label="A"):
    """Split a U(N) ten-form into the part proportional to s_1 = c_1 (cancellable by an axion with X_2 = s_1) and the rest, and factor each:
    the number of Green–Schwarz fields a spectrum needs is the number of nonzero factorised pieces (plus the irreducible quintic, which
    none can remove)."""
    I = sp.expand(I); s1 = s(label, 1); s5 = s(label, 5)
    quintic = I.coeff(s5).subs(s1, 0) if I.has(s5) else sp.Integer(0)
    rest = sp.expand(I - quintic * s5)
    axion = sp.expand(sum(cf * m for m, cf in rest.as_coefficients_dict().items() if m.has(s1)))
    two_form = sp.expand(rest - axion)
    def fac(e): return [(sp.expand(f), k) for f, k in sp.factor_list(e)[1]] if e != 0 else []
    return {"irreducible_quintic": quintic, "axion_piece": axion, "axion_X8": sp.expand(axion / s1) if axion != 0 else 0, "axion_factors": fac(axion),
            "two_form_piece": two_form, "two_form_factors": fac(two_form),
            "green_schwarz_fields_needed": int(axion != 0) + int(two_form != 0), "cancellable": quintic == 0}


# ------------------------------------------------------------------ B1: complex representations reduced on the split backgrounds
STACK_MODELS = {"S1": (("c", 3, (0, 0)), ("L", 2, (-3, -1)), ("1", 1, (1, -3))),
                "S2": (("c", 3, (0, 0)), ("L", 2, (3, -1)), ("1", 1, (1, 3)))}
SM_NAMES = {("c", "L"): "(3,2)", ("c", "1"): "(3,1)", ("L", "1"): "(1,2)", ("c", "c"): "Lambda^2 3 = 3bar / Sym^2 3 = 6", ("L", "L"): "Lambda^2 2 = 1 / Sym^2 2 = 3", ("1", "1"): "singlet"}

def spin_twisted_index(bidegree):
    """chi(S, K^{1/2} (x) L) for L of bidegree (a, b) on X0(143) x 143a1: chi(X, S0(a)) chi(E, N_b) = a b (VI.1)."""
    a, b = bidegree; return a * b

def stack_blocks(rep, model="S2", stacks=None):
    """Decompose one complex representation of U(N), N = sum n_i, under the split background (+) L_i (x) C^{n_i}:
    rep in {'fund', 'antifund', 'antisym2', 'sym2', 'adjoint'}.  Each block is (gauge content, line bundle bidegree, multiplicity)
    with its spin-twisted index a b and parity.  Blocks of Lambda^2 / Sym^2 of a single stack carry L_i^2, index 4 a_i b_i: EVEN."""
    stacks = stacks if stacks is not None else STACK_MODELS[model]
    add = lambda d, e: (d[0] + e[0], d[1] + e[1]); neg = lambda d: (-d[0], -d[1])
    rows = []
    def row(content, bid, mult, origin):
        idx = spin_twisted_index(bid); rows.append({"content": content, "bidegree": bid, "multiplicity": mult, "index": idx, "odd": idx % 2 == 1, "origin": origin})
    if rep == "fund":
        for n, r, d in stacks: row(f"{n}: fundamental of U({r})", d, 1, "L_i")
    elif rep == "antifund":
        for n, r, d in stacks: row(f"{n}: antifundamental of U({r})", neg(d), 1, "L_i^-1")
    elif rep in ("antisym2", "sym2"):
        for (ni, ri, di), (nj, rj, dj) in itertools.combinations(stacks, 2):
            row(f"({ni},{nj}) bifundamental", add(di, dj), 1, "L_i L_j")
        for n, r, d in stacks:
            if rep == "antisym2" and r >= 2: row(f"{n}: Lambda^2 of U({r}) (dim {r * (r - 1) // 2})", add(d, d), 1, "L_i^2")
            if rep == "sym2": row(f"{n}: Sym^2 of U({r}) (dim {r * (r + 1) // 2})", add(d, d), 1, "L_i^2")
    elif rep == "adjoint":
        for (ni, ri, di), (nj, rj, dj) in itertools.permutations(stacks, 2):
            row(f"({ni},{nj}bar) bifundamental", add(di, neg(dj)), 1, "L_i L_j^-1")
    else: raise ValueError(rep)
    return rows

def b1_index_table(model="S2"):
    """Every SM-type block available from fundamentals, antifundamentals, Lambda^2 and Sym^2 of U(6) on the S1/S2 split background
    with its spin-twisted index; the odd-index blocks are the only candidates for a family count of 3."""
    out = {}
    for rep in ("fund", "antifund", "antisym2", "sym2", "adjoint"): out[rep] = stack_blocks(rep, model)
    odd = [(rep, r["content"], r["bidegree"], r["index"]) for rep, rows in out.items() for r in rows if r["odd"]]
    three = [t for t in odd if abs(t[3]) == 3]
    return {"model": model, "blocks": out, "odd_index_blocks": odd, "index_pm3_blocks": three}

def split_stack_parity_lemma(box=4):
    """LEMMA (exact).  On a split background (line bundles times trivial rank), every block of Lambda^2 or Sym^2 of a single stack is
    L_i^2 with index 4 a_i b_i (even), every adjoint block L_i L_j^-1 has net chirality 0 (parity theorem), and the only blocks that can
    carry an odd index are the fundamental/antifundamental blocks L_i^{+-1} (index +-a_i b_i) and the mixed blocks L_i L_j of Lambda^2/Sym^2
    (index (a_i + a_j)(b_i + b_j)).  Consequence for a U(N) parent with SM stacks c, L, 1: u^c and e^c cannot come from Lambda^2 of the
    colour or doublet stack (the SU(5)-like 10 = (3,2) + 3bar + 1 assignment fails on a split background), so colour antitriplets
    must come from antifundamental-type blocks, and u^c, d^c then need a further U(1) (or a non-split colour bundle) to differ in
    hypercharge (ledger row 3)."""
    rng = range(-box, box + 1)
    squares_even = all((2 * a * 2 * b) % 2 == 0 for a in rng for b in rng)
    return {"squares_even": squares_even, "adjoint_net_zero": True, "odd_candidates": ["L_i^{+-1} (fund/antifund)", "L_i L_j, i != j (Lambda^2/Sym^2 mixed blocks)"],
            "statement": "index(L_i^2 block) = 4 a_i b_i is even; index(L_i L_j) = (a_i + a_j)(b_i + b_j); index(L_i^{+-1}) = +-a_i b_i"}


def su6_three_stack_scan(index=3, degree_box=3):
    """Scan bidegrees (a_i, b_i) in a box for the three stacks c, L, 1 (U(6) = U(3) x U(2) x U(1)) and list the assignments in which the
    (3,2) block of Lambda^2 (Q from L_c L_L), the colour antitriplet block of an antifundamental (L_c^-1) and the doublet block of an
    antifundamental (L_L^-1) all have |index| = `index` with the sign pattern of one chirality (Q: +index, antitriplet: -index for the
    conjugate field, doublet: -index) — the minimal SU(5)-like B1 content Lambda^2(6) + 6bar without the Lambda^2 colour block.  EXACT scan."""
    rng = range(-degree_box, degree_box + 1); out = []
    for ac, bc, aL, bL in itertools.product(rng, repeat=4):
        iQ = (ac + aL) * (bc + bL); iu = -(ac * bc); iL = -(aL * bL)
        if iQ == index and iu == -index and iL == -index: out.append({"d_c": (ac, bc), "d_L": (aL, bL), "index_Q": iQ, "index_3bar": iu, "index_2bar": iL})
    return {"solutions": out, "count": len(out), "note": "e^c from Lambda^2 2 = L_L^2 has even index (lemma); it must come from a singlet-stack block"}


# ------------------------------------------------------------------ B1 finite search: Lambda^2(N) + (16 - N) fundamentals on split stacks
def b1_content_blocks(stacks, n_fund, n_anti=1, n_sym=0):
    """Net 4D left-handed content of the U(N) spectrum n_anti Lambda^2(N) (+) n_sym Sym^2(N) (+) n_fund N (net signed multiplicities; a
    conjugate representation of the opposite 8D chirality is the same net content) on a split background.  stacks: tuple of
    (name, rank, (a, b)); N = sum of ranks.  Each block: (colour rep, weak rep, charge vector over the stack U(1)s) -> net LH count
    = spin-twisted index times multiplicity; a negative index on B is recorded on the conjugate block B*.  EXACT (VI.1 index a b per
    block; Kuenneth)."""
    names = [n for n, _, _ in stacks]; idx = {n: i for i, n in enumerate(names)}; k = len(stacks)
    ranks = {n: r for n, r, _ in stacks}
    def colour(n): return "3" if ranks[n] == 3 else "1"
    def weak(n): return "2" if ranks[n] == 2 else "1"
    blocks = {}
    conj = {"3": "3bar", "3bar": "3", "1": "1", "6": "6bar", "6bar": "6", "2": "2", "3w": "3w"}
    def put(col, wk, charge, count):
        if count == 0: return
        if count < 0: col = conj[col]; charge = tuple(-c for c in charge); count = -count
        key = (col, wk, charge); blocks[key] = blocks.get(key, 0) + count
    for (ni, ri, di), (nj, rj, dj) in itertools.combinations(stacks, 2):          # mixed blocks L_i L_j of Lambda^2 and Sym^2 alike
        ch = [0] * k; ch[idx[ni]] += 1; ch[idx[nj]] += 1
        col = "3" if 3 in (ri, rj) else "1"; wk = "2" if 2 in (ri, rj) else "1"
        put(col, wk, tuple(ch), (n_anti + n_sym) * (di[0] + dj[0]) * (di[1] + dj[1]))
    for n, r, d in stacks:                                                       # L_i^2 blocks: Lambda^2 (3bar, singlet) and Sym^2 (6, weak triplet, singlet)
        ch = [0] * k; ch[idx[n]] = 2; ch = tuple(ch); i2 = 4 * d[0] * d[1]
        if r == 3: put("3bar", "1", ch, n_anti * i2); put("6", "1", ch, n_sym * i2)
        elif r == 2: put("1", "1", ch, n_anti * i2); put("1", "3w", ch, n_sym * i2)
        else: put("1", "1", ch, n_sym * i2)
    for n, r, d in stacks:                                                       # fundamentals L_i
        ch = [0] * k; ch[idx[n]] = 1; put(colour(n), weak(n), tuple(ch), n_fund * d[0] * d[1])
    return blocks

def _solve_rational(rows, rhs, n):
    """Exact Gaussian elimination over Q: rows (list of coefficient lists), rhs; returns (solution list with None for free unknowns,
    consistent flag).  Free unknowns are reported as None; the caller pins them."""
    M = [[Fraction(c) for c in r] + [Fraction(b)] for r, b in zip(rows, rhs)]; piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; M[r] = [v / M[r][c] for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0: M[i] = [a - M[i][c] * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    if any(all(v == 0 for v in row[:-1]) and row[-1] != 0 for row in M): return None, False
    sol = [None] * n
    for i, c in enumerate(piv):
        if all(M[i][j] == 0 for j in range(n) if j != c): sol[c] = M[i][-1]
    return sol, True

def hypercharge_assignment(blocks, targets):
    """Solve Y = sum_i y_i q_i exactly (rationals) for the stack hypercharges y from `targets` = {block key: hypercharge}.  Returns None if
    inconsistent; otherwise {"y": [...], "pinned": bool, "block_hypercharges": {...}} — with unpinned stacks reported as None and the
    dependent block hypercharges only when every stack is pinned (row 3 of the ledger)."""
    k = len(next(iter(blocks))[2])
    sol, ok = _solve_rational([list(key[2]) for key in targets], [Fraction(Y) for Y in targets.values()], k)
    if not ok: return None
    pinned = all(v is not None for v in sol)
    hy = {key: sum(c * yy for c, yy in zip(key[2], sol)) for key in blocks} if pinned else None
    return {"y": sol, "pinned": pinned, "block_hypercharges": hy}

def sm_net_content(blocks, hy):
    """Net left-handed counts by SM representation after identifying a block with its conjugate (2 ~ 2bar): returns the dictionary
    {(colour, weak, Y): net} with every entry nonzero, hypercharges as given in `hy`."""
    net = {}
    for key, v in blocks.items():
        col, wk, _ = key; Y = hy[key]
        if col == "3bar": key2 = ("3", wk, -Y); v2 = -v            # record antitriplets as negative triplets
        elif col == "1" and wk == "2": key2 = ("1", "2", Y if Y <= 0 else -Y); v2 = v if Y <= 0 else -v
        elif col == "1": key2 = ("1", "1", abs(Y)); v2 = v if Y >= 0 else -v
        else: key2 = (col, wk, Y); v2 = v
        net[key2] = net.get(key2, 0) + v2
    return {k: v for k, v in net.items() if v != 0}

SM_NET = lambda f: {("3", "2", Fraction(1, 6)): f, ("3", "1", Fraction(2, 3)): -f, ("3", "1", Fraction(-1, 3)): -f, ("1", "2", Fraction(-1, 2)): f, ("1", "1", Fraction(1)): f}

def b1_scan(k_singlets=2, box=3, families=3, n_anti=1, n_sym=0):
    """Exact finite search over split backgrounds for the quintic-free content Lambda^2(N) + (16 - N) N, N = 5 + k_singlets, with stacks
    c (U(3)), L (U(2)) and k singlet stacks, bidegrees in [-box, box]^2 (colour stack with a_c b_c = 0, forced by the Lambda^2 colour
    block 4 a_c b_c).  Accepts an assignment when, for some rational stack hypercharges y, the net left-handed content (conjugates
    identified, 2 ~ 2bar) is EXACTLY `families` Standard-Model families plus hypercharge-neutral singlets: no chiral exotic of any kind.
    The hypercharge y is pinned by Q, the two antitriplets, and the first doublet and singlet blocks (each tried at its SM values).
    Returns the accepted assignments and counters at each stage.  EXACT scan."""
    N = 5 + k_singlets; n_fund = -(N - 16) * n_anti - (N + 16) * n_sym            # quintic-free (sun_content_conditions)
    rng = range(-box, box + 1); pairs = [(a, b) for a in rng for b in rng]
    c_pairs = [(a, b) for a, b in pairs if a * b == 0]
    names = ["c", "L"] + [f"s{j}" for j in range(1, k_singlets + 1)]
    accepted = []; counters = {"Q_ok": 0, "two_antitriplets": 0, "coloured_content_is_SM": 0, "hypercharge_pinned": 0, "exact_SM": 0}
    seen = set()
    for dc in c_pairs:
        for dL in pairs:
            if (dc[0] + dL[0]) * (dc[1] + dL[1]) != families: continue
            counters["Q_ok"] += 1
            for ds in itertools.product(pairs, repeat=k_singlets):
                stacks = tuple(zip(names, (3, 2) + (1,) * k_singlets, (dc, dL) + ds))
                anti = [((dc[0] + d[0]) * (dc[1] + d[1])) for d in ds]
                if sum(1 for v in anti if v == -families) < 2: continue
                counters["two_antitriplets"] += 1
                blocks = b1_content_blocks(stacks, n_fund, n_anti, n_sym)
                col = {key: v for key, v in blocks.items() if key[0] != "1" or key[1] == "3w"}
                Qs = {key: v for key, v in col.items() if key[0] == "3" and key[1] == "2"}; anti3 = {key: v for key, v in col.items() if key[0] == "3bar" and key[1] == "1"}
                if sum(Qs.values()) != families or len(Qs) != 1 or len(col) != 3 or sorted(anti3.values()) != [families, families]: continue
                counters["coloured_content_is_SM"] += 1
                (q_key,) = Qs; a1, a2 = list(anti3)
                dbl = [key for key in blocks if key[0] == "1" and key[1] == "2"]; sng = [key for key in blocks if key[0] == "1" and key[1] == "1"]
                found = None
                pin_menu = [(d, (Fraction(-1, 2), Fraction(1, 2))) for d in dbl] + [(s_, (Fraction(1), Fraction(-1), Fraction(0))) for s_ in sng]
                def pinned_fits(targets, menu):
                    """Exhaustive: every acceptable assignment gives the next unpinned block one of its SM values."""
                    fit = hypercharge_assignment(blocks, targets)
                    if fit is None: return
                    if fit["pinned"]: yield fit; return
                    for i, (key, values) in enumerate(menu):
                        if key in targets: continue
                        for val in values: yield from pinned_fits({**targets, key: val}, menu[i + 1:])
                        return
                for Yu, Yd in ((Fraction(-2, 3), Fraction(1, 3)), (Fraction(1, 3), Fraction(-2, 3))):
                    for fit in pinned_fits({q_key: Fraction(1, 6), a1: Yu, a2: Yd}, pin_menu):
                        counters["hypercharge_pinned"] += 1
                        net = sm_net_content(blocks, fit["block_hypercharges"])
                        core = {kk: v for kk, v in net.items() if not (kk[0] == "1" and kk[1] == "1" and kk[2] == 0)}
                        if core == SM_NET(families):
                            found = {"stacks": stacks, "y": fit["y"], "blocks": blocks, "hypercharges": fit["block_hypercharges"], "net": net}; break
                    if found: break
                if found:
                    counters["exact_SM"] += 1
                    sig = (tuple(d for _, _, d in found["stacks"]), tuple(found["y"]))
                    if sig not in seen: seen.add(sig); accepted.append(found)
    return {"N": N, "n_fund": n_fund, "n_anti": n_anti, "n_sym": n_sym, "k_singlets": k_singlets, "box": box, "counters": counters, "accepted": accepted, "n_accepted": len(accepted)}


def b1_trilemma(scan):
    """The three-stack trilemma transported to the B1 class: for every accepted assignment of `b1_scan`, (Y) the bidegree Yukawa rule for
    (Q, u^c) and (Q, d^c) with the package convention (the block, not its conjugate), (H) the slope-free ratio of each Higgs block
    L_L L_s^-1, (C) the colour-slope ratio of the rank-5 rest (w = 2 d_L + sum d_s - 5 d_c), and the counts: Y+H (both Higgs loci exist and
    coincide), Y+C, and Y+H+C at ONE common r.  Box 3: 60 accepted, Y+H and Y+C disjoint (no assignment has both loci); box 4: 156 accepted,
    Y+H+C never at a common r — the loci that do coexist sit at r_H = 2 with r_c = 5/2 (and the mirror 1/2, 2/5).  EXACT (finite)."""
    from . import product_surface as PS
    rows = []
    for acc in scan["accepted"]:
        d = {n: dd for n, _, dd in acc["stacks"]}; hy = acc["hypercharges"]; singlets = [n for n, r, _ in acc["stacks"] if r == 1]
        Q = (d["c"][0] + d["L"][0], d["c"][1] + d["L"][1]); su = sd = None
        for key in acc["blocks"]:
            if key[0] == "3bar" and key[1] == "1":
                s_ = singlets[key[2][2:].index(-1)]
                if hy[key] == Fraction(-2, 3): su = s_
                else: sd = s_
        Y = all(PS.yukawa_bidegree_rule(Q, (d["c"][0] + d[s_][0], d["c"][1] + d[s_][1]))["allowed"] for s_ in (su, sd))
        H = {lab: (d["L"][0] - d[s_][0], d["L"][1] - d[s_][1]) for lab, s_ in (("u", su), ("d", sd))}
        rH = {lab: (Fraction(-h[0], h[1]) if h[1] and Fraction(-h[0], h[1]) > 0 else None) for lab, h in H.items()}
        w = (2 * d["L"][0] + sum(d[s_][0] for s_ in singlets) - 5 * d["c"][0], 2 * d["L"][1] + sum(d[s_][1] for s_ in singlets) - 5 * d["c"][1])
        rc = Fraction(-w[0], w[1]) if w[1] and Fraction(-w[0], w[1]) > 0 else None
        rows.append({"stacks": acc["stacks"], "Q": Q, "u_stack": su, "d_stack": sd, "Y": Y, "H_blocks": H, "r_H": rH, "w": w, "r_c": rc, "y": acc["y"]})
    YH = [r for r in rows if r["Y"] and r["r_H"]["u"] is not None and r["r_H"]["u"] == r["r_H"]["d"]]
    YC = [r for r in rows if r["Y"] and r["r_c"] is not None]
    YHC = [r for r in YH if r["r_c"] == r["r_H"]["u"]]
    both = [r for r in YH if r["r_c"] is not None]
    return {"rows": rows, "n": len(rows), "Y+H": len(YH), "Y+C": len(YC), "Y+H+C_common_r": len(YHC), "H_and_C_at_different_r": [(r["r_H"]["u"], r["r_c"]) for r in both],
            "conclusion": "no accepted assignment has a Yukawa, a slope-free Higgs and a colour-matched slope at one ratio" if not YHC else "Y+H+C found"}
