"""R2C-07 (2026-09-19): the product surface S = X0(143) x E (E = 143a1) — index and slope decoupled (exact arithmetic).

On the curve, chirality (index = deg) and the Higgs mass (slope = deg) are one number (R2C-06).  On a surface they are not: for
L = M (x) N (M on X, N on E, degrees dM, dN):  c1(L)^2 = 2 dM dN,  fermion index = dM dN (K_E trivial; S0 on X),  slope mu(L) = dM A_E + dN A_X.
Family bundle S0(sum P) (x) N with (dM, dN) = (3, -1): index -3 — three chiral families of one chirality (h^0(X, S0(sum P)) = 3 pure, h^1(E, N) = 1).
Gauge-vertex triangle (R2C-06 in two degrees): L_H = -(L_1 + L_2) = degrees (-6, 2), the same ray as the family bundle, so all slopes vanish
together on the locus A_X = 3 A_E — a polystable point: chiral families, massless Higgs, no tachyon.  Higgs zero modes there: the vector
modes H^1(S, L_H) = h^1(X, O(-2 sum P)) h^0(E, N^2) = 18 x 2 = 36.  Off the locus, m^2 ~ slope: the Higgs slope 2(A_X - 3 A_E) and the family
blocks' vector-mode slope -(A_X - 3 A_E) have opposite signs, so for A_X < 3 A_E the Higgs is tachyonic (electroweak breaking) while the
family-block scalars are massive: the electroweak scale is a Kaehler-modulus deviation, not a flux number.  CC-31 (R2C-08): the earlier claim of a factorised
Yukawa I x T for two (3,-1) families was WRONG — two families of the same Kuenneth parity have no gauge-vertex Yukawa (see below).
OPEN (new gates): the parent is an 8D gauge theory (four internal real dimensions): 8D chirality rules and anomalies replace R2C-01/03;
stabilisation of (A_X, A_E) near the locus; the other blocks' massless moduli at the locus; gravity with two Kaehler moduli."""
import itertools
import sympy as sp

A_X, A_E = sp.symbols("A_X A_E", positive=True)

def curve_cohomology(degree, spin_twisted):
    """(h^0, h^1) on X0(143): spin-twisted (S0 (x) O(D), index = degree, pure for our divisors) or untwisted line bundle of that degree
    (h^0 - h^1 = degree - 12; h^1 = 12 - degree + h^0(O(-D)) for negative degree D on the CM divisors: h^0 = 0)."""
    d = int(degree)
    if spin_twisted: return (d, 0) if d > 0 else ((0, -d) if d < 0 else (2, 2))
    if d < 0: return (0, 12 - d)
    if d == 0: return (1, 13)
    return (d - 12 + (1 if d in (3, 6) else 0), (1 if d in (3, 6) else 0))       # h^0(O(sum P)) = 1 (hecke.gonality_lower_bound) and h^0(O(2 sum P)) = 1 (canonical.gates.gate_petri_w13_quotient): both proved in v0.33.0

def torus_cohomology(degree):
    n = int(degree); return (n, 0) if n > 0 else ((0, -n) if n < 0 else (1, 1))

def product_block(dM, dN):
    hMf = curve_cohomology(dM, True); hMv = curve_cohomology(dM, False); hN = torus_cohomology(dN)
    ferm = (hMf[0] * hN[0] + hMf[1] * hN[1], hMf[0] * hN[1] + hMf[1] * hN[0]); vec_h1 = hMv[0] * hN[1] + hMv[1] * hN[0]
    return {"degrees": (dM, dN), "c1_squared": 2 * dM * dN, "fermion_index": dM * dN, "fermion_zero_modes": ferm, "vector_zero_modes_H1": vec_h1, "slope": sp.expand(dM * A_E + dN * A_X)}

def triangle(L1, L2):
    LH = (-(L1[0] + L2[0]), -(L1[1] + L2[1])); return {"L_H": LH, "block": product_block(*LH), "slope_zero_locus": sp.solve(sp.Eq(product_block(*LH)["slope"], 0), A_X)}

def m4_record():
    fam = product_block(3, -1); tri = triangle((3, -1), (3, -1)); H = tri["block"]
    dev = sp.Symbol("delta")                                                        # A_X = 3 A_E - delta
    return {"model_id": "M4", "surface": "X0(143) x 143a1", "family_block": fam, "higgs_block": H, "polystable_locus": tri["slope_zero_locus"],
            "higgs_zero_modes_at_locus": H["vector_zero_modes_H1"], "family_zero_modes": fam["fermion_zero_modes"],
            "off_locus": {"higgs_slope": sp.expand(H["slope"].subs(A_X, 3 * A_E - dev)), "family_vector_slope": sp.expand(fam["slope"].subs(A_X, 3 * A_E - dev)),
                          "electroweak_side": "delta > 0 (A_X < 3 A_E): Higgs tachyonic, family-block scalars massive"},
            "yukawa": "NONE for two (3,-1) families (same Kuenneth parity; CC-31).  Mixed-origin solutions (R2C-08): Q (3,+-1) from the curve, u^c (+-1,-3) from the torus, Higgs (-4,2) at A_X = 2 A_E or (-2,4) at A_X = A_E/2; M = A v B^T, rank <= rank v (top-only at one VEV)",
            "open_gates": ["8D chirality and anomaly rules", "Kaehler moduli stabilisation near A_X = 3 A_E", "masses of the other blocks' moduli", "gravity with two Kaehler moduli"], "status": "research record (R2C-07)"}


# ------------------------------------------------ R2C-08 (2026-09-19): 8D chirality selection rule on the product surface, exhaustive scan
def kunneth_parity(dM, dN):
    """Internal chirality of the family mode: H^0 (positive degree) counts 0, H^1 (negative degree) counts 1, on each factor; parity mod 2."""
    return ((dM < 0) + (dN < 0)) % 2

def yukawa_selection_rule(fam1, fam2):
    """The gauge-vertex Yukawa is the (0,2)-form int_S psi_1 ^ psi_2 ^ A_H valued in K_S (spinor bundles K^{1/2} K^{1/2} = K; the block
    bundles multiply to O).  It needs (q1, q2, q_H) with q1 + q2 + q_H = 2 and q_H = 1 for the (0,1) Higgs polarisation, i.e. OPPOSITE Kuenneth
    parities of the two families; with the (1,0) polarisation (a section of L_H (x) K) the form would be valued in K^2 and cannot be integrated.
    The 8D chiralities of the two blocks are then equal (the 4D-L/4D-R requirement), as in R2C-01."""
    p1, p2 = kunneth_parity(*fam1), kunneth_parity(*fam2)
    return {"parities": (p1, p2), "yukawa_allowed": p1 != p2, "reason": "opposite Kuenneth parity required" if p1 != p2 else "same parity: no (0,2)-form of the right degree (K vs K^2)"}

def family_types(index=3): return [(a, b) for a in (index, -index, 1, -1) for b in (index, -index, 1, -1) if abs(a * b) == index]

def scan_family_pairs(index=3):
    """All ordered pairs of family types with |index| = 3 and opposite parity; the Higgs block and its slope-free locus (if positive)."""
    out = []
    for f1 in family_types(index):
        for f2 in family_types(index):
            if not yukawa_selection_rule(f1, f2)["yukawa_allowed"]: continue
            H = (-(f1[0] + f2[0]), -(f1[1] + f2[1])); mu = H[0] * A_E + H[1] * A_X; loc = sp.solve(sp.Eq(mu, 0), A_X)
            ok = bool(loc) and loc[0] != 0 and (loc[0] / A_E).is_positive
            origins = tuple("curve" if abs(f[0]) == index else "torus" for f in (f1, f2))
            out.append({"Q": f1, "u": f2, "H": H, "origins": origins, "slope_free_locus": loc[0] if ok else None})
    return out

def curve_curve_no_go(index=3):
    rows = [r for r in scan_family_pairs(index) if r["origins"] == ("curve", "curve")]
    return {"pairs": len(rows), "slope_free": [r for r in rows if r["slope_free_locus"] is not None], "no_go": all(r["slope_free_locus"] is None for r in rows)}

def mixed_origin_solutions(index=3):
    return [r for r in scan_family_pairs(index) if r["slope_free_locus"] is not None]

def rank_structure():
    """With Q's family index on the curve and u^c's on the torus, M_ij(v) = sum_{m,alpha} v_{m alpha} a_i^{(m)} b_j^{(alpha)} = (A v B^T)_ij:
    rank M <= rank v.  One Higgs VEV gives rank one (only the top massive at leading order); the lighter masses are set by the subleading
    singular values of the VEV matrix — a hierarchy mechanism, not a prediction of its size.  The torus factors are triple products of
    theta functions on 143a1 (degrees 1, 3, 4 for the (3,-1)/(-1,-3)/(-2,4) solution)."""
    return {"mass_matrix": "M = A v B^T", "rank_bound": "rank(v)", "one_VEV": "rank 1: top only", "torus_factor": "theta triple products on 143a1",
            "v0330_note": "with the u^c curve factor 2-dimensional the general bound is rank M <= min(3, 2 * kunneth_rank(v)); see theta_torus.up_mass_rank_bound"}


# ------------------------------------------------ CC-32 / R2C-09: the FULL selection rule (bidegrees sum to (1,1)) and the two solutions
def bidegree(dM, dN): return (1 if dM < 0 else 0, 1 if dN < 0 else 0)

def yukawa_bidegree_rule(fam1, fam2):
    """Lambda^{0,2}(S) = Lambda^{0,1}(X) (x) Lambda^{0,1}(E): the bidegrees (q_X, q_E) of the two families and the (0,1) Higgs must sum to (1,1).
    Parity alone (R2C-08) was necessary, not sufficient — CC-32.  The Higgs bidegree is the complement; a (0,0) Higgs would be an elementary scalar."""
    s = (bidegree(*fam1)[0] + bidegree(*fam2)[0], bidegree(*fam1)[1] + bidegree(*fam2)[1])
    return {"family_bidegree_sum": s, "higgs_bidegree": (1 - s[0], 1 - s[1]) if s in ((1, 0), (0, 1)) else None, "allowed": s in ((1, 0), (0, 1))}

def _curve_untwisted(d):
    if d < 0: return (0, 12 - d)
    if d == 0: return (1, 13)
    return (1, 13 - d)                                                   # CM-point divisors of degree 1..12: one section, h^1 = 13 - d

def scan_family_pairs_full(index=3):
    out = []
    for f1 in family_types(index):
        for f2 in family_types(index):
            r = yukawa_bidegree_rule(f1, f2)
            if not r["allowed"]: continue
            H = (-(f1[0] + f2[0]), -(f1[1] + f2[1])); hb = r["higgs_bidegree"]; modes = _curve_untwisted(H[0])[hb[0]] * torus_cohomology(H[1])[hb[1]]
            mu = H[0] * A_E + H[1] * A_X; loc = sp.solve(sp.Eq(mu, 0), A_X); ok = bool(loc) and loc[0] != 0
            out.append({"Q": f1, "u": f2, "H": H, "higgs_bidegree": hb, "higgs_modes": modes, "slope_free_locus": loc[0] if ok else None})
    return out

def surface_solutions(index=3):
    """The two mixed-origin solutions (mirror pairs): S1 = ((3,1),(1,-3)) with Higgs (-4,2), 32 modes, A_X = 2 A_E; S2 = ((-3,1),(1,3)) with Higgs
    (2,-4), 4 modes, A_X = A_E/2.  Every curve x curve pair fails, and so does every torus x torus pair (a zero degree on one factor).
    Conventions made explicit in v0.33.0 (compendium VI.3): a family's bidegree is that of its net chiral zero modes (`bidegree`); the
    degree-1 curve factor is S0(P) with P one of P1, P2, P3 (for P = P4 the S2 Higgs factor O(sum P - P) has no sections); the S2 count
    h^0(O(P2 + P3)) = 1 uses that X0(143) is not hyperelliptic (hecke.gonality_lower_bound)."""
    return [r for r in scan_family_pairs_full(index) if r["slope_free_locus"] is not None and r["higgs_modes"] > 0]


# ------------------------------------------------ v0.33.0 (2026-09-23): the Higgs level on the surface from Chapter V (compendium VI.3, remark)
def higgs_slope_mass(H, bidegree_H):
    """Lowest level of the Higgs polarisation of a block of bidegree H = (H_X, H_E) whose (0,1)-leg sits on the factor with bidegree_H = (1,0)
    (leg on X) or (0,1) (leg on E), for the product connection and product metric: the leg factor contributes the vector level
    -2 pi |deg|/A (V.2, curvature cancels pointwise), the other factor its scalar Landau level +2 pi deg/A (V.1).  The result is
    2 pi mu(L_H)/(A_X A_E): massless exactly on the slope-free locus, tachyonic for mu < 0, massive for mu > 0.  S1: H = (-4, 2), leg on X;
    S2: H = (2, -4), leg on E."""
    HX, HE = H
    if bidegree_H == (1, 0): m2 = -2 * sp.pi * abs(HX) / A_X + 2 * sp.pi * HE / A_E
    elif bidegree_H == (0, 1): m2 = 2 * sp.pi * HX / A_X - 2 * sp.pi * abs(HE) / A_E
    else: raise ValueError("bidegree_H must be (1,0) or (0,1)")
    mu = HX * A_E + HE * A_X
    return {"m2_lowest": sp.simplify(m2), "slope": mu, "equals_2pi_slope_over_areas": sp.simplify(m2 - 2 * sp.pi * mu / (A_X * A_E)) == 0,
            "massless_locus": sp.solve(sp.Eq(mu, 0), A_X)}


# ================================================================== v0.33.1 (2026-09-27): exact block census of the split surface backgrounds
# For a split background V = (+)_i L_i^{(+) n_i} on S = X x E (product metric, areas A_X = r A_E; product HYM connection) every
# off-diagonal gauge block L_i (x) L_j^{-1} is a line bundle of bidegree (a, b).  Its internal (0,1)-forms have an X-leg and an E-leg,
# and the Jacobi operator of T29 separates on product modes: the leg carrying the form obeys V.2/V.3 on its factor (tachyonic
# polarisation 2 d*d - |B|, massive 2 dbar*dbar + |B|), the other factor contributes its scalar Bochner spectrum (V.1/V.4).  On the flat
# torus every level is a Landau level, (2n+1)|b| in units of 2 pi/A_E; on the curve only the ground level is exact (attained iff the
# relevant h^0 is positive, Lemma E.1 and RR supply the multiplicities).  So the census is exact for ground x ladder products and
# gives (i) the harmonic (Dolbeault-closed) negative directions, i.e. the Kuenneth H^1 classes with negative level, and (ii) a lower
# bound on the Morse index from all exactly known product modes.  Units: 2 pi/A_E; x = 1/r = A_E/A_X; 2 pi mu(L)/(A_X A_E) = a x + b.
S2_STACKS = (("c", 3, (0, 0)), ("L", 2, (3, -1)), ("1", 1, (1, 3)))      # Q = L_c L_L^-1 = (-3, 1), u^c = L_1 L_c^-1 = (1, 3), H = L_L L_1^-1 = (2, -4)
S1_STACKS = (("c", 3, (0, 0)), ("L", 2, (-3, -1)), ("1", 1, (1, -3)))    # Q = (3, 1), u^c = (1, -3), H = (-4, 2)
STACK_MODELS = {"S1": S1_STACKS, "S2": S2_STACKS}

def h0_cm_divisor(m, genus_Y=6):
    """h^0(X, O(D)) and h^0(K_X(-D)) for D = m1 P1 + m2 P2 + m3 P3 (+ m4 P4), the P_i the W13-fixed CM points.  EXACT where returned.
    W13 fixes each P_i, so H^0(K_X(-D)) splits into W13-eigenspaces.  Invariant differentials are odd in a local coordinate z with
    W13: z -> -z, descend to Y = X/W13 (genus 6) and satisfy ord_{P_i} >= m_i iff the descended form vanishes to order floor(m_i/2) at
    Q_i = pi(P_i); anti-invariant ones are even and satisfy ord >= m_i iff ord >= 2 ceil(m_i/2).  With gonality(Y) >= 4 (Petri gate) an
    effective divisor E on Y of degree <= 3 has h^0(K_Y(-E)) = 6 - deg E, and with gonality(X) >= 4 the value conditions at <= 3 of the
    P_i on the 7-dimensional anti-invariant space are independent.  So for ceil(m_i/2) <= 1 at at most three points and sum floor(m_i/2)
    <= 3:  h^0(K_X(-D)) = (6 - sum floor(m_i/2)) + (7 - #{i : m_i >= 1}),  and Riemann-Roch gives h^0(O(D)).  Examples: sum P -> 1,
    2 sum P -> 1 (Lemma E.2), sum P + P1 -> 1 (closes the S1 degree-4 caveat), 2 sum P + P1 -> 1."""
    m = tuple(m) + (0,) * (4 - len(m))
    pts = [i for i in range(4) if m[i] >= 1]
    if len(pts) > 3 or any(-(-mi // 2) > 1 for mi in m) or sum(mi // 2 for mi in m) > 3:
        return None
    h0K = (genus_Y - sum(mi // 2 for mi in m)) + (7 - len(pts))
    deg = sum(m); return {"h0_K_minus_D": h0K, "h0_O_D": h0K + deg + 1 - 13, "degree": deg}

def curve_h0_flux(a):
    """h^0(X, O(D)) for the flux line bundle of degree a supported on the marked CM points, D = sum P (a = 3), P1 (1), P2 + P3 (2),
    sum P + P1 (4), 2 sum P (6), and their negatives.  All exact: gonality >= 4 (hecke.gonality_lower_bound) for |a| <= 3, and
    `h0_cm_divisor` (the W13 eigenspace decomposition with the Petri gate on the quotient) for a = 4 and a = 6."""
    if a < 0: return 0, None
    if a == 0: return 1, None
    table = {1: (1,), 2: (0, 1, 1), 3: (1, 1, 1), 4: (2, 1, 1), 5: (2, 2, 1), 6: (2, 2, 2)}
    if a in table: return h0_cm_divisor(table[a])["h0_O_D"], None
    return 1, f"degree {a}: not covered by h0_cm_divisor"

def curve_h1_flux(a):
    """h^1 = h^0 - a + 12 (Riemann-Roch, genus 13), with the same caveat as curve_h0_flux."""
    h0, note = curve_h0_flux(a); return h0 - a + 12, note

def block_census(r, model="S2", stacks=None, torus_levels=6):
    """Exact census of every gauge block of a split background at A_X/A_E = r.

    Returns {"r", "x", "blocks": [...], "summary": {...}}.  Each block lists, for its X-leg and E-leg, the ground level (units
    2 pi/A_E), its multiplicity per component, whether the ground mode is Dolbeault-closed (a Kuenneth H^1 class), the exactly
    known torus ladder above it, and the number of components n_comp = n_i n_j.  The summary counts negative exact-product modes
    (a lower bound on the Morse index of the split background), the harmonic negative directions among them, the harmonic
    massless directions, and the critical ratios where a ground level crosses zero.  Compendium V.1–V.3 and VI.3 (remark).
    At r = 1/2 for S2 this reproduces the coupled-relaxation study's 279 negative directions with 90 Dolbeault-closed, and the
    4 + 56 massless Higgs coefficients per component."""
    stacks = stacks if stacks is not None else STACK_MODELS[model]
    x = sp.Rational(1) / sp.nsimplify(r)
    blocks, notes = [], []
    neg = neg_closed = zero_closed = zero_total = 0; crit = set()
    for (ni, ranki, di), (nj, rankj, dj) in itertools.permutations(stacks, 2):
        a, b = di[0] - dj[0], di[1] - dj[1]; ncomp = ranki * rankj
        h1a, note1 = curve_h1_flux(a); h0abs, note0 = curve_h0_flux(abs(a))
        for note in (note1, note0):
            if note: notes.append(f"block {ni}{nj} ({a},{b}): {note}")
        # X-leg: (0,1)-form on X in L_a  (x)  scalar on E in N_b
        gX = a * x + abs(b); mX = h1a * (abs(b) if b else 1); closedX = b >= 0
        ladderX = [(a * x + abs(b) * (2 * n + 1), h1a * abs(b)) for n in range(1, torus_levels)] if b else []
        # E-leg: scalar on X in L_a  (x)  (0,1)-form on E in N_b
        gE = abs(a) * x + (-abs(b) if b < 0 else (3 * b if b > 0 else 0)); mE = h0abs * (abs(b) if b else 1); closedE = (a >= 0) and (b <= 0)
        ladderE = [(abs(a) * x + (abs(b) * (2 * n - 1) if b < 0 else b * (2 * n + 3)), h0abs * abs(b)) for n in range(1, torus_levels)] if b else []
        legs = {"X_leg": {"ground": gX, "multiplicity": mX, "closed": closedX, "ladder": ladderX},
                "E_leg": {"ground": gE, "multiplicity": mE, "closed": closedE, "ladder": ladderE}}
        for leg in legs.values():
            g, m = leg["ground"], leg["multiplicity"]
            if g < 0:
                neg += m * ncomp
                if leg["closed"]: neg_closed += m * ncomp
            elif g == 0:
                zero_total += m * ncomp
                if leg["closed"]: zero_closed += m * ncomp
            for lv, lm in leg["ladder"]:
                if lv < 0: neg += lm * ncomp
        if a < 0 and b: crit.add(sp.Rational(-a, abs(b)))               # X-leg ground a x + |b| = 0  ->  r = |a|/|b|
        if b < 0 and a: crit.add(sp.Rational(abs(a), abs(b)))          # E-leg ground |a| x - |b| = 0  ->  r = |a|/|b|
        blocks.append({"block": f"L_{ni} (x) L_{nj}^-1", "bidegree": (a, b), "n_comp": ncomp, "slope_level": a * x + b, **legs})
    # diagonal blocks: Wilson-line moduli H^1(S, O) = 13 (X-leg) + 1 (E-leg) per component
    ndiag = sum(rk * rk for _, rk, _ in stacks); zero_closed += 14 * ndiag; zero_total += 14 * ndiag
    return {"model": model if stacks is STACK_MODELS.get(model) else "custom", "r": sp.nsimplify(r), "x": x, "units": "2 pi / A_E",
            "blocks": blocks, "notes": sorted(set(notes)),
            "summary": {"negative_exact_product_modes": neg, "harmonic_negative": neg_closed, "harmonic_massless": zero_closed,
                        "massless_total_exact": zero_total, "diagonal_components": ndiag, "critical_r": sorted(crit),
                        "polystable": all(b["slope_level"] == 0 for b in blocks)}}

def census_scan(model="S2", ratios=(sp.Rational(1, 4), sp.Rational(1, 3), sp.Rational(1, 2), 1, 2, 3, 4)):
    """Summary rows of `block_census` across ratios: (r, negative, harmonic negative, harmonic massless)."""
    return [(sp.nsimplify(r), *(lambda s: (s["negative_exact_product_modes"], s["harmonic_negative"], s["harmonic_massless"]))(block_census(r, model)["summary"])) for r in ratios]


# ================================================================== v0.33.1 (2026-09-27): the S1/S2 extension graphs (holomorphic recombination data)
def extension_graph(model="S2", r=None):
    """For every ordered pair of stacks (i, j): dim Ext^1(L_j, L_i) = h^1(S, L_i L_j^-1) and dim Ext^2 = h^2, from Kuenneth, with the
    level of the corresponding harmonic classes at r (units 2 pi/A_E; negative = tachyonic, zero = massless, positive = massive) —
    the exact data every holomorphic recombination of the split background starts from.  Iterated extensions are holomorphic bundles by
    construction (no obstruction theory), so the remaining question for each candidate is mu-stability at the chosen r."""
    stacks = STACK_MODELS[model]; x = None if r is None else sp.Rational(1) / sp.nsimplify(r)
    def curve(a):   # (h0, h1) of the untwisted flux bundle of degree a on X
        h0 = curve_h0_flux(a)[0]; return h0, h0 - a + 12
    out = []
    for (ni, ranki, di), (nj, rankj, dj) in itertools.permutations(stacks, 2):
        a, b = di[0] - dj[0], di[1] - dj[1]; (h0X, h1X), (h0E, h1E) = curve(a), torus_cohomology(b)
        ext1_X, ext1_E = h1X * h0E, h0X * h1E                          # X-leg classes H^1(X) (x) H^0(E) and E-leg classes H^0(X) (x) H^1(E)
        ext2 = h1X * h1E                                                # H^2(S, L) = H^1(X) (x) H^1(E)
        lev_X = None if x is None else a * x + abs(b); lev_E = None if x is None else abs(a) * x + (-abs(b) if b < 0 else (3 * b if b > 0 else 0))
        out.append({"extension": f"0 -> L_{ni} -> ? -> L_{nj} -> 0", "class_bundle": (a, b), "n_comp": ranki * rankj,
                    "ext1": ext1_X + ext1_E, "ext1_X_leg": ext1_X, "ext1_E_leg": ext1_E, "ext2": ext2,
                    "level_X_leg": lev_X if ext1_X else None, "level_E_leg": lev_E if ext1_E else None})
    return {"model": model, "r": None if r is None else sp.nsimplify(r), "edges": out}

def s1_colour_recombination_dims(k_colours=3):
    """S1 at its wall r = 2: the only tachyonic harmonic classes are the 9 of Ext^1(L_c, L_1) = C^3 (x) H^1(E, N_-3).  Recombining
    k colour components with L_1 by classes xi_1..xi_k in H^1(E, N_-3) that are linearly independent gives a rank-(k+1) extension E_k;
    gluing the doublet stack L_L on top uses Ext^1(E_k, L_L) = ker( Ext^1(L_1, L_L) -> Ext^2(L_c, L_L)^k ), the cup product with the
    xi's.  Ext^1(L_1, L_L) = H^1(X, O(-sum P - P)) (x) H^0(E, N_2) (16 x 2 = 32, the Higgs classes) and the cup product factorises:
    on X, multiplication by the constant section of O(P), H^1(O(-sum P - P)) -> H^1(O(-sum P)), surjective with 1-dimensional kernel;
    on E, the pairing H^1(N_-3) (x) H^0(N_2) -> H^1(N_-1), which is the S1 torus factor B (3 x 2, rank 2 — VI.5).  Hence the classes
    eta_{r alpha} must satisfy sum_alpha eta_{r alpha} B(xi_i)_alpha = 0 for the 15 non-kernel r and each i: with rank(B restricted to
    the xi's) = min(k, 2) the surviving Higgs-type gluing space has dimension 15 * (2 - min(k, 2)) + 2.  EXACT (Kuenneth + VI.5 rank)."""
    rank = min(k_colours, 2); dim = 15 * (2 - rank) + 2
    return {"k_colours": k_colours, "tachyonic_classes_used": 3 * k_colours, "rank_E_k": k_colours + 1,
            "slope_E_k_at_r2_units_A_E": sp.Rational(-5, k_colours + 1), "dim_Ext1_E_k_L_L": dim,
            "note": "Ext^1(L_c, L_L) = 0 in S1: the doublet stack can only be reached through L_1 (Higgs classes)"}

def cm_point_cp_structure():
    """CP on X0(143) (complex conjugation tau -> -conj(tau), the real structure of the curve over Q) pairs the W13-fixed points as
    {P1, P2}, {P3, P4} (they share j = -82306.31.. and j = 6896962306.31.., the two roots of the class polynomial of discriminant -52),
    while W11 pairs them as {P1, P3}, {P2, P4} (u -> -u; Part A, Lemma II.3.B).  The family divisor sum P = P1 + P2 + P3 ("omit P4")
    is mapped by CP to "omit P3", by W11 to "omit P2" and by CP W11 to "omit P1": the four choices form one orbit of V4 = <CP, W11>
    and none of them is CP-invariant.  Since Aut(X0(143)) is the Atkin–Lehner group, "omit P4" ~ "omit P2" and "omit P3" ~ "omit P1"
    are the two holomorphically inequivalent three-family models, and they are CP conjugates of each other: the flux choice breaks CP
    explicitly and geometrically.  The torus factor does not: 143a1 has rational j, its tau lies on the CP-symmetric line Re tau = 1/2,
    and every rephasing-invariant quartet phase of the closed-form B (S1) and B' (S2) is 0 or pi (`torus_cp_test`)."""
    return {"CP_pairs": (("P1", "P2"), ("P3", "P4")), "W11_pairs": (("P1", "P3"), ("P2", "P4")),
            "u_values": {"P1": "-i/sqrt(13)", "P2": "+i/sqrt(13)", "P3": "+i/sqrt(13)", "P4": "-i/sqrt(13)"},
            "orbit_of_family_divisor": {"identity": "omit P4", "CP": "omit P3", "W11": "omit P2", "CP W11": "omit P1"},
            "conclusion": "the three-family flux breaks CP explicitly; the torus factor is CP-conserving"}

def torus_cp_test(tol=1e-12):
    """Rephasing-invariant quartet phases arg(B_jb B_j'b' conj(B_jb') conj(B_j'b)) of the closed-form torus factors: all 0 or pi means
    B is equivalent to a real matrix, i.e. the torus contributes no CP-violating phase."""
    from . import theta_torus as TT
    import mpmath as mp
    out = {}
    for name, k2 in (("S1", 2), ("S2", 3)):
        B = TT.torus_factor_closed_form(1, k2)["B"]; phases = []
        for (j, jj) in itertools.combinations(range(B.rows), 2):
            for (b, bb) in itertools.combinations(range(B.cols), 2):
                q = B[j, b] * B[jj, bb] * mp.conj(B[j, bb]) * mp.conj(B[jj, b])
                if abs(q) > 1e-20: phases.append(float(mp.arg(q)))
        out[name] = {"quartet_phases": phases, "real_up_to_rephasing": all(min(abs(p), abs(abs(p) - float(mp.pi))) < tol for p in phases)}
    return out


# ================================================================== v0.33.1 (2026-09-27): the S1 iterated extensions are mu-unstable for r < 5
def s1_iterated_extension_instability(r=2):
    """THEOREM (exact).  Let V be any rank-6 holomorphic bundle on S = X0(143) x 143a1 admitting a filtration whose graded pieces are
    the six S1 line bundles L_c^3 = O^3, L_L^2 = (-3,-1)^2, L_1 = (1,-3) — i.e. any iterated extension of the S1 split background.
    Then V is mu-unstable for the product polarisation with r = A_X/A_E < 5; in particular at the S1 wall r = 2.

    Proof.  mu(V) = -(5 + 5r)/6 A_E.  (i) If some colour line L_c is a subsheaf of V it destabilises (slope 0).  Since
    Ext^1(L_c, L_L) = H^1((-3,-1)) = 0 and Ext^1(L_c, F) injects into Ext^1(L_c, L_1) = H^0(X, O(P)) (x) H^1(E, N_-3) = C^3 for every
    F built from the L_L's and L_1, the colour stack can only sit on top of L_1, attached by classes s (x) zeta(v), zeta: C^3 -> H^1(E, N_-3);
    if zeta has a kernel, that colour line splits off — so zeta must be an isomorphism.  (ii) For a point e of E the class delta_e in
    H^1(E, N_-3) (the coboundary of the point, dual to evaluation at e on H^0(N_3)) spans the kernel of H^1(E, N_-3) -> H^1(E, N_-3(e)),
    so the twisted colour line L_c(-X x {e}) v with zeta(v) = delta_e lifts to the colour–L_1 extension; the further lift to V is
    obstructed in H^1(L_L(X x {e}))^2 = (H^1(X, O(-sum P)) (x) H^0(E, N_-1(e)))^2, which vanishes for every e except the single point
    e_0 with N_1 = O(e_0).  Since zeta is onto, every delta_e is hit; pick e != e_0.  (iii) mu(L_c(-X x {e})) = -A_X = -r A_E exceeds
    mu(V) exactly when r < 5.  QED.  So no mu-stable bundle at the S1 wall is an iterated extension of the split stacks, for any choice
    of the 9 tachyonic classes and of the Higgs-type gluing (32 + 2 classes) — the S1 analogue of the S2 relaxation studies, obtained
    without obstruction theory.  For r >= 5 the argument gives nothing; stable bundles of this topology at large r are expected from
    fiberwise stability (E-degree -5 coprime to the rank), and are not iterated extensions of these six line bundles."""
    r = sp.nsimplify(r)
    mu_V = -(5 + 5 * r) / 6; mu_twisted = -r; mu_split = sp.Integer(0)
    return {"r": r, "mu_V_units_A_E": mu_V, "destabilisers": {"split colour line L_c": mu_split, "twisted colour line L_c(-X x {e})": mu_twisted},
            "unstable": bool(mu_twisted > mu_V), "threshold_r": 5,
            "ingredients": {"Ext1(L_c, L_L)": 0, "Ext1(L_c, L_1) per colour": 3, "kernel of H1(N_-3) -> H1(N_-3(e))": "C delta_e",
                            "H1(L_L(X x e)) for e != e_0": 0, "exception": "e = e_0 only (N_1 = O(e_0)); irrelevant since zeta is onto"}}


# ================================================================== v0.33.1 (2026-09-27): the colour slope condition and the three-stack trilemma
def colour_slope_condition(q, u):
    """An SU(3)-preserving HYM vacuum needs the colour summand L_c (x) C^3 to be a direct summand of a polystable bundle, hence
    mu(L_c) = mu(W) for W the (possibly recombined) rank-3 rest, whose c_1 relative to colour is w = 2 (L_L - L_c) + (L_1 - L_c) = -2q + u
    in terms of the family types q = Q, u = u^c (VI.3 conventions).  mu(W (x) L_c^-1) = (w_X A_E + w_E A_X)/3 vanishes for some positive
    r = A_X/A_E iff w_X w_E < 0.  Returns (w, matchable, r_colour)."""
    w = (-2 * q[0] + u[0], -2 * q[1] + u[1]); ok = w[0] * w[1] < 0
    return {"w": w, "colour_matchable": ok, "r_colour": sp.Rational(-w[0], w[1]) if ok else None}

def three_stack_trilemma():
    """THEOREM (exact, finite check over the 64 ordered family-type pairs).  In the three-stack framework (colour^3, doublet^2, singlet)
    on X0(143) x 143a1 with family types of index +-3, no pair (Q, u^c) satisfies all three of
      (Y) a gauge-vertex Yukawa (VI.2 admissibility),
      (H) a Higgs massless at leading order (slope-free Higgs block, VI.3),
      (C) an SU(3)-preserving polystable vacuum at some r (`colour_slope_condition`).
    (Y)+(H) are exactly the four VI.3 solutions, all with w of equal signs ((-5,-5), (7,1), (-5,-5), (1,7)): no r matches colour, whatever
    W recombines into — which is why both split walls close.  (H)+(C) forces q = u: same Kuenneth type, no Yukawa (CC-31).  (Y)+(C) is
    possible (eight pairs) but their Higgs block has same-sign entries, i.e. the Higgs is tachyonic at every r (electroweak breaking at the
    compactification scale).  Any two of the three can be had; not all three."""
    types = [(a, b) for a in (3, -3, 1, -1) for b in (3, -3, 1, -1) if abs(a * b) == 3]
    out = {"Y+H": [], "H+C": [], "Y+C": [], "Y+H+C": []}
    for q, u in itertools.product(types, types):
        flags = (int(q[0] < 0) + int(u[0] < 0), int(q[1] < 0) + int(u[1] < 0)); Y = flags in ((1, 0), (0, 1))
        H = (-(q[0] + u[0]), -(q[1] + u[1])); Hf = H[0] * H[1] < 0
        C = colour_slope_condition(q, u)["colour_matchable"]
        rec = (q, u, H, colour_slope_condition(q, u)["w"])
        if Y and Hf: out["Y+H"].append(rec)
        if Hf and C: out["H+C"].append(rec)
        if Y and C: out["Y+C"].append(rec)
        if Y and Hf and C: out["Y+H+C"].append(rec)
    out["conclusion"] = "no family-type pair is simultaneously Yukawa-admissible, Higgs-slope-free and colour-slope-matchable"
    return out


# ================================================================== v0.33.1 (2026-09-27): chirality parity on the spin surface, and the rank-2 net-index lemma
SPIN_TWIST = (12, 0)                     # K_S^{1/2} = S0 (x) O_E has bidegree (12, 0); K_S = (24, 0)

def adjoint_net_chirality(block, twist=SPIN_TWIST):
    """Net number of left-handed fermions in the bifundamental (N_i, N_j-bar) of an 8D gauge theory whose fermions transform in the
    ADJOINT, reduced on X0(143) x 143a1 with spinors twisted by a line bundle R (default the spin structure K^{1/2}).

    In 8D no reality condition relates the (i,j) and (j,i) blocks at fixed 8D chirality (charge conjugation flips 8D chirality), so both
    blocks are independent and net = index(R (x) L) - index(R (x) L^-1) = chi(R (x) L) - chi(R (x) L^-1) = c_1(L).(2R - K_S)  (Riemann–Roch).
    With chi(X, O(a)) = a - 12 and chi(E, N_b) = b:  net = 2 [a rho_E + b (rho_X - 12)]  for L = (a, b), R = (rho_X, rho_E).  EXACT."""
    a, b = block; rx, re = twist
    chi = lambda d: (d[0] - 12) * d[1]
    net = chi((rx + a, re + b)) - chi((rx - a, re - b))
    assert net == 2 * (a * re + b * (rx - 12))
    return {"block": (a, b), "twist": (rx, re), "index_block": chi((rx + a, re + b)), "index_conjugate": chi((rx - a, re - b)),
            "net_chirality": net, "vector_like": net == 0}

def chirality_parity_theorem(box=6):
    """THEOREM (exact).  On a spin Kaehler surface S, for fermions in the adjoint of a gauge group broken by a background bundle, twisted by any
    spin^c structure with determinant class c (c characteristic: c.D = D.D mod 2), the net chirality of every bifundamental block F is
    c_1(F).c = c_1(F)^2 (mod 2).  NS(X0(143) x 143a1) is an even lattice (K_S = 2 K^{1/2}; Wu), so the net chirality is EVEN: three
    generations cannot come from adjoint matter on this surface, for any flux, any recombination and any twist.  With the untwisted spin
    structure it is ZERO: the S1/S2 family blocks (index +-3) are cancelled by their conjugate blocks (same index, Serre duality preserving
    degree parity in complex dimension 2) — they are vector-like pairs.  Returns the exact check over a box of blocks and twists and the
    S1/S2 bookkeeping."""
    rng = range(-box, box + 1)
    odd = [(a, b, rx, re) for a in rng for b in rng for rx in range(0, 25) for re in rng
           if adjoint_net_chirality((a, b), (rx, re))["net_chirality"] % 2]
    untwisted = {adjoint_net_chirality((a, b))["net_chirality"] for a in rng for b in rng}
    fam = {n: adjoint_net_chirality(blk) for n, blk in (("S1_Q", (3, 1)), ("S1_u", (1, -3)), ("S2_Q", (-3, 1)), ("S2_u", (1, 3)))}
    return {"odd_cases": odd, "untwisted_nets": untwisted, "families": fam,
            "statement": "net chirality = c1(F).c = c1(F)^2 mod 2 = even on X0(143) x 143a1; zero for the spin twist"}

def rank2_net_index_lemma(summand1, summand2, r):
    """LEMMA (exact; parallelogram law).  For the index form Q(a, b) = a b (the spin-twisted index of a block, VI.1): let W1, W2 be equal-slope
    summands, each a line bundle m or a rank-2 extension with constituents m +- e whose extension block 2e satisfies Q(2e) <= 0 (necessary for
    Ext^1 != 0 in the destabilising-free direction: an extension block is never in the (-,-) quadrant, and a (+,+) block has positive slope).
    Then the net index of the block Hom(W2, W1) is  n1 n2 Q(m1 - m2) + n2 * sum Q(e1) + n1 * sum Q(e2) <= 0,  because m1 - m2 lies on the
    slope-zero ray where Q = -r b^2 <= 0.  Consequence (under a chiral parent, where one block per pair carries the families and VI.3's sign
    rule holds): Q and u^c can never have the opposite net indices -3, +3 in an SU(3) x SU(2) x U(1)_Y-preserving polystable vacuum whose
    stable summands have rank <= 2 — the recombination search over options (u^c partner, colour, doublet) finds no candidate, and this is why.
    Returns the net index and the decomposition."""
    import fractions
    Qf = lambda v: v[0] * v[1]
    r = fractions.Fraction(r)
    mu = lambda W: sum(fractions.Fraction(c[0]) + c[1] * r for c in W) / len(W)
    assert mu(summand1) == mu(summand2), "summands must have equal slope"
    net = sum(Qf((c1[0] - c2[0], c1[1] - c2[1])) for c1 in summand1 for c2 in summand2)
    return {"net_index": net, "bound": "<= 0", "holds": net <= 0}
