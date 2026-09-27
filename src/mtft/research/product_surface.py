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

def curve_h0_flux(a):
    """h^0(X, O(D)) for the flux line bundle of degree a supported on the marked CM points (P in {P1, P2, P3}).
    Exact for |a| <= 3 (gonality >= 4: hecke.gonality_lower_bound) and for a = 6 = 2 sum P (canonical.gates.gate_petri_w13_quotient);
    degree 4 (sum P + P) is 1 unless X carries a g^1_4 through those points, which is not excluded — returned as (1, 'open: 1 or 2')."""
    if a < 0: return 0, None
    if a == 0: return 1, None
    if a <= 3 or a == 6: return 1, None
    return 1, "open: 1 unless a g^1_4 passes through the marked points"

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
