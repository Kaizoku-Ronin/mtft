"""PV-01 (2026-09-29): the flux (Yang–Mills) energy of a split background on S = X0(143) x 143a1 as an exact function of the shape modulus
r = A_X/A_E, its topological floor, the shape minimiser, and the Einstein-frame scaling that decides what the flux energy stabilises —
the `cosmology.vacuum_energy(r)` tool of the parent-action ledger (docs/SM/PARENT_ACTION_REQUIREMENTS.md, C).

Exact facts (product metric, product HYM connections on the line bundles).  For L of bidegree (a, b):  F = f_X omega_X + f_E omega_E with
f_X = 2 pi a / A_X, f_E = 2 pi b / A_E, |omega_X|^2 = |omega_E|^2 = 1 pointwise, orthogonal; integrating over A_X A_E = r A_E^2:
    E(L; r) := int_S |F|^2 = 4 pi^2 ( a^2 / r + b^2 r )      — scale-invariant (Yang–Mills is conformal in four real dimensions)
             = 4 pi^2 [ (a + b r)^2 / r  -  2 a b ]           — the |Lambda F|^2 term (slope squared) minus c_1(L)^2 = 2 a b
    E >= 8 pi^2 |a b| = 4 pi^2 |c_1(L)^2|  with equality iff r = |a/b|  (AM >= GM on the two legs: (anti-)self-duality of F).
For a split background V = (+) L_i (x) C^{n_i}:  E(r) = 4 pi^2 (A/r + B r),  A = sum n_i a_i^2,  B = sum n_i b_i^2,
    minimiser r* = sqrt(A/B),  E(r*) = 8 pi^2 sqrt(A B) >= 8 pi^2 sum n_i |a_i b_i|  (Cauchy–Schwarz; equality iff every |a_i/b_i| is the same).
The split connection is a Yang–Mills critical point at every r; it is the energy minimum of its topological class only where it is
polystable (all slopes equal), which for S1/S2 never happens (block census).  The floor of the class is the Hermitian–Yang–Mills energy
    E_HYM(r) = 4 pi^2 [ Delta + A_1^2 / r + B_1^2 r ] / N,   c_1(V) = (A_1, B_1),  Delta = 2 N c_2 - (N - 1) c_1^2  (Bogomolov discriminant),
attained by a stable bundle of that topology if one exists (Delta >= 0 necessary); E_split(r) = E_HYM(r) exactly at a polystable r.
Einstein frame (8D gauge theory + gravity reduced on S; M_P^2 = M_8^6 Vol, Vol = A_X A_E): the flux term V = E(r) / (2 g_8^2) is
volume-independent in the 8D frame and scales as Vol^-2 in the Einstein frame; the curvature term of the hyperbolic curve,
-(M_8^6 / 2) int_S R = + 48 pi M_8^6 A_E = 48 pi M_8^6 (Vol / r)^{1/2}, scales as r^{-1/2} Vol^{-3/2}.  Both are positive and decreasing
in Vol: no stationary point in the volume (decompactification runaway, the surface analogue of compactification.einstein_frame_potential),
while at fixed volume the shape r is stabilised (E is convex, unbounded at 0 and infinity when A, B > 0).  A cosmological constant is not
selected by these two terms; they do not contain the papers' delta^-6 e^{-2/alpha} (CC-34).

Status: EXACT for the energies and their minimisers; the Einstein-frame statement is EXACT within the stated reduction (no fermion
Casimir energy, no moduli potential from the parent's extra terms, no warping)."""
import itertools
import sympy as sp

R = sp.Symbol("r", positive=True)
STACK_MODELS = {"S1": (("c", 3, (0, 0)), ("L", 2, (-3, -1)), ("1", 1, (1, -3))),
                "S2": (("c", 3, (0, 0)), ("L", 2, (3, -1)), ("1", 1, (1, 3)))}

def intersection(d, e):
    """Intersection form of NS(X0(143) x 143a1) on product classes: (a, b).(a', b') = a b' + a' b."""
    return d[0] * e[1] + e[0] * d[1]

def line_bundle_energy(a, b, r=R):
    """int_S |F|^2 / (4 pi^2) for the HYM product connection on the line bundle of bidegree (a, b) at A_X/A_E = r, with its decomposition."""
    r = sp.nsimplify(r); E = sp.Integer(a) ** 2 / r + sp.Integer(b) ** 2 * r
    return {"E_over_4pi2": E, "slope_term": (a + b * r) ** 2 / r, "c1_squared": 2 * a * b, "floor_over_4pi2": 2 * abs(a * b),
            "self_dual_ratio": (sp.Rational(abs(a), abs(b)) if b else None), "identity": sp.simplify(E - ((a + b * r) ** 2 / r - 2 * a * b)) == 0}

def split_topology(stacks):
    """Rank, c_1 = (A_1, B_1), c_2, c_1^2 and the Bogomolov discriminant Delta = 2 N c_2 - (N - 1) c_1^2 of V = (+) L_i (x) C^{n_i}."""
    N = sum(n for _, n, _ in stacks)
    A1 = sum(n * d[0] for _, n, d in stacks); B1 = sum(n * d[1] for _, n, d in stacks)
    c2 = sum(ni * nj * intersection(di, dj) for (_, ni, di), (_, nj, dj) in itertools.combinations(stacks, 2)) + sum(n * (n - 1) // 2 * intersection(d, d) for _, n, d in stacks)
    c1sq = 2 * A1 * B1; Delta = 2 * N * c2 - (N - 1) * c1sq
    return {"rank": N, "c1": (A1, B1), "c2": c2, "c1_squared": c1sq, "Delta": Delta, "bogomolov_ok": Delta >= 0}

def flux_energy(r=R, model="S2", stacks=None):
    """E(r)/(4 pi^2) of the split background, its exact minimiser and value, the topological floor and the Cauchy–Schwarz gap."""
    stacks = stacks if stacks is not None else STACK_MODELS[model]
    A = sum(n * d[0] ** 2 for _, n, d in stacks); B = sum(n * d[1] ** 2 for _, n, d in stacks)
    E = A / R + B * R
    r_star = sp.sqrt(sp.Rational(A, B)) if B else None; E_star = 2 * sp.sqrt(A * B) if B else None
    floor = 2 * sum(n * abs(d[0] * d[1]) for _, n, d in stacks)
    slopes = [d[0] + d[1] * R for _, _, d in stacks]
    return {"model": model if stacks is STACK_MODELS.get(model) else "custom", "A": A, "B": B, "E_over_4pi2": E, "E_at_r": sp.nsimplify(E.subs(R, sp.nsimplify(r))),
            "r_star": r_star, "E_star_over_4pi2": E_star, "floor_over_4pi2": floor, "cauchy_schwarz_gap": sp.simplify(E_star - floor) if E_star is not None else None,
            "slopes_units_A_E": slopes, "polystable_r": _common_root(slopes), "convex": bool(A > 0 and B > 0)}

def _common_root(slopes):
    """Positive r at which every slope equals the first (polystability of the split sum); "all" when there is a single slope."""
    if len(slopes) < 2: return "all"
    roots = None
    for s in slopes[1:]:
        diff = sp.expand(s - slopes[0])
        sol = "all" if diff == 0 else set(x for x in sp.solve(diff, R) if x.is_positive)
        if sol == "all": continue
        roots = sol if roots is None else roots & sol
    return "all" if roots is None else sorted(roots)

def hym_floor(r=R, model="S2", stacks=None):
    """E_HYM(r)/(4 pi^2) = [Delta + A_1^2/r + B_1^2 r]/N for the topological class of the split background, its minimiser r = |A_1/B_1|,
    and the excess E_split - E_HYM at r (zero iff polystable there)."""
    stacks = stacks if stacks is not None else STACK_MODELS[model]; top = split_topology(stacks); N = top["rank"]; A1, B1 = top["c1"]
    E_hym = (top["Delta"] + sp.Integer(A1) ** 2 / R + sp.Integer(B1) ** 2 * R) / N
    E_split = flux_energy(model=model, stacks=stacks)["E_over_4pi2"]
    rr = sp.nsimplify(r)
    return {**top, "E_HYM_over_4pi2": sp.simplify(E_hym), "r_hym": (sp.Rational(abs(A1), abs(B1)) if B1 else None), "E_HYM_min_over_4pi2": sp.Rational(top["Delta"] + 2 * abs(A1 * B1), N) if B1 else None,
            "excess_split_minus_HYM_at_r": sp.simplify((E_split - E_hym).subs(R, rr)), "excess_polynomial": sp.factor(sp.simplify((E_split - E_hym) * R * N))}

def einstein_frame_scaling():
    """Exponents of the two classical terms in the 4D Einstein-frame potential V(r, Vol) after reducing an 8D gauge theory with gravity on S:
    flux  ~ E(r) Vol^-2  (8D-frame value volume-independent, Weyl factor (Vol_0/Vol)^2); curvature  ~ +48 pi M_8^6 r^-1/2 Vol^-3/2  (hyperbolic
    curve, int_X R = 4 pi chi = -96 pi).  Both positive, both decreasing in Vol: no volume stationary point; shape stabilised at fixed volume."""
    vol, c_flux, c_curv = sp.symbols("Vol c_flux c_curv", positive=True)
    V = c_flux * (sp.Symbol("A", positive=True) / R + sp.Symbol("B", positive=True) * R) * vol ** -2 + c_curv * R ** sp.Rational(-1, 2) * vol ** sp.Rational(-3, 2)
    dV = sp.diff(V, vol)
    return {"V": V, "dV_dVol": sp.simplify(dV), "volume_stationary_point": False, "reason": "every term positive with a negative power of Vol",
            "shape_stationary_at_fixed_volume": True, "flux_exponent": -2, "curvature_exponent": sp.Rational(-3, 2), "curvature_r_exponent": sp.Rational(-1, 2),
            "chi_X0_143": -24, "int_R_X": -96 * sp.pi, "statement": "decompactification runaway; no cosmological constant selected by the classical flux + curvature terms"}

def landscape_report(model="S2", ratios=(sp.Rational(1, 3), sp.Rational(1, 2), 1, 2, 3)):
    """E_split(r), E_HYM(r) and the excess at the census ratios, with the exact minimisers."""
    fe = flux_energy(model=model); hf = hym_floor(model=model)
    rows = [(rr, fe["E_over_4pi2"].subs(R, rr), hf["E_HYM_over_4pi2"].subs(R, rr), sp.simplify((fe["E_over_4pi2"] - hf["E_HYM_over_4pi2"]).subs(R, rr))) for rr in ratios]
    return {"model": model, "rows_r_Esplit_EHYM_excess": rows, "r_star_split": fe["r_star"], "E_star_split": fe["E_star_over_4pi2"], "floor_sum_blocks": fe["floor_over_4pi2"],
            "r_hym": hf["r_hym"], "E_HYM_min": hf["E_HYM_min_over_4pi2"], "Delta": hf["Delta"], "polystable_r": fe["polystable_r"]}
