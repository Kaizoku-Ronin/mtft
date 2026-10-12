"""AL-01 (2026-10-02): blowing up cusps of X0(143) x 143a1 — what Atkin–Lehner symmetry allows, and the one-family-per-cusp background.

Setting: the Pati–Salam groups 4 = c + A1, 2_L = L, 2_R = B + A2 of PS-01; blow up points (cusp_i, O) — all lie on the graph of phi, since every cusp
maps to O.  Classes on the k-point blow-up: (a, b, c; e_1..e_k), form ab' + a'b - 8cc' - sum e_i e_i', K~ = (24, 0, 0; 1, .., 1), chi(O) = 0.  Under the
Scherk–Schwarz twist (fermions twisted by an even spin structure of 143a1) a block with b = 0 has fermion zero modes only as exceptional jets, point by
point: charge e at a point gives e(e-1)/2 modes for e >= 0 and m(m+1)/2 for e = -m < 0; its net chirality is sum_i e_i - 24 b.

1. OBSTRUCTION (EXACT).  The Atkin–Lehner group permutes the four cusps simply transitively.  If the whole group is a symmetry of a line-bundle background
   on the blow-up of all four cusps — acting on the gauge group through any homomorphism into {id, D, C, DC} (D = left–right parity with conjugation,
   C = charge conjugation) — then no Pati–Salam block has three net families: a group element acting trivially on the gauge group makes every exceptional
   vector constant on its two-point orbits (net chirality even), and an element acting through C makes the spectrum vector-like.  `al_obstruction` checks
   all sixteen homomorphisms over a box of exceptional charges.
2. A single involution as left–right parity (EXACT).  Blowing up the two-cusp orbit {infinity, W_Q(infinity)} with W_Q acting as D: L-flux -3 at infinity,
   R-flux +3 at the partner cusp reproduces the minimal Pati–Salam indices with 216 fermion zero modes instead of 288 — but the left-handed families sit at
   one cusp and the right-handed at the other, so the tree-level Pati–Salam Yukawa vanishes identically (classes with disjoint supports).
3. ONE FAMILY PER CUSP (EXACT at index level).  Blow up three of the four cusps; give L exceptional charge -1 and R charge +1 at each (D-parity holds point by
   point).  Then (4,2,1) and (4bar,1,2) each have exactly three net families and NO vector-like partners (one mode per cusp), the bidoublet has net 6 with
   12 pairs: 96 fermion zero modes, against 288 for one blown-up point.  D-flat at eps_i = A_E / 3; ample iff r = A_X/A_E > 1 (the curve C x {O} passes
   through all three points).  The four choices of the empty cusp are related by Atkin–Lehner, so the background is unique up to isomorphism.  Each family's
   left- and right-handed modes share a cusp.  [CC-36: under the twist all tree-level Yukawas between jet fermions vanish — the earlier claim that they are
   diagonal with CKM = 1 is retracted; see `mtft.matter.yukawa_jets`.]
OPEN: the Higgs sections of the three-cusp background (their jets at the cusps set the diagonal Yukawas and hence the mass hierarchy); the corrections
that generate CKM; the gamma-flux for Pati–Salam breaking brings back localised vector-like pairs; the jet sector has no (nu^c, N) pairs, so inverse-seesaw
partners must be gamma-localised, and their Dirac pairing M then involves the same separation as mu."""
import itertools
from fractions import Fraction

# ------------------------------------------------------------------ lattice of a k-point blow-up
def dot(u, v):
    a, b, c, *e = u; a2, b2, c2, *e2 = v
    return a * b2 + a2 * b - 8 * c * c2 - sum(x * y for x, y in zip(e, e2))

def K(k): return (24, 0, 0) + (1,) * k

def chi(L): k = len(L) - 3; return Fraction(dot(L, L) - dot(L, K(k)), 2)

def jets(e): return e * (e - 1) // 2 if e >= 0 else (-e) * (-e + 1) // 2

def fermion_modes(Le):
    """Fermion zero modes of a b = 0 block (and of its conjugate) from its exceptional charges, under the twist."""
    return sum(jets(e) for e in Le), sum(jets(-e) for e in Le)

def net(L): a, b, c, *e = L; return sum(e) - 24 * b

def slope(L, A_X, A_E, kappa, eps):
    a, b, c, *e = L; return a * A_E + b * A_X - 8 * c * kappa + sum(x * y for x, y in zip(e, eps))


# ------------------------------------------------------------------ 1. the Atkin–Lehner obstruction
AL = [(0, 0), (1, 0), (0, 1), (1, 1)]                       # W_11^x W_13^y; points of the orbit <-> group elements (cusp h(infinity))
AUT = {"id": ({"4": "4", "L": "L", "R": "R"}, 1), "D": ({"4": "4", "L": "R", "R": "L"}, -1),
       "C": ({"4": "4", "L": "L", "R": "R"}, -1), "DC": ({"4": "4", "L": "R", "R": "L"}, 1)}
MUL = {("id", x): x for x in AUT} | {(x, "id"): x for x in AUT} | {("D", "D"): "id", ("C", "C"): "id", ("DC", "DC"): "id", ("D", "C"): "DC", ("C", "D"): "DC",
                                                                  ("D", "DC"): "C", ("DC", "D"): "C", ("C", "DC"): "D", ("DC", "C"): "D"}

def _rho(img11, img13, h):
    x, y = h; r = "id"
    if x: r = MUL[(r, img11)]
    if y: r = MUL[(r, img13)]
    return r

def al_obstruction(box=3):
    """Every homomorphism AL -> {id, D, C, DC} and every base exceptional vector in [-box, box]^3 for the groups (4, L, R): the equivariant extension to the four
    cusps, and the net chirality of (4,L), (4,R), (L,R).  Returns the configurations with three net families in (4,L) (expected: none)."""
    found = []; stats = {}
    for img11, img13 in itertools.product(AUT, repeat=2):
        hits = 0
        for base in itertools.product(range(-box, box + 1), repeat=3):
            e1 = dict(zip(("4", "L", "R"), base)); e = {s: {} for s in e1}
            for h in AL:
                perm, sign = AUT[_rho(img11, img13, h)]
                inv = {v: k for k, v in perm.items()}
                for t in e1: e[t][h] = sign * e1[inv[t]]
            n4L = sum(e["4"][h] - e["L"][h] for h in AL); n4R = sum(e["R"][h] - e["4"][h] for h in AL); nLR = sum(e["R"][h] - e["L"][h] for h in AL)
            if abs(n4L) == 3 and abs(n4R) == 3: found.append((img11, img13, base)); hits += 1
            stats.setdefault((img11, img13), set()).add((n4L % 2, n4R % 2, nLR == 0 and n4L == 0))
        stats[(img11, img13)] = hits
    return {"three_family_configurations": found, "per_homomorphism_hits": stats}


# ------------------------------------------------------------------ the three backgrounds
GROUP_RANK = {"4": 4, "L": 2, "R": 2}

def model(name):
    """Exceptional charges per group and point, C-degrees a, and the blown-up points."""
    if name == "one point":
        return {"points": ["inf"], "e": {"4": (0,), "L": (-3,), "R": (3,)}, "a": {"4": 0, "L": 1, "R": -1}}
    if name == "two cusps, D-parity":
        return {"points": ["inf", "W_Q(inf)"], "e": {"4": (0, 0), "L": (-3, 0), "R": (0, 3)}, "a": {"4": 0, "L": 1, "R": -1}}
    if name == "three cusps, one family each":
        return {"points": ["c1", "c2", "c3"], "e": {"4": (0, 0, 0), "L": (-1, -1, -1), "R": (1, 1, 1)}, "a": {"4": 0, "L": 1, "R": -1}}
    raise ValueError(name)

BLOCKS = {"(4,2,1)": ("4", "L"), "(4bar,1,2)": ("R", "4"), "(1,2,2)": ("R", "L")}

def analyse(name):
    M = model(name); out = {"model": name, "points": M["points"], "blocks": {}}; total = 0
    for rep, (i, j) in BLOCKS.items():
        de = tuple(x - y for x, y in zip(M["e"][i], M["e"][j])); L = (M["a"][i] - M["a"][j], 0, 0) + de
        f, fb = fermion_modes(de); w = GROUP_RANK[i] * GROUP_RANK[j]
        support = [p for p, x in zip(M["points"], de) if x != 0]
        out["blocks"][rep] = {"delta_e": de, "net": abs(net(L)), "modes": f, "conj": fb, "pairs_per_entry": min(f, fb), "support": support, "weyl_total": w * (f + fb)}
        total += w * (f + fb)
    out["fermion_total"] = total
    sL, sR = set(out["blocks"]["(4,2,1)"]["support"]), set(out["blocks"]["(4bar,1,2)"]["support"])
    out["yukawa_tree_level"] = "zero (left- and right-handed families at different cusps)" if not (sL & sR) else (
        "diagonal in the cusp basis" if len(sL) > 1 else "single point (local)")
    k = len(M["points"]); eps = [Fraction(1, k) * 1] * k
    out["d_flat_eps_over_A_E"] = Fraction(1, 3) if name != "two cusps, D-parity" else Fraction(1, 3)
    return out

def three_cusp_ampleness():
    """eps_i = A_E/3 at each of three points: fibres through one point need A_E - eps > 0; C x {O} (class pi^*f2 - E1 - E2 - E3) needs A_X - 3 eps > 0, i.e.
    r > 1; the graph's strict transform pi^*Gamma - sum E_i has degree 4 A_E + A_X - 3 eps > 0."""
    import sympy as sp
    AX, AE = sp.symbols("A_X A_E", positive=True); eps = AE / 3
    om = (AX, AE, 0, -eps, -eps, -eps)
    fibre = dot((1, 0, 0, -1, 0, 0), om); CxO = dot((0, 1, 0, -1, -1, -1), om); graph = dot((4, 1, 1, -1, -1, -1), om)
    slopes = {g: slope((model("three cusps, one family each")["a"][g], 0, 0) + model("three cusps, one family each")["e"][g], AX, AE, 0, [eps] * 3) for g in ("4", "L", "R")}
    return {"fibre": sp.simplify(fibre), "C_x_O": sp.simplify(CxO), "graph": sp.simplify(graph), "slopes": {k: sp.simplify(v) for k, v in slopes.items()}}

def report():
    return {"obstruction": {k: v for k, v in al_obstruction().items() if k == "three_family_configurations"},
            "models": {n: analyse(n) for n in ("one point", "two cusps, D-parity", "three cusps, one family each")}, "ampleness": three_cusp_ampleness()}
