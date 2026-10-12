"""FS-01 (2026-10-02): scanning SU(2)_R-breaking fluxes at the single blown-up cusp — the toolkit's stack-level extension.

Stacks c (3), A1 (1), L (2), B (1), A2 (1) keep their Matter-04 exceptional charges (0, 0, -3, 3, 3), so every index and every fermion jet count is fixed.  Two
chirality-neutral knobs remain per stack:
  * the curve degree a_s (degree on X0(143)): under the twist, fermion zero modes of b = 0 blocks are jets whatever a_s is (H^*(S, M (x) eta) = 0 for every M);
    a_s only moves slopes and bosonic curve cohomology;
  * the gamma-charge c_s: moves slopes, localises 4 Delta c^2 extra pairs per entry orientation, and makes a Higgs entry with Delta c != 0 gamma-localised.
Slopes (A_E = 1): mu_s = a_s + e_s eps - 8 c_s kappa; stack FI xi_s = n_s mu_s.  Pati–Salam -> SM by ONE Standard-Model-singlet condensate requires xi_L = xi_B = 0
and xi supported on {A1, A2} (or on {colour + A1, A2} when SU(4) is unbroken by the flux), with 0 < eps < A_E.

1. Curve degrees alone suffice (EXACT).  With every c_s = 0 there are D-flat single-singlet backgrounds, e.g. a = (c 0, A1 1, L 1, B -1, A2 -2) at eps = A_E/3, or the
   SU(4)-preserving a = (1, 1, 1, -1, -5): Pati–Salam breaking needs no gamma-flux, no gamma-localised extras (288 fermion zero modes), and keeps every Higgs entry
   global.  This supersedes PS-01's 'the gamma direction is needed' (CC-38).
2. TK-1 lifted (EXACT).  With curve degrees free, gamma-charge +-1 on nu^c, an SU(4)-preserving flux, a global (L, B) Higgs and a D-flat single-singlet condensate are
   compatible; `scan` lists them.
3. A gamma-localised Higgs feeds torus charge linearly (STRUCTURAL).  EPS-01 item 2 removes torus dependence from DENSITIES; a Higgs entry with gamma-charge Delta
   enters the Yukawa linearly as a plane wave e^{i k(x) w} with k proportional to Delta phi', so its sector gets eps_sector ~ (2.1 ... 6.7) |Delta| / sqrt(r)
   (the EPS-01 estimate with c -> Delta).  The up and down sectors use different Higgs entries once SU(2)_R is broken: (L, A2) for up, (L, B) for down, so
   eps_d / eps_u = |Delta_d / Delta_u| when both are gamma-localised — a ratio of integers to compare with 2.30.  A global (Delta = 0) Higgs leaves its sector
   with the condensate-misalignment eps of EPS-01 only.
4. The Higgs-charge rule (EXACT).  With SU(4) kept by the flux (c_colour = c_A1), tracelessness gives 4 c_4 = -(4 c_L - Delta_u - Delta_d), so
   Delta_u + Delta_d = 0 (mod 4).  The admissible eps_d/eps_u = |Delta_d/Delta_u| near the data are 3 (charges 1, 3) and 5/3 (3, 5); a ratio of 2 needs (4, 8),
   7/3 needs (6, 14).  The minimal-charge prediction is eps_d/eps_u = 3 against 2.30 (the proportionality carries order-one sector factors).
   With gamma-charge +-1 on nu^c as well, tracelessness becomes 8 c_4 -+ 4 + 3 Delta_u - Delta_d = 0: equal-sign ratio 3 is impossible, opposite-sign ratio 3 needs
   charges (2, 6) (e.g. c = (1, 1, 0, -6, 2), 3104 fermion zero modes).  The cheap options therefore separate: ratio 3 without a charged nu^c, or a charged nu^c
   with ratio 1, 5/3 or a global up-type Higgs.
`scan` enumerates traceless (a, c) in a box and returns the D-flat single-singlet backgrounds with their charges, eps-ratio prediction and cost."""
from dataclasses import dataclass
from fractions import Fraction
import itertools
import math

STACKS = ("c", "A1", "L", "B", "A2")
N = {"c": 3, "A1": 1, "L": 2, "B": 1, "A2": 1}
EXC = {"c": 0, "A1": 0, "L": -3, "B": 3, "A2": 3}
DATA_RATIO = 2.3049


@dataclass(frozen=True)
class StackFlux:
    a: tuple            # curve degrees (c, A1, L, B, A2)
    g: tuple            # gamma-charges (c, A1, L, B, A2)

    def A(self, s): return self.a[STACKS.index(s)]
    def G(self, s): return self.g[STACKS.index(s)]

    def traceless(self):
        return sum(N[s] * self.A(s) for s in STACKS) == 0 and sum(N[s] * self.G(s) for s in STACKS) == 0

    def su4(self): return self.A("c") == self.A("A1") and self.G("c") == self.G("A1")
    def su2r(self): return self.A("B") == self.A("A2") and self.G("B") == self.G("A2")

    def slope(self, s, eps, kappa): return self.A(s) + EXC[s] * eps - 8 * self.G(s) * kappa

    def solve_kahler(self):
        """Exact (eps, kappa) with xi_L = xi_B = 0 (and xi_c = 0 when SU(4) is broken); returns None if inconsistent or out of range.  A one-parameter family is
        reported with kappa = 0 when allowed, else the smallest |kappa| making eps fall in (0, 1)."""
        aL, aB, gL, gB = Fraction(self.A("L")), Fraction(self.A("B")), self.G("L"), self.G("B")
        # aL - 3 eps - 8 gL k = 0 ;  aB + 3 eps - 8 gB k = 0  ->  (aL + aB) = 8 k (gL + gB)
        sols = []
        if gL + gB != 0:
            k = (aL + aB) / (8 * (gL + gB)); sols = [k]
        elif aL + aB == 0:
            sols = [Fraction(0)] + [Fraction(m, 16) for m in range(-16, 17) if m]
        else:
            return None
        for k in sols:
            eps = (aL - 8 * gL * k) / 3
            if not (0 < eps < 1): continue
            if not self.su4() and self.A("c") - 8 * self.G("c") * k != 0: continue
            mu = {s: self.slope(s, eps, k) for s in STACKS}
            pair = ("4", "A2") if self.su4() else ("A1", "A2")
            if self.su4() and mu["c"] != mu["A1"]: continue
            if mu["A2"] == 0: continue
            return {"eps": eps, "kappa": k, "pair": pair, "mu": mu}
        return None

    def charges(self):
        G = self.G
        return {"nu_c": G("A1") - G("A2"), "delta_u": G("L") - G("A2"), "delta_d": G("L") - G("B"), "Q": G("c") - G("L"), "u_c": G("c") - G("A2"),
                "d_c": G("c") - G("B"), "lepton_doublet": G("A1") - G("L"), "e_c": G("A1") - G("B"), "leptoquark": G("c") - G("A1")}

    def fermion_total(self):
        extra = 0
        for s, t in itertools.combinations(STACKS, 2):
            d = self.G(s) - self.G(t); extra += N[s] * N[t] * 2 * 4 * d * d
        return 288 + extra

    def eps_prediction(self):
        ch = self.charges(); du, dd = abs(ch["delta_u"]), abs(ch["delta_d"])
        if du and dd: return {"ratio": dd / du, "mode": "both Higgs entries gamma-localised"}
        if dd and not du: return {"ratio": math.inf, "mode": "down gamma-localised, up from condensate misalignment"}
        if du and not dd: return {"ratio": 0.0, "mode": "up gamma-localised, down from condensate misalignment"}
        return {"ratio": 1.0, "mode": "both global: eps_u = eps_d (EPS-01)"}

    def summary(self):
        k = self.solve_kahler()
        if k is None: return None
        p = self.eps_prediction()
        return {"a": self.a, "g": self.g, "eps": k["eps"], "kappa": k["kappa"], "pair": k["pair"], "su4": self.su4(), "su2r": self.su2r(),
                "charges": self.charges(), "fermions": self.fermion_total(), "eps_ratio": p["ratio"], "eps_mode": p["mode"]}


def scan(abox=3, gbox=2):
    """All traceless (a, c) with entries in [-abox, abox] x [-gbox, gbox] that are D-flat with one Standard-Model-singlet condensate (A2 determined by tracelessness)."""
    out = []
    rng_a = range(-abox, abox + 1); rng_g = range(-gbox, gbox + 1)
    for ac, aA1, aL, aB in itertools.product(rng_a, rng_a, range(1, abox + 1), rng_a):
        aA2 = -(3 * ac + aA1 + 2 * aL + aB)
        if abs(aA2) > 3 * abox: continue
        for gc, gA1, gL, gB in itertools.product(rng_g, rng_g, rng_g, rng_g):
            gA2 = -(3 * gc + gA1 + 2 * gL + gB)
            if abs(gA2) > 2 * gbox: continue
            s = StackFlux((ac, aA1, aL, aB, aA2), (gc, gA1, gL, gB, gA2)).summary()
            if s: out.append(s)
    return out

def best(rows, nu_one=None, su4=None, global_down=None, both_local=None, key="cost"):
    sel = [r for r in rows if (nu_one is None or (abs(r["charges"]["nu_c"]) == 1) == nu_one) and (su4 is None or r["su4"] == su4)
           and (global_down is None or (r["charges"]["delta_d"] == 0) == global_down)
           and (both_local is None or (r["charges"]["delta_u"] != 0 and r["charges"]["delta_d"] != 0) == both_local)]
    if key == "ratio":
        sel.sort(key=lambda r: (abs(math.log(r["eps_ratio"] / DATA_RATIO)) if 0 < r["eps_ratio"] < math.inf else 99, r["fermions"]))
    else:
        sel.sort(key=lambda r: (r["fermions"], sum(abs(x) for x in r["a"])))
    return sel

def report(abox=3, gbox=2):
    rows = scan(abox, gbox)
    minimal = best([r for r in rows if all(x == 0 for x in r["g"])])
    tk1 = best(rows, nu_one=True, su4=True, global_down=True)
    ratio = best(rows, su4=True, both_local=True, key="ratio")
    return {"n_solutions": len(rows), "curve_degree_only": minimal[:3], "tk1_lifted": tk1[:3], "ratio_best": ratio[:5]}


def higgs_charge_rule(max_charge=8):
    """Admissible (Delta_u, Delta_d) with SU(4) kept: Delta_u + Delta_d = 0 mod 4, both nonzero; the smallest-charge representative of each ratio."""
    reps = {}
    for du in range(1, max_charge + 1):
        for dd in range(1, 2 * max_charge + 1):
            if (du + dd) % 4: continue
            r = Fraction(dd, du)
            if r not in reps or du + dd < sum(reps[r]): reps[r] = (du, dd)
    return reps

def families(abox=3, gbox=3, rows=None):
    """The four background families of the scan, cheapest representative each (pass `rows` to reuse a scan)."""
    rows = scan(abox, gbox) if rows is None else rows
    pick = lambda sel: min(sel, key=lambda r: (r["fermions"], sum(abs(x) for x in r["a"]))) if sel else None
    return {"minimal (curve degrees only)": pick([r for r in rows if all(x == 0 for x in r["g"])]),
            "ratio 3 (SU(4) kept, Higgs charges 1 and 3)": pick([r for r in rows if r["su4"] and abs(r["charges"]["delta_u"]) == 1 and abs(r["charges"]["delta_d"]) == 3]),
            "neutrino-optimal, up Higgs global": pick([r for r in rows if r["su4"] and r["charges"]["delta_u"] == 0 and r["charges"]["delta_d"] != 0 and abs(r["charges"]["nu_c"]) == 1]),
            "neutrino-optimal, both Higgs localised": pick([r for r in rows if r["su4"] and r["charges"]["delta_u"] and r["charges"]["delta_d"] and abs(r["charges"]["nu_c"]) == 1])}


def ratio_three_with_charged_nu_c():
    """Explicit witness: SU(4) kept, nu^c charge -1, Higgs charges (-2, 6): ratio 3, and its D-flat completion by curve degrees at kappa = 0."""
    g = (1, 1, 0, -6, 2)
    for a in itertools.product(range(-3, 4), repeat=4):
        ac, aA1, aL, aB = a
        if ac != aA1: continue
        aA2 = -(3 * ac + aA1 + 2 * aL + aB)
        s = StackFlux((ac, aA1, aL, aB, aA2), g).summary()
        if s: return s
    return None
