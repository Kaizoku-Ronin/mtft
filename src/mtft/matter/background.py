"""mtft.matter.background — one object for a Pati–Salam background on blow-ups of X0(143) x 143a1 (toolkit, 2026-10-02).

A Background records where points (cusp, O) are blown up, the exceptional charges e, curve degrees a, gamma-charges c and Z/2 Wilson lines (= the spin
structure eta) of the three Pati–Salam groups 4 = colour + A1, L, R = B + A2, and the even spin structure of the Scherk–Schwarz twist.  `analyse()` assembles the
results of Matter-04 ... EPS-01 for it:

  indices (net chirality sum e - 24 b), the twist pattern of every block (YK-01), fermion zero modes (exceptional jets, gamma-localised 4c^2, and for untwisted
  matter the curve copies of YK-01 where modelled), D-flatness of the exceptional data at eps_i = A_E/3, ampleness conditions, the anomaly table (PS-01),
  the tree-level Yukawa verdict (YK-01 jet theorem and classification), the one-loop stability of the relative-U(1) Wilson line (PA-01), the messenger scale
  (PA-01), quark-mixing admissibility (CK-01: families at different cusps mix by at most 1/143) and the family structure (CK-01 / EPS-01).

Statements inherit the classes of the waves they come from (see `mtft.matter.ledger`).  Presets: `PRESETS`."""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
import math

GROUP_RANK = {"4": 4, "L": 2, "R": 2}
CUSPS = ("inf", "0", "1/11", "1/13")
BLOCK_ORIENT = {"(4,2,1)": ("4", "L"), "(4bar,1,2)": ("R", "4"), "(1,2,2)": ("R", "L")}


def _jets(e):
    return e * (e - 1) // 2 if e >= 0 else (-e) * (-e + 1) // 2


@dataclass
class Background:
    name: str
    points: tuple
    e: dict
    a: dict = field(default_factory=lambda: {"4": 0, "L": 1, "R": -1})
    gamma: dict = field(default_factory=lambda: {"4": 0, "L": 0, "R": 0})
    wilson: dict = field(default_factory=lambda: {"4": 0, "L": 0, "R": 0})
    spin: str = "complex pair"
    note: str = ""

    # ---------------------------------------------------------------- per-block data
    def block(self, i, j):
        de = tuple(x - y for x, y in zip(self.e[i], self.e[j]))
        return {"a": self.a[i] - self.a[j], "c": self.gamma[i] - self.gamma[j], "e": de, "wilson": (self.wilson[i] - self.wilson[j]) % 2}

    def twist(self, i, j):
        w = self.block(i, j)["wilson"]
        return {"bosons": "untwisted" if w == 0 else "twisted", "fermions": "untwisted" if (1 + w) % 2 == 0 else "twisted"}

    def net(self, i, j):
        return sum(self.block(i, j)["e"])

    def fermion_modes(self, i, j):
        """(modes, conjugate modes) per rank-one entry; None where the curve copies are not modelled."""
        b = self.block(i, j); tw = self.twist(i, j)["fermions"]; c = b["c"]
        if tw == "twisted" or c != 0:
            loc = 4 * c * c
            return (sum(_jets(x) for x in b["e"]) + loc, sum(_jets(-x) for x in b["e"]) + loc)
        if set(abs(x) for x in b["e"]) <= {1} and abs(b["a"]) == 1 and len(self.points) == 3 and len(set(b["e"])) == 1:
            from mtft.matter.yukawa_jets import untwisted_matter_entry
            u = untwisted_matter_entry(); return (u["R"], u["Rbar"])
        return None

    # ---------------------------------------------------------------- assembled analyses
    def indices(self):
        return {rep: abs(self.net(i, j)) for rep, (i, j) in BLOCK_ORIENT.items()}

    def zero_modes(self):
        out, total = {}, 0
        for rep, (i, j) in BLOCK_ORIENT.items():
            fm = self.fermion_modes(i, j); w = GROUP_RANK[i] * GROUP_RANK[j]
            out[rep] = None if fm is None else w * (fm[0] + fm[1])
            total = None if (total is None or fm is None) else total + out[rep]
        return {"by_block": out, "total": total}

    def d_flat(self, A_E=Fraction(3)):
        eps = A_E / 3
        slopes = {g: self.a[g] * A_E + sum(x * eps for x in self.e[g]) for g in GROUP_RANK}
        return {"slopes_at_eps_AE_over_3": slopes, "d_flat": all(v == 0 for v in slopes.values())}

    def ampleness(self):
        k = len(self.points)
        return {"fibre through one point": "A_E - eps > 0", f"X0(143) x O through {k} point(s)": f"A_X > {k} eps = {Fraction(k, 3)} A_E",
                "graph of phi through the points": f"4 A_E + A_X - {k} eps > 0"}

    def anomalies(self):
        n = self.indices()
        if n != {"(4,2,1)": 3, "(4bar,1,2)": 3, "(1,2,2)": 6}:
            return {"note": "anomaly table implemented for the minimal Pati–Salam content (3, 3, 6) only"}
        from mtft.matter.ps_breaking import anomalies
        A = anomalies(); return {"nonabelian_cancel": A["SU4^3"] == 0 and A["SU4^2 Q4"] == 0, "witten_even": A["Witten_SU2L_doublets"] % 2 == 0 and A["Witten_SU2R_doublets"] % 2 == 0,
                                 "stack_U1s_anomalous": (A["SU2L^2 QL"], A["SU2R^2 QR"], A["SU2L^2 Q4"]) != (0, 0, 0)}

    def yukawa(self):
        higgs = self.twist("R", "L")["bosons"] == "untwisted"; matter = self.twist("4", "L")["fermions"] == "untwisted" or self.twist("R", "4")["fermions"] == "untwisted"
        if higgs and matter:
            return "tree-level Yukawas possible (matter untwisted; YK-01)"
        if not matter:
            return "tree-level Yukawas vanish: all matter fermions are exceptional jets (YK-01 theorem); masses need messenger exchange (PB-01, EPS-01)"
        return "no Phi-type Higgs boson (bidoublet twisted): no superpotential Yukawa"

    def wilson_line_stability(self):
        """Relative U(1) between the quark–lepton group and L, R: only (4,2,1) and (4bar,1,2) are charged (PA-01)."""
        signs = []
        for i, j in (("4", "L"), ("R", "4")):
            t = self.twist(i, j); signs.append(+1 if (t["fermions"], t["bosons"]) == ("twisted", "untwisted") else -1 if (t["fermions"], t["bosons"]) == ("untwisted", "twisted") else 0)
        s = sum(signs)
        return "minimum (stable)" if s > 0 else "maximum (rolls to the jet vacuum)" if s < 0 else "flat at one loop"

    def messenger_scale(self):
        from mtft.matter.path_a import messenger_masses
        m = messenger_masses(); return m["real (P = 1/2)"] if self.spin == "real" else m["complex pair (P = tau/2)"]

    def _supports(self):
        left = [p for p, x in zip(self.points, self.block("4", "L")["e"]) if x != 0]
        right = [p for p, x in zip(self.points, self.block("R", "4")["e"]) if x != 0]
        return left, right

    def ckm_admissible(self):
        left, right = self._supports()
        if set(left) != set(right):
            return {"family_points": (left, right), "admissible": False,
                    "reason": "left- and right-handed families at different cusps: their Yukawa vanishes identically (AL-01)"}
        ok = len(left) == 1
        return {"family_points": left, "admissible": ok,
                "reason": "all families at one point (local mixing)" if ok else "families at different cusps mix by at most 1/143 < |V_cb| (CK-01, CC-37)"}

    def family_structure(self):
        left, right = self._supports(); de = self.block("4", "L")["e"]
        if set(left) != set(right):
            return {"type": "left-handed and right-handed families at different cusps", "texture": "Yukawa zero (disjoint supports, AL-01)"}
        if len(left) == 1:
            m = max(abs(x) for x in de)
            return {"type": "jet-degree triplet at one point", "jet_degree": m - 1, "torus_charges": list(range(m)) if m == 3 else None,
                    "texture": "rank one at leading order; [[1, eps, eps^2], [eps, eps^2, 0], [eps^2, 0, 0]] with eps the condensate misalignment (EPS-01)"}
        yuk = "tree-level, local per cusp" if self.yukawa().startswith("tree-level Yukawas possible") else "tree Yukawas vanish (CC-36)"
        return {"type": "one family per blown-up cusp", "texture": f"{yuk}; inter-cusp mixing at most 1/143 (CC-37)"}

    def higgs(self):
        b = self.block("R", "L"); tw = self.twist("R", "L")["bosons"]
        if tw == "untwisted" and b["c"] == 0:
            return "global torus-constant Phi-type Higgs available (PB-01 mechanism applies)"
        if tw == "untwisted":
            return f"bidoublet carries gamma-charge {b['c']}: the Higgs zero modes are gamma-localised, no global torus-constant Higgs (PB-01 needs revisiting)"
        return "bidoublet bosons twisted: no Phi-type Higgs"

    def localisation(self):
        out = {}
        for rep, (i, j) in BLOCK_ORIENT.items():
            c = self.block(i, j)["c"]
            if c: out[rep] = {"gamma_charge": c, "bosons": f"phi^-1(E[{abs(c)}]) ({4 * c * c} points{', the cusps' if abs(c) == 1 else ''})",
                              "fermions": f"phi^-1({abs(c)}Q = P_eta) ({4 * c * c} points)", "disjoint": True}
        return out

    def analyse(self):
        return {"name": self.name, "points": self.points, "spin": self.spin, "indices": self.indices(), "zero_modes": self.zero_modes(), "d_flat": self.d_flat(),
                "ampleness": self.ampleness(), "anomalies": self.anomalies(), "twists": {rep: self.twist(i, j) for rep, (i, j) in BLOCK_ORIENT.items()},
                "yukawa": self.yukawa(), "higgs": self.higgs(), "wilson_line_stability": self.wilson_line_stability(), "messenger_M_sqrtA_E": self.messenger_scale(),
                "ckm": self.ckm_admissible(), "families": self.family_structure(), "localisation": self.localisation(), "note": self.note}

    def text(self):
        r = self.analyse(); z = r["zero_modes"]
        lines = [f"Background {r['name']}  ({self.note})", f"  blown-up points: {', '.join(r['points'])}   spin structure: {r['spin']}",
                 f"  indices: {r['indices']}", f"  fermion zero modes: {z['by_block']}  total {z['total']}",
                 f"  D-flat at eps = A_E/3: {r['d_flat']['d_flat']}   ampleness: {'; '.join(f'{k}: {v}' for k, v in r['ampleness'].items())}",
                 f"  twists: " + "; ".join(f"{k} bosons {v['bosons']}, fermions {v['fermions']}" for k, v in r["twists"].items()),
                 f"  Yukawa: {r['yukawa']}", f"  Higgs: {r['higgs']}", f"  relative-U(1) Wilson line: {r['wilson_line_stability']}", f"  messenger scale M sqrt(A_E) = {r['messenger_M_sqrtA_E']:.3f}",
                 f"  quark mixing: {'admissible' if r['ckm']['admissible'] else 'EXCLUDED'} — {r['ckm']['reason']}", f"  families: {r['families']['type']} — {r['families']['texture']}"]
        if r["localisation"]:
            lines.append("  gamma-localisation: " + "; ".join(f"{k}: c = {v['gamma_charge']}, bosons {v['bosons']}, fermions {v['fermions']}" for k, v in r["localisation"].items()))
        if isinstance(r["anomalies"], dict) and "nonabelian_cancel" in r["anomalies"]:
            a = r["anomalies"]; lines.append(f"  anomalies: nonabelian cancel {a['nonabelian_cancel']}, Witten even {a['witten_even']}, stack U(1)s anomalous {a['stack_U1s_anomalous']}")
        return "\n".join(lines)


def _ps(name, gamma, note):
    return Background(name, ("inf",), {"4": (0,), "L": (-3,), "R": (3,)}, gamma=gamma, note=note)

PRESETS = {
    "matter04": Background("matter04", ("inf",), {"4": (0,), "L": (-3,), "R": (3,)}, note="canonical jet vacuum: one blown-up cusp (Matter-04, restored by CC-37)"),
    "ps-cheap": _ps("ps-cheap", {"4": 1, "L": -1, "R": -1}, "Pati–Salam breaking by gamma-flux (1, 1, -1, -1, -1), NU-01"),
    "ps-nu-optimal": _ps("ps-nu-optimal", {"4": 1, "L": -2, "R": 0}, "gamma-flux with charge one on nu^c, NU-02"),
    "ps-ps01": _ps("ps-ps01", {"4": 1, "L": 0, "R": -2}, "gamma-flux of PS-01"),
    "two-cusps-dparity": Background("two-cusps-dparity", ("inf", "W_Q(inf)"), {"4": (0, 0), "L": (-3, 0), "R": (0, 3)}, note="AL-01: an AL involution as left–right parity"),
    "three-cusps": Background("three-cusps", ("inf", "0", "1/11"), {"4": (0, 0, 0), "L": (-1, -1, -1), "R": (1, 1, 1)}, note="AL-01 one family per cusp (retired by CC-37)"),
    "matter-untwisted": Background("matter-untwisted", ("inf", "0", "1/11"), {"4": (0, 0, 0), "L": (-1, -1, -1), "R": (1, 1, 1)}, wilson={"4": 1, "L": 0, "R": 0},
                                   note="YK-01 matter untwisting on three cusps (a one-loop maximum, PA-01)"),
}

def preset(name):
    return PRESETS[name]


def ps_symmetric_gamma_tradeoff(box=3):
    """Pati–Salam-symmetric, traceless gamma-fluxes (c_4, c_L, c_R) with 4 c_4 + 2 c_L + 2 c_R = 0 and c_4 != 0.  Returns the gamma-charge of the nu^c entry for
    those that keep the bidoublet free of gamma-flux (a global Higgs, as PB-01 needs).  EXACT: c_L = c_R forces c(nu^c) = 2 c_4, never +-1."""
    rows = []
    for c4 in range(-box, box + 1):
        for cL in range(-2 * box, 2 * box + 1):
            for cR in range(-2 * box, 2 * box + 1):
                if c4 == 0 or 4 * c4 + 2 * cL + 2 * cR != 0: continue
                rows.append({"c": (c4, cL, cR), "nu_c_charge": c4 - cR, "bidoublet_charge": cL - cR})
    global_higgs = [r for r in rows if r["bidoublet_charge"] == 0]
    return {"with_global_higgs": global_higgs, "min_abs_nu_c_charge_with_global_higgs": min(abs(r["nu_c_charge"]) for r in global_higgs),
            "nu_c_charge_one_possible_with_global_higgs": any(abs(r["nu_c_charge"]) == 1 for r in global_higgs)}
