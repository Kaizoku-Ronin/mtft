"""mtft.parent — candidate parent actions as data, and the parent-action ledger's gates as verdicts (PV-01, 2026-09-29).

A ParentCandidate is a finite data set: spacetime dimension, gauge group, fermion content in the parent group (representation kinds
with signed multiplicities), a supersymmetry flag, the internal background (split stacks on X0(143) x 143a1 with bidegrees), and the
family/Higgs identification when the candidate declares one.  `gates(candidate)` runs every ledger row
(docs/SM/PARENT_ACTION_REQUIREMENTS.md, A.1-A.10 and C) that the data allow and returns a verdict per row: True/False when the row is
decided EXACTLY, None when it is OPEN for that candidate, each with its witness.  There is no "viable" boolean: `ledger_table` prints
the scoreboard.  Gate 4, route 2 (one parent, three internal modes) is the target the catalogue is measured against.

Catalogue (research records, not endorsements): the three-stack adjoint 8D gauge theory on S1 and S2 (the ledger's status column),
the E7 6D record M3, the C3X 6D candidate, and three B1 exemplars found by `anomaly8d.b1_scan` — U(8) with Lambda^2(8) + 8 x 8 on a
five-stack split background: Yukawa-admissible with both Higgs blocks slope-free at r = 2 but colour never matchable (YH); colour
matchable at r = 3/2 but Higgs massive at every r (YC); and, from the box-4 scan, both loci present at DIFFERENT shapes, r_H = 2 and
r_c = 5/2 (HC).  The B1 trilemma these exhibit — never Y, H and C at one ratio — is a finite exact statement of the box-3/box-4 scans
(`anomaly8d.b1_trilemma`), not a theorem."""
from dataclasses import dataclass, field
from typing import Optional, Tuple
import sympy as sp
from .research import product_surface as PS, anomaly8d as A8, vacuum_energy as VE, unified_parent as UP, gravitational_anomaly as GA, parents as PR

Stack = Tuple[str, int, Tuple[int, int]]

@dataclass(frozen=True)
class ParentCandidate:
    model_id: str
    dimension: int
    group: str
    content: Tuple[Tuple[str, int], ...]                 # (representation kind, signed multiplicity) — 'adjoint', 'fund', 'antisym2', 'sym2', 'u1(q)'
    supersymmetric: Optional[bool] = None                # None = undeclared (the pending programme decision)
    stacks: Optional[Tuple[Stack, ...]] = None           # split background on X0(143) x 143a1
    family_blocks: Optional[dict] = None                 # {'Q': (a, b), 'u': (a, b), 'd': (a, b)} block bidegrees, package convention (the block, not its conjugate)
    higgs_blocks: Optional[dict] = None                  # {'u': (a, b), 'd': (a, b)} adjoint blocks carrying the Higgs
    notes: str = ""

def _verdict(passed, witness, status="EXACT"): return {"pass": passed, "witness": witness, "status": status}

def gates(c: ParentCandidate):
    """Every ledger row the candidate's data decide.  Rows: 1 odd chirality; 2 Yukawa between chiral families; 3 u^c/d^c in different
    summands; 4 SM unbroken at M_KK (colour slope); 5 no destabilising twisted colour line; 6 tachyon-free background; 7 light Higgs;
    8 strong CP; 9 anomaly integrality (6D: tensor p2; 8D: quintic-free + Green–Schwarz fields); 10 hypercharge normalisation;
    C vacuum energy of the background (flux landscape)."""
    g = {}
    kinds = dict(c.content); adjoint_only = set(kinds) <= {"adjoint"}
    # 1. odd chirality
    if c.dimension == 8 and adjoint_only:
        g["1_odd_chirality"] = _verdict(False, "adjoint matter on the spin surface: net chirality even, zero for the spin twist (chirality_parity_theorem)")
    elif c.dimension == 8 and c.stacks:
        odd = [(rep, r["content"], r["bidegree"], r["index"]) for rep in ("fund", "antisym2", "sym2") if kinds.get(rep) for r in A8.stack_blocks(rep, stacks=c.stacks) if r["odd"]]
        g["1_odd_chirality"] = _verdict(bool(odd), {"odd_index_blocks": odd, "lemma": A8.split_stack_parity_lemma()["statement"]})
    elif c.dimension == 6:
        g["1_odd_chirality"] = _verdict(True, "index on a curve is a degree (R2C-06): odd degrees exist", "EXACT (curve)")
    else: g["1_odd_chirality"] = _verdict(None, "no background declared", "OPEN")
    # 2. Yukawa between chiral families (bidegree rule)
    fb = c.family_blocks or {}
    if fb.get("Q") and fb.get("u") and fb.get("d"):
        yu = PS.yukawa_bidegree_rule(fb["Q"], fb["u"]); yd = PS.yukawa_bidegree_rule(fb["Q"], fb["d"])
        g["2_yukawa"] = _verdict(yu["allowed"] and yd["allowed"], {"up": yu, "down": yd})
    elif fb.get("Q") and fb.get("u") and adjoint_only:
        yu = PS.yukawa_bidegree_rule(fb["Q"], fb["u"])
        g["2_yukawa"] = _verdict(False, {"up_rule": yu, "reason": "adjoint blocks are vector-like pairs (parity theorem): no chiral family to couple; d^c absent"})
    elif adjoint_only and c.stacks: g["2_yukawa"] = _verdict(False, "split polystable vacua have no gauge-vertex Yukawa between chiral families (one-ray lemma); families vector-like anyway")
    else: g["2_yukawa"] = _verdict(None, "no family blocks declared", "OPEN")
    # 3. u^c, d^c in different summands
    if fb.get("u") and fb.get("d"):
        g["3_uc_dc_summands"] = _verdict(fb["u"] != fb["d"] or "distinct singlet stacks" in c.notes, {"u": fb["u"], "d": fb["d"], "rule": A8.split_stack_parity_lemma()["statement"]})
    elif adjoint_only and c.stacks: g["3_uc_dc_summands"] = _verdict(False, "three stacks: one antitriplet species only (hypercharge lemma forces >= 2 singlet classes)")
    else: g["3_uc_dc_summands"] = _verdict(None, "no assignment", "OPEN")
    # 4. SM unbroken at M_KK: colour slope matchable
    if c.stacks:
        dc = next(d for n, r, d in c.stacks if r == 3); rest = [(n, r, d) for n, r, d in c.stacks if r != 3]
        w = (sum(r * (d[0] - dc[0]) for _, r, d in rest), sum(r * (d[1] - dc[1]) for _, r, d in rest))
        rc = sp.Rational(-w[0], w[1]) if w[1] and sp.Rational(-w[0], w[1]) > 0 else None
        g["4_colour_slope"] = _verdict(rc is not None, {"w": w, "r_colour": rc, "note": "necessary for an SU(3)-preserving polystable vacuum"})
    else: g["4_colour_slope"] = _verdict(None, "no background", "OPEN")
    # 5. destabilising twisted colour line (proved for S1; other backgrounds open)
    if c.stacks == PS.STACK_MODELS["S1"]: g["5_no_twisted_colour_destabiliser"] = _verdict(False, PS.s1_iterated_extension_instability(2))
    else: g["5_no_twisted_colour_destabiliser"] = _verdict(None, "theorem proved for the S1 stacks only", "OPEN")
    # 6. tachyon-free: polystable at some r
    if c.stacks:
        fe = VE.flux_energy(stacks=c.stacks); poly = fe["polystable_r"]
        g["6_tachyon_free"] = _verdict(bool(poly) and poly != [], {"polystable_r": poly, "slopes": fe["slopes_units_A_E"], "census": "block_census(r, stacks=...) for the Morse index"})
    else: g["6_tachyon_free"] = _verdict(None, "no background", "OPEN")
    # 7. light Higgs at leading order
    if c.higgs_blocks:
        loci = {k: (sp.Rational(-H[0], H[1]) if H[1] and sp.Rational(-H[0], H[1]) > 0 else None) for k, H in c.higgs_blocks.items()}
        common = set(v for v in loci.values()) if loci else set()
        g["7_light_higgs"] = _verdict(all(v is not None for v in loci.values()) and len(common) == 1, {"slope_free_loci": loci})
    else: g["7_light_higgs"] = _verdict(None, "no Higgs block declared", "OPEN")
    # 8. strong CP
    g["8_strong_cp"] = _verdict(None, {"curve": PS.cm_point_cp_structure()["conclusion"], "torus": "CP-conserving (torus_cp_test)"}, "OPEN (lift must be phase-aligned or an axion exists)")
    # 9. anomalies
    if c.dimension == 6:
        n = sum(m * UP.ADJOINT_DIM.get(c.group, 0) for k, m in c.content if k == "adjoint")
        g["9_anomaly"] = _verdict(GA.tensor_integrality(n)["integral"] if n else None, GA.tensor_integrality(n) if n else "content not in adjoint-dimension table", "EXACT (CC-33)" if n else "OPEN")
    elif c.dimension == 8:
        N = sum(r for _, r, _ in c.stacks) if c.stacks else None
        spec = []
        for k, m in c.content:
            if k == "adjoint": spec.append((A8.adjoint("A", N), m))
            elif k == "fund": spec.append((A8.fund("A", N), m))
            elif k == "antisym2": spec.append((A8.antisym2("A", N), m))
            elif k == "sym2": spec.append((A8.sym2("A", N), m))
        I = A8.I10(spec); d = A8.gs_decomposition(I)
        g["9_anomaly"] = _verdict(d["cancellable"], {"irreducible_quintic": d["irreducible_quintic"], "green_schwarz_fields_needed": d["green_schwarz_fields_needed"], "axion_X8": d["axion_X8"], "two_form_factors": d["two_form_factors"]},
                                  "EXACT (quintic-free; Green–Schwarz completion assumed, global anomalies open)")
    else: g["9_anomaly"] = _verdict(None, "dimension not covered", "OPEN")
    # 10. hypercharge normalisation
    g["10_hypercharge_normalisation"] = _verdict(True if adjoint_only and c.stacks else None, "CC-26 (three-stack)" if adjoint_only and c.stacks else "to be recomputed for the candidate's U(1)s", "EXACT" if adjoint_only and c.stacks else "OPEN")
    # C. vacuum energy of the background
    if c.stacks:
        fe = VE.flux_energy(stacks=c.stacks); hf = VE.hym_floor(stacks=c.stacks)
        g["C_vacuum_energy"] = _verdict(None, {"E_split_over_4pi2": fe["E_over_4pi2"], "r_star": fe["r_star"], "E_star": fe["E_star_over_4pi2"], "Delta": hf["Delta"], "r_hym": hf["r_hym"], "E_HYM_min": hf["E_HYM_min_over_4pi2"],
                                               "volume": VE.einstein_frame_scaling()["statement"]}, "EXACT (flux landscape); the vacuum itself OPEN")
    else: g["C_vacuum_energy"] = _verdict(None, "no background", "OPEN")
    g["susy"] = _verdict(None, "undeclared" if c.supersymmetric is None else c.supersymmetric, "programme decision")
    return g

# ---------------------------------------------------------------- catalogue
THREE_STACK_8D_S1 = ParentCandidate("THREE_STACK_8D_S1", 8, "U(6) = U(3) x U(2) x U(1)", (("adjoint", 1),), None, PS.STACK_MODELS["S1"],
                                    {"Q": (3, 1), "u": (1, -3), "d": None}, {"u": (-4, 2)}, "the ledger's status column; families vector-like (parity theorem)")
THREE_STACK_8D_S2 = ParentCandidate("THREE_STACK_8D_S2", 8, "U(6) = U(3) x U(2) x U(1)", (("adjoint", 1),), None, PS.STACK_MODELS["S2"],
                                    {"Q": (-3, 1), "u": (1, 3), "d": None}, {"u": (2, -4)}, "mirror of S1")
_B1_STACKS_YH = (("c", 3, (0, 1)), ("L", 2, (3, 0)), ("s1", 1, (-1, 2)), ("s2", 1, (-1, 2)), ("s3", 1, (0, 1)))
_B1_STACKS_YC = (("c", 3, (0, 1)), ("L", 2, (3, 0)), ("s1", 1, (3, -2)), ("s2", 1, (3, -2)), ("s3", 1, (0, 1)))
B1_U8_YH = ParentCandidate("B1_U8_LAMBDA2_YH", 8, "U(8) = U(3) x U(2) x U(1)^3", (("antisym2", 1), ("fund", 8)), None, _B1_STACKS_YH,
                           {"Q": (3, 1), "u": (-1, 3), "d": (-1, 3)}, {"u": (4, -2), "d": (4, -2)},
                           "anomaly8d.b1_scan exemplar: exact SM content, Yukawa-admissible, both Higgs slope-free at r = 2, colour never matchable; u^c and d^c on distinct singlet stacks of equal bidegree")
B1_U8_YC = ParentCandidate("B1_U8_LAMBDA2_YC", 8, "U(8) = U(3) x U(2) x U(1)^3", (("antisym2", 1), ("fund", 8)), None, _B1_STACKS_YC,
                           {"Q": (3, 1), "u": (3, -1), "d": (3, -1)}, {"u": (0, 2), "d": (0, 2)},
                           "anomaly8d.b1_scan exemplar: exact SM content, Yukawa-admissible, colour matchable at r = 3/2, Higgs massive at every r; distinct singlet stacks")
_B1_STACKS_HC = (("c", 3, (3, 0)), ("L", 2, (0, 1)), ("s1", 1, (-4, 3)), ("s2", 1, (-4, 3)), ("s3", 1, (3, 0)))
B1_U8_HC = ParentCandidate("B1_U8_LAMBDA2_HC", 8, "U(8) = U(3) x U(2) x U(1)^3", (("antisym2", 1), ("fund", 8)), None, _B1_STACKS_HC,
                           {"Q": (3, 1), "u": (-1, 3), "d": (-1, 3)}, {"u": (4, -2), "d": (4, -2)},
                           "anomaly8d.b1_scan box-4 exemplar: exact SM content, Yukawa-admissible, both Higgs slope-free at r = 2 AND colour matchable at r = 5/2 — two different shapes; distinct singlet stacks")
M3_E7_6D = ParentCandidate("M3_E7_6D", 6, "E7", (("adjoint", 1),), None, None, None, None, UP.M3_RECORD["kind"])
C3X_6D = ParentCandidate("C3X_6D", 6, PR.C3X["global_group"], (), False, None, None, None, "three INPUT copies; " + PR.C3X["status"]["vacuum"])
CATALOGUE = {c.model_id: c for c in (THREE_STACK_8D_S1, THREE_STACK_8D_S2, B1_U8_YH, B1_U8_YC, B1_U8_HC, M3_E7_6D, C3X_6D)}

def ledger_table(candidates=None):
    """Scoreboard: rows = ledger gates, columns = candidates; entries T/F/- (open)."""
    cands = candidates or list(CATALOGUE.values()); res = {c.model_id: gates(c) for c in cands}
    rows = sorted({k for g in res.values() for k in g}, key=lambda s: (s[0].isdigit() and int(s.split("_")[0]) or 99, s))
    sym = lambda v: "T" if v is True else ("F" if v is False else "-")
    lines = ["gate".ljust(34) + "".join(c.model_id[:18].ljust(20) for c in cands)]
    for r in rows: lines.append(r.ljust(34) + "".join(sym(res[c.model_id].get(r, {}).get("pass")).ljust(20) for c in cands))
    return {"table": "\n".join(lines), "verdicts": res}
