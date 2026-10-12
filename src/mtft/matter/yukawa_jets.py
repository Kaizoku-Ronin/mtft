"""YK-01 (2026-10-02): Yukawa couplings of exceptional-jet families — a vanishing theorem, the only tree-level way out, and its price.

1. THEOREM (EXACT).  Under the Scherk–Schwarz twist every fermion zero mode of a b = 0 block is an exceptional jet: a class delta(s) in H^1(pi^*F(-mE)) from
   0 -> pi^*F(-mE) -> pi^*F -> pi^*F|_{mE} -> 0 with H^*(S, F) = 0.  Its Dolbeault representative is (dbar chi) h / sigma^m — h a local holomorphic extension
   of the jet, chi one cutoff equal to 1 near the exceptional curves, sigma the section of O(E) — so two such classes multiply to (dbar chi ^ dbar chi)(...) = 0
   pointwise (equivalently: both are Cech classes of the two-set cover {S~ - E, tube(E)}, whose Cech H^2 vanishes).  The zero-mode Yukawa
   int psi_1 ^ psi_2 ^ Phi_H is a cup product (it depends only on Dolbeault classes), and with the twist the only fermion zero modes are A-type, so EVERY
   tree-level Yukawa between two jet fermions vanishes — at one blown-up point or at different ones, whatever the Higgs section.  This covers every background
   since Matter-04 (one point, two cusps, three cusps).  Masses must come from non-jet modes, from higher (Massey) products through the Pati–Salam condensate,
   or from loops.
2. CLASSIFICATION (EXACT).  A superpotential Yukawa needs a Phi-type Higgs boson; twisted blocks have none (H^0(K~ x L x N) = 0 for nontrivial flat N), so the
   bidoublet must be untwisted for bosons: N_L = N_R.  Then (4,2,1) and (4bar,1,2) carry the same fermion twist eta N_4 N_L^-1, so the left- and right-handed
   matter are either both jets (Yukawas zero) or both untwisted.  Tree-level Yukawas exist iff N_4 = eta N_L: a Z/2 Wilson line equal to the spin structure on
   the quark–lepton group.  `classify` enumerates every Z/2 assignment and finds exactly this one (up to an overall shift).
3. THE MATTER-UNTWISTED SPECTRUM (EXACT twist pattern).  With that Wilson line: quark and lepton fermions untwisted, their scalar partners twisted (jets only);
   Higgs bosons untwisted, Higgsinos twisted (jets only); gauge bosons untwisted, gauginos twisted (removed) — the Standard-Model pattern of fermionic matter,
   bosonic Higgs and gauge fields.  SU(4) and left–right parity are kept (eta^-1 = eta), indices and slopes do not change.
4. PRICE (EXACT counts on the three-cusp background; genus 13, X0(143) non-hyperelliptic, M_L = O(P) with P a cusp that is not blown up, M_R = M_L^-1,
   evaluation of the constant section at the three cusps of rank 1).  Each matter entry gets h^*(X0(143)) (x) h^*(143a1) copies: 25 vector-like pairs and the
   net 3, so (4,2,1) and (4bar,1,2) carry 424 Weyl fermions each, against 24 in the jet-only background; with the 48 bidoublet jets, 896 in all.
   Global-mode Yukawas are curve-level cup products (the Matter-02 coupling slice is their untwisted ancestor); the light families appear only after the
   vector-like pairs are given masses.
Corrections: CC-36 (`docs/SM/CC36_JET_YUKAWA_RETRACTION.md`)."""
import itertools

GROUPS = ("4", "L", "R")
BLOCKS = {"adj 4": ("4", "4"), "adj L": ("L", "L"), "adj R": ("R", "R"), "(4,2,1)": ("4", "L"), "(4bar,1,2)": ("4", "R"), "(1,2,2)": ("L", "R")}

def twist_table(N):
    """N: group -> 0 (trivial) or 1 (eta) as Z/2 Wilson lines on 143a1.  A block's bosons see N_i N_j^-1; its fermions see eta N_i N_j^-1.
    Returns, per block, whether bosons / fermions are untwisted (keep the product cohomology) or twisted (exceptional jets only)."""
    out = {}
    for name, (i, j) in BLOCKS.items():
        b = (N[i] - N[j]) % 2; f = (1 + N[i] - N[j]) % 2
        out[name] = {"bosons": "untwisted" if b == 0 else "jets", "fermions": "untwisted" if f == 0 else "jets"}
    return out

def tree_yukawa_possible(N):
    t = twist_table(N)
    higgs_ok = t["(1,2,2)"]["bosons"] == "untwisted"
    matter_global = t["(4,2,1)"]["fermions"] == "untwisted" or t["(4bar,1,2)"]["fermions"] == "untwisted"
    return higgs_ok and matter_global

def classify():
    sols = [N for N in (dict(zip(GROUPS, v)) for v in itertools.product((0, 1), repeat=3)) if tree_yukawa_possible(N)]
    shared_twist = all(twist_table(dict(zip(GROUPS, v)))["(4,2,1)"]["fermions"] == twist_table(dict(zip(GROUPS, v)))["(4bar,1,2)"]["fermions"]
                       for v in itertools.product((0, 1), repeat=3) if twist_table(dict(zip(GROUPS, v)))["(1,2,2)"]["bosons"] == "untwisted")
    return {"solutions": sols, "matter_blocks_share_fermion_twist_when_higgs_untwisted": shared_twist}

MATTER_UNTWISTED = {"4": 1, "L": 0, "R": 0}

def sm_pattern():
    t = twist_table(MATTER_UNTWISTED)
    return {"quark_lepton_fermions": t["(4,2,1)"]["fermions"], "squarks_sleptons": t["(4,2,1)"]["bosons"], "higgs_bosons": t["(1,2,2)"]["bosons"],
            "higgsinos": t["(1,2,2)"]["fermions"], "gauge_bosons": t["adj 4"]["bosons"], "gauginos": t["adj 4"]["fermions"]}


# ------------------------------------------------------------------ counts on the three-cusp background
G = 13

def curve_h(deg, effective_degree_two=True):
    """(h0, h1) of a line bundle on X0(143) for the degrees used here: O(-P), O(P) (P not a Weierstrass gap issue for degree one), O(+-2P)."""
    table = {-2: (0, G + 1), -1: (0, G), 0: (1, G), 1: (1, G - 1), 2: (1, G - 2)}
    return table[deg]

def product_h(deg):
    """h^0, h^1, h^2 of M (x) O_E on C x E (E-cohomology (1, 1))."""
    h0, h1 = curve_h(deg); return (h0, h0 + h1, h1)

def blown_up_h(deg, e_per_point, k=3, eval_rank=1):
    """Cohomology of pi^*(M (x) O_E) (x) O(e sum E_i) on the k-point blow-up, untwisted, e in {-1, 0, +1}."""
    h0, h1, h2 = product_h(deg)
    if e_per_point >= 0: return (h0, h1, h2)                                   # pi_* O(E) = O, R^1 pi_* O(E) = 0
    # e = -1: 0 -> L(-E) -> L -> + L|_{E_i} -> 0; evaluation H^0(S, L) -> C^k of rank eval_rank
    return (h0 - eval_rank, h1 + (k - eval_rank), h2)

def untwisted_matter_entry():
    """Hom(L -> 4) = pi^*(M_L^-1 (x) O)(+E) and Hom(4 -> L) = pi^*(M_L (x) O)(-E) with M_L = O(P): chiral multiplets in R: h1(L) + h2(L^-1);
    in Rbar: h1(L^-1) + h2(L)."""
    a = blown_up_h(-1, +1); b = blown_up_h(+1, -1)
    R = a[1] + b[2]; Rbar = b[1] + a[2]
    return {"Hom(L->4)": a, "Hom(4->L)": b, "R": R, "Rbar": Rbar, "net": Rbar - R, "pairs": min(R, Rbar), "per_entry_total": R + Rbar}

def spectrum_comparison():
    e = untwisted_matter_entry()
    untw = {"(4,2,1)": 8 * e["per_entry_total"], "(4bar,1,2)": 8 * e["per_entry_total"], "(1,2,2)": 48}
    jet = {"(4,2,1)": 24, "(4bar,1,2)": 24, "(1,2,2)": 48}
    return {"jet_only": jet, "jet_only_total": sum(jet.values()), "matter_untwisted": untw, "matter_untwisted_total": sum(untw.values())}

def report():
    return {"classification": classify(), "sm_pattern": sm_pattern(), "entry": untwisted_matter_entry(), "comparison": spectrum_comparison()}
