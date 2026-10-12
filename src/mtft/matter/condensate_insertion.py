"""CI-01 (2026-10-06): the condensate insertion — the lepton pattern lives in the tails, where the condensate cannot reach.

1. The projector is rank one (EXACT).  The condensate's internal direction is fixed by its bundle (the nu^c entry of the (4bar, 1, 2) block), so Sigma Sigma^dagger =
   |Sigma|^2 P_lepton everywhere.  Every lepton-specific amplitude carries the same density as the SU(4)-singlet condensate channel: the lepton/quark ratio of an entry is
   1 + rho w_ij, with rho the universal ratio of the projector to the singlet coefficient and w_ij the share of the quark entry carried by the condensate channel.
2. Tails are universal (EXACT).  An entry carried by a family's tail involves only zero-mode profiles fixed by flux that commutes with SU(4); its lepton factor is 1
   (times a dressing weight <= 1).  The programme's 2–3 entries are tail-carried (TO-01, FL-01), so the muon's 3.09 forces the 3–3 lepton factor <= 1/1.758 = 0.569:
   y_tau/y_b <= 0.57 at matching against 1.40–1.53 — a factor >= 2.46 (1.76 even with b–tau unification).
3. The entry pattern (NUMERICAL, running-independent).  The 2–3 entries need 1.758 relative to the 3–3 entry and the 1–2 entries 0.959: the lepton tails 1.76 times the
   quark tails with the local entries untouched.  That predicts f_1 f_2 = f_3^2; the data give 0.919, inside m_d/m_s's 8% uncertainty.
4. The condensate channel cannot carry the 2–3 entries (ESTIMATE: gaussian jets of width l_1 = 0.268, holomorphic 1-jet of the condensate at P1, singlet coefficient
   set to zero, the most favourable case).  The entries reach P1 through the third family's tail either way, so the channel's excess over the 1–2 entries is bounded by
   the tail's lever arm in jet units, Lambda = 2d/(n l_1): s[1/2 + sqrt(1/4 + Lambda^2/8)] with s the 1–2 insertion.  The data need Lambda >= 5.9 (f_3 = 1.53) to 7.7
   (1.40), and no Lambda at all below f_3 = 1.043; the programme's n = 3 tail over d = 0.95–1.21 gives 2.4–3.0.
5. Verdict (STRUCTURAL).  The single condensate produces the lepton pattern neither by dressing (QL-01) nor by insertion.  The pattern asks for SU(4) breaking in the zero
   modes themselves — lepton tails 1.76 times quark tails — e.g. a discrete Wilson line along B–L displacing the lepton profiles.  CC-53 supersedes QL-3.
OPEN: the B–L Wilson line; the overall b–tau normalisation stays GC-01's threshold problem."""
import math

L1 = math.sqrt(0.225 / math.pi)          # blow-up radius at P1 (TA-01: eps_1 = 0.225 A_E)
SEPARATION = (0.95, 1.21)                 # P1–P2 separation that V_cb needs (TA-01, DD-01)


def _ql():
    from mtft.matter import quark_lepton_split as Q
    return Q

def r_B():
    return math.sqrt(_ql().z_e())

def r_A():
    Q = _ql(); return math.sqrt(Q.first_over_third() * Q.z_e())

def product_relation():
    return r_A() ** 2

def tail_shortfall(mu):
    """Factor by which tail-carried 2–3 entries miss: the 3–3 lepton factor allowed (1/r_B) against the running y_tau/y_b."""
    return _ql().tau_over_b(mu) * r_B()

def max_excess(s, lam):
    return s * (0.5 + math.sqrt(0.25 + lam * lam / 8))

def lever_needed(f3):
    s = r_A() * f3 - 1
    if s <= 0: return math.inf
    need = (r_B() / r_A() - 1) * (1 + s) / s
    return math.sqrt(8 * ((need - 0.5) ** 2 - 0.25))

def lever_arm(d, n=3, l1=L1):
    return 2 * d / (n * l1)

def threshold_f3():
    return 1 / r_A()

def report():
    Q = _ql(); f3 = [Q.tau_over_b(m) for m in (8.9e11, 2.1e14)]
    return {"r_B": r_B(), "r_A": r_A(), "product": product_relation(), "tail_shortfall": [tail_shortfall(m) for m in (8.9e11, 2.1e14)],
            "lever_needed": [lever_needed(f) for f in f3], "lever_arm_n3": [lever_arm(d) for d in SEPARATION], "f3_floor": threshold_f3()}
