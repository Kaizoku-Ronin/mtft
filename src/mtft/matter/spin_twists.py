"""Matter-04 (2026-10-02): the adjoint question answered through spinors — spin structures on 143a1, a Scherk–Schwarz twist of the
BHV bulk theory, and what it does to adjoint and charged zero modes.

The problem (exact within the BHV dictionary used by Matter-02/03).  With the topological twist the 8D gaugino becomes (0,p)-forms and
the 8D scalar a (2,0)-form, so each unbroken gauge generator carries H^{0,*}(S, O) = (1, 14, 13): one gaugino and 27 adjoint chiral
multiplets on S = X0(143) x 143a1 (h^{0,1} = 13 + 1, h^{0,2} = 13 x 1).  These are birational invariants: every blow-up has the same
numbers.  For SU(3), N = 1 running gives b0 = 3*3 - 27*3 - n_f = -78 at n_f = 6: no asymptotic freedom while they are light.

The spinor lever.  Fermions are sections of the spinor bundle; a Kahler surface's spinors are (0,*)-forms twisted by a square root of
K.  Choices of that square root (with the R-symmetry bundle held fixed) form a torsor over H^1(S, Z2) = Z2^28 = J0(143)[2] x E[2]:
fermions live in Omega^{0,*}(zeta (x) eta) with zeta a 2-torsion line bundle on X0(143) and eta one on 143a1.  BHV takes zeta = eta = O
(supersymmetry preserved).  Any other choice makes the would-be preserved supercharge antiperiodic around a cycle: Scherk–Schwarz.
Bosons are tensors; they do not see zeta or eta.  Adjoint fermion zero modes are H^*(C, zeta) (x) H^*(E, eta) (Kunneth; pullback to the
blow-up changes nothing):
    eta = O, zeta = O    : (1, 14, 13)  — 28 per generator (BHV)
    eta = O, zeta != O   : (0, 12, 12)  — 24 (a nontrivial flat bundle on a genus-13 curve keeps g - 1 = 12)
    eta != O (any zeta)  : (0, 0, 0)    — none: on an elliptic curve a nontrivial flat bundle has no cohomology at all.
So only the genus-one factor can remove every adjoint fermion, and it does so for each of its three EVEN spin structures — the nontrivial
2-torsion points of 143a1.  Those are the roots of the 2-division polynomial 4x^3 - 4x^2 - 4x - 7 (y = -1/2), irreducible over Q with
discriminant 16 Delta = -2^4 * 11 * 13^2: they are permuted by S3 over the field with quadratic resolvent Q(sqrt(-11)).  The fourth spin
structure (eta = O, periodic-periodic) is the odd one and is what BHV uses.

What survives (EXACT, cohomology).  (i) Family chirality: eta is torsion, so c1(L).(2R - K) is unchanged — e - 24b on the blow-up; the
three families from the exceptional flux remain.  (ii) Every charged FERMION block with torsion-free E-degree b = 0 (all Matter-02 blocks)
sees eta (x) N nontrivial and keeps only the jets at the blown-up point: h^1 = (9q^2 + 3q)/2, h^0 = h^2 = 0 — the curve (genus-13) part
disappears for quarks, leptons and higgsinos alike, while the bosonic Higgs keeps every section.  (iii) Bosonic adjoints are untouched:
1 gauge field + 14 Wilson-line + 13 Phi complex scalars per generator at tree level, now without supersymmetric protection.
One loop (STRUCTURAL; computed only in the flat-circle form here): Scherk–Schwarz along a cycle is a thermal boundary condition.  For the
Wilson line theta along that cycle the potential is V = -sum_n [B_n - (-1)^n F_n] cos(n theta); with matched spectra (B_n = F_n > 0) only odd
n survive, V = -2 sum_{n odd} B_n cos(n theta), minimised at theta = 0 with V''(0) = 2 sum_{n odd} n^2 B_n > 0.  The other 26 bosonic
adjoints are expected to get positive thermal-type masses; computing them on the blow-up of X0(143) x 143a1 is the next calculation.
Below that scale QCD is SU(3) with six quark flavours: b0 = 11 - 2*6/3 = 7 > 0.

Status: zero-mode counts EXACT (Kunneth, Riemann–Roch, flat bundles on curves); 2-torsion arithmetic EXACT; one-loop sign STRUCTURAL in the
flat-circle form; the full one-loop spectrum, the Higgs mass competition (tree-level slope term vs positive one-loop term) and the index-forced
jet pairs are OPEN.  Physical scales are not assigned."""
from fractions import Fraction
import math
import sympy as sp

G_CURVE = 13
A_INVARIANTS = (0, -1, 1, -1, -2)          # 143a1: y^2 + y = x^3 - x^2 - x - 2


# ------------------------------------------------------------------ 143a1: its four spin structures
def two_division_polynomial():
    """4x^3 + b2 x^2 + 2 b4 x + b6 for 143a1, its discriminant (= 16 Delta), irreducibility and the 2-torsion points."""
    a1, a2, a3, a4, a6 = A_INVARIANTS; x = sp.Symbol("x")
    b2 = a1 * a1 + 4 * a2; b4 = 2 * a4 + a1 * a3; b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    delta = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    f = sp.expand(4 * x ** 3 + b2 * x ** 2 + 2 * b4 * x + b6); P = sp.Poly(f, x)
    disc = sp.discriminant(f, x)
    roots = [complex(r) for r in P.nroots(n=20)]
    return {"polynomial": f, "Delta": delta, "discriminant": disc, "disc_equals_16_Delta": disc == 16 * delta,
            "irreducible_over_Q": P.is_irreducible, "galois_group": "S3" if (P.is_irreducible and not sp.sqrt(disc).is_rational) else "C3 or reducible",
            "quadratic_resolvent": f"Q(sqrt({_squarefree_part(disc)}))",
            "two_torsion_points": [(r, -0.5) for r in roots], "y_of_2_torsion": Fraction(-1, 2)}

def _squarefree_part(n):
    n = int(n); s = -1 if n < 0 else 1
    for p, k in sp.factorint(abs(n)).items():
        if k % 2: s *= p
    return s

def elliptic_spin_structures():
    """The four spin structures of 143a1 as theta characteristics: O (odd, h^0 = h^1 = 1) and the three nontrivial 2-torsion bundles (even,
    no cohomology)."""
    pts = two_division_polynomial()["two_torsion_points"]
    rows = [{"label": "eta_0 = O", "parity": "odd", "boundary_conditions": "periodic-periodic", "h0": 1, "h1": 1, "point": None}]
    for i, p in enumerate(pts, 1):
        rows.append({"label": f"eta_{i}", "parity": "even", "boundary_conditions": "antiperiodic on at least one cycle", "h0": 0, "h1": 0, "point": p})
    return rows


# ------------------------------------------------------------------ zero modes
def flat_cohomology_curve(genus, trivial):
    """(h0, h1) of a degree-zero line bundle on a curve: (1, g) if trivial, (0, g - 1) otherwise (no sections; Riemann–Roch)."""
    return (1, genus) if trivial else (0, genus - 1)

def kunneth(hc, he):
    """Cohomology of a box product on C x E from (h0, h1) on each factor: degrees 0, 1, 2."""
    return (hc[0] * he[0], hc[0] * he[1] + hc[1] * he[0], hc[1] * he[1])

def adjoint_fermion_zero_modes(zeta_trivial=True, eta_trivial=True, genus=G_CURVE):
    """H^*(S, zeta (x) eta) for the fermion twist zeta on X0(143), eta on 143a1 (BHV: both trivial).  Unchanged by blowing up."""
    h = kunneth(flat_cohomology_curve(genus, zeta_trivial), flat_cohomology_curve(1, eta_trivial))
    return {"h": h, "total": sum(h), "gaugino": h[0], "adjoint_chiral_fermions": h[1] + h[2]}

def bosonic_adjoint_zero_modes(genus=G_CURVE):
    """Per unbroken generator, independent of the spin structure: gauge field H^0(O) = 1; Wilson lines H^{0,1} = g + 1; Phi in H^{2,0} = g * 1."""
    return {"vector": 1, "wilson_line_scalars": genus + 1, "phi_scalars": genus, "complex_scalars": 2 * genus + 1}

def charged_fermion_block(q, eta_trivial=True):
    """Fermion zero modes of a Matter-02 block L^q (b = 0).  eta trivial: the BHV cohomology (h1 = a-modes, phi = h0(K L^q)).  eta nontrivial:
    eta (x) N is a nontrivial flat bundle on E, the product cohomology vanishes and only exceptional jets survive: h1 = (9q^2 + 3q)/2."""
    if eta_trivial:
        from mtft.matter.collinear import cohomology
        c = cohomology(q); return {"a": c.h1, "phi": c.phi_dim, "chiral": c.h1 + c.phi_dim, "gauge": c.h0}
    jets = (9 * q * q + 3 * q) // 2
    return {"a": jets, "phi": 0, "chiral": jets, "gauge": 0}

def fermion_census(eta_trivial=True):
    """|Q| = 1 fermion census for the Matter-02/03 background (split level), by block, for BHV (eta trivial) or Scherk–Schwarz (eta even)."""
    from mtft.matter.charged_modes import _COMPONENTS
    rows = {"Qminus": [], "Qplus": []}
    for i, ni, qi in _COMPONENTS:
        for j, nj, qj in _COMPONENTS:
            ch = qi - qj
            if ch not in (-1, 1): continue
            rows["Qminus" if ch == -1 else "Qplus"].append({"target": i, "source": j, "q": ni - nj, **charged_fermion_block(ni - nj, eta_trivial)})
    totals = {k: sum(r["chiral"] for r in v) for k, v in rows.items()}
    return {"rows": rows, "totals": totals}

def family_chirality(c1, twist_torsion=True):
    """Net chirality c1(L).(2R - K~) with R = zeta (x) eta torsion: e - 24 b on the blow-up (identical to BHV)."""
    from mtft.surface.blowup import net_chirality
    return net_chirality(c1)


# ------------------------------------------------------------------ one loop and running
def ss_wilson_potential(theta, n_max=400, power=5, matched=True):
    """Flat-circle (Hosotani/Scherk–Schwarz) form: bosons -sum B_n cos(n theta), antiperiodic fermions +sum (-1)^n F_n cos(n theta),
    B_n = F_n = n^-power.  Returns (V_B, V_F, V)."""
    VB = -sum(math.cos(n * theta) / n ** power for n in range(1, n_max + 1))
    VF = sum((-1) ** n * math.cos(n * theta) / n ** power for n in range(1, n_max + 1)) if matched else 0.0
    return VB, VF, VB + VF

def ss_mass_squared_at_symmetric_point(n_max=400, power=5):
    """V''(0) = 2 sum_{n odd} n^2 B_n with B_n = n^-power (> 0)."""
    return 2 * sum(n * n / n ** power for n in range(1, n_max + 1, 2))

def qcd_b0(n_flavours=6, adjoint_chirals=0, supersymmetric=False):
    """One-loop SU(3) coefficient b0 (positive = asymptotically free).  N = 1: 3*3 - 3*n_adj - n_f; non-supersymmetric: 11 - 2 n_f / 3."""
    return 9 - 3 * adjoint_chirals - n_flavours if supersymmetric else Fraction(11) - Fraction(2 * n_flavours, 3)

def ss_fermion_total(stack_e, ranks):
    """Total charged fermion zero modes under the Scherk–Schwarz twist for split stacks with E-degree b = 0, exceptional fluxes e_i and ranks
    r_i (stacks separated by gauge Wilson lines where classes coincide): sum_{i<j} r_i r_j (e_i - e_j)^2 (each unordered pair: T(m) + T(m-1) = m^2
    jets), with net chiral part sum r_i r_j |e_i - e_j|.  EXACT given the twist."""
    names = list(stack_e); tot = 0; net = 0
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            d = stack_e[names[a]] - stack_e[names[b]]; w = ranks[names[a]] * ranks[names[b]]
            tot += w * d * d; net += w * abs(d)
    return {"fermion_zero_modes": tot, "block_net_chiral": net}

RANKS = {"c": 3, "L": 2, "B": 1, "A1": 1, "A2": 1}
MATTER02_E = {"c": 0, "L": -3, "B": 3, "A1": -6, "A2": 9}            # collinear family, charges (0, 1, -1, 2, -3) x (1, 0, -3)
MINIMAL_E = {"c": 0, "L": -3, "B": 3, "A1": 0, "A2": 3}              # charges (0, 1, -1, 0, -1) x (1, 0, -3): A1 ~ c and A2 ~ B in class

def minimal_flux_assignment():
    """Among D-flat b = 0 assignments at eps/A_E = 1/3 (exceptional fluxes all multiples of 3) with exact SM net indices (phi_L = -3, phi_B = 3,
    phi_A1 + phi_A2 = 3, phi_c = 0), the Scherk–Schwarz fermion total 288 is minimal (phi_A1 in {0, 3}); the Matter-02 collinear choice gives 1152.
    The minimal classes coincide pairwise (A1 with c, A2 with B): gauge Wilson lines on 143a1 must separate them, or the gauge group enhances to
    U(4) x U(2) x U(2) (Pati–Salam type) — untested."""
    from mtft.surface.blowup import sm_phi_conditions, slope, ample
    out = []
    for a in range(-9, 13, 3):
        e = {"c": 0, "L": -3, "B": 3, "A1": a, "A2": 3 - a}
        rec = sm_phi_conditions(e); tot = ss_fermion_total(e, RANKS)
        classes = {k: (-v // 3, 0, v) for k, v in e.items()}                  # a_i = -e_i/3 makes every slope vanish at eps = A_E/3
        flat = all(slope(c, 1, 1, sp.Rational(1, 3)) == 0 for c in classes.values()) and ample(1, 1, sp.Rational(1, 3))
        out.append({"phi_A1": a, "exact_SM": rec["exact_SM"], "d_flat": flat, **tot,
                    "coinciding_classes": sorted({tuple(sorted(k for k in classes if classes[k] == cl)) for cl in classes.values() if list(classes.values()).count(cl) > 1})})
    best = min(out, key=lambda r: r["fermion_zero_modes"])
    return {"scan": out, "minimum": best["fermion_zero_modes"], "matter02": ss_fermion_total(MATTER02_E, RANKS)["fermion_zero_modes"]}

def report():
    cen_b, cen_s = fermion_census(True), fermion_census(False)
    return {"spin_structures": elliptic_spin_structures(), "two_division": {k: v for k, v in two_division_polynomial().items() if k != "two_torsion_points"},
            "adjoint_fermions": {"BHV": adjoint_fermion_zero_modes(True, True), "curve_twist_only": adjoint_fermion_zero_modes(False, True),
                                 "scherk_schwarz_on_143a1": adjoint_fermion_zero_modes(True, False)},
            "bosonic_adjoints": bosonic_adjoint_zero_modes(), "fermion_census_BHV": cen_b["totals"], "fermion_census_SS": cen_s["totals"],
            "b0": {"SUSY_with_27_adjoints": qcd_b0(6, 27, True), "SUSY_no_adjoints": qcd_b0(6, 0, True), "QCD_below_SS_scale": qcd_b0(6)},
            "ss_mass_squared_flat_circle": ss_mass_squared_at_symmetric_point(), "family_chirality_example": family_chirality((-1, 0, 3))}
