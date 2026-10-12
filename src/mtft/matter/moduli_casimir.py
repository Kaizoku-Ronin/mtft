"""MS-01 (2026-10-05): can the Scherk–Schwarz Casimir energy fix the curve's size?  Only if one sign comes out right — and a first count says it doesn't.

1. The structure (EXACT, one loop).  With the twist on the torus only, bosons and fermions differ only on odd windings w (chi_F(w) = -1), where each degree of freedom
   counts once with weight 2.  Heat-kernel expansion along the curve, Tr e^{-t Delta_C} = A_X/(4 pi t) + a_1 + ..., and Poisson resummation on the torus give
   V_J = -K [ (384/pi) N A_X A_E^-3 S0 + 128 a_1 A_E^-2 S1 ],  S_k = sum over odd windings of the unit-area lattice of |w|^-(2k+6),
   with N the twisted degrees of freedom and a_1 their summed curvature coefficient.  In the Einstein frame (divide by Vol^2 = (A_X A_E)^2) this is
   V_E = alpha / A_X + beta / A_X^2 with alpha < 0 (the thermal-like leading term) and beta of sign -a_1.
2. The condition (EXACT).  A minimum in A_X exists iff beta > 0, i.e. a_1 < 0, at A_X / A_E = (2 pi / 3)(|a_1| / N)(S1 / S0).  The curve's Euler characteristic
   chi(X_0(143)) = -24 enters a_1 with the sign that matters: for scalars a_1 = N chi / 6 < 0 and the minimum sits at (8 pi / 3)(S1 / S0) A_E.
3. The numbers (NUMERICAL).  For 143a1's unit-area lattice, S1 = 2.107219 (P = 1/2) and 3.238510 (the complex pair) — exactly Matter-05's one-loop sums, a cross-check
   of the twist convention — and S1/S0 = 1.326 and 1.081.  The scalar count would put the curve at A_X = 11.1 A_E (real) or 9.06 A_E (pair): the same order as
   V_cb's 18–36 (SP-2), a factor two to four short.
4. The sign (STRUCTURAL, first count).  The 8D vector multiplet is not eight scalar pairs: per real component, scalars give chi/6, the curve 1-form -2chi/3 in total,
   spinors (D^2 = Laplacian + R/4) -chi/12.  Summing bosons and fermions gives a_1 = -chi/3 = +8 > 0: the wrong sign, no minimum, the curve collapses.  The matter
   blocks' flux terms enter a_1 at the same order (through the curve degrees), so the verdict waits on the full twisted spectrum.
5. The FI direction (STRUCTURAL).  The blow-up sizes enter only through local curvature at the next order (a_2), so at this order the Casimir energy does not fix
   the FI against GS-1's relaxation either.
MS-02 (2026-10-05): the full twisted spectrum.
6. Dolbeault pairs (EXACT, heat kernel + Hodge theory).  With the topological twist every field of a block is a form on the curve valued in the block's bundle L
   (degree d), and the fields come in pairs: (0,0) with (0,1).  The two Laplacians share their nonzero spectrum, so their traces differ by the index d + 1 - g:
   a_1(0,0) = chi/6 + d/2, a_1(0,1) = -chi/3 - d/2, sum -chi/6.  The flux cancels pair by pair.  Bosons and fermions of a block share these operators (the twist acts
   on the torus only), each block has four complex pairs: a_1 = -4 chi/3 = +32 for N = 16, a_1/N = -chi/12 = +2 — whatever the block's flux, twist or Wilson line.
7. No minimum (EXACT at this order).  The curve-size-independent term follows the total curvature c_1(S~).[omega] = chi(C) A_E - sum eps_k, negative for genus >= 2
   and made more negative by every blow-up.  So beta/alpha > 0 for every block and every twist: V = alpha/A_X + beta/A_X^2 has no minimum.  In the jet vacuum
   (all blocks thermal-signed) the curve collapses.  A sphere (chi > 0) would allow a minimum; a genus-13 curve does not.  MS-4's first count is confirmed.
8. What this leaves (STRUCTURAL).  Nothing in the model stabilises the Kaehler moduli: the D-terms relax (GS-1) and the one-loop Casimir energy runs.  The
   consistent reading is a fixed geometry — the areas as parameters of the background, as in MTFT's arithmetic picture.  Then GS-1's relaxation does not occur,
   the Casimir energy is a constant vacuum energy, and A_X/A_E stays a parameter that V_cb measures (SP-2) unless a canonical prescription fixes it.
OPEN: a canonical prescription for A_X/A_E; non-perturbative or gravitational stabilisation if the moduli are dynamical."""
import math


def odd_sums(spin="1/2", R=60):
    from mtft.research.theta_torus import tau_143a1
    tau = complex(tau_143a1()); s = 1 / math.sqrt(tau.imag)
    phase = {"1/2": lambda m, n: n % 2, "tau/2": lambda m, n: m % 2, "(1+tau)/2": lambda m, n: (m + n) % 2}[spin]
    S0 = S1 = 0.0
    for m in range(-R, R + 1):
        for n in range(-R, R + 1):
            if (m, n) == (0, 0) or phase(m, n) == 0: continue
            r2 = abs(s * (m + n * tau)) ** 2
            S0 += r2 ** -4; S1 += r2 ** -3
    return S0, S1

def minimum_ratio(a1_over_N, spin="1/2"):
    """A_X / A_E at the minimum, or None when a_1 >= 0 (no minimum)."""
    if a1_over_N >= 0: return None
    S0, S1 = odd_sums(spin)
    return (2 * math.pi / 3) * abs(a1_over_N) * S1 / S0

def a1_counts(chi=-24):
    bosons = 6 * (chi / 6) + (-2 * chi / 3)                 # six scalar-like real components on the curve, one 1-form (two components)
    fermions = 8 * (-chi / 12)                               # eight real spinor components, D^2 = Laplacian + R/4
    return {"scalar_toy_per_dof": chi / 6, "bosons": bosons, "fermions": fermions, "vector_multiplet": bosons + fermions}

def dolbeault_a1(chi, d):
    g = 1 - chi // 2
    a0 = chi / 6 + d / 2
    a1 = a0 - (d + 1 - g)
    return a0, a1, a0 + a1

def block_a1(chi):
    """Summed curvature coefficient (real degrees of freedom) of one block: four complex Dolbeault pairs, bosons and fermions."""
    return 4 * 2 * (-chi / 6), 16

def total_curvature(chi, A_E, eps):
    return chi * A_E - sum(eps)

def minimum_possible(chi):
    a1, N = block_a1(chi); return a1 / N < 0

def report():
    out = {sp: odd_sums(sp) for sp in ("1/2", "tau/2", "(1+tau)/2")}
    return {"sums": out, "ratio": {sp: out[sp][1] / out[sp][0] for sp in out},
            "scalar_minimum": {sp: minimum_ratio(-24 / 6, sp) for sp in out}, "sym_a1": -(-24) / 3, "sym_minimum": minimum_ratio(8 / 16)}
