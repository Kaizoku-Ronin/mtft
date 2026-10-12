"""Matter-05 (2026-10-02): one-loop masses of the 27 bosonic adjoints per gauge generator under the Scherk–Schwarz twist of Matter-04 —
13 curve Wilson lines (H^{0,1}(X0(143))), 1 torus Wilson line (H^{0,1}(143a1)), 13 Phi scalars (H^{2,0}) — closing the adjoint question at one loop.

Setting.  BHV-twisted 8D Yang–Mills on R^4 x C x E (C = X0(143), E = 143a1, any product Kahler metric), fermions twisted by an even spin structure
eta on E.  One loop: V = -(1/32 pi^2) int dt t^-3 Str e^{-t M^2}.  Twisted supersymmetry on C pairs the C-spectra of bosons and fermions for every flat
connection and every holomorphic Phi, so Str = n Tr_C(t) P(t) with n real bosonic degrees of freedom per charged direction and
    P(t) = Tr_E(t) - Tr_{E,eta}(t) = (A_E / 2 pi t) sum_{w in Lambda_E, eta(w) = -1} exp(-|w|^2 / 4t) > 0      (Poisson resummation on the torus).
SIGN THEOREM (EXACT on the product; any product metric; gauge-multiplet loops).  V(a, theta, phi) >= V(0) with strictly positive Hessian:
  * curve Wilson lines a: every bosonic species on C has, up to a-independent index constants, the nonzero spectrum of the scalar Laplacian twisted by the
    flat bundle L_a (isospectrality on curves); by the method of images Tr e^{-t Delta_a} = sum_h cos(a.h) N_h(t) with N_h(t) > 0 (the heat kernel of the
    scalar Laplacian on the universal cover is positive), so Tr decreases away from a = 0 and V increases.  The Hessian is a sum over closed geodesics
    gamma of X0(143) (hyperbolic classes of Gamma_0(143)) weighted by (alpha(a).[gamma])^2 > 0.
  * Phi: a holomorphic vev adds the potential |[Phi, X]|^2 = |alpha(phi)|^2 |omega(z)|^2 >= 0 to every charged species; by min-max every eigenvalue rises,
    Tr_C decreases, V increases.  Holds for matter loops too (every field is a U(8)-adjoint block).
  * torus Wilson line theta: V = -c n sum_w (1 - eta(w)) cos(alpha(theta).w) W_w with W_w > 0; only odd windings survive; Hessian 2 sum (alpha.w)^2 W_w > 0.
    Holds for matter loops too when no background gauge Wilson lines are present.
  OPEN: the matter-loop (flux-block) contribution to the curve Wilson-line masses — magnetic Laplacians carry phases in their image sums, so the sign is not
  fixed by this argument.  The blow-up breaks the product near the exceptional curve; corrections are bounded by the area fraction eps/A_C = 1/(3r).
CLOSED FORMS (leading order in A_E / R_C^2, small-t heat kernel on C; unit-area lattice sums S_eta = sum_{eta(w)=-1} |w|^-6, T_eta = sum w w^T |w|^-8):
    m_Phi^2   = (n / 2 pi^4) g_4^2 (A_C / A_E^2) S_eta  sum_alpha |alpha(e)|^2          (complex Phi, V = m^2 |phi|^2)
    M_theta^2 = (6 n / pi^4) g_4^2 (A_C / A_E^2) T_eta  sum_alpha |alpha(e)|^2          (real 2-vector, V = theta^T M^2 theta / 2)
  The curve Wilson lines are suppressed relative to these: their image sums need a closed geodesic of C; the shortest has trace 4 in Gamma_0(143),
  length 2 arccosh 2 = 2.6339 (curvature -1), e.g. [[19, -2], [143, -15]]; hyperbolic area 56 pi.  They are the lightest adjoints, the more so the
  smaller 143a1 is compared with the curvature radius of X0(143).
143a1 (tau = 1/2 + 1.0232745927 i, rhombic) — the even spin structures split 1 + 2, as the 2-torsion does over Q:
  antiperiodic on the tau-cycle  (real 2-torsion point P = 1/2):       S = 2.107219,  T eigenvalues 0.473605, 1.633614
  antiperiodic on the 1- or (1+tau)-cycle (conjugate pair tau/2, (1+tau)/2): S = 3.238510, T eigenvalues 0.770478, 2.468032
  The conjugate pair puts antiperiodicity on the shortest lattice vector: adjoints 1.537 times heavier (Phi).
QCD (SU(3), six flavours): b0 = 11 - 4 - (1/2) N_real_adjoint_scalars: 7 below every adjoint mass; -6 with the 26 real curve Wilson-line octets light;
-20 with all 54.  Asymptotic freedom holds below the lightest adjoint; above it QCD runs the other way up to the compactification scale, so that window
must stay short — a constraint on the sizes: 143a1 should not be much smaller than X0(143)'s curvature radius (r = A_C / A_E of order the hyperbolic
area 56 pi is natural; the blow-up correction 1/(3r) is then about 0.2 per cent).

Status: sign theorem EXACT on the product (gauge loops); closed forms EXACT at leading order in the stated regime; lattice sums NUMERICAL (converged to
1e-8); curve Wilson-line magnitudes and matter-loop signs OPEN; physical scales unassigned."""
import math
import numpy as np


def tau_143a1():
    from mtft.research.theta_torus import tau_143a1 as _t
    return complex(_t())

SPIN_STRUCTURES = {(0, 1): {"two_torsion_point": "1/2 (real)", "antiperiodic_cycle": "tau"},
                   (1, 0): {"two_torsion_point": "tau/2 (complex)", "antiperiodic_cycle": "1"},
                   (1, 1): {"two_torsion_point": "(1+tau)/2 (complex)", "antiperiodic_cycle": "1+tau"}}

def lattice_sums(eps, N=300, tau=None):
    """Unit-area lattice Lambda = (Z + tau Z)/sqrt(Im tau); odd windings eta(m + n tau) = (-1)^(eps1 m + eps2 n) = -1.  Returns S = sum |w|^-6 and
    T = sum w w^T |w|^-8 (2x2) — the arithmetic factors of m_Phi^2 and M_theta^2."""
    tau = tau if tau is not None else tau_143a1(); s = math.sqrt(tau.imag); S = 0.0; T = np.zeros((2, 2))
    for m in range(-N, N + 1):
        for n in range(-N, N + 1):
            if (eps[0] * m + eps[1] * n) % 2 == 0: continue
            w = (m + n * tau) / s; r2 = w.real * w.real + w.imag * w.imag
            S += r2 ** -3; v = np.array([w.real, w.imag]); T += np.outer(v, v) / r2 ** 4
    return S, T

def spin_structure_table(N=300):
    out = {}
    for eps, info in SPIN_STRUCTURES.items():
        S, T = lattice_sums(eps, N); ev = np.linalg.eigvalsh(T)
        out[eps] = {**info, "S": S, "T_eigenvalues": tuple(ev), "theta_over_phi": tuple(12 * e / S for e in ev)}
    return out

def mass_coefficients(n=8):
    """K_Phi, K_theta in m^2 = K g_4^2 (A_C / A_E^2) (lattice factor) (group factor)."""
    return {"K_Phi": n / (2 * math.pi ** 4), "K_theta": 6 * n / math.pi ** 4}

def systole_X0_143():
    """Smallest |trace| > 2 in Gamma_0(143): a(t - a) = 1 mod 143 solvable.  Trace 3 has no solution; trace 4 does."""
    N = 143
    for t in range(3, 10):
        sols = [a for a in range(N) if (a * (t - a) - 1) % N == 0]
        if sols:
            a = sols[0]; d = t - a; b = (a * d - 1) // N
            return {"trace": t, "length": 2 * math.acosh(t / 2), "example": ((a, b), (N, d)), "area": math.pi / 3 * 168}
    return None


# ------------------------------------------------------------------ the ingredients of the sign theorem, checked on solvable models
def circle_twisted_heat_trace(t, a, ell=1.0, kmax=200):
    """Tr e^{-t Delta_a} on a circle of length ell with flat twist a: sum_k exp(-t (2 pi (k + a)/ell)^2) (maximal at a = 0: method of images)."""
    return sum(math.exp(-t * (2 * math.pi * (k + a) / ell) ** 2) for k in range(-kmax, kmax + 1))

def torus_P(t, eps=(0, 1), N=None, tau=None):
    """P(t) = (1/2 pi t) sum_{eta(w)=-1} exp(-|w|^2/4t) on the unit-area 143a1 lattice (A_E = 1): positive; -> 1 as t -> infinity (the boson zero mode,
    absent for antiperiodic fermions).  The winding cutoff grows with sqrt(t)."""
    tau = tau if tau is not None else tau_143a1(); s = math.sqrt(tau.imag); tot = 0.0
    N = N if N is not None else max(40, int(14 * math.sqrt(t)) + 10)
    for m in range(-N, N + 1):
        for n in range(-N, N + 1):
            if (eps[0] * m + eps[1] * n) % 2 == 0: continue
            w = (m + n * tau) / s; tot += math.exp(-(w.real ** 2 + w.imag ** 2) / (4 * t))
    return tot / (2 * math.pi * t)

def minmax_shift(potential, n=200, L=1.0):
    """Lowest eigenvalues of -d^2/dx^2 on a circle with and without a nonnegative potential (finite differences): every eigenvalue rises."""
    h = L / n; D = np.zeros((n, n))
    for i in range(n): D[i, i] = 2 / h ** 2; D[i, (i + 1) % n] = -1 / h ** 2; D[i, (i - 1) % n] = -1 / h ** 2
    x = np.arange(n) * h; V = np.diag([potential(xx) for xx in x])
    return np.linalg.eigvalsh(D)[:20], np.linalg.eigvalsh(D + V)[:20]

def qcd_b0_windows():
    """SU(3), six flavours: 11 - 4 - N_real_adjoint_scalars * (1/6) * 3."""
    b = lambda nreal: 11 - 4 - nreal * 0.5
    return {"below_all_adjoints": b(0), "curve_Wilson_lines_light_only": b(26), "all_27_complex_light": b(54)}

def report(N=300):
    return {"spin_structures": spin_structure_table(N), "mass_coefficients": mass_coefficients(), "systole": systole_X0_143(), "b0": qcd_b0_windows(),
            "blowup_correction_at_r_56pi": 1 / (3 * 56 * math.pi)}
