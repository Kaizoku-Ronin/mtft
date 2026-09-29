"""CC-35 (2026-09-29): reflection positivity of the MTFT lattice action — the Polyakov-loop term is NOT covered by the Osterwalder–Seiler
argument, the lemma Papers 13 (§4.2) and 24 (Theorem 2.4, Step 4) use for it is false, and at strong coupling the reflection positivity the
papers define fails.  Tools, the counterexample, and the repair.

The action (mtft.lattice.MTFTAction):  S = S_Wilson + kappa sum_x sum_n a_n(y) (1 - Re Tr P(x)^n / N),  a_n = w_n e^{-2 pi y n},  w_1 = 0.
Across a time-reflection plane the Polyakov factor couples the positive and negative halves of each temporal line through the class function
    k(U) = exp( (kappa/N) sum_n a_n Re Tr U^n ),
and the Osterwalder–Seiler mechanism needs k to be POSITIVE-DEFINITE on the group: every coefficient khat(rho) in k = sum khat(rho) chi_rho
non-negative (Peter–Weyl / Bochner).

1. The lemma is false (EXACT).  "c >= 0 makes exp(c Re Tr P^n) positive-definite" holds for n = 1 only.  Tr U^n is the power sum p_n, a VIRTUAL
   character: p_n = sum_{k=0}^{n-1} (-1)^k s_{(n-k, 1^k)} (Murnaghan–Nakayama), e.g. Tr U^2 = chi_Sym2 - chi_Lambda2 (the Adams operation
   psi^2 of anomaly8d).  On SU(3): Tr U^2 = chi_6 - chi_3bar, Tr U^3 = chi_10 - chi_8 + 1.  With w_1 = 0 the tower starts at n = 2 and the
   first-order coefficients are, exactly,
        SU(3) 3bar:     (kappa/6)(a_4 - a_2)   (< 0 for y > log(sqrt 2)/(2 pi) = 0.0552)
        SU(3) adjoint:  -(kappa/3) a_3         (< 0 for every y > 0)
        SU(2) spin 1/2: -(kappa/2) a_3         (< 0 for every y > 0)
   (`first_order_coefficients`).  The cancellation that rescues exp(sum h^n p_n / n) = sum h^k chi_Sym^k (a heavy-boson determinant) needs the
   n = 1 term, which log 1 = 0 removes.
2. The counterexample (NUMERICAL, Weyl integration; EXACT at first order).  Site reflection Theta(x0, x) = (-x0, x) — the papers' definition —
   with F a matrix element of the positive-half Polyakov line V_+ in the representation rho: at beta = 0 the measure factorises over spatial
   sites and  <(Theta F)* F> = khat(rho-bar)/(d_rho^2 khat(1)).  At the package defaults (kappa = 1, y = 0.18174) this is negative in the
   3bar and adjoint channels of SU(3) and the spin-1/2 channel of SU(2); by continuity it stays negative for small beta > 0.  The kernel becomes
   positive-definite (over the representations tested) only for kappa between 500 and 700 (SU(3)) or 30 and 60 (SU(2)) at y_c.
3. What is NOT decided.  For LINK reflections (plane between time slices — the ones that build the transfer matrix) the two crossing temporal
   links enter P as free Haar variables at beta = 0 and integrate the Polyakov coupling away (the kernel is the constant khat(1)); at beta > 0 the
   status is open.  Positivity restricted to gauge-invariant functions is also open.  So a positive transfer matrix is not ruled out — it is
   unproven.  Paper 24's own Step 3 argues for n x 1 rectangles, a different action; rectangles crossing the plane off-centre are not of the
   form A Theta(A)^dagger and are not covered by the argument as written either.
4. Two structural points recorded with the correction.  The Polyakov term singles out the time direction and needs a compact time circle, so
   the action is not Euclidean (O(4) / hypercubic) invariant — OS1 can hold in the continuum only if the term becomes irrelevant.  And
   mu_N(y) = min_m sum n^2 a_n (1 - cos 2 pi n m / N) > 0 is the curvature of the classical Polyakov potential at the centre-symmetric point (a
   sum of non-negative terms): a stiffness, not a spectral gap of any transfer matrix.
5. The repair (EXACT).  Attach the weights to genuine characters: k_sym(U) = exp((kappa/N) sum_n a_n Re chi_Sym^n(U)) is positive-definite
   for every kappa >= 0 (sums, products and exponentials of positive-type class functions are positive-type).  Placing the tower on plaquettes,
   sum_n a_n Re chi_Sym^n(U_plaquette), restores hypercubic symmetry and falls under Osterwalder–Seiler; the eigenphase potential changes, so
   mu_N and the spectral lock at y_c must be recomputed for the repaired action.  Universality then predicts the same continuum Yang–Mills:
   the arithmetic tower can shape the lattice road but cannot by itself supply the continuum gap.

Status: items 1 and 5 EXACT; item 2 EXACT at first order in kappa and NUMERICAL (Weyl integration, representations with l1 + l2 <= 10 for
SU(3), spins <= 8 for SU(2)) beyond it; item 3 OPEN; item 4 structural."""
import math
import itertools
import numpy as np
import sympy as sp

from mtft.arithmetic import weight_array

Y_DEFAULT = 0.18174
KAPPA_DEFAULT = 1.0


def a_n(y=Y_DEFAULT, n_max=40):
    w = weight_array(n_max)
    return [float(w[n - 1]) * math.exp(-2 * math.pi * y * n) for n in range(1, n_max + 1)]


# ------------------------------------------------------------------ exact: power sums as virtual characters
def hook_expansion(n):
    """Murnaghan–Nakayama: p_n = sum_{k=0}^{n-1} (-1)^k s_{(n-k, 1^k)} as a list of (partition, sign)."""
    return [(tuple([n - k] + [1] * k), (-1) ** k) for k in range(n)]

def reduce_to_su(partition, N):
    """The SU(N) irrep of a partition: None if more than N rows; otherwise subtract full columns and return (l1, ..., l_{N-1})."""
    lam = list(partition) + [0] * (N - len(partition))
    if len(partition) > N: return None
    m = lam[N - 1]
    return tuple(l - m for l in lam[:N - 1])

def conjugate_rep(lam, N):
    """Highest weight of the dual: (l1, ..., l_{N-1}, 0) -> (l1 - l_{N-1}, ..., l1 - l2, l1) truncated, i.e. conj_i = l1 - l_{N+1-i}."""
    full = list(lam) + [0]
    return tuple(full[0] - full[N - i] for i in range(1, N))

def dim_su(lam, N):
    full = list(lam) + [0]
    num = 1; den = 1
    for i in range(N):
        for j in range(i + 1, N):
            num *= full[i] - full[j] + j - i; den *= j - i
    return num // den

def power_sum_characters(n, N):
    """Tr U^n on SU(N) as {irrep: integer coefficient} (EXACT)."""
    out = {}
    for part, sgn in hook_expansion(n):
        lam = reduce_to_su(part, N)
        if lam is not None: out[lam] = out.get(lam, 0) + sgn
    return {k: v for k, v in out.items() if v}

def first_order_coefficients(N, y=Y_DEFAULT, n_max=40, kappa=1):
    """Leading-order (in kappa) character coefficients of k(U) = exp((kappa/N) sum_n a_n Re Tr U^n): (kappa/N) sum_n a_n (c_{n,rho} + c_{n,rho-bar})/2.
    EXACT given the a_n."""
    a = a_n(y, n_max); out = {}
    for n in range(1, n_max + 1):
        for lam, c in power_sum_characters(n, N).items():
            for rep in (lam, conjugate_rep(lam, N)):
                out[rep] = out.get(rep, 0.0) + kappa / N * a[n - 1] * c / 2
    return {k: v for k, v in out.items() if abs(v) > 0}


# ------------------------------------------------------------------ numerical: Weyl integration on SU(2), SU(3)
class _Torus:
    def __init__(self, N, M):
        self.N = N
        if N == 2:
            th = 2 * np.pi * (np.arange(M) + 0.5) / M; z = np.exp(1j * th); self.Z = [z, np.conj(z)]
            mu = np.sin(th) ** 2
        elif N == 3:
            th = 2 * np.pi * (np.arange(M) + 0.5) / M; T1, T2 = np.meshgrid(th, th, indexing="ij")
            z1, z2 = np.exp(1j * T1), np.exp(1j * T2); self.Z = [z1, z2, np.conj(z1 * z2)]
            mu = np.abs(z1 - z2) ** 2 * np.abs(z1 - self.Z[2]) ** 2 * np.abs(z2 - self.Z[2]) ** 2
        else: raise ValueError("N = 2 or 3")
        self.mu = mu / np.sum(mu)
        e = [np.ones_like(self.Z[0])] + [sum(np.prod(c, axis=0) for c in itertools.combinations(self.Z, k)) for k in range(1, N + 1)]
        self.h = [np.ones_like(self.Z[0])]
        for k in range(1, 40):                                     # complete homogeneous (Sym^k characters) by Newton–Girard in e_k
            self.h.append(sum((-1) ** (i + 1) * e[i] * self.h[k - i] for i in range(1, min(k, N) + 1)))
    def p(self, n): return sum(z ** n for z in self.Z)
    def H(self, k): return self.h[k] if k >= 0 else np.zeros_like(self.Z[0])
    def schur(self, lam):
        """Jacobi–Trudi s_lam = det(h_{lam_i - i + j})."""
        lam = list(lam) + [0] * (self.N - 1 - len(lam)); n = len(lam)
        if n == 1: return self.H(lam[0])
        return self.H(lam[0]) * self.H(lam[1]) - self.H(lam[0] + 1) * self.H(lam[1] - 1)

def default_reps(N):
    return [(l1, l2) for l1 in range(0, 11) for l2 in range(0, l1 + 1) if l1 + l2 <= 10] if N == 3 else [(t,) for t in range(0, 17)]

def kernel_coefficients(N, kappa=KAPPA_DEFAULT, y=Y_DEFAULT, weights="power", reps=None, M=None, n_max=30):
    """Normalised character coefficients khat(rho)/khat(1) of the Polyakov kernel.  weights = "power" (the MTFT action, Re Tr U^n) or "sym"
    (the repair, Re chi_Sym^n).  NUMERICAL (midpoint Weyl integration; the integrands are smooth and periodic)."""
    M = M or (180 if N == 3 else 2000); T = _Torus(N, M); a = a_n(y, n_max)
    f = sum(a[n - 1] * (T.p(n).real if weights == "power" else T.H(n).real) for n in range(1, n_max + 1))
    k = np.exp(kappa / N * f); norm = float(np.sum(k * T.mu))
    reps = reps or default_reps(N)
    return {lam: float(np.real(np.sum(k * np.conj(T.schur(lam)) * T.mu))) / norm for lam in reps}

def site_reflection_pairing(N, rep, kappa=KAPPA_DEFAULT, y=Y_DEFAULT):
    """<(Theta F)* F> / Z at beta = 0 for F = a diagonal matrix element of the positive-half Polyakov line in rep: khat(rep-bar)/(d_rep^2 khat(1))."""
    c = kernel_coefficients(N, kappa, y, reps=[conjugate_rep(rep, N)])
    return c[conjugate_rep(rep, N)] / dim_su(rep, N) ** 2

def min_coefficient(N, kappa, y=Y_DEFAULT, weights="power"):
    c = kernel_coefficients(N, kappa, y, weights=weights); lam = min(c, key=c.get)
    return c[lam], lam

def positivity_bracket(N, y=Y_DEFAULT, lo=None, hi=None, iters=12):
    """Bisection for the kappa at which the power-weight kernel's lowest tested coefficient turns non-negative (NUMERICAL, truncated rep set)."""
    lo = lo if lo is not None else (30.0 if N == 2 else 500.0); hi = hi if hi is not None else (60.0 if N == 2 else 700.0)
    assert min_coefficient(N, lo, y)[0] < 0 <= min_coefficient(N, hi, y)[0]
    for _ in range(iters):
        mid = math.sqrt(lo * hi)
        if min_coefficient(N, mid, y)[0] < 0: lo = mid
        else: hi = mid
    return lo, hi


def link_reflection_note():
    return ("Link reflection (plane between t = 0 and t = 1): P = C1 V+ C2 V- with C1, C2 the crossing links; at beta = 0, "
            "int k(C1 X) dC1 = khat(1) for every X (Haar invariance), so the Polyakov coupling across the plane is erased.  Status at beta > 0: OPEN.")

def cc35_report():
    fo3 = first_order_coefficients(3); fo2 = first_order_coefficients(2)
    c3 = kernel_coefficients(3, reps=[(1, 1), (2, 1), (1, 0)]); c2 = kernel_coefficients(2, reps=[(1,)])
    return {"lemma_false": {"TrU2_SU3": power_sum_characters(2, 3), "TrU3_SU3": power_sum_characters(3, 3), "TrU3_SU2": power_sum_characters(3, 2)},
            "first_order_per_kappa": {"SU3_3bar": fo3.get((1, 1)), "SU3_adjoint": fo3.get((2, 1)), "SU2_spin_half": fo2.get((1,))},
            "defaults": {"kappa": KAPPA_DEFAULT, "y": Y_DEFAULT, "SU3_3bar": c3[(1, 1)], "SU3_adjoint": c3[(2, 1)], "SU2_spin_half": c2[(1,)]},
            "site_reflection_pairing_SU3_adjoint": site_reflection_pairing(3, (2, 1)),
            "link_reflection": link_reflection_note(),
            "repair_min_coefficient_kappa_1_10_100": [min_coefficient(3, kp, weights="sym")[0] for kp in (1, 10, 100)]}
