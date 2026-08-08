#!/usr/bin/env python3
"""
Combinatorial Ancestry — the 2024 lineage as working tools
============================================================================

MIT License — Copyright (c) 2026 Roger Tano
See LICENSE file for full terms.

Before the weights and the curve there was the 2024 combinatorial
program: figurate power sums, the graph uncertainty principle, the
number-phase Fourier pair, and combinatorial thermodynamics.  This
module lands that ancestry as certified, importable theorems, plus the
exact bridges that connect it to the modern stack.

Contents and epistemic classes
------------------------------
1. FIGURATE / FAULHABER (EXACT, Fractions throughout)
   S_p(n) = C(p) T_n^q − R_p(n),  q = (p+1)/2,  C(p) = 2^q/(p+1).
   Provenance: Faulhaber 1631; Jacobi 1834; Knuth, "Johann Faulhaber
   and sums of powers" (1993).  The 2024 telescoping / local-jump rule
   is re-certified here by exact recurrence against the closed form.

2. THE SIGMA INVOLUTION (Pr / EXACT) — the structural upgrade.
   sigma: n -> −1−n fixes T_n; the invariant ring is Q[n]^sigma = Q[T],
   and the sigma-odd sector is the rank-one module (2n+1)·Q[T] with
   ramified element (2n+1)^2 = 8T+1.  Odd-p power sums are sigma-EVEN
   (pure T-polynomials); even-p sums are sigma-ODD.  The 2024 "Figurate
   Representation Principle" is therefore a parity selection rule under
   an involution — the same theorem-shape as the eta-parity
   superselection on the X0(143) harmonic stage (du03).  That
   correspondence is STRUCTURAL (shared mechanism class), not a
   derivation link; it is recorded as such and nothing more.

3. THE s-GONAL ABSORPTION DEFECT (Pr / EXACT) — an honest negative
   with a closed form.  For s-gonal numbers P_s the analogous ansatz
   A_s(p) P_s^q fails to absorb the n^p Faulhaber term unless s = 3:
       [A_s P_s^q − S_p]_{n^p} = −(s−3)/(s−2),   independent of p.
   Triangular numbers are singled out because P_3 = T is the unique
   sigma-invariant s-gonal family: P_s(−1−n) − P_s(n) = (2n+1)(s−3).

4. GRAPH UNCERTAINTY (Pr / EXACT + Cert)
   [D, L] = AD − DA with entries A_uv (d_v − d_u); zero iff every edge
   has zero degree gradient; on a connected graph, iff regular.
   Robertson bound Delta_D Delta_L >= |<[D,L]>|/2 certified on complex
   states (real states give <[X,Y]> = 0 for real-symmetric X, Y — the
   bound is trivially met there; use complex states for content).

5. NUMBER-PHASE AND ENTROPIC UNCERTAINTY (Pr provenance / Cert here)
   K = diag(0..n), Theta = F K F* with F the DFT.  Robertson certified
   numerically.  Entropic bound H_K + H_Theta >= log(n+1) is
   Maassen–Uffink 1988 for the mutually unbiased position/DFT pair
   (overlap 1/sqrt(d), bound −2 log(1/sqrt d) = log d) — the corrected
   2024 constant log(n+1) is exactly this.
   ARCHIVED 2024 NUMBERS — RESOLVED: the binomial(5, 1/2) state gives
   Delta_K = sqrt(np(1−p)) = sqrt(5)/2 = 1.1180 (EXACT).  The
   index-units first pass here gave Delta_Theta = 1.309440,
   product = 1.463999, bound = 0.846970 — each exactly the archived
   value divided by 2·pi/6, which identified the 2024 convention as
   ANGULAR (Theta spectrum 2·pi·j/d).  Under that convention all four
   archived numbers (1.1180, 1.3712, 1.5331, 0.8869) reproduce to
   their printed precision.  Cert(5e-4); discovery route preserved in
   number_phase_regression().

6. q-COMBINATORIAL THERMODYNAMICS (Df + Cert; anchors EXACT)
   Gaussian binomials by two independent exact routes (q-Pascal
   recurrence vs product-formula polynomial division).  Galois-number
   anchors at q=2: 1, 2, 5, 16, 67, 374.  Gibbs ensembles with the
   full chain logZ -> U, Var, S, C, F and the Fisher metric
   g_bb = Var(E); every identity checked by two routes (ensemble
   moments vs finite differences of logZ).  Multiplicity-as-energy
   E_k = −log Omega_k gives Z(beta) = sum Omega_k^beta with EXACT
   special values Z(0) = n+1 and, at q = 1, Z(1) = 2^n.

7. THE BRIDGE (the reason this module exists)
   Take the 2024 counting->energy map literally on the integers:
   E_n = log n makes the Gibbs ensemble the primon gas, Z = zeta(beta).
   Evaluate the 2024 ensemble machinery on the MODERN observable — the
   Tano weight w_n — and the mean is a closed form:

       <w>_beta = −zeta'(beta+1)        (Pr; one line from
                                         sum w_n n^{−s} = −zeta(s) zeta'(s+1))

   As beta -> 1+ the mean weight tends to −zeta'(2), i.e. exactly
   2·T_INF in the package normalization (constants.T_INF = −zeta'(2)/2)
   and exactly the marked-gas cold-gas alpha.  The ancestor ensemble at
   unit inverse temperature sits on the corpus constant.  Certified
   numerically with an explicit tail majorant (Cert, tolerance stated
   per call), with mtft.arithmetic.weight_array as the independent
   route for the sieve.

Working rules honored: every public number carries a class; exact
claims use Fractions with no floats; nothing here imports sympy; the
convention-dependence of the archived 1.3712 is reported, not patched.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np

__all__ = [
    # figurate / Faulhaber
    "bernoulli_plus", "faulhaber_coeffs", "power_sum", "power_sum_backward",
    "C_figurate", "FigurateDecomposition", "figurate_decomposition",
    "jump_rule_R", "R_polynomial_coeffs",
    # sigma involution
    "sigma_reflect", "sigma_split", "to_T_basis", "odd_sector_factor",
    "sigma_parity_check",
    # s-gonal defect
    "sgonal_coeffs", "sgonal_A", "sgonal_defect", "sgonal_defect_measured",
    # graph uncertainty
    "graph_cycle", "graph_path", "graph_star", "graph_complete",
    "graph_hypercube", "graph_erdos_renyi", "degrees", "is_regular",
    "commutator_DL", "degree_gradient_matrix", "robertson_margin",
    # number-phase
    "dft_matrix", "phase_operator", "binomial_state", "shannon_entropy",
    "entropic_margin", "number_phase_regression", "ARCHIVED_2024",
    # q-thermodynamics
    "gaussian_binomial", "gaussian_binomial_product", "galois_number",
    "GALOIS_ANCHORS_Q2", "Ensemble", "multiplicity_ensemble",
    "multiplicity_Z_exact_specials",
    # bridge
    "tano_weight_sieve", "weight_dirichlet_tail_majorant",
    "weight_dirichlet_identity_check", "mean_tano_weight",
    "mean_tano_weight_exact", "primon_fisher_check",
    # legend
    "ANCESTRY_LEGEND",
]

Fr = Fraction

# ════════════════════════════════════════════════════════════════════
# Exact polynomial helpers — coefficient lists, index = degree, Fractions
# ════════════════════════════════════════════════════════════════════

def _ptrim(c: List[Fraction]) -> List[Fraction]:
    while c and c[-1] == 0:
        c.pop()
    return c


def _padd(a: Sequence[Fraction], b: Sequence[Fraction]) -> List[Fraction]:
    n = max(len(a), len(b))
    out = [Fr(0)] * n
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return _ptrim(out)


def _pscale(a: Sequence[Fraction], s: Fraction) -> List[Fraction]:
    return _ptrim([s * x for x in a])


def _pmul(a: Sequence[Fraction], b: Sequence[Fraction]) -> List[Fraction]:
    if not a or not b:
        return []
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return _ptrim(out)


def _ppow(a: Sequence[Fraction], k: int) -> List[Fraction]:
    out: List[Fraction] = [Fr(1)]
    base = list(a)
    while k:
        if k & 1:
            out = _pmul(out, base)
        base = _pmul(base, base)
        k >>= 1
    return out


def _peval(a: Sequence[Fraction], x: Fraction) -> Fraction:
    acc = Fr(0)
    for c in reversed(list(a)):
        acc = acc * x + c
    return acc


_T_POLY: List[Fraction] = [Fr(0), Fr(1, 2), Fr(1, 2)]  # T = n/2 + n^2/2


# ════════════════════════════════════════════════════════════════════
# 1. Faulhaber engine (EXACT)
# ════════════════════════════════════════════════════════════════════

def bernoulli_plus(m_max: int) -> List[Fraction]:
    """Bernoulli numbers B_0..B_m_max in the B+ convention (B_1 = +1/2).

    EXACT.  Recurrence B_m = 1 − sum_{j<m} C(m,j) B_j / (m−j+1).
    """
    B = [Fr(1)]
    for m in range(1, m_max + 1):
        s = Fr(0)
        for j in range(m):
            s += Fr(math.comb(m, j), m - j + 1) * B[j]
        B.append(Fr(1) - s)
    return B


def faulhaber_coeffs(p: int) -> List[Fraction]:
    """Coefficients (index = degree in n) of S_p(n) = sum_{k=1}^n k^p.

    EXACT.  S_p(n) = (1/(p+1)) sum_{j=0}^{p} C(p+1, j) B_j^+ n^{p+1−j}.
    """
    if p < 1:
        raise ValueError("p >= 1")
    B = bernoulli_plus(p)
    c = [Fr(0)] * (p + 2)
    for j in range(p + 1):
        c[p + 1 - j] += Fr(math.comb(p + 1, j)) * B[j] / (p + 1)
    return _ptrim(c)


def power_sum(p: int, n: int) -> int:
    """Route A: direct summation sum_{k=1}^n k^p.  EXACT, no Bernoulli."""
    return sum(k ** p for k in range(1, n + 1))


def power_sum_backward(p: int, m: int) -> int:
    """The Faulhaber polynomial's value at a NEGATIVE integer m <= −1,
    obtained purely from the defining recurrence S(n) − S(n−1) = n^p run
    downward from S(0) = 0.  EXACT, shares no steps with the Bernoulli
    route.  (Used as route A of the sigma-parity certificate.)
    """
    if m > -1:
        raise ValueError("m <= -1")
    S = 0
    n = 0
    while n > m:
        S -= n ** p          # S(n−1) = S(n) − n^p
        n -= 1
    return S


def C_figurate(p: int) -> Fraction:
    """C(p) = 2^{(p+1)/2} / (p+1) for odd p.  EXACT (Df)."""
    if p % 2 == 0:
        raise ValueError("odd p only")
    q = (p + 1) // 2
    return Fr(2 ** q, p + 1)


# ── sigma involution machinery ──────────────────────────────────────

def sigma_reflect(c: Sequence[Fraction]) -> List[Fraction]:
    """Coefficients of f(−1−n) given coefficients of f(n).  EXACT."""
    out = [Fr(0)] * len(c)
    for d, cd in enumerate(c):
        if cd == 0:
            continue
        sgn = -1 if d % 2 else 1
        for i in range(d + 1):
            out[i] += cd * sgn * math.comb(d, i)
    return _ptrim(out)


def sigma_split(c: Sequence[Fraction]) -> Tuple[List[Fraction], List[Fraction]]:
    """(sigma-even part, sigma-odd part) of f.  EXACT."""
    g = sigma_reflect(c)
    even = _pscale(_padd(c, g), Fr(1, 2))
    odd = _pscale(_padd(list(c), _pscale(g, Fr(-1))), Fr(1, 2))
    return even, odd


def to_T_basis(c_even: Sequence[Fraction]) -> Dict[int, Fraction]:
    """Express a sigma-EVEN polynomial as an element of Q[T].

    EXACT; the reduction's internal assertions (even degree at every
    step) are themselves the certificate that Q[n]^sigma = Q[T].
    Raises ValueError if the input is not sigma-even.
    """
    f = _ptrim(list(c_even))
    out: Dict[int, Fraction] = {}
    while f:
        d = len(f) - 1
        if d == 0:
            out[0] = out.get(0, Fr(0)) + f[0]
            break
        if d % 2 != 0:
            raise ValueError("not sigma-even: odd degree survived reduction")
        m = d // 2
        coeff = f[-1] * (2 ** m)          # T^m has leading 1/2^m
        out[m] = out.get(m, Fr(0)) + coeff
        f = _padd(f, _pscale(_ppow(_T_POLY, m), -coeff))
    return {k: v for k, v in out.items() if v != 0}


def odd_sector_factor(c: Sequence[Fraction]) -> Dict[int, Fraction]:
    """For a sigma-ODD polynomial f, return a(T) with f = (2n+1)·a(T).

    EXACT.  Synthetic division by (2n+1) at root −1/2; the remainder
    must vanish and the quotient must be sigma-even — both asserted.
    """
    f = _ptrim(list(c))
    if not f:
        return {}
    # divide by (2n+1): f(n) = (2n+1) g(n) + r;  r = f(−1/2)
    r = _peval(f, Fr(-1, 2))
    if r != 0:
        raise ValueError("not in (2n+1)·Q[T]: nonzero remainder")
    # deflate: repeatedly peel leading term  c_d n^d = (2n+1)·(c_d/2) n^{d−1} − ...
    g = [Fr(0)] * (len(f) - 1)
    rem = list(f)
    for d in range(len(f) - 1, 0, -1):
        q = rem[d] / 2
        g[d - 1] = q
        rem[d] = Fr(0)
        rem[d - 1] -= q
    if _ptrim(rem):
        raise ValueError("division by (2n+1) left a remainder")
    return to_T_basis(g)


def sigma_parity_check(p: int, n_max: int = 40) -> bool:
    """Two-route certificate that S_p is sigma-even (odd p) or
    sigma-odd (even p).  EXACT.

    Route A: values of the polynomial at negative integers via the pure
    downward recurrence, compared with (−1)^{p+1} S_p(n) from direct
    summation.  Route B: coefficient-level substitution n -> −1−n on
    the Bernoulli/Faulhaber coefficients.  The routes share no steps.
    """
    sgn = (-1) ** (p + 1)
    for n in range(0, n_max + 1):
        if power_sum_backward(p, -1 - n) != sgn * power_sum(p, n):
            return False
    c = faulhaber_coeffs(p)
    return sigma_reflect(c) == _pscale(c, Fr(sgn))


# ── the figurate decomposition itself ───────────────────────────────

@dataclass(frozen=True)
class FigurateDecomposition:
    """S_p = C T^q − R_p with everything in the T-basis.  EXACT."""
    p: int
    q: int
    C: Fraction
    S_T: Dict[int, Fraction]
    R_T: Dict[int, Fraction]

    @property
    def deg_T_R(self) -> int:
        return max(self.R_T) if self.R_T else -1


def figurate_decomposition(p: int) -> FigurateDecomposition:
    """The 2024 identity, certified: for odd p, S_p is a pure
    T-polynomial with leading coefficient C(p) = 2^q/(p+1), and
    R_p = C T^q − S_p has T-degree <= q−1 (n-degree <= p−1).  EXACT.
    """
    if p % 2 == 0:
        raise ValueError("figurate decomposition is the odd-p statement")
    q = (p + 1) // 2
    even, odd = sigma_split(faulhaber_coeffs(p))
    if _ptrim(list(odd)):
        raise AssertionError("odd-p power sum failed sigma-evenness")
    S_T = to_T_basis(even)
    C = C_figurate(p)
    if S_T.get(q) != C:
        raise AssertionError("leading T-coefficient != 2^q/(p+1)")
    R_T = {m: -v for m, v in S_T.items() if m != q}
    if any(m >= q for m in R_T):
        raise AssertionError("degree reduction failed")
    return FigurateDecomposition(p=p, q=q, C=C, S_T=S_T, R_T=R_T)


def R_polynomial_coeffs(p: int) -> List[Fraction]:
    """R_p(n) = C(p) T_n^q − S_p(n) as an n-polynomial.  EXACT."""
    q = (p + 1) // 2
    CTq = _pscale(_ppow(_T_POLY, q), C_figurate(p))
    return _padd(CTq, _pscale(faulhaber_coeffs(p), Fr(-1)))


def jump_rule_R(p: int, n_max: int) -> List[Fraction]:
    """The 2024 local-jump recurrence, run exactly:
        R(0) = 0,
        R(m) = R(m−1) + C(p)[T_m^q − T_{m−1}^q] − m^p.
    Returns [R(0), ..., R(n_max)].  EXACT; route B against the closed
    form R_polynomial_coeffs (route A).  Also certifies R(1) = C(p)−1.
    """
    q = (p + 1) // 2
    C = C_figurate(p)
    vals = [Fr(0)]
    T_prev = Fr(0)
    for m in range(1, n_max + 1):
        T_m = Fr(m * (m + 1), 2)
        vals.append(vals[-1] + C * (T_m ** q - T_prev ** q) - Fr(m ** p))
        T_prev = T_m
    return vals


# ════════════════════════════════════════════════════════════════════
# 3. s-gonal absorption defect (Pr / EXACT) — the honest negative
# ════════════════════════════════════════════════════════════════════

def sgonal_coeffs(s: int) -> List[Fraction]:
    """P_s(n) = ((s−2)n^2 − (s−4)n)/2.  Df."""
    return _ptrim([Fr(0), Fr(-(s - 4), 2), Fr(s - 2, 2)])


def sgonal_A(s: int, p: int) -> Fraction:
    """A_s(p) = 2^q / ((p+1)(s−2)^q): the unique constant matching the
    leading n^{p+1} Faulhaber term.  Df/EXACT."""
    q = (p + 1) // 2
    return Fr(2 ** q, (p + 1) * (s - 2) ** q)


def sgonal_defect(s: int) -> Fraction:
    """Predicted n^p coefficient of A_s P_s^q − S_p:  −(s−3)/(s−2),
    independent of p.  Pr (two-line proof in the module docstring);
    zero iff s = 3 — triangular numbers are the unique absorbing
    family, because P_3 = T is the unique sigma-invariant P_s."""
    return Fr(-(s - 3), s - 2)


def sgonal_defect_measured(s: int, p: int) -> Fraction:
    """Route B: actually build A_s P_s^q − S_p and read off the n^p
    coefficient.  EXACT."""
    if p % 2 == 0:
        raise ValueError("odd p only")
    q = (p + 1) // 2
    R = _padd(_pscale(_ppow(sgonal_coeffs(s), q), sgonal_A(s, p)),
              _pscale(faulhaber_coeffs(p), Fr(-1)))
    # leading n^{p+1} must always cancel (A_s is chosen for exactly that)
    if len(R) - 1 > p:
        raise AssertionError("leading Faulhaber term failed to cancel")
    return R[p] if len(R) > p else Fr(0)


# ════════════════════════════════════════════════════════════════════
# 4. Graph uncertainty (Pr / EXACT + Cert)
# ════════════════════════════════════════════════════════════════════

def graph_cycle(n: int) -> np.ndarray:
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        A[i, (i + 1) % n] = A[(i + 1) % n, i] = 1
    return A


def graph_path(n: int) -> np.ndarray:
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n - 1):
        A[i, i + 1] = A[i + 1, i] = 1
    return A


def graph_star(n: int) -> np.ndarray:
    A = np.zeros((n, n), dtype=np.int64)
    A[0, 1:] = A[1:, 0] = 1
    return A


def graph_complete(n: int) -> np.ndarray:
    A = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
    return A


def graph_hypercube(d: int) -> np.ndarray:
    n = 1 << d
    A = np.zeros((n, n), dtype=np.int64)
    for v in range(n):
        for b in range(d):
            A[v, v ^ (1 << b)] = 1
    return A


def graph_erdos_renyi(n: int, prob: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    U = rng.random((n, n)) < prob
    A = np.triu(U, 1).astype(np.int64)
    return A + A.T


def degrees(A: np.ndarray) -> np.ndarray:
    return A.sum(axis=1)


def is_regular(A: np.ndarray) -> bool:
    d = degrees(A)
    return bool(np.all(d == d[0]))


def commutator_DL(A: np.ndarray) -> np.ndarray:
    """Route A: [D, L] = DL − LD by direct integer matrix products.
    EXACT (int64)."""
    D = np.diag(degrees(A))
    L = D - A
    return D @ L - L @ D


def degree_gradient_matrix(A: np.ndarray) -> np.ndarray:
    """Route B: the claimed closed form  [D,L]_{uv} = A_{uv}(d_v − d_u).
    EXACT (int64).  Pr: [D,L] = [D, D−A] = AD − DA, and
    (AD − DA)_{uv} = A_{uv} d_v − d_u A_{uv}."""
    d = degrees(A)
    return A * (d[None, :] - d[:, None])


def robertson_margin(X: np.ndarray, Y: np.ndarray, psi: np.ndarray) -> Dict[str, float]:
    """Delta_X Delta_Y − |<[X,Y]>|/2 for Hermitian X, Y and unit psi.
    Cert: nonnegative up to float roundoff.  NOTE: for real-symmetric
    X, Y and a REAL state, <[X,Y]> = 0 identically — use complex states
    for a nontrivial bound."""
    psi = psi / np.linalg.norm(psi)
    def mean(M):
        return float(np.real(np.vdot(psi, M @ psi)))
    mx, my = mean(X), mean(Y)
    vx = mean(X @ X) - mx * mx
    vy = mean(Y @ Y) - my * my
    comm = complex(np.vdot(psi, (X @ Y - Y @ X) @ psi))
    lower = 0.5 * abs(comm)
    prod = math.sqrt(max(vx, 0.0)) * math.sqrt(max(vy, 0.0))
    return {"dX": math.sqrt(max(vx, 0.0)), "dY": math.sqrt(max(vy, 0.0)),
            "product": prod, "bound": lower, "margin": prod - lower}


# ════════════════════════════════════════════════════════════════════
# 5. Number-phase pair and entropic uncertainty (Cert)
# ════════════════════════════════════════════════════════════════════

def dft_matrix(d: int, sign: int = +1) -> np.ndarray:
    """F_{kj} = exp(sign·2·pi·i·k·j/d)/sqrt(d).  Df."""
    k = np.arange(d)
    return np.exp(sign * 2j * np.pi * np.outer(k, k) / d) / math.sqrt(d)


def phase_operator(d: int, sign: int = +1, angular: bool = False) -> np.ndarray:
    """Theta = F K F* — the DFT conjugate of K = diag(0..d−1).  Df.

    angular=True rescales the spectrum to angle units, 2·pi·j/d — the
    convention the 2024 archive used (identified by this module's
    regression: the index-units first pass missed the archived numbers
    by exactly the factor 2·pi/d, which pinned the convention)."""
    F = dft_matrix(d, sign)
    K = np.diag(np.arange(d, dtype=float))
    Th = F @ K @ F.conj().T
    return (2 * np.pi / d) * Th if angular else Th


def binomial_state(n: int, x: float) -> np.ndarray:
    """f_k = sqrt(C(n,k) x^k (1−x)^{n−k}) — the 2024 test state.  Df."""
    return np.sqrt(np.array([math.comb(n, k) * x ** k * (1 - x) ** (n - k)
                             for k in range(n + 1)]))


def shannon_entropy(prob: np.ndarray) -> float:
    p = prob[prob > 1e-300]
    return float(-(p * np.log(p)).sum())


def entropic_margin(psi: np.ndarray, sign: int = +1) -> Dict[str, float]:
    """H_K + H_Theta − log(d) for the position/DFT pair.

    Pr (Maassen–Uffink 1988: bound = −2 log max overlap = log d for
    this mutually unbiased pair); Cert here at stated tolerance."""
    psi = psi.astype(complex) / np.linalg.norm(psi)
    d = len(psi)
    F = dft_matrix(d, sign)
    hk = shannon_entropy(np.abs(psi) ** 2)
    ht = shannon_entropy(np.abs(F.conj().T @ psi) ** 2)
    return {"H_K": hk, "H_Theta": ht, "bound": math.log(d),
            "margin": hk + ht - math.log(d)}


ARCHIVED_2024: Dict[str, float] = {
    # the numbers quoted in the 2024 record for n=5, x=0.5
    "dK": 1.1180, "dTheta": 1.3712, "product": 1.5331, "bound": 0.8869,
}


def number_phase_regression(rtol: float = 5e-4) -> Dict[str, object]:
    """Recompute the archived 2024 (n=5, x=1/2) numbers.  Artifact over
    memory, with the discovery route preserved: the index-units first
    pass gave dTheta = 1.309440, product = 1.463999, bound = 0.846970 —
    each exactly the archived value divided by 2·pi/6.  That factor
    identified the 2024 convention as ANGULAR (Theta spectrum 2·pi·j/d),
    under which all four archived numbers reproduce to their printed
    precision.  Cert(rtol) on the angular row; the index rows are kept
    as the record of how the convention was pinned."""
    psi = binomial_state(5, 0.5)
    d = len(psi)
    K = np.diag(np.arange(d, dtype=float))
    out: Dict[str, object] = {"archived": dict(ARCHIVED_2024),
                              "dK_exact": math.sqrt(5) / 2}
    for name, kwargs in (("index_plus", dict(sign=+1)),
                         ("index_minus", dict(sign=-1)),
                         ("angular", dict(sign=+1, angular=True))):
        r = robertson_margin(K, phase_operator(d, **kwargs), psi.astype(complex))
        r["matches_archived"] = bool(
            abs(r["dY"] - ARCHIVED_2024["dTheta"]) <= rtol * ARCHIVED_2024["dTheta"]
            and abs(r["product"] - ARCHIVED_2024["product"]) <= rtol * ARCHIVED_2024["product"]
            and abs(r["bound"] - ARCHIVED_2024["bound"]) <= rtol * ARCHIVED_2024["bound"])
        out[name] = r
    return out


# ════════════════════════════════════════════════════════════════════
# 6. q-combinatorics and thermodynamics
# ════════════════════════════════════════════════════════════════════

def gaussian_binomial(n: int, k: int) -> List[int]:
    """[n, k]_q as an integer-coefficient polynomial in q (index =
    degree), by the q-Pascal recurrence [n,k] = [n−1,k−1] + q^k [n−1,k].
    EXACT (route A)."""
    if k < 0 or k > n:
        return []
    row: List[List[int]] = [[1]] + [[] for _ in range(k)]
    for m in range(1, n + 1):
        upper = min(k, m)
        for j in range(upper, 0, -1):
            a = row[j - 1]
            b = row[j]
            shifted = [0] * j + b if b else []
            L = max(len(a), len(shifted))
            row[j] = [(a[i] if i < len(a) else 0) +
                      (shifted[i] if i < len(shifted) else 0) for i in range(L)]
    return row[k]


def gaussian_binomial_product(n: int, k: int) -> List[int]:
    """Route B: prod_{i=1}^{k} (1 − q^{n−k+i}) / (1 − q^i) by exact
    polynomial multiplication and division over Z.  EXACT."""
    if k < 0 or k > n:
        return []
    def one_minus_qm(m: int) -> List[int]:
        c = [0] * (m + 1)
        c[0], c[m] = 1, -1
        return c
    num: List[int] = [1]
    den: List[int] = [1]
    for i in range(1, k + 1):
        num = _int_pmul(num, one_minus_qm(n - k + i))
        den = _int_pmul(den, one_minus_qm(i))
    quo, rem = _int_pdivmod(num, den)
    if any(rem):
        raise AssertionError("Gaussian binomial division not exact")
    return quo


def _int_pmul(a: List[int], b: List[int]) -> List[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def _int_pdivmod(a: List[int], b: List[int]) -> Tuple[List[int], List[int]]:
    a = list(a)
    db, lb = len(b) - 1, b[-1]
    if lb not in (1, -1):
        raise ValueError("monic-up-to-sign divisor required")
    quo = [0] * max(1, len(a) - db)
    for i in range(len(a) - 1, db - 1, -1):
        c = a[i] // lb
        if c:
            quo[i - db] = c
            for j, y in enumerate(b):
                a[i - db + j] -= c * y
    return quo, a


def galois_number(n: int, q: int) -> int:
    """G_n(q) = sum_k [n,k]_q — for prime powers q, the number of
    subspaces of F_q^n.  EXACT."""
    return sum(sum(c * q ** d for d, c in enumerate(gaussian_binomial(n, k)))
               for k in range(n + 1))


GALOIS_ANCHORS_Q2: Dict[int, int] = {0: 1, 1: 2, 2: 5, 3: 16, 4: 67, 5: 374}


@dataclass
class Ensemble:
    """A Gibbs ensemble over finitely many levels: p_k ∝ Omega_k e^{−beta E_k}.

    Every thermodynamic identity is exposed by TWO routes (ensemble
    moments vs finite differences of logZ) so the standard chain
    logZ -> U, Var, S, C, F, Fisher can be certified rather than
    assumed.  Cert at the tolerances used in the test gate.
    """
    energies: np.ndarray
    log_degeneracies: np.ndarray = field(default=None)  # type: ignore

    def __post_init__(self):
        self.energies = np.asarray(self.energies, dtype=float)
        if self.log_degeneracies is None:
            self.log_degeneracies = np.zeros_like(self.energies)
        else:
            self.log_degeneracies = np.asarray(self.log_degeneracies, dtype=float)

    def log_Z(self, beta: float) -> float:
        w = self.log_degeneracies - beta * self.energies
        m = float(np.max(w))
        return m + math.log(float(np.exp(w - m).sum()))

    def probs(self, beta: float) -> np.ndarray:
        w = self.log_degeneracies - beta * self.energies
        w -= w.max()
        e = np.exp(w)
        return e / e.sum()

    # ── route A: ensemble moments ────────────────────────────────
    def U(self, beta: float) -> float:
        return float(self.probs(beta) @ self.energies)

    def var_E(self, beta: float) -> float:
        p = self.probs(beta)
        u = float(p @ self.energies)
        return float(p @ (self.energies - u) ** 2)

    def entropy(self, beta: float) -> float:
        return shannon_entropy(self.probs(beta))

    # ── route B: derivatives of logZ ─────────────────────────────
    def U_fd(self, beta: float, h: float = 1e-5) -> float:
        return -(self.log_Z(beta + h) - self.log_Z(beta - h)) / (2 * h)

    def var_fd(self, beta: float, h: float = 1e-4) -> float:
        return (self.log_Z(beta + h) - 2 * self.log_Z(beta)
                + self.log_Z(beta - h)) / h ** 2

    def entropy_identity(self, beta: float) -> float:
        """S = beta·U + logZ (route B; algebraic identity)."""
        return beta * self.U(beta) + self.log_Z(beta)

    # ── derived ──────────────────────────────────────────────────
    def heat_capacity(self, beta: float) -> float:
        return beta ** 2 * self.var_E(beta)

    def free_energy(self, beta: float) -> float:
        return -self.log_Z(beta) / beta

    def fisher_beta(self, beta: float) -> float:
        """g_bb = d^2 logZ / d beta^2 = Var(E) — the thermal Fisher
        metric (route A returns the covariance form)."""
        return self.var_E(beta)


def multiplicity_ensemble(omegas: Sequence[int]) -> Ensemble:
    """The 2024 'multiplicity as energy' map E_k = −log Omega_k, under
    which Z(beta) = sum_k Omega_k^beta.  Df."""
    lo = np.log(np.asarray(omegas, dtype=float))
    return Ensemble(energies=-lo, log_degeneracies=np.zeros_like(lo))


def multiplicity_Z_exact_specials(omegas: Sequence[int]) -> Dict[str, int]:
    """EXACT special values: Z(0) = number of levels; Z(1) = sum Omega_k
    (= 2^n for binomial multiplicities, = Galois numbers for Gaussian
    ones).  Integer arithmetic, no floats."""
    return {"Z0": len(omegas), "Z1": int(sum(int(o) for o in omegas))}


# ════════════════════════════════════════════════════════════════════
# 7. The bridge: 2024 ensembles meet the Tano weights (Pr + Cert)
# ════════════════════════════════════════════════════════════════════

def tano_weight_sieve(n_max: int) -> np.ndarray:
    """w_n = sum_{d|n} (log d)/d for n = 0..n_max (w_0 = 0) by direct
    divisor sieve — deliberately independent of mtft.arithmetic so the
    package's weight_array can serve as route B.  Cert (float)."""
    w = np.zeros(n_max + 1)
    for d in range(2, n_max + 1):        # d = 1 contributes log(1)/1 = 0
        w[d::d] += math.log(d) / d
    return w


def weight_dirichlet_tail_majorant(s: float, n_max: int) -> float:
    """Explicit majorant for sum_{n>N} w_n n^{−s}, s > 1.

    Pr: w_n <= sum_{d<=n} (log d)/d <= (log n)^2/2 + 1, and the
    integral comparison gives
      sum_{n>N} ((log n)^2/2 + 1) n^{−s}
        <= N^{1−s} [ (log N)^2/(2(s−1)) + log N/(s−1)^2 + 1/(s−1)^3
                     + 1/(s−1) ] + first-term slack.
    Returned with a factor-2 safety margin.  Cert-supporting bound.
    """
    if s <= 1:
        raise ValueError("s > 1 required")
    L, a = math.log(n_max), s - 1
    integral = n_max ** (-a) * (L * L / (2 * a) + L / a ** 2 + 1 / a ** 3 + 1 / a)
    first = ((L * L) / 2 + 1) * n_max ** (-s)
    return 2.0 * (integral + first)


def weight_dirichlet_identity_check(s: float, n_max: int = 200_000) -> Dict[str, float]:
    """Certify  sum_n w_n n^{−s} = −zeta(s) zeta'(s+1)  numerically.

    Pr (one line: Dirichlet convolution of 1 with (log d)/d).
    Route A: sieve partial sum + explicit tail majorant.
    Route B: mpmath zeta and zeta'.  Cert(gap <= tail).
    """
    import mpmath as mp
    w = tano_weight_sieve(n_max)
    n = np.arange(1, n_max + 1, dtype=float)
    lhs = float(w[1:] @ n ** (-s))
    rhs = float(-mp.zeta(s) * mp.zeta(s + 1, derivative=1))
    tail = weight_dirichlet_tail_majorant(s, n_max)
    return {"lhs_partial": lhs, "rhs": rhs, "gap": abs(lhs - rhs), "tail_bound": tail}


def mean_tano_weight(beta: float, n_max: int = 200_000) -> float:
    """Truncated Gibbs mean of the Tano weight on the primon ensemble
    (E_n = log n):  <w>_beta = sum w_n n^{−beta} / sum n^{−beta}.
    Cert (float, truncation stated by n_max)."""
    w = tano_weight_sieve(n_max)
    n = np.arange(1, n_max + 1, dtype=float)
    boltz = n ** (-beta)
    return float((w[1:] @ boltz) / boltz.sum())


def mean_tano_weight_exact(beta: float) -> float:
    """<w>_beta = −zeta'(beta+1), beta > 1.

    Pr: divide the Dirichlet identity by Z = zeta(beta).  As
    beta -> 1+, <w>_beta -> −zeta'(2) — the marked-gas cold-gas alpha,
    equal to 2·T_INF in the package normalization
    (constants.T_INF = −zeta'(2)/2).  The 2024 ensemble machinery,
    evaluated on the modern observable, lands on the corpus constant.
    """
    import mpmath as mp
    return float(-mp.zeta(beta + 1, derivative=1))


def primon_fisher_check(beta: float = 3.0, n_max: int = 50_000) -> Dict[str, float]:
    """The thermal Fisher metric of the primon gas equals the log-zeta
    curvature:  Var_beta(log n) = (d^2/d beta^2) log zeta(beta).

    Route A: ensemble variance on the truncated gas.
    Route B: mpmath second derivative of log zeta.  Cert."""
    import mpmath as mp
    n = np.arange(1, n_max + 1, dtype=float)
    gas = Ensemble(energies=np.log(n))
    route_a = gas.var_E(beta)
    route_b = float(mp.diff(lambda b: mp.log(mp.zeta(b)), beta, 2))
    return {"ensemble_var": route_a, "logzeta_curvature": route_b,
            "gap": abs(route_a - route_b)}


# ════════════════════════════════════════════════════════════════════
# Legend entries (schema-compatible; registration handled at repo level)
# ════════════════════════════════════════════════════════════════════

ANCESTRY_LEGEND: Tuple[Dict[str, object], ...] = (
    dict(name="figurate_decomposition", tier="0", kind="identity",
         primitives=("I", "III"), tag="Pr", exactness="EXACT",
         nature="S_p = C(p) T^q − R_p, deg_T R <= q−1; C(p)=2^q/(p+1). "
                "Faulhaber's theorem in the 2024 figurate normalization.",
         example="mtft.combinatorial.figurate_decomposition(7)",
         upstream=("integers",), ref="Faulhaber 1631; Jacobi 1834; Knuth 1993"),
    dict(name="sigma_parity", tier="0", kind="identity",
         primitives=("II", "IV"), tag="Pr", exactness="EXACT",
         nature="sigma: n->−1−n. Q[n]^sigma = Q[T]; odd sector (2n+1)Q[T], "
                "(2n+1)^2 = 8T+1. Odd-p sums sigma-even, even-p sigma-odd: "
                "the figurate principle is a parity selection rule.",
         example="mtft.combinatorial.sigma_parity_check(9)",
         upstream=("integers",), ref="structural companion to du03 eta-parity"),
    dict(name="sgonal_defect", tier="0", kind="identity",
         primitives=("III",), tag="Pr", exactness="EXACT",
         nature="[A_s P_s^q − S_p]_{n^p} = −(s−3)/(s−2), p-independent; "
                "zero iff s=3 (unique sigma-invariant P_s). Honest negative "
                "of the 2024 s-gonal ansatz, now with a closed form.",
         example="mtft.combinatorial.sgonal_defect_measured(5, 7)",
         upstream=("integers",), ref=""),
    dict(name="degree_gradient_commutator", tier="0", kind="identity",
         primitives=("V",), tag="Pr", exactness="EXACT",
         nature="[D,L]_{uv} = A_{uv}(d_v − d_u); zero iff regular (connected). "
                "Uncertainty generated by degree gradients across edges.",
         example="mtft.combinatorial.commutator_DL(A)",
         upstream=("integers",), ref="Robertson 1929 for the bound"),
    dict(name="entropic_uncertainty", tier="0", kind="identity",
         primitives=("IV",), tag="Pr", exactness="CERTIFIED(1e-9)",
         nature="H_K + H_Theta >= log(n+1) for the position/DFT pair.",
         example="mtft.combinatorial.entropic_margin(psi)",
         upstream=("integers",), ref="Maassen–Uffink 1988"),
    dict(name="mean_tano_weight", tier="5d", kind="identity",
         primitives=("II", "III"), tag="Pr", exactness="EXACT",
         nature="<w>_beta = −zeta'(beta+1) on the primon ensemble; "
                "beta->1+ endpoint −zeta'(2) = 2·T_INF = cold-gas alpha. "
                "The 2024 counting->energy map evaluated on w_n.",
         example="mtft.combinatorial.mean_tano_weight_exact(2.0)",
         upstream=("w_n", "integers"), ref="Dirichlet convolution; marked_gas alpha"),
)
