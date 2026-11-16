r"""mtft.gravity.tetrad — a tetrad for the squared hopfion, frozen into its null flow; the two helicities (MODEL; v1.1.0, 2026-11-16).

From the lab's record MP23 (14 November 2026; derived there, tested in space, unreviewed).  To first order a
tetrad is e^μ = dx^μ + B^μ, with metric perturbation h_μν = B_μν + B_νμ.  With a = x − iy, b = t − i − z,
d = r² − (t − i)², β = 2a/d, γ = 1 − α = 2ib/d and the constant null covectors ∂a = (0, 1, −i, 0), ∂b = (1, 0, 0, −1):

    B_μ = Re[(i a² (∂_μ b) dγ − b² (∂_μ a) dβ) / d²],
    𝓕_μ = dB_μ = Re[G_μ dα ∧ dβ],   G_μ = (i/d)(a ∂_μ b − b ∂_μ a).

Each of the four field strengths is the light's own 2-form dα ∧ dβ times a weight that is constant along every
ray, so the tetrad is frozen into the null flow k: i_k 𝓕_μ = 0.  The curvature of h is the squared hopfion
W = n ⊗ n/d⁵ (E = Re W, B_ij = ½ε_jkl R_0ikl = −Im W), and its Ricci tensor vanishes.

The helicity 3-form of a frame is e_a ∧ de^a, the Nieh–Yan 3-form of a teleparallel frame; its part of second
order is B_a ∧ dB^a.  For this frame it lives on the label sphere, so it does not change in time:

    H_g = ∫ B_a ∧ dB^a = (2/3)π²,      H_em = ∫ A ∧ dA = π²   for the light field F = n/d³.

The number belongs to the frame.  In transverse-traceless gauge Huang, Long, Zhao and Zhou (arXiv:2608.30154, §6.4)
find a helicity of the same kind for the type-N hopfion that changes in time; the two normalisations have not been
compared.  The lab's fit found this frame to be the only one of its holomorphic class up to degree 4; frames
outside that class are not excluded.  The ratio 2/3 comes out of one integral and has no mechanism (AG-D5).

Recorded misses of MP23, kept: the forecast that the square of the light's potential, d·σ ⊗ σ, is a potential for
the squared hopfion (it is not: its own Ricci tensor is 5.9 times |W|); and a prediction of 3/4 for the ratio, from
a slip of algebra.

Derivatives here are exact: every quantity is a rational function of (t, x, y, z), carried with its first and
second derivatives (``Jet``).  ``exhibit`` runs in double precision, at any number of digits (mpmath), or over the
Gaussian rationals, where the differences are exactly zero.  ``certificate`` turns exact zeros at random rational
events into a statement about every event, with a stated chance of error (Schwartz–Zippel).
"""
from __future__ import annotations

import math
from fractions import Fraction

import mpmath as mp
import numpy as np

H_EM = math.pi ** 2
H_G = 2 * math.pi ** 2 / 3
ETA = (-1, 1, 1, 1)
DEGREE_BOUND = 22                                    # of the numerators of the differences that ``exhibit`` tests (see ``certificate``)
_EPS3 = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}


class GQ:
    """An exact Gaussian rational, re + i·im with Fractions: just enough arithmetic for ``Jet``."""
    __slots__ = ("im", "re")

    def __init__(self, re=0, im=0):
        self.re, self.im = Fraction(re), Fraction(im)

    @staticmethod
    def of(v):
        return v if isinstance(v, GQ) else GQ(v)

    real = property(lambda self: self.re)
    imag = property(lambda self: self.im)

    def conjugate(self):
        return GQ(self.re, -self.im)

    def __add__(self, o):
        o = GQ.of(o)
        return GQ(self.re + o.re, self.im + o.im)
    __radd__ = __add__

    def __neg__(self):
        return GQ(-self.re, -self.im)

    def __sub__(self, o):
        o = GQ.of(o)
        return GQ(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return GQ.of(o) - self

    def __mul__(self, o):
        o = GQ.of(o)
        return GQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = GQ.of(o)
        n = o.re * o.re + o.im * o.im
        return GQ((self.re * o.re + self.im * o.im) / n, (self.im * o.re - self.re * o.im) / n)

    def __rtruediv__(self, o):
        return GQ.of(o) / self

    def __eq__(self, o):
        o = GQ.of(o)
        return self.re == o.re and self.im == o.im

    def __hash__(self):
        return hash((self.re, self.im))

    def __repr__(self):
        return f"GQ({self.re}, {self.im})"


_SYSTEMS = {"float": (complex, 1j), "mp": (mp.mpc, mp.mpc(0, 1)), "exact": (GQ, GQ(0, 1))}


class Jet:
    """A value with its gradient and its Hessian in the four coordinates: exact first and second derivatives of any
    expression built with +, −, × and ÷."""
    __slots__ = ("g", "h", "v")

    def __init__(self, v, g, h):
        self.v, self.g, self.h = v, g, h

    @staticmethod
    def variable(v, index, zero, one):
        return Jet(v, [one if i == index else zero for i in range(4)], [[zero] * 4 for _ in range(4)])

    def _scale(self, s):
        return Jet(self.v * s, [x * s for x in self.g], [[x * s for x in row] for row in self.h])

    def __add__(self, o):
        if isinstance(o, Jet):
            return Jet(self.v + o.v, [x + y for x, y in zip(self.g, o.g)], [[x + y for x, y in zip(r, s)] for r, s in zip(self.h, o.h)])
        return Jet(self.v + o, self.g, self.h)
    __radd__ = __add__

    def __neg__(self):
        return self._scale(-1)

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Jet):
            return self._scale(o)
        a, b = self, o
        return Jet(a.v * b.v, [a.g[i] * b.v + a.v * b.g[i] for i in range(4)],
                   [[a.h[i][j] * b.v + a.g[i] * b.g[j] + a.g[j] * b.g[i] + a.v * b.h[i][j] for j in range(4)] for i in range(4)])
    __rmul__ = __mul__

    def reciprocal(self):
        r = 1 / self.v
        r2, r3 = r * r, r * r * r
        return Jet(r, [-x * r2 for x in self.g], [[-self.h[i][j] * r2 + 2 * self.g[i] * self.g[j] * r3 for j in range(4)] for i in range(4)])

    def __truediv__(self, o):
        return self * o.reciprocal() if isinstance(o, Jet) else self._scale(1 / o)

    def __rtruediv__(self, o):
        return self.reciprocal() * o


def _setup(x, scalar):
    """The jets of a, b, d and of the gradients of β and γ at the event x = (t, x, y, z), in the number system asked for."""
    conv, I = _SYSTEMS[scalar]
    if scalar == "exact":
        vals = [GQ(Fraction(v)) for v in x]
    elif scalar == "mp":
        vals = [mp.mpc(v) for v in x]
    else:
        vals = [complex(float(v)) for v in x]
    zero, one = conv(0), conv(1)
    t, X, Y, Z = (Jet.variable(vals[i], i, zero, one) for i in range(4))
    a, b = X - Y * I, t - I - Z
    tm = t - I
    d = X * X + Y * Y + Z * Z - tm * tm
    da, db = [zero, one, -I, zero], [one, zero, zero, -one]
    dd_ = [tm * (-2), X * 2, Y * 2, Z * 2]
    inv = d.reciprocal()
    inv2 = inv * inv
    dal = [inv * (-2 * I * db[m]) + b * dd_[m] * inv2 * (2 * I) for m in range(4)]
    dbe = [inv * (2 * da[m]) - a * dd_[m] * inv2 * 2 for m in range(4)]
    return {"I": I, "zero": zero, "a": a, "b": b, "d": d, "inv": inv, "inv2": inv2, "da": da, "db": db, "dal": dal, "dbe": dbe}


def _complex_potential(x, scalar="float"):
    """C_μν with B_μν = Re C_μν, as jets: (i a² (∂_μ b)(∂_ν γ) − b² (∂_μ a)(∂_ν β))/d², with ∂γ = −∂α."""
    s = _setup(x, scalar)
    a2, b2, I = s["a"] * s["a"], s["b"] * s["b"], s["I"]
    return [[(a2 * (I * s["db"][m]) * (-s["dal"][n]) - b2 * s["da"][m] * s["dbe"][n]) * s["inv2"] for n in range(4)] for m in range(4)], s


def _parts(x, scalar="float"):
    """B, its first and its second derivatives, as real numbers of the chosen system: B[μ][ν], dB[ρ][μ][ν], d2B[ρ][σ][μ][ν]."""
    C, s = _complex_potential(x, scalar)
    R4 = range(4)
    B = [[C[m][n].v.real for n in R4] for m in R4]
    dB = [[[C[m][n].g[r].real for n in R4] for m in R4] for r in R4]
    d2B = [[[[C[m][n].h[r][q].real for n in R4] for m in R4] for q in R4] for r in R4]
    return B, dB, d2B, s


def _fields(s):
    """n, W = n ⊗ n/d⁵ and the null flow k = (1, v) from a setup, in its number system."""
    a, b, d, I = s["a"].v, s["b"].v, s["d"].v, s["I"]
    n = [4 * (b * b - a * a), -4 * I * (a * a + b * b), 8 * a * b]
    d5 = d * d * d * d * d
    W = [[n[i] * n[j] / d5 for j in range(3)] for i in range(3)]
    nc = [v.conjugate() for v in n]
    cr = [n[1] * nc[2] - n[2] * nc[1], n[2] * nc[0] - n[0] * nc[2], n[0] * nc[1] - n[1] * nc[0]]
    nn = (n[0] * nc[0] + n[1] * nc[1] + n[2] * nc[2]).real
    return n, W, [1] + [(I * cr[i]).real / nn for i in range(3)]


def _riemann(d2B):
    """The linearised Riemann tensor of h = B + Bᵀ: R_mnrs = ½(∂_n∂_r h_ms + ∂_m∂_s h_nr − ∂_m∂_r h_ns − ∂_n∂_s h_mr)."""
    R4 = range(4)

    def D(r, q, m, n):
        return d2B[r][q][m][n] + d2B[r][q][n][m]
    return [[[[(D(n, r, m, s_) + D(m, s_, n, r) - D(m, r, n, s_) - D(n, s_, m, r)) / 2 for s_ in R4] for r in R4] for n in R4] for m in R4]


def exhibit(x, scalar: str = "float") -> dict:
    """T-G1's criterion at one event x = (t, x, y, z).

    Returns the sizes, relative to the field, of: i_k 𝓕 (``frozen``), E − Re W and B + Im W (``E``, ``B``), the Ricci
    tensor (``ricci``) and R·k (``Rk``).  ``scalar`` is "float", "mp" (at the current mpmath precision) or "exact"
    (x given as fractions or integers; then the differences are exact and ``exact`` says whether all are zero)."""
    _, dB, d2B, s = _parts(x, scalar)
    R4, R3 = range(4), range(3)
    _, W, k = _fields(s)
    F = [[[dB[r][m][n] - dB[n][m][r] for n in R4] for m in R4] for r in R4]                 # F[ρ][μ][ν] = ∂_ρ B_μν − ∂_ν B_μρ
    froz = [sum(k[r] * F[r][m][n] for r in R4) for m in R4 for n in R4]
    R = _riemann(d2B)
    E = [[R[0][i + 1][0][j + 1] for j in R3] for i in R3]
    Bm = [[sum(e * R[0][i + 1][kk + 1][ll + 1] for (jj, kk, ll), e in _EPS3.items() if jj == j) / 2 for j in R3] for i in R3]
    ric = [sum(ETA[r] * R[r][m][r][n] for r in R4) for m in R4 for n in R4]
    Rk = [sum(R[m][n][r][q] * k[q] for q in R4) for m in R4 for n in R4 for r in R4]
    dE = [E[i][j] - W[i][j].real for i in R3 for j in R3]
    dBm = [Bm[i][j] + W[i][j].imag for i in R3 for j in R3]
    if scalar == "exact":
        zero = all(v == 0 for v in froz + dE + dBm + ric + Rk)
        return {"exact": zero, "frozen": 0.0 if zero else float("nan"), "E": 0.0 if zero else float("nan"), "B": 0.0 if zero else float("nan"),
                "ricci": 0.0 if zero else float("nan"), "Rk": 0.0 if zero else float("nan")}

    def big(vals):
        return max(abs(v) for v in vals)
    sF, sW, sR = big([v for a in F for b in a for v in b]), big([abs(v) for row in W for v in row]), big([v for a in R for b in a for c in b for v in c])
    return {"frozen": float(big(froz) / sF), "E": float(big(dE) / sW), "B": float(big(dBm) / sW), "ricci": float(big(ric) / sW), "Rk": float(big(Rk) / sR)}


def certificate(points: int = 10, grid: int = 2000, seed: int = 143) -> dict:
    """A randomized certificate that the identities of ``exhibit`` hold at every event, not only where they were tried.

    Each of the 114 differences (i_k 𝓕: 16; E − Re W and B + Im W: 9 each; Ricci: 16; R·k: 64) is a polynomial of total
    degree at most ``DEGREE_BOUND`` in (t, x, y, z) over a denominator that has no real zero.  The count: B_μν = Re(N/d⁴)
    with deg N ≤ 4; each derivative raises the numerator's degree by one and the power of d by one; a real part of
    N/d^k has a numerator of degree deg N + 2k; W's real and imaginary parts have degree 14 over |d|¹⁰; and the flow k
    brings a numerator and a denominator of degree 4.  That gives 19 for i_k 𝓕, 18 for the curvature and 22 for R·k.

    By the Schwartz–Zippel lemma a non-zero polynomial of degree D vanishes at a point whose coordinates are drawn
    independently and uniformly from a set of S values with probability at most D/S.  Here every coordinate is drawn
    from the 2·grid + 1 values k/1000, |k| ≤ grid, and the differences are evaluated exactly over the Gaussian
    rationals.  ``false_pass_below`` is the bound, over all 114 differences, on the chance that an identity which is not
    true passes at all the points.  The points come from a seeded generator; the bound is the lemma's, for a uniform draw."""
    import random
    rng = random.Random(seed)
    values = 2 * grid + 1
    holds = all(exhibit(tuple(Fraction(rng.randint(-grid, grid), 1000) for _ in range(4)), "exact")["exact"] for _ in range(points))
    return {"holds": holds, "points": points, "values_per_coordinate": values, "degree_bound": DEGREE_BOUND, "differences": 114,
            "false_pass_below": 114 * (DEGREE_BOUND / values) ** points}


def potential(x) -> np.ndarray:
    """B_μν at the event x (first index the leg, second the form index)."""
    return np.array(_parts(x)[0], float)


def field_strength(x) -> np.ndarray:
    """𝓕[ρ, μ, ν] = ∂_ρ B_μν − ∂_ν B_μρ: the 2-form (in ρ, ν) of the leg μ."""
    dB = np.array(_parts(x)[1], float)
    return dB - np.transpose(dB, (2, 1, 0))


def weight(x) -> np.ndarray:
    """G_μ = (i/d)(a ∂_μ b − b ∂_μ a), the weight of each leg: constant along every ray and linear in the labels."""
    s = _setup(x, "float")
    a, b, d = s["a"].v, s["b"].v, s["d"].v
    return np.array([1j / d * (a * s["db"][m] - b * s["da"][m]) for m in range(4)])


def strength_from_weight(x) -> np.ndarray:
    """Re[G_μ (dα ∧ dβ)_ρν], in the layout of ``field_strength``: the second route to the field strength."""
    s = _setup(x, "float")
    dal, dbe = np.array([j.v for j in s["dal"]]), np.array([j.v for j in s["dbe"]])
    two = dal[:, None] * dbe[None, :] - dbe[:, None] * dal[None, :]
    return np.real(weight(x)[None, :, None] * two[:, None, :])


def metric(x) -> np.ndarray:
    """h_μν = B_μν + B_νμ."""
    B = potential(x)
    return B + B.T


def curvature(x) -> dict:
    """E_ij = R_0i0j, B_ij = ½ε_jkl R_0ikl and the Ricci tensor of h at x, with the squared hopfion W = n ⊗ n/d⁵."""
    _, _, d2B, s = _parts(x)
    R = np.array(_riemann(d2B), float)
    eps = np.zeros((3, 3, 3))
    for key, e in _EPS3.items():
        eps[key] = e
    return {"E": R[0, 1:, 0, 1:], "B": 0.5 * np.einsum("jkl,ikl->ij", eps, R[0, 1:, 1:, 1:]), "ricci": np.einsum("r,rmrn->mn", np.array(ETA, float), R),
            "W": np.array(_fields(s)[1], complex), "riemann": R}


def null_flow(x) -> np.ndarray:
    """k = (1, v): the null direction of the light rays at x."""
    return np.array(_fields(_setup(x, "float"))[2], float)


# ----------------------------------------------------------------------------------------------------------------- helicities
_DA, _DB = np.array([0, 1, -1j, 0]), np.array([1, 0, 0, -1], complex)


def _three_forms(al, be, dal, dbe):
    """B_a ∧ dB^a and A ∧ dA on three tangent vectors, from the derivatives of the labels along them (arrays [3, …])."""
    gam, dgam = 1 - al, -dal
    Bv = 0.25 * np.real(gam ** 2 * _DA[:, None, None] * dbe[None] + 1j * be ** 2 * _DB[:, None, None] * dgam[None])

    def pair(u, v, j, k):
        return u[j] * v[k] - u[k] * v[j]

    def dBv(j, k):
        return 0.5 * np.real((gam * _DA[:, None] - 1j * be * _DB[:, None]) * pair(dgam, dbe, j, k)[None])
    Kg = sum(ETA[m] * (Bv[m, 0] * dBv(1, 2)[m] - Bv[m, 1] * dBv(0, 2)[m] + Bv[m, 2] * dBv(0, 1)[m]) for m in range(4))
    Av = 0.5 * np.real(al * dbe - be * dal)

    def dAv(j, k):
        return np.real(pair(dal, dbe, j, k))
    return Kg, Av[0] * dAv(1, 2) - Av[1] * dAv(0, 2) + Av[2] * dAv(0, 1)


def helicities_on_sphere(n_eta: int = 40, n_xi: int = 48) -> dict:
    """H_g and H_em as integrals over the label sphere (Gauss–Legendre in η, the trapezoid rule in ξ₁ and ξ₂)."""
    g, w = np.polynomial.legendre.leggauss(n_eta)
    eta, w = (g + 1) * math.pi / 4, w * math.pi / 4
    xi = np.arange(n_xi) * 2 * math.pi / n_xi
    Et, X1, X2 = np.meshgrid(eta, xi, xi, indexing="ij")
    al, be = np.cos(Et) * np.exp(1j * X1), np.sin(Et) * np.exp(1j * X2)
    dal = np.array([1j * al, -np.sin(Et) * np.exp(1j * X1), 0 * al])                      # along ξ₁, η, ξ₂: the orientation space gives the sphere
    dbe = np.array([0 * be, np.cos(Et) * np.exp(1j * X2), 1j * be])
    Kg, Ke = _three_forms(al.reshape(-1), be.reshape(-1), dal.reshape(3, -1), dbe.reshape(3, -1))
    wgt = (w[:, None, None] * np.ones_like(Et) * (2 * math.pi / n_xi) ** 2).reshape(-1)
    Hg, He = float(np.sum(Kg * wgt)), float(np.sum(Ke * wgt))
    return {"H_g": Hg, "H_em": He, "ratio": Hg / He}


def helicities_on_slice(t: float, n: int = 88, L: float = 1.2) -> dict:
    """The same two integrals over the slice of time t in space, mapped to a cube: they do not depend on t."""
    g, w = np.polynomial.legendre.leggauss(n)
    xs, ws = L * g / (1 - g * g), w * L * (1 + g * g) / (1 - g * g) ** 2
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    Wt = (ws[:, None, None] * ws[None, :, None] * ws[None, None, :]).reshape(-1)
    X, Y, Z = X.reshape(-1), Y.reshape(-1), Z.reshape(-1)
    a, b = X - 1j * Y, t - 1j - Z
    d = X * X + Y * Y + Z * Z - (t - 1j) ** 2
    al, be = 1 - 2j * b / d, 2 * a / d
    da, db, dd_ = [1, -1j, 0], [0, 0, -1], [2 * X, 2 * Y, 2 * Z]
    dal = np.array([-2j * db[i] / d + 2j * b * dd_[i] / d ** 2 for i in range(3)])
    dbe = np.array([2 * da[i] / d - 2 * a * dd_[i] / d ** 2 for i in range(3)])
    Kg, Ke = _three_forms(al, be, dal, dbe)
    Hg, He = float(np.sum(Kg * Wt)), float(np.sum(Ke * Wt))
    return {"H_g": Hg, "H_em": He, "ratio": Hg / He}
