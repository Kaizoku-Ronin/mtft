r"""mtft.matter.twist_costs — what a twist buys and costs at the Pati–Salam point (EXACT at one loop; TC-01, 2026-11-19).

The setting is VE-01's.  Eight Wilson-line eigenvalues θ_i ∈ (ℝ/ℤ)² sit on the unit-area torus of
143a1 (stacks c, A1, L, B, A2 of ranks 3, 1, 2, 1, 1), and a twist puts a weight W(w) ≥ 0 on the
windings w = (m, n):

    V(θ) = −Σ_w W(w) |Tr U^w|²,      Tr U^w = Σ_i exp(2πi w·θ_i)      (in units of c = 1.5 r/π⁸ M_KK⁴).

The weights.  K(w) = |w|⁻⁸ with |w|² = |m + nτ|²/Im τ (VE-01's leading term; ``a1=True`` adds
(2π/3r)|w|⁻⁶).  A change of spin structure P ∈ E[2] (MATTER-04) gives W = K on the windings with
χ_P(w) = −1 and 0 on the rest; a continuous twist P = (0, (1 − δ)/2), possible only outside the
Calabi–Yau setting, gives W = K sin⁴(π⟨w, P⟩) (FF-9's f_vec/16); a spin structure ζ of the curve adds
χ_P(w) κ_ζ(w)/2 (FF-14).  Any non-negative weight can be passed in.

The Pati–Salam pattern on an offset h ∈ E[2]∖O puts c and A1 at h and L, B, A2 at 0.  With s_i = −1
on c, A1 and +1 elsewhere, and the integer windings w = (m, n) as vectors,

    M_A = Σ_{χ_h(w)=+1} W(w) w wᵀ,      M_B = Σ_{χ_h(w)=−1} W(w) w wᵀ,
    Hess V = −8π² M_B ⊗ s sᵀ + 8π² M_A ⊗ (8·1 − J)          (16 × 16, ``pattern_hessian``),

and V = −64 Σ_{χ_h(w)=+1} W(w) there.  The relative U(1) T_S = s ⊗ e has curvature 64π² eᵀ(M_A − M_B)e;
the twelve flat directions of the balanced point have 64π² eig(M_A), six times each; two directions
(the overall U(1)) are exactly flat.  Route 2 is the direct sum and its Richardson second differences
(``hessian_by_differences``), compared on the Hessian's own scale: a direction whose curvature
nearly cancels makes a relative error per direction meaningless (FF-13's recorded metric fault).

The twist table (EXACT, MATTER-04 and YK-01 §2).  For P ∈ E[2] the quark–lepton fermions in the
offset-h blocks are untwisted exactly when P = h, and Higgs bosons never see P.  So at the pattern

    tree Yukawas  ⟺  P = h  ⟺  V = 0  ⟺  M_A = 0  ⟺  the twelve are flat at one loop,

on all nine pairs (``twist_row``, ``complementarity``): among the spin structures of 143a1 the pin of
the twelve and the tree Yukawas exclude each other.  The block of the condensate's boson (A1–A2) is
twisted in every row.

The continuous twist's spectrum (EXACT, KK masses on the torus; SC-01's M_SS = 1/(2√Im τ) M_KK).  At
P = (0, (1 − δ)/2) on the ledger's offset h = p₂: quark–lepton fermions δ·M_SS (vector-like pairs),
Φ-type Higgs bosons 2δ·M_SS, Φ in an offset-h block (1 − 2δ)·M_SS (``continuous_masses``, with a
brute-force minimum over the dual lattice as the second route).

The pull on τ (ESTIMATE; FF-4, FF-5).  The γ-flux holds τ at 143a1 with the well V_γ0 (cosh d − 1),
V_γ0 = 8/(π² g₄² r) M_KK⁴.  A one-loop energy V(τ) (ledger units) pulls τ off by sinh d = c |∇V|_hyp/V_γ0
(``hold_pull``), over the window r ∈ {230, 280}, g₄² ∈ {0.31, 0.37}.  The jet vacuum's pull is
8.7–10.6 M_KK⁴ (FF-4), the gate of this function.

Provenance: records FF-12 and FF-13 (19 November 2026), project experiments FF12_Hold_Plan.md and
FF13_Twist_Plan.md.  What this is not: the Calabi–Yau setting allows only the four points of E[2];
a continuous twist is a model setting, and nothing here fixes the twist dynamically.
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction as Fr
from functools import lru_cache

import numpy as np

STACKS = ("c", "A1", "L", "B", "A2")
RANK = {"c": 3, "A1": 1, "L": 2, "B": 1, "A2": 1}
STACK_OF = tuple(s for s in STACKS for _ in range(RANK[s]))
LEFT = ("c", "A1")
"""The stacks at the offset h in the Pati–Salam pattern."""

S_SIGN = np.array([-1.0 if s in LEFT else 1.0 for s in STACK_OF])

E2 = {"O": (Fr(0), Fr(0)), "p1": (Fr(1, 2), Fr(0)), "p2": (Fr(0), Fr(1, 2)), "p3": (Fr(1, 2), Fr(1, 2))}
"""The 2-torsion points of 143a1 in torus coordinates (1-cycle, τ-cycle); the ledger's twist is p2."""

MATTER = ("c-L", "c-B", "c-A2", "A1-L", "A1-B", "A1-A2")
"""The quark–lepton blocks: Q, d^c, u^c, L, e^c, ν^c."""

HIGGS = ("L-B", "L-A2")
"""The blocks of the Yukawa (Φ-type) Higgs bosons."""

WINDOW = tuple((r, g2) for r in (230.0, 280.0) for g2 in (0.31, 0.37))
"""(r, g₄²) over which the hold is judged (FH-3, FF-4)."""


def tau() -> complex:
    """τ of 143a1 (the ledger's torus)."""
    from mtft.matter.vacuum_energy import tau as _tau
    return _tau()


def m_ss(t: complex | None = None) -> float:
    """SC-01's Scherk–Schwarz scale 1/(2√Im τ), in units of M_KK."""
    t = tau() if t is None else t
    return 1 / (2 * math.sqrt(t.imag))


def c_factor(r: float) -> float:
    """c = 1.5 r/π⁸ (M_KK⁴): ledger units to M_KK⁴ (VE-01)."""
    return 1.5 * r / math.pi ** 8


def well(r: float, g2: float) -> float:
    """V_γ0 = 8/(π² g₄² r) (M_KK⁴): the depth scale of the γ-flux's well (FF-4)."""
    return 8 / (math.pi ** 2 * g2 * r)


# ------------------------------------------------------------------ windings and weights
@lru_cache(maxsize=16)
def windings(N: int = 40, t: complex | None = None):
    """(m, n, |w|²) for 0 < max(|m|, |n|) ≤ N on the unit-area torus of τ (default 143a1)."""
    if N < 1:
        raise ValueError("N must be at least 1")
    t = tau() if t is None else complex(t)
    ax = np.arange(-N, N + 1)
    m, n = np.meshgrid(ax, ax, indexing="ij")
    m, n = m.ravel(), n.ravel()
    keep = ~((m == 0) & (n == 0))
    m, n = m[keep].astype(float), n[keep].astype(float)
    w2 = np.abs(m + n * t) ** 2 / t.imag
    for a in (m, n, w2):
        a.setflags(write=False)
    return m, n, w2


def K(w2: np.ndarray, r: float | None = None) -> np.ndarray:
    """VE-01's weight |w|⁻⁸, plus (2π/3r)|w|⁻⁶ when r is given."""
    out = w2 ** -4
    if r is not None:
        out = out + (2 * math.pi / (3 * r)) * w2 ** -3
    return out


def _point(P):
    if len(P) != 2:
        raise ValueError("a twist is a point (P1, P2) of the torus")
    return Fr(P[0]) % 1, Fr(P[1]) % 1


def is_two_torsion(P) -> bool:
    p = _point(P)
    return all((2 * x) % 1 == 0 for x in p)


def chi(P, m: np.ndarray, n: np.ndarray) -> np.ndarray:
    """χ_P(w) = exp(2πi⟨w, P⟩) = ±1 for P in E[2], exactly."""
    p = _point(P)
    if not is_two_torsion(p):
        raise ValueError("χ_P is a sign only for P in E[2]")
    par = (np.rint(2 * float(p[0]) * m) + np.rint(2 * float(p[1]) * n)).astype(int) % 2
    return 1 - 2 * par


def weights_spin(P, N: int = 40, t: complex | None = None, r: float | None = None) -> np.ndarray:
    """W = K on the windings with χ_P(w) = −1, else 0: a change of spin structure P ∈ E[2]."""
    m, n, w2 = windings(N, t)
    return np.where(chi(P, m, n) == -1, K(w2, r), 0.0)


def weights_continuous(delta: float, N: int = 40, t: complex | None = None, r: float | None = None) -> np.ndarray:
    """W = K sin⁴(π n (1 − δ)/2): the continuous twist P = (0, (1 − δ)/2), a distance δ from p2 along the τ-cycle."""
    if not 0 <= delta <= 1:
        raise ValueError("δ must lie in [0, 1]")
    m, n, w2 = windings(N, t)
    return K(w2, r) * np.sin(math.pi * n * (1 - delta) / 2) ** 4


# ------------------------------------------------------------------ the pattern, route 1: closed forms
def pattern(h) -> np.ndarray:
    """θ (8 × 2) of the Pati–Salam pattern on the offset h: c and A1 at h, L, B and A2 at 0."""
    p = _point(h)
    if p == (0, 0):
        raise ValueError("the pattern needs a nonzero offset h")
    return np.array([[float(p[0]), float(p[1])] if s in LEFT else [0.0, 0.0] for s in STACK_OF])


def pattern_matrices(W: np.ndarray, h, N: int = 40, t: complex | None = None):
    """(M_A, M_B): Σ W w wᵀ over the windings with χ_h = +1 and with χ_h = −1."""
    m, n, _ = windings(N, t)
    W = np.asarray(W, float)
    if W.shape != m.shape:
        raise ValueError("the weight must have one entry per winding of the grid")
    if (W < 0).any():
        raise ValueError("weights must be non-negative")
    ch = chi(h, m, n)
    w = np.stack([m, n], 1)
    A, B = ch == 1, ch == -1
    MA = (W[A, None, None] * w[A][:, :, None] * w[A][:, None, :]).sum(0)
    MB = (W[B, None, None] * w[B][:, :, None] * w[B][:, None, :]).sum(0)
    return MA, MB


def pattern_hessian(W: np.ndarray, h, N: int = 40, t: complex | None = None) -> np.ndarray:
    """The 16 × 16 Hessian of V at the pattern, flattened as (eigenvalue, cycle)."""
    MA, MB = pattern_matrices(W, h, N, t)
    return -8 * math.pi ** 2 * np.kron(np.outer(S_SIGN, S_SIGN), MB) + 8 * math.pi ** 2 * np.kron(8 * np.eye(8) - np.ones((8, 8)), MA)


def pattern_energy(W: np.ndarray, h, N: int = 40, t: complex | None = None) -> float:
    """V at the pattern: −64 Σ_{χ_h(w)=+1} W(w)."""
    m, n, _ = windings(N, t)
    return float(-64 * np.asarray(W, float)[chi(h, m, n) == 1].sum())


def pattern_curvatures(W: np.ndarray, h, N: int = 40, t: complex | None = None) -> dict:
    """V, the T_S curvatures 64π² eig(M_A − M_B), the twelve's 64π² eig(M_A) and the Hessian's signature."""
    MA, MB = pattern_matrices(W, h, N, t)
    H = pattern_hessian(W, h, N, t)
    ev = np.linalg.eigvalsh(H)
    sc = max(float(np.abs(ev).max()), 1e-300)
    sig = (int((ev < -1e-9 * sc).sum()), int((np.abs(ev) <= 1e-9 * sc).sum()), int((ev > 1e-9 * sc).sum()))
    return {"V": pattern_energy(W, h, N, t), "T_S": tuple(64 * math.pi ** 2 * np.linalg.eigvalsh(MA - MB)),
            "twelve": tuple(64 * math.pi ** 2 * np.linalg.eigvalsh(MA)), "signature": sig, "MA": MA, "MB": MB}


# ------------------------------------------------------------------ route 2: the direct sum
def V_direct(theta: np.ndarray, W: np.ndarray, N: int = 40, t: complex | None = None) -> float:
    """V(θ) = −Σ_w W(w)|Tr U^w|² summed directly."""
    m, n, _ = windings(N, t)
    th = np.asarray(theta, float)
    ph = 2 * math.pi * (np.outer(th[:, 0], m) + np.outer(th[:, 1], n))
    tr2 = np.cos(ph).sum(0) ** 2 + np.sin(ph).sum(0) ** 2
    return float(-(np.asarray(W, float) * tr2).sum())


def hessian_by_differences(W: np.ndarray, h, directions, step: float = 1e-3, N: int = 40, t: complex | None = None) -> list:
    """vᵀ(Hess V)v along unit directions v (8 × 2 arrays) by Richardson second differences of the direct sum."""
    th = pattern(h)
    out = []
    for d in directions:
        d = np.asarray(d, float)
        nrm = np.linalg.norm(d)
        if nrm == 0:
            raise ValueError("a direction must be nonzero")
        d = d / nrm

        def g(u, d=d):
            return (V_direct(th + u * d, W, N, t) - 2 * V_direct(th, W, N, t) + V_direct(th - u * d, W, N, t)) / u ** 2
        out.append((4 * g(step / 2) - g(step)) / 3)
    return out


def standard_directions(seed: int = 1313, n_random: int = 4) -> list:
    """T_S on both cycles, eight relative directions and a few random relative ones."""
    rng = np.random.default_rng(seed)
    dirs = [np.outer(S_SIGN, e) for e in ([1.0, 0.0], [0.0, 1.0])]
    rel = [np.array([1, -1, 0, 0, 0, 0, 0, 0.0]), np.array([1, 1, -2, 0, 0, 0, 0, 0.0]),
           np.array([0, 0, 0, 0, 1, -1, 0, 0.0]), np.array([0, 0, 0, 1, -1, -1, 1, 0.0])]
    for v in rel:
        for e in ([1.0, 0.0], [0.0, 1.0]):
            dirs.append(np.outer(v, e))
    for _ in range(n_random):
        d = rng.normal(size=(8, 2))
        dirs.append(d - d.mean(0))
    return dirs


def compare_routes(W: np.ndarray, h, N: int = 40, t: complex | None = None, step: float = 1e-3) -> float:
    """max |closed − differences| over the standard directions, relative to the Hessian's largest eigenvalue."""
    H = pattern_hessian(W, h, N, t)
    sc = max(float(np.abs(np.linalg.eigvalsh(H)).max()), 1.0)
    dirs = standard_directions()
    fd = hessian_by_differences(W, h, dirs, step, N, t)
    worst = 0.0
    for d, f in zip(dirs, fd):
        v = (d / np.linalg.norm(d)).reshape(-1)
        worst = max(worst, abs(float(v @ H @ v) - f) / sc)
    return worst


# ------------------------------------------------------------------ the twist table
def block_status(P, h) -> dict:
    """For each block of stacks: its offset, and whether fermions (which see P) and bosons (which do not) are untwisted."""
    p, hh = _point(P), _point(h)
    out = {}
    for i, j in itertools.combinations(STACKS, 2):
        off = hh if ((i in LEFT) != (j in LEFT)) else (Fr(0), Fr(0))
        out[f"{i}-{j}"] = {"offset": off, "fermions untwisted": tuple((off[k] + p[k]) % 1 for k in (0, 1)) == (0, 0),
                           "bosons untwisted": off == (0, 0)}
    return out


def twist_row(P, h, N: int = 40, t: complex | None = None) -> dict:
    """One row of the twist table: spin structure P ∈ E[2] of 143a1 on the pattern with offset h ∈ E[2]∖O."""
    p, hh = _point(P), _point(h)
    if not (is_two_torsion(p) and is_two_torsion(hh)):
        raise ValueError("the twist table is for P and h in E[2]; use weights_continuous for a continuous twist")
    W = weights_spin(p, N, t)
    cv = pattern_curvatures(W, hh, N, t)
    st = block_status(p, hh)
    matter = all(st[b]["fermions untwisted"] for b in MATTER)
    higgs = all(st[b]["bosons untwisted"] for b in HIGGS)
    ma_zero = bool(np.abs(cv["MA"]).max() == 0)
    row = {"P": tuple(str(x) for x in p), "h": tuple(str(x) for x in hh), "V": cv["V"], "T_S": cv["T_S"], "twelve": cv["twelve"],
           "signature": cv["signature"], "matter fermions untwisted": matter, "Higgs bosons untwisted": higgs,
           "condensate boson untwisted": st["A1-A2"]["bosons untwisted"], "tree Yukawas": matter and higgs,
           "twelve pinned": bool(min(cv["twelve"]) > 0), "MA = 0": ma_zero}
    row["complementary"] = (matter == (p == hh)) and (matter == (cv["V"] == 0)) and (matter == ma_zero) and (ma_zero != row["twelve pinned"])
    return row


def twist_rows(N: int = 40, t: complex | None = None) -> dict:
    """The nine rows P, h ∈ {p1, p2, p3} (P = O is the untwisted, supersymmetric case and has W = 0)."""
    keys = ("p1", "p2", "p3")
    return {f"P={a}, h={b}": twist_row(E2[a], E2[b], N, t) for a in keys for b in keys}


def complementarity(N: int = 40, t: complex | None = None) -> bool:
    """Whether tree Yukawas, P = h, V = 0, M_A = 0 and flat twelve go together on all nine pairs."""
    return all(r["complementary"] for r in twist_rows(N, t).values())


# ------------------------------------------------------------------ the continuous twist's spectrum
def kk_min(x, kmax: int = 6, t: complex | None = None) -> float:
    """The lightest KK mass (M_KK units) of a field whose torus momenta are shifted by x = (x1, x2): brute force."""
    t = tau() if t is None else t
    y, best = t.imag, math.inf
    for k1 in range(-kmax, kmax + 1):
        for k2 in range(-kmax, kmax + 1):
            a, b = k1 + x[0], k2 + x[1]
            best = min(best, math.sqrt(y * a * a + (b - a * t.real) ** 2 / y))
    return best


def continuous_masses(delta: float, t: complex | None = None, brute: bool = False) -> dict:
    """The continuous twist P = (0, (1 − δ)/2) on the offset h = p2: the three lightest masses, closed form (and brute force)."""
    if not 0 < delta < 0.5:
        raise ValueError("δ must lie in (0, 1/2)")
    M = m_ss(t)
    out = {"matter fermion": delta * M, "Phi-type Higgs": 2 * delta * M, "Phi in an offset-h block": (1 - 2 * delta) * M}
    if brute:
        P = (0.0, (1 - delta) / 2)
        out = {"closed": out, "brute force": {
            "matter fermion": min(kk_min((P[0], 0.5 + P[1]), t=t), kk_min((-P[0], 0.5 - P[1]), t=t)),
            "Phi-type Higgs": kk_min((2 * P[0], 2 * P[1]), t=t),
            "Phi in an offset-h block": kk_min((2 * P[0], 0.5 + 2 * P[1]), t=t)}}
    return out


# ------------------------------------------------------------------ the pull on τ
def hyperbolic_gradient(V, t: complex, eps: float = 1e-5) -> float:
    """|∇V|_hyp = Im τ |∇V| at τ, by central differences of a function of τ."""
    fx = (V(t + eps) - V(t - eps)) / (2 * eps)
    fy = (V(t + 1j * eps) - V(t - 1j * eps)) / (2 * eps)
    return t.imag * math.hypot(fx, fy)


def hold_pull(V, t: complex | None = None, window=WINDOW, eps: float = 1e-5) -> dict:
    """sinh d = c(r)|∇V|_hyp / V_γ0(r, g₄²) over the window, for V(τ) in ledger units (ESTIMATE)."""
    if not window:
        raise ValueError("the window is empty")
    t = tau() if t is None else t
    g = hyperbolic_gradient(V, t, eps)
    sinhd = {f"{r},{g2}": c_factor(r) * g / well(r, g2) for r, g2 in window}
    return {"|grad V|_hyp": g, "V": V(t), "sinh d": sinhd, "max sinh d": max(sinhd.values()), "min sinh d": min(sinhd.values()),
            "held": max(sinhd.values()) < 1.0}


def pattern_energy_of_tau(P, h, N: int = 40):
    """τ ↦ V at the pattern for the spin structure P ∈ E[2] (ledger units), for ``hold_pull``."""
    def V(t):
        return pattern_energy(weights_spin(P, N, complex(t)), h, N, complex(t))
    return V


def jet_vacuum_pull(N: int = 40) -> tuple:
    """c|∇V|_hyp of the jet vacuum (P = p2, every offset 0) at r = 230 and 280: FF-4's 8.7 to 10.6 M_KK⁴."""
    def V(t):
        m, n, w2 = windings(N, complex(t))
        return float(-64 * K(w2)[chi(E2["p2"], m, n) == -1].sum())
    g = hyperbolic_gradient(V, tau())
    return c_factor(230.0) * g, c_factor(280.0) * g
