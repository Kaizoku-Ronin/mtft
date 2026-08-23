"""mtft.boundary — the cusp/elliptic boundary layer of the canonical ring.

The canonical quadratic differentials sit strictly inside the weight-4
Hecke module:

    0 -> S_4^(2)/H^0(2K) -> S_4/H^0(2K) -> S_4/S_4^(2) -> 0
         elliptic quotient    boundary B_N    cusp quotient

with dim B_N = nu_infinity + nu_2 + nu_3.  For squarefree N = pq the cusp
quotient is 4-dimensional (four cusps, on which the Atkin-Lehner group
acts freely transitively, so C[cusps] is the regular representation of
(Z/2)^2 and each character appears exactly once).

Leakage.  A Hecke operator applied to H^0(2K) generally leaves it.  The
composite

    l_A : H^0(2K) --A--> S_4 -->> B_N

has a *signature* (r_cusp, r_ell).  The arc's finding is that good primes
saturate both quotients while level primes degenerate along explicitly
identifiable channels.

NOT A DIRECT SUM.  B_N is not canonically B_cusp (+) B_ell; only the
short exact sequence above is canonical.  The two ranks are well defined
as successive quotients without choosing a splitting.

Provenance: `_data/PROVENANCE.txt`.  The GP scripts that produced the
census live in `studies/rt1b_2026aug/`, with `FREEZE_HASHES.txt`.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path

__all__ = [
    "LEVEL", "CENSUS", "RELATIONS_143", "THEOREM_TARGETS", "HYPOTHESES",
    "nu2", "nu3", "boundary_dim", "s4_qexpansions", "cusp_functionals",
    "operator", "data_path",
]

LEVEL = 143

#: Discovery census, all odd squarefree N = pq <= 160 with genus >= 3 that
#: pass canonical-multiplication surjectivity.  Signatures are
#: (r_cusp, r_ell).  `good` is the common signature of every tested good
#: prime at that level.  Source: studies/rt1b_2026aug/rt1b4_boundary_census.gp
CENSUS = {
    #  N: (p, q, genus, b_N, nu2, nu3, good, U_p, U_q)
    51: (3, 17, 5, 4, 0, 0, (4, 0), (2, 0), (0, 0)),
    55: (5, 11, 5, 4, 0, 0, (4, 0), (1, 0), (1, 0)),
    57: (3, 19, 5, 6, 0, 2, (4, 2), (2, 1), (0, 0)),
    65: (5, 13, 5, 8, 4, 0, (4, 4), (2, 1), (1, 0)),
    69: (3, 23, 7, 4, 0, 0, (4, 0), (2, 0), (0, 0)),
    77: (7, 11, 7, 4, 0, 0, (4, 0), (1, 0), (1, 0)),
    85: (5, 17, 7, 8, 4, 0, (4, 4), (2, 1), (1, 0)),
    87: (3, 29, 9, 4, 0, 0, (4, 0), (2, 0), (0, 0)),
    91: (7, 13, 7, 8, 0, 4, (4, 4), (2, 1), (1, 0)),
    93: (3, 31, 9, 6, 0, 2, (4, 2), (2, 2), (0, 0)),
    95: (5, 19, 9, 4, 0, 0, (4, 0), (2, 0), (1, 0)),
    111: (3, 37, 11, 6, 0, 2, (4, 2), (2, 2), (0, 0)),
    115: (5, 23, 11, 4, 0, 0, (4, 0), (2, 0), (1, 0)),
    119: (7, 17, 11, 4, 0, 0, (4, 0), (2, 0), (1, 0)),
    123: (3, 41, 13, 4, 0, 0, (4, 0), (2, 0), (0, 0)),
    129: (3, 43, 13, 6, 0, 2, (4, 2), (2, 2), (0, 0)),
    133: (7, 19, 11, 8, 0, 4, (4, 4), (2, 1), (1, 0)),
    141: (3, 47, 15, 4, 0, 0, (4, 0), (2, 0), (0, 0)),
    143: (11, 13, 13, 4, 0, 0, (4, 0), (2, 0), (1, 0)),
    145: (5, 29, 13, 8, 4, 0, (4, 4), (2, 2), (1, 0)),
    155: (5, 31, 15, 4, 0, 0, (4, 0), (2, 0), (1, 0)),
    159: (3, 53, 17, 4, 0, 0, (4, 0), (2, 0), (0, 0)),
}

#: Exact cusp-leakage description at N = 143, in the coordinates
#: c_d(F) = a_1(F | W_d) with W_d from PARI's mfatkininit.
#: The vanishing support and the ranks are normalization-invariant; the
#: coefficient 1 in c_1 = c_11 is pinned to that normalization.
RELATIONS_143 = {
    "T_l (l = 2, 3, 5, 7)": {"rank": 4, "zero": (), "extra": None},
    "U_11": {"rank": 2, "zero": (11, 143), "extra": None},
    "U_13": {"rank": 1, "zero": (13, 143), "extra": "c_1 = c_11"},
}

#: The two separated theorem targets left open by the arc.
THEOREM_TARGETS = {
    "vanishing_support": (
        "im(cusp leakage of U_p) is contained in span{ e_d : d | N, p does "
        "not divide d }. Verified at all 22 discovery levels. For N = pq "
        "this forces r_cusp(U_p) <= 2 with no census. Governs WHICH "
        "channels can exist."
    ),
    "surviving_pair_equality": (
        "c_p . U_q = c_1 . U_q on H^0(2K), for the larger level prime. "
        "Verified with constant exactly 1 at all 12 discovery levels with "
        "p >= 5. Governs WHY the surviving pair collapses to a line."
    ),
}

#: Pre-registered hypotheses and their status.
HYPOTHESES = {
    "Q1": ("every good prime gives (4, nu_2 + nu_3)", "HOLDS 66/66 in discovery"),
    "Q2": ("r_cusp(U_q) = 1 for the larger level prime", "FALSIFIED: (0,0) at all ten N = 3q levels"),
    "Q3": ("smaller level prime reaches the elliptic quotient", "HOLDS 9/9 in discovery"),
    "Q4": ("r_ell <= 1 for level primes", "FALSIFIED: r_ell = 2 at 93, 111, 129, 145"),
    "H1": ("N = 3q, mu_2 surjective => U_q H^0(2K) inside H^0(2K)",
           "PRE-REGISTERED, holdout 160 < N <= 400 not completed"),
}

_DATA = Path(__file__).resolve().parent / "_data"


def data_path(name: str) -> Path:
    p = _DATA / name
    if not p.exists():
        raise FileNotFoundError(f"{name} missing from mtft/boundary/_data")
    return p


def _read(name):
    """Return (rows, scale) for a shipped labelled-CSV matrix."""
    scale, rows, header = 1, [], False
    for line in data_path(name).read_text().splitlines():
        s = line.strip()
        if s.split()[:2] == ["#", "SCALE"]:
            scale = int(s.split()[2])
        if not s or s.startswith("#"):
            continue
        parts = s.split(",")
        if not header:
            header = True
            continue
        rows.append([int(x) for x in parts[1:]])
    return rows, scale


def s4_qexpansions():
    """Integral basis of S_4(Gamma_0(143)), 141 x 40, rows = coeff of q^n."""
    return _read("s4_143_qexpansions.txt")[0]


def cusp_functionals():
    """The four c_d as exact Fractions: 4 x 40, rows d = 1, 11, 13, 143."""
    rows, scale = _read("s4_143_cusp_functionals.txt")
    return [[Fraction(v, scale) for v in r] for r in rows]


def operator(name: str):
    """T_2, U_11 or U_13 on the S_4 basis, exact Fractions, 40 x 40."""
    if name not in ("T2", "U11", "U13"):
        raise ValueError("operator must be one of T2, U11, U13")
    rows, scale = _read(f"s4_143_op_{name}.txt")
    return [[Fraction(v, scale) for v in r] for r in rows]


def _kron(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def nu2(N: int) -> int:
    """Number of order-2 elliptic points of Gamma_0(N), N odd squarefree."""
    out = 1
    for p in _primes_dividing(N):
        out *= 1 + _kron(-1, p)
    return out


def nu3(N: int) -> int:
    """Number of order-3 elliptic points of Gamma_0(N), N odd squarefree."""
    out = 1
    for p in _primes_dividing(N):
        out *= 1 if p == 3 else 1 + _kron(-3, p)
    return out


def _primes_dividing(N):
    out, n, d = [], N, 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def boundary_dim(N: int) -> int:
    """dim B_N = nu_infinity + nu_2 + nu_3, with nu_infinity = 4 for N = pq."""
    return 4 + nu2(N) + nu3(N)
