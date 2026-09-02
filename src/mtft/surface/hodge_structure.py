"""mtft.surface.hodge_structure — the true Hodge structure on H_1(X0(N), R) in the cycle basis.

J_true = Q^-1 S Q with Q = [Re P; Im P] the real period matrix of the
rational cuspidal basis over the tree/cotree cycles and S = multiplication
by i.  Non-tautological gates: the two Riemann bilinear relations against the
recovered intersection form (symplectic; Jint J symmetric positive definite)
and commutation with the exact integer Hecke/AL transport (unrelated
computation).  Overall sign fixed by positivity.

E2 on an elliptic block: for a rational newform block (kernel of T_p - a_p,
rank 2), the saturated lattice with J restricted is an elliptic curve; its
j-invariant must lie in the j-list of the isogeny class (PARI ellisomat).

The Whitney family of :mod:`hodge` must converge to J_true (conformal
invariance of the 1-form star + Voevodsky-Shabat); the Siegel distance is
reported as DIAGNOSTIC_CONVERGING until a rate is certified.
"""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from math import gcd
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np
import sympy as sp

from ..gprun import find_gp
from .hodge import siegel_distance


def _rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


@dataclass
class HodgeStructure:
    N: int
    J: np.ndarray                    # cycle basis, sign fixed by Riemann positivity
    G: np.ndarray                    # Jint J, symmetric positive definite
    gates: Dict[str, float]
    positivity_min_eigenvalue: float
    period_singular_values: Tuple[float, float]

    @property
    def claim_class(self) -> str:
        return "CERTIFIED_NUMERICAL"


def from_periods(N: int, Q: np.ndarray, Jint: np.ndarray, ops: Dict[str, np.ndarray]) -> HodgeStructure:
    n = Q.shape[0]
    g = n // 2
    S = np.block([[np.zeros((g, g)), -np.eye(g)], [np.eye(g), np.zeros((g, g))]])
    sv = np.linalg.svd(Q, compute_uv=False)
    J = np.linalg.solve(Q, S @ Q)
    G = Jint.astype(float) @ J
    Gs = (G + G.T) / 2
    ev = np.linalg.eigvalsh(Gs)
    if ev[0] < 0:
        J, G, Gs, ev = -J, -G, -Gs, -ev[::-1]
    gates = {
        "J_squared_plus_identity_rel": _rel(J @ J, -np.eye(n)),
        "riemann_I_symplectic_rel": _rel(J.T @ Jint @ J, Jint),
        "riemann_II_symmetry_rel": _rel(G, G.T),
    }
    for k, A in ops.items():
        gates[f"commutator_{k}_rel"] = _rel(J @ A, A @ J)
    return HodgeStructure(N, J, Gs, gates, float(ev[0]), (float(sv[-1]), float(sv[0])))


def gates_pass(hs: HodgeStructure, tol: float = 1e-10) -> bool:
    return all(v < tol for v in hs.gates.values()) and hs.positivity_min_eigenvalue > 0


def family_distances(hs: HodgeStructure, family: List[Dict]) -> List[Dict]:
    out = []
    for lv in family:
        d, dmax = siegel_distance(hs.G, lv["G"])
        out.append({"level": lv["level"], "siegel_distance": d, "siegel_max_log": dmax,
                    "rel_frobenius": _rel(lv["J"], hs.J),
                    "wedge_sv_min_max": (float(lv["wedge_singular_values"][0]), float(lv["wedge_singular_values"][-1]))})
    return out


# ------------------------------------------------- elliptic block j-invariant
def _integer_kernel_basis(M: np.ndarray) -> np.ndarray:
    """Saturated Z-basis of ker(M) in Z^n (columns), via sympy nullspace + Smith."""
    ns = sp.Matrix(M.astype(int).tolist()).nullspace()
    cols = []
    for v in ns:
        den = sp.ilcm(*[sp.fraction(sp.nsimplify(x))[1] for x in v]) if len(v) else 1
        cols.append([int(x * den) for x in v])
    K = sp.Matrix(cols).T                      # n x r
    r = K.shape[1]
    # saturate: content of r x r minors
    from itertools import combinations
    content = 0
    for rows in combinations(range(K.shape[0]), r):
        content = gcd(content, int(K.extract(list(rows), list(range(r))).det()))
        if content == 1:
            break
    if content != 1:
        from ..integral_lattice import saturate
        primes = [p for p in range(2, content + 1) if content % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]
        Ksat, _ = saturate(np.array(K.tolist(), dtype=object), primes)
        K = sp.Matrix(np.array(Ksat, dtype=object).tolist())
    return np.array(K.tolist(), dtype=float)


def elliptic_block_tau(hs: HodgeStructure, T: np.ndarray, eigenvalue: int) -> complex:
    """tau of the elliptic curve (saturated rank-2 lattice ker(T - a) with J restricted)."""
    K = _integer_kernel_basis(T - eigenvalue * np.eye(T.shape[0], dtype=np.int64))
    if K.shape[1] != 2:
        raise ValueError(f"kernel rank {K.shape[1]} != 2")
    v1, v2 = K[:, 0], K[:, 1]
    coef, *_ = np.linalg.lstsq(K, hs.J @ v1, rcond=None)
    a, b = coef
    if b < 0:                                   # orient so that Im tau > 0
        coef, *_ = np.linalg.lstsq(K[:, ::-1], hs.J @ v2, rcond=None)
        a, b = coef
    return complex((1j - a) / b)


def klein_j(tau: complex, dps: int = 30) -> complex:
    with mp.workdps(dps):
        return complex(1728 * mp.kleinj(mp.mpc(tau.real, tau.imag)))


def isogeny_class_j_list(ainvs: List[int], timeout: int = 60) -> List[complex]:
    gp = find_gp()
    if gp is None:
        raise RuntimeError("PARI/GP not found")
    import tempfile
    from pathlib import Path
    script = f"E=ellinit({list(ainvs)}); L=ellisomat(E); print(\"J=\", vector(#L[1], i, ellinit(L[1][i][1]).j)); print(\"COND=\", ellglobalred(E)[1]);\n"
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "isog.gp"
        path.write_text(script)
        proc = subprocess.run([gp, "-q", str(path)], capture_output=True, text=True, timeout=timeout)
    m = re.search(r"J=\s*\[(.*?)\]", proc.stdout, re.S)
    if not m:
        raise RuntimeError(proc.stdout + proc.stderr)
    return [complex(float(sp.Rational(x.strip()))) for x in m.group(1).split(",")]


def elliptic_block_check(hs: HodgeStructure, T: np.ndarray, eigenvalue: int, ainvs: List[int],
                         rtol: float = 1e-8) -> Dict:
    tau = elliptic_block_tau(hs, T, eigenvalue)
    j = klein_j(tau)
    jl = isogeny_class_j_list(ainvs)
    best = min(jl, key=lambda z: abs(z - j))
    rel = abs(best - j) / max(abs(best), 1.0)
    return {"tau": tau, "j": j, "isogeny_class_j": jl, "closest_rel_error": float(rel),
            "status": "PASS" if rel < rtol else "FAIL", "class": "CERTIFIED_NUMERICAL"}
