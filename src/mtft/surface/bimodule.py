"""mtft.surface.bimodule — doubled-space (real spectral triple) census on H_1(X0(N), R).

Instrument for AF-08/AF-09-type questions.  Given an operator alphabet A on V
(matrices in a common frame, with an optional positive metric so that adjoints
are transposes after orthonormalisation) and a twist W, build H = V ⊕ V with

    left action   π(a) = diag(a, W a W⁻¹)
    real structure J_F = swap,   opposite action b° = J_F π(b)ᵀ J_F
    grading Γ = diag(I, −I),     Dirac D_M = [[0, M], [Mᵀ, 0]]

and measure, with ABSOLUTE scales (per unit generator / unit M):

    order_zero        max ‖[π(a), b°]‖ / (‖a‖‖b‖)
    first_order_space {M : [[D_M, π(a)], b°] = 0 for generators a, b}
                      (generators suffice given order zero)
    one_form_census   sizes of [D_M, π(a)] for M in that space; a "one-form
                      dimension" counts only forms above an absolute threshold
    sector_support    fraction of one-form norm inside a projector

Exact facts recorded at N = 143 (frozen data): the untwisted doubling fails
order zero because U13 is not normal; twisting by W13 or W143 restores it
through U13* = W13 U13 W13; every configuration has identically vanishing
one-forms (involutive twists: [D, a] ∝ λ_{σ²(i)} − λ_i = 0).  A closure
dimension or a one-form count that a random alphabet also produces is not
evidence (falsifier F8); run ``random_alphabet`` as a control.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

import numpy as np


def orthonormal_frame(metric: np.ndarray):
    """(R, R⁻¹) with R = G^{1/2}; conjugating by R makes the metric the identity."""
    w, U = np.linalg.eigh((metric + metric.T) / 2)
    if w[0] <= 0:
        raise ValueError("metric must be positive definite")
    R = U @ np.diag(np.sqrt(w)) @ U.T
    Ri = U @ np.diag(1 / np.sqrt(w)) @ U.T
    return R, Ri


def sector_dimensions(W1: np.ndarray, W2: np.ndarray) -> Dict[str, int]:
    """Joint (±,±) eigen-multiplicities of two commuting involutions."""
    n = W1.shape[0]
    I = np.eye(n)
    out = {}
    for a in (1, -1):
        for b in (1, -1):
            out[f"({'+' if a > 0 else '-'},{'+' if b > 0 else '-'})"] = int(round(np.trace((I + a * W1) @ (I + b * W2) / 4)))
    return out


def tensor_sector_dimensions(W1: np.ndarray, W2: np.ndarray) -> Dict[str, int]:
    """Sectors of the synchronised action W⊗W on V⊗V (Hom-space dimensions)."""
    return sector_dimensions(np.kron(W1, W1), np.kron(W2, W2))


def adjoint_identity(U: np.ndarray, W: np.ndarray, metric: Optional[np.ndarray] = None) -> float:
    """Relative residual of U* = W U W with U* the metric adjoint (transpose if metric is None)."""
    Ustar = U.T if metric is None else np.linalg.solve(metric, U.T @ metric)
    return float(np.linalg.norm(Ustar - W @ U @ W) / np.linalg.norm(U))


@dataclass
class Doubling:
    alphabet: Dict[str, np.ndarray]
    twist: np.ndarray
    n: int

    @classmethod
    def build(cls, alphabet: Dict[str, np.ndarray], metric: Optional[np.ndarray] = None,
              twist: Optional[np.ndarray] = None) -> "Doubling":
        names = list(alphabet)
        mats = [np.asarray(alphabet[k], dtype=float) for k in names]
        n = mats[0].shape[0]
        W = np.eye(n) if twist is None else np.asarray(twist, dtype=float)
        if metric is not None:
            R, Ri = orthonormal_frame(np.asarray(metric, dtype=float))
            mats = [R @ A @ Ri for A in mats]
            W = R @ W @ Ri
        return cls(dict(zip(names, mats)), W, n)

    # representations on V ⊕ V
    def left(self, a: np.ndarray) -> np.ndarray:
        Z = np.zeros((self.n, self.n))
        return np.block([[a, Z], [Z, self.twist @ a @ self.twist.T]])

    def opposite(self, b: np.ndarray) -> np.ndarray:
        Z = np.zeros((self.n, self.n))
        return np.block([[self.twist @ b.T @ self.twist.T, Z], [Z, b.T]])

    def dirac(self, M: np.ndarray) -> np.ndarray:
        Z = np.zeros((self.n, self.n))
        return np.block([[Z, M], [M.T, Z]])

    def order_zero(self) -> Dict:
        worst = 0.0
        for a in self.alphabet.values():
            for b in self.alphabet.values():
                worst = max(worst, np.linalg.norm(self.left(a) @ self.opposite(b) - self.opposite(b) @ self.left(a))
                            / (np.linalg.norm(a) * np.linalg.norm(b)))
        return {"max_rel_residual": float(worst), "status": "PASS" if worst < 1e-10 else "FAIL"}

    def first_order_space(self, tol: float = 1e-9) -> np.ndarray:
        """Columns: M (as n²-vectors) with [[D_M, π(a)], b°] = 0; residual < tol per unit M."""
        n = self.n
        rows = []
        gens = list(self.alphabet.values())
        norms = [np.linalg.norm(a) for a in gens]
        for i in range(n):
            for j in range(n):
                M = np.zeros((n, n))
                M[i, j] = 1
                D = self.dirac(M)
                parts = []
                for a, na in zip(gens, norms):
                    X = D @ self.left(a) - self.left(a) @ D
                    for b, nb in zip(gens, norms):
                        parts.append(((X @ self.opposite(b) - self.opposite(b) @ X) / (na * nb)).reshape(-1))
                rows.append(np.concatenate(parts))
        L = np.array(rows)
        g, V = np.linalg.eigh(L @ L.T)
        return V[:, g < tol * tol]

    def one_form_census(self, F: np.ndarray, threshold: float = 1e-6) -> Dict:
        n = self.n
        forms, maxsize = [], 0.0
        for v in F.T:
            M = v.reshape(n, n)
            D = self.dirac(M)
            for a in self.alphabet.values():
                f = (D @ self.left(a) - self.left(a) @ D) / np.linalg.norm(a)
                maxsize = max(maxsize, float(np.linalg.norm(f)))
                forms.append(f.reshape(-1))
        forms = np.array(forms) if forms else np.zeros((0, 4 * n * n))
        sv = np.linalg.svd(forms, compute_uv=False) if len(forms) else np.zeros(1)
        dim = int(np.sum(sv > threshold))
        return {"first_order_dimension": int(F.shape[1]), "max_one_form_size": maxsize,
                "one_form_dimension": dim, "forms": forms, "threshold": threshold}

    def sector_support(self, forms: np.ndarray, projector: np.ndarray, threshold: float = 1e-6) -> Optional[float]:
        n = self.n
        Z = np.zeros((n, n))
        Q = np.block([[projector, Z], [Z, self.twist @ projector @ self.twist.T]])
        big = [f.reshape(2 * n, 2 * n) for f in forms if np.linalg.norm(f) > threshold]
        if not big:
            return None
        return float(sum(np.linalg.norm(Q @ f @ Q) for f in big) / sum(np.linalg.norm(f) for f in big))


def census(alphabet: Dict[str, np.ndarray], metric: Optional[np.ndarray], twists: Dict[str, Optional[np.ndarray]],
           threshold: float = 1e-6) -> List[Dict]:
    """One row per twist: order-zero, first-order dimension, one-form size/dimension."""
    rows = []
    for label, W in twists.items():
        dbl = Doubling.build(alphabet, metric, W)
        oz = dbl.order_zero()
        F = dbl.first_order_space()
        of = dbl.one_form_census(F, threshold)
        rows.append({"twist": label, "order_zero": oz["max_rel_residual"], "order_zero_status": oz["status"],
                     "first_order_dimension": of["first_order_dimension"],
                     "max_one_form_size": of["max_one_form_size"], "one_form_dimension": of["one_form_dimension"]})
    return rows


def random_alphabet(n: int, k: int = 2, symmetric: bool = True, seed: int = 0) -> Dict[str, np.ndarray]:
    """Control alphabet (falsifier F8): generic operators, symmetric if requested."""
    rng = np.random.default_rng(seed)
    out = {}
    for i in range(k):
        A = rng.standard_normal((n, n))
        out[f"R{i}"] = (A + A.T) / 2 if symmetric else A
    return out


def x0143_census(threshold: float = 1e-6) -> Dict:
    """The AF-09 experiment on the frozen N=143 data: A = R[T2, T3, U13], twists 1, W11, W13, W143."""
    from .frozen import x0143
    d = x0143()
    Jint = d["intersection_cycles"].astype(float)
    G = Jint @ d["J_true"]
    G = (G + G.T) / 2
    if np.linalg.eigvalsh(G)[0] < 0:
        G = -G
    alphabet = {k: d[k].astype(float) for k in ("T2", "T3", "U13")}
    twists = {"untwisted": None, "W11": d["W11"].astype(float), "W13": d["W13"].astype(float), "W143": d["W143"].astype(float)}
    rows = census(alphabet, G, twists, threshold)
    return {"rows": rows,
            "sectors_V": sector_dimensions(d["W11"].astype(float), d["W13"].astype(float)),
            "sectors_VxV": tensor_sector_dimensions(d["W11"].astype(float), d["W13"].astype(float)),
            "AL_adjoint_identity_U13": adjoint_identity(d["U13"].astype(float), d["W13"].astype(float), G),
            "AL_adjoint_identity_U11": adjoint_identity(d["U11"].astype(float), d["W11"].astype(float), G)}


# ------------------------------------------------ v0.27.1: full real-spectral-triple gate set
def cyclic_permutation(dim: int, k: int) -> np.ndarray:
    """Permutation matrix on (C^dim)^{⊗k} cycling the tensor factors (slot i -> slot i+1)."""
    import itertools
    idx = list(itertools.product(range(dim), repeat=k))
    pos = {s: i for i, s in enumerate(idx)}
    P = np.zeros((dim ** k, dim ** k))
    for i, s in enumerate(idx):
        P[pos[(s[-1],) + s[:-1]], i] = 1.0
    return P


def kron_power(A: np.ndarray, k: int) -> np.ndarray:
    out = np.array(A, dtype=float)
    for _ in range(k - 1):
        out = np.kron(out, A)
    return out


def tensor_alphabet(alphabet: Dict[str, np.ndarray], k: int, max_dim: int = 4096) -> Dict[str, np.ndarray]:
    """Slot-wise generators a⊗I⊗…, I⊗a⊗I, … of A^{⊗k} (dense; refuses beyond max_dim)."""
    n = next(iter(alphabet.values())).shape[0]
    if n ** k > max_dim:
        raise RuntimeError(f"dense tensor alphabet dimension {n ** k} exceeds max_dim={max_dim}; use character orbits")
    I = np.eye(n)
    out = {}
    for name, a in alphabet.items():
        for slot in range(k):
            mats = [np.asarray(a, float) if s == slot else I for s in range(k)]
            M = mats[0]
            for m in mats[1:]:
                M = np.kron(M, m)
            out[f"{name}@{slot}"] = M
    return out


def symmetric_dirac_block(W: np.ndarray, h: np.ndarray) -> np.ndarray:
    """K_h = W h + h W^{-1}: self-adjoint (h self-adjoint, W orthogonal) and real (J D = D J)."""
    return W @ h + h @ np.linalg.inv(W)


class RealTriple(Doubling):
    """Doubling with the complete finite real-spectral-triple gate set."""

    def axiom_gates(self, D: np.ndarray, tol: float = 1e-10) -> Dict:
        n = self.n
        Z = np.zeros((n, n))
        I = np.eye(n)
        Jr = np.block([[Z, I], [I, Z]])            # real structure: swap (composed with conjugation over C)
        G = np.block([[I, Z], [Z, -I]])
        nrm = lambda A: float(np.linalg.norm(A))
        scale = max(nrm(D), 1e-300)
        gates = {
            "selfadjoint": nrm(D - D.T) / scale,
            "reality_JD_equals_DJ": nrm(Jr @ D @ Jr - D) / scale,
            "grading_D_odd": nrm(G @ D + D @ G) / scale,
            "grading_algebra_even": max(nrm(G @ self.left(a) - self.left(a) @ G) / nrm(a) for a in self.alphabet.values()),
            "J_Gamma_anticommute": nrm(Jr @ G + G @ Jr) / nrm(G),
            "order_zero": self.order_zero()["max_rel_residual"],
        }
        fo = 0.0
        for a in self.alphabet.values():
            X = D @ self.left(a) - self.left(a) @ D
            for b in self.alphabet.values():
                fo = max(fo, nrm(X @ self.opposite(b) - self.opposite(b) @ X) / (scale * nrm(a) * nrm(b)))
        gates["first_order"] = fo
        # first-order is scaled by ||D|| ||a|| ||b|| (absolute), never by the one-form norm
        gates["max_one_form_size"] = max(nrm(D @ self.left(a) - self.left(a) @ D) / nrm(a) for a in self.alphabet.values())
        gates["status"] = "PASS" if all(v < tol for k, v in gates.items() if k not in ("status", "max_one_form_size")) else "FAIL"
        return gates

    def pairing(self, basis: Optional[Dict[str, np.ndarray]] = None) -> Dict:
        """Q_ij = Tr(Γ π(e_i) e_j°) on a basis of the algebra (defaults to alphabet ∪ {I}); the finite
        Poincaré-duality gate is nondegeneracy on the FULL algebra basis (pass primitive idempotents)."""
        n = self.n
        Z = np.zeros((n, n))
        I = np.eye(n)
        G = np.block([[I, Z], [Z, -I]])
        B = dict(basis) if basis else {**self.alphabet, "I": I}
        names = list(B)
        Q = np.array([[np.trace(G @ self.left(B[i]) @ self.opposite(B[j])) for j in names] for i in names])
        rank = int(np.linalg.matrix_rank(Q, tol=1e-9 * max(1.0, np.abs(Q).max())))
        return {"Q": Q, "basis": names, "rank": rank, "nullity": len(names) - rank,
                "unit_in_radical": bool("I" in B and np.abs(Q[names.index("I")]).max() < 1e-9),
                "nondegenerate_on_supplied_basis": rank == len(names)}

    def differential_rank(self, D: np.ndarray) -> int:
        forms = np.array([(D @ self.left(a) - self.left(a) @ D).reshape(-1) for a in self.alphabet.values()])
        return int(np.linalg.matrix_rank(forms, tol=1e-9 * max(1.0, np.abs(forms).max())))


def two_state_cyclic_control(k: int = 3) -> Dict:
    """Astra's exact synthetic control (two-state factor algebra, cyclic twist), all gates."""
    n = 2 ** k
    P = cyclic_permutation(2, k)
    E = np.eye(n)
    alphabet = {f"e{i}": np.diag(E[i]) for i in range(n)}
    T = RealTriple.build(alphabet, None, P)
    h = np.diag(np.arange(2, n + 2, dtype=float))
    rows = {}
    for label, M in (("[[0,Wh],[Wh^T,0]]", P @ h), ("[[0,Wh],[Wh,0]]", None), ("K_h=Wh+hW^-1", symmetric_dirac_block(P, h)), ("K=W+W^-1", P + P.T)):
        if M is None:
            Z = np.zeros((n, n)); D = np.block([[Z, P @ h], [P @ h, Z]])
        else:
            Z = np.zeros((n, n)); D = np.block([[Z, M], [M.T, Z]])
        g = T.axiom_gates(D)
        rows[label] = {k_: (round(v, 12) if isinstance(v, float) else v) for k_, v in g.items()}
        rows[label]["rank_d"] = T.differential_rank(D)
    rows["pairing"] = {k_: v for k_, v in T.pairing().items() if k_ != "Q"}
    return rows
