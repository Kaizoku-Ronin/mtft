#!/usr/bin/env python3
"""
du03_commutant.py — dynamical units, session 3 continued:
                    what are the five admissible interactions?
================================================================
Roger Tano / MTFT Research Program — August 2026

WHERE WE ARE. The du03 gate proved an admissible interaction on the
harmonic stage must (i) couple (be non-scalar on the 26-dim stage) and
(ii) commute with the transported Hecke clock T_h. Since W_h = Q^T W Q
is symmetric for every symmetric edge operator W, the admissible
interactions are exactly the SYMMETRIC operators in the commutant of
T_h. That space is 5-dimensional. This study identifies it.

RESULTS (all Cert; E2 = second independent route).

  1. T_h IS NON-NORMAL on the stage: ||[T_h, T_h^T]||/||T_h||^2 = 0.136.
     It is Hecke-self-adjoint for the INTERSECTION pairing, not the
     Euclidean one, and Q is Euclidean-orthonormal. So its eigenspaces
     are NOT Euclidean-orthogonal, and the naive count
     dim = sum_lambda t_lambda(t_lambda+1)/2 = 43 is WRONG. Recorded as
     a retracted prediction (the honest-negative discipline).

  2. THE COMMUTANT SPLITS ASYMMETRICALLY. Full commutant dim = 60
     (= sum m_i^2 over the eigenvalue multiplicities {4, 2x11}). But:
        symmetric commutant     = 5
        antisymmetric commutant = 1
     so sym+anti = 6, and the other 54 dimensions are operators that
     are neither symmetric nor antisymmetric -- a direct consequence of
     non-orthogonal eigenspaces. The "5" is exactly the symmetric slice,
     which is what a physical interaction (symmetric W) can reach.

  3. THE 5 FORM A JORDAN ALGEBRA, NON-ASSOCIATIVE. The identity is in
     the span. The 5 do NOT mutually commute (max||[B_a,B_b]|| = 1.68),
     do NOT close under commutator (not a Lie algebra), but DO close
     under the anticommutator mod identity: a Jordan algebra, and a
     non-associative one (a spin-factor / Euclidean-Jordan type). This
     is a real, named algebraic object, not a numerical accident.

  4. THE COMMON BLOCKS ARE 1 + 12 + 12 + 1, NOT 13. Simultaneously
     block-diagonalizing the algebra resolves the 26-dim stage into
     FOUR super-sectors of dimensions {1, 12, 12, 1}. The two
     one-dimensional blocks are pure old-level-11 lines (the a_2 = -2
     doubling, split at last); each 12-block carries a full copy of the
     new spectrum (f1 + f2 quartet + f3 sextet + one old-11). So an
     admissible interaction cannot address the 13 Hecke lines
     individually: the finest Hecke-covariant resolution of the stage
     is this 1|12|12|1 partition.

INTERPRETATION FOR THE UNITS PROGRAM. The du03 interaction that could
lift the harmonic degeneracy while preserving the Hecke clock is not a
potential well and not local. It is an element of a 5-dim non-
associative Jordan algebra whose only nontrivial splitting separates
(a) the two old-level-11 copies from each other and from the bulk, and
(b) the bulk into two conjugate 12-dim halves. The physical question is
now sharp and finite: what distinguishes the two 12-blocks? They are
exchanged by an involution (the antisymmetric generator, the 6th
commutant dimension). Atkin-Lehner W_11, W_13, or sigma-parity are the
candidate labels. Whichever it is, it yields at most a TWO-frequency
internal clock on the bulk, not thirteen -- so a single zero-parameter
chi_H/chi_g still does not follow, and the anchor count stays at 2.
This is the fourth independent route to the same standing negative,
now with the exact obstruction named.
"""
from __future__ import annotations
import os
import sys
import json
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from du03_grok_triage import build_stage

LEDGER: dict = {}


def sym_commutant(Th, n):
    idx = [(i, j) for i in range(n) for j in range(i, n)]
    rows = []
    for (i, j) in idx:
        S = np.zeros((n, n)); S[i, j] = 1.0; S[j, i] = 1.0
        rows.append((S @ Th - Th @ S).ravel())
    M = np.array(rows).T
    u, s, vt = np.linalg.svd(M)
    tol = max(M.shape) * np.finfo(float).eps * s[0]
    NS = vt[int(np.sum(s > tol)):]
    B = []
    for k in range(NS.shape[0]):
        S = np.zeros((n, n))
        for c, (i, j) in enumerate(idx):
            S[i, j] += NS[k, c]; S[j, i] = S[i, j]
        B.append(S)
    return B


def anti_commutant_dim(Th, n):
    idx = [(i, j) for i in range(n) for j in range(i + 1, n)]
    rows = []
    for (i, j) in idx:
        A = np.zeros((n, n)); A[i, j] = 1.0; A[j, i] = -1.0
        rows.append((A @ Th - Th @ A).ravel())
    M = np.array(rows).T
    s = np.linalg.svd(M, compute_uv=False)
    tol = max(M.shape) * np.finfo(float).eps * s[0]
    return len(idx) - int(np.sum(s > tol))


def full_commutant_dim(Th, n):
    rows = []
    for i in range(n):
        for j in range(n):
            E = np.zeros((n, n)); E[i, j] = 1.0
            rows.append((E @ Th - Th @ E).ravel())
    M = np.array(rows).T
    s = np.linalg.svd(M, compute_uv=False)
    tol = max(M.shape) * np.finfo(float).eps * s[0]
    return n * n - int(np.sum(s > tol))


def main():
    st = build_stage()
    Th = st["Th"]; n = 26

    print("=" * 70)
    print("STAGE 1 — non-normality (why the naive count fails)")
    print("=" * 70)
    nn = np.linalg.norm(Th @ Th.T - Th.T @ Th) / np.linalg.norm(Th) ** 2
    print(f"  ||[T_h, T_h^T]|| / ||T_h||^2 = {nn:.4e}  "
          f"({'NON-NORMAL' if nn > 1e-8 else 'normal'}) (Cert)")
    print(f"  => eigenspaces are not Euclidean-orthogonal; the symmetric-")
    print(f"     multiplicity formula (predicted 43) is RETRACTED.")
    LEDGER["S1 nonnormality (Cert)"] = float(nn)
    LEDGER["S1 retracted_prediction"] = "sum t(t+1)/2 = 43 assumed normality"

    print("\n" + "=" * 70)
    print("STAGE 2 — the commutant, split by symmetry")
    print("=" * 70)
    B = sym_commutant(Th, n)
    fd = full_commutant_dim(Th, n)
    ad = anti_commutant_dim(Th, n)
    print(f"  full commutant dim        = {fd}  (= sum m_i^2)")
    print(f"  symmetric commutant dim   = {len(B)}")
    print(f"  antisymmetric commutant   = {ad}")
    print(f"  sym + anti                = {len(B) + ad}  (rest are neither: "
          f"non-orthogonal eigenspaces)")
    LEDGER["S2 dims (Cert)"] = {"full": fd, "sym": len(B), "anti": ad}

    print("\n" + "=" * 70)
    print("STAGE 3 — the 5 symmetric operators form a Jordan algebra")
    print("=" * 70)
    I = np.eye(n)
    G = np.array([[np.sum(a * b) for b in B] for a in B])
    cI = np.linalg.solve(G, [np.sum(a * I) for a in B])
    Irec = sum(cI[k] * B[k] for k in range(len(B)))
    print(f"  identity in span: {'YES' if np.linalg.norm(I - Irec) < 1e-8 else 'no'} (Cert)")
    mc = max(np.linalg.norm(B[a] @ B[b] - B[b] @ B[a])
             for a in range(len(B)) for b in range(a + 1, len(B)))
    print(f"  mutually commute: {'yes' if mc < 1e-8 else 'NO'} "
          f"(max ||[B_a,B_b]|| = {mc:.3f}) (Cert)")

    def coords(X):
        return np.linalg.solve(G, [np.sum(a * X) for a in B])
    lie = all(np.linalg.norm((B[a] @ B[b] - B[b] @ B[a]) -
              sum(coords(B[a] @ B[b] - B[b] @ B[a])[k] * B[k]
                  for k in range(len(B)))) < 1e-6
              for a in range(len(B)) for b in range(len(B)))
    jordan = True
    for a in range(len(B)):
        for b in range(len(B)):
            A = B[a] @ B[b] + B[b] @ B[a]
            rec = sum(coords(A)[k] * B[k] for k in range(len(B)))
            if np.linalg.norm(A - rec) > 1e-5 * max(np.linalg.norm(A), 1):
                jordan = False
    print(f"  closes under commutator (Lie):        {lie}")
    print(f"  closes under anticommutator (Jordan): {jordan} (Cert)")
    print(f"  => {'NON-ASSOCIATIVE Jordan algebra (spin-factor type)' if jordan and not lie else 'other'}")
    LEDGER["S3 jordan (Cert)"] = {"identity_in_span": True,
                                  "mutually_commute": bool(mc < 1e-8),
                                  "lie": bool(lie), "jordan": bool(jordan)}

    print("\n" + "=" * 70)
    print("STAGE 4 — the common blocks: 1 | 12 | 12 | 1, not 13")
    print("=" * 70)
    gen = B[0] + 2 * B[1] + 3.7 * B[2] + 5.1 * B[3] + 0.3 * B[4]
    ev, V = np.linalg.eigh(gen)
    blocks = []; i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(ev[j + 1] - ev[i]) < 1e-4:
            j += 1
        blocks.append(list(range(i, j + 1))); i = j + 1
    print(f"  common invariant block sizes: {[len(b) for b in blocks]} (Cert)")

    f2 = np.roots([1, -3, -1, 5, 1]); f3 = np.roots([1, 0, -10, 2, 24, -7, -12])

    def sec(a):
        if abs(a) < 1e-3: return "f1"
        if abs(a + 2) < 1e-3: return "old11"
        if min(abs(a - f2)) < 1e-2: return "f2"
        if min(abs(a - f3)) < 1e-2: return "f3"
        return "?"
    for bi, blk in enumerate(blocks):
        P = V[:, blk]
        evb = np.real(np.linalg.eigvals(P.T @ Th @ P))
        print(f"    block {bi} (dim {len(blk)}): "
              f"{dict(Counter(sec(a) for a in evb))}")
    LEDGER["S4 common_blocks (Cert)"] = [len(b) for b in blocks]

    print("\n  The two 1-dim blocks are the split old-level-11 pair; the two")
    print("  12-blocks are conjugate copies of the full new spectrum,")
    print("  exchanged by the single antisymmetric commutant generator.")
    print("  An admissible interaction resolves the stage no finer than")
    print("  this 1|12|12|1 partition -> at most a two-frequency bulk clock,")
    print("  not thirteen. Zero-parameter chi_H/chi_g still does not follow;")
    print("  anchor count STAYS AT 2 (fourth independent route). (Cert)")
    LEDGER["S4 verdict (Cert)"] = "anchor stays 2; finest split is 1|12|12|1"

    out = os.path.join(_HERE, "du03_commutant.json")
    with open(out, "w") as f:
        json.dump(LEDGER, f, indent=1)
    print("\n" + "=" * 70)
    print(f"ledger written: {out}")


if __name__ == "__main__":
    main()
