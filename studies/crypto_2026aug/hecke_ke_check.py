"""
hecke_ke_check — is Hecke-CDH (eprint 2025/1681, Def. 11) actually hard?

The paper's Sec 3.3 dismisses the eigenbasis attack with:

  "For general forms that are not eigenforms or do not lie in a common
   Hecke eigenspace, the action of H cannot be expressed as a simple
   product of eigenvalues, this attack is inapplicable."

That sentence is the thing to test.  If {f_j} is a simultaneous
eigenbasis for the whole Hecke algebra, then EVERY f decomposes as
f = sum_j c_j f_j and

        H(f) = sum_j ( prod_i lambda_j(N_i) ) c_j f_j

is diagonal in that basis for every f, eigenform or not.  If so, an
eavesdropper never needs to recover N_1..N_s at all: she reads the
eigenvalue tuple off the public value and multiplies componentwise.

Run: build H_A, H_B as products of Hecke operators exactly as in the
protocol (Fig. 1), form the true shared secret H_A(H_B(f)), then have
Eve reconstruct it from public data alone.  Compare.

Hecke matrices come from the mtft engine's own Manin-symbol machinery
(level 143, weight 2, 29-dim) -- integer matrices, no floating point
in the protocol itself.
"""
from __future__ import annotations
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import sympy as sp
from x0143_particle_box import P1, tessellation, ModularSymbols

N = 143
random.seed(20260805)


def hecke(ms, p):
    """Exact integer matrix of T_p on the 29-dim space."""
    M = ms.hecke_on_quotient(p)
    return sp.Matrix(np.rint(np.array(M, dtype=float)).astype(int).tolist())


def main():
    p1 = P1(N)
    tessellation(N)
    ms = ModularSymbols(p1)
    d = ms.dim
    print(f"space dimension d = {d}   (paper's experiments use d ~ 4..10)")

    # Hecke operators at primes not dividing N.  3-square-free N_i means
    # T_{N_i} = T_pa T_pb T_pc, so the whole protocol is products of T_p.
    S_A = [2, 3, 5, 7, 17, 19]        # Alice's prime pool
    S_B = [23, 29, 31, 37, 41, 43]    # Bob's prime pool  (disjoint)
    T = {p: hecke(ms, p) for p in S_A + S_B}
    print(f"built T_p for p in {S_A + S_B}")

    # commutativity sanity
    allp = S_A + S_B
    comm = all((T[a]*T[b] - T[b]*T[a]).is_zero_matrix
               for a in allp for b in allp)
    print(f"all T_p commute: {comm}")

    # ---- the protocol, run honestly -----------------------------------
    # Alice: s 3-square-free integers over S_A
    s = 3
    A_ints = [tuple(random.sample(S_A, 3)) for _ in range(s)]
    B_ints = [tuple(random.sample(S_B, 3)) for _ in range(s)]
    print(f"\nAlice's N_i (3-square-free): "
          f"{[a*b*c for a, b, c in A_ints]}")
    print(f"Bob's   M_i (3-square-free): "
          f"{[a*b*c for a, b, c in B_ints]}")

    HA = sp.eye(d)
    for (a, b, c) in A_ints:
        HA = HA * T[a] * T[b] * T[c]
    HB = sp.eye(d)
    for (a, b, c) in B_ints:
        HB = HB * T[a] * T[b] * T[c]

    # public form f: a generic (NON-eigenform) element, integer coords
    f = sp.Matrix([random.randint(1, 9) for _ in range(d)])

    pub_A = HA * f            # Alice sends
    pub_B = HB * f            # Bob sends
    shared = HA * (HB * f)    # the true shared secret
    print(f"\nshared secret (first 5 coords): {list(shared)[:5]}")

    # ---- Eve, using ONLY f, pub_A, pub_B -------------------------------
    # Step 1: simultaneous eigen-decomposition of the Hecke algebra.
    # Public: depends only on N and k, not on any secret.
    Tgen = T[2] + 2*T[3] + 3*T[5]          # generic element, splits the algebra
    Tn = np.array(Tgen.tolist(), dtype=float)
    w, P = np.linalg.eig(Tn)
    Pinv = np.linalg.inv(P)

    fv = np.array(f.tolist(), dtype=float).ravel()
    av = np.array(pub_A.tolist(), dtype=float).ravel()
    bv = np.array(pub_B.tolist(), dtype=float).ravel()

    fh, ah, bh = Pinv @ fv, Pinv @ av, Pinv @ bv

    # Step 2: on each joint eigenspace H acts as a scalar.  Read it off.
    lamA = np.where(np.abs(fh) > 1e-9, ah / np.where(np.abs(fh) > 1e-9, fh, 1), 0)
    lamB = np.where(np.abs(fh) > 1e-9, bh / np.where(np.abs(fh) > 1e-9, fh, 1), 0)

    # Step 3: the shared secret is the componentwise product.
    eve_hat = lamA * lamB * fh
    eve = np.real(P @ eve_hat)

    truth = np.array(shared.tolist(), dtype=float).ravel()
    rel = np.linalg.norm(eve - truth) / np.linalg.norm(truth)

    print(f"Eve's reconstruction (first 5): "
          f"{[round(v, 3) for v in eve[:5]]}")
    print(f"\nrelative error |Eve - truth| / |truth| = {rel:.3e}")
    print(f"EVE RECOVERS THE SHARED SECRET: {rel < 1e-6}")
    print("\nEve never recovered N_1..N_s.  She did one eigendecomposition")
    print(f"of a {d}x{d} matrix and {d} divisions.  The NP-hardness of")
    print("recovering the factorization is irrelevant to CDH.")


if __name__ == "__main__":
    main()
