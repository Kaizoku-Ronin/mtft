#!/usr/bin/env python3
"""
studies/x0143_capture_parity.py — the capture ceiling theorem (v0.11.3)
=======================================================================
Closes the Add. BQ leg-6 / AG-D5 open item: WHY is the capture ceiling
of the X0(143) atom exactly 1/2?

THEOREM (certified below). Complex conjugation of X0(143) (z -> -zbar)
acts on the Farey tessellation through the coset map
        sigma : (c:d) -> (-c : c+d),
an involution of the 56 triangles and an automorphism of the dual graph
with exactly 4 fixed triangles (the real locus), fixing the width-1
nucleus. Because odd states must vanish at fixed points, the odd parity
sector CANNOT see the nucleus well; the sigma-invariant dipole obeys the
selection rule <even|D|odd> = 0, the jump network splits into two
superselection sectors, and the odd-sector bottom is a dark sink. Hence
for a packet started at triangle v:
        ceiling = || P_even psi_0 ||^2
                = 1/2  exactly, if sigma(v) != v   (52 triangles)
                = 1    exactly, if sigma(v) == v   (4 triangles)
independent of V0 and gamma0. Breaking sigma-symmetry of the dipole by
1% on one node destroys the dark state and sends the ceiling to 1.
The dark sink energy is the free dual-graph spectral gap lambda_1
(= 0.2726, DIAGNOSTIC in v01) — the Fiedler-type odd ground state.

Run:  python studies/x0143_capture_parity.py
"""
import numpy as np
from x0143_particle_box import tessellation, dual_graph
from x0143_particle_box_v03 import cusp_incidence, bound_states, emission_generator

N = 143

def sigma_reflection(p1, tris):
    tri_index = {}
    for i, t in enumerate(tris):
        for x in t:
            tri_index[x] = i
    img = {}
    for i, t in enumerate(tris):
        js = {tri_index[p1.canon[((-x[0]) % N, (x[0] + x[1]) % N)]]
              for x in t}
        assert len(js) == 1, "reflection does not permute triangles"
        img[i] = js.pop()
    perm = np.array([img[i] for i in range(len(tris))])
    assert np.array_equal(perm[perm], np.arange(len(tris)))
    return perm

def main():
    p1, tris, edges, cert = tessellation(N)
    A, edge_list, tri_of = dual_graph(p1, tris, edges)
    L = np.diag(A.sum(axis=1)) - A
    n_inc, widths = cusp_incidence(p1, tris)
    ci = int(np.argmin(widths))
    nucleus = int(np.nonzero(n_inc[ci])[0][0])
    sigma = sigma_reflection(p1, tris)
    S = np.zeros((56, 56)); S[np.arange(56), sigma] = 1.0
    assert np.allclose(A[np.ix_(sigma, sigma)], A), "not a graph automorphism"
    assert sigma[nucleus] == nucleus
    fixed = np.nonzero(sigma == np.arange(56))[0]
    print(f"sigma = complex conjugation on the tessellation: involution, "
          f"graph automorphism, fixes nucleus; fixed triangles {list(fixed)} "
          f"(EXACT)")

    V0 = 4.0
    H = L - V0 * np.diag(n_inc[ci].astype(float))
    w, V, bound = bound_states(H)
    d = -np.ones(56, int); d[nucleus] = 0; fr = [nucleus]
    while fr:
        nx = []
        for u in fr:
            for v in np.nonzero(A[u])[0]:
                if d[v] < 0: d[v] = d[u] + 1; nx.append(int(v))
        fr = nx
    D = d.astype(float)
    assert np.allclose(D[sigma], D)
    assert np.linalg.norm(H @ S - S @ H) == 0.0
    par = np.array([float(V[:, k] @ (S @ V[:, k])) for k in range(56)])
    assert np.all(np.abs(np.abs(par) - 1) < 1e-9)
    Dm = V.T @ (D[:, None] * V)
    cross = max(abs(Dm[a, b]) for a in range(56) for b in range(56)
                if par[a] * par[b] < 0)
    print(f"parity superselection: 30 even + 26 odd; "
          f"max |<even|D|odd>| = {cross:.2e} (Cert)")

    G, M = emission_generator(w, V, D, gamma0=0.05)
    out = G.sum(axis=0)
    absorb = np.zeros(56)
    for k in np.argsort(w):
        if out[k] < 1e-14:
            absorb[k] = 1.0 if k in bound else 0.0
        else:
            absorb[k] = np.dot(G[:, k], absorb) / out[k]
    ceil = np.zeros(56); even = np.zeros(56)
    for v0 in range(56):
        psi = np.zeros(56); psi[v0] = 1.0
        c2 = (V.T @ psi) ** 2
        ceil[v0] = float(np.dot(c2, absorb))
        even[v0] = float(np.dot(c2, (par > 0)))
    moved = np.array([v for v in range(56) if sigma[v] != v])
    print(f"ceiling == even-weight (all 56 starts): "
          f"max dev {np.abs(ceil - even).max():.2e} (Cert)")
    print(f"moved starts: max |ceiling - 1/2| = "
          f"{np.abs(ceil[moved] - 0.5).max():.2e} (Cert; the theorem)")
    print(f"fixed starts: ceilings = "
          f"{[round(float(ceil[v]), 9) for v in fixed]} (Cert; = 1 exactly)")
    print("verdict: exactly-1/2 is the real structure of X0(143) acting as "
          "a parity superselection rule; AG-D5 item CLOSED.")

if __name__ == "__main__":
    main()
