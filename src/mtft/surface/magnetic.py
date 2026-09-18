"""KK-TOWER-01: magnetic Bochner spectra of line bundles O(D) on X0(143) in the compact (Liouville) metric.

P1 finite elements with Peierls phases on edges.  The Chern connection of O(D) with its HYM metric is, in the unitary
frame, the coexact 1-form (1/2)*d log h_D; on the mesh it is obtained by solving the dual-graph Laplacian for the edge
phases with prescribed triangle fluxes: the uniform smooth curvature -2 pi deg(D)/A_c per unit area, plus invisible 2 pi n
lumps in the triangles carrying the divisor points (cusp points: a triangle adjacent to the truncation circle).  The
stiffness matrix is conformally invariant; the mass matrix carries e^{2u}.  Validation: the holomorphic sections appear
as a Landau level at 2 pi deg/A_c with multiplicity h^0 (degree 30: 18 modes, +0.7%; degree 15: 3 modes; degree 12: 2).
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from .hym import assemble, solve_liouville, locate, Evaluator
from . import rrspace as RR

class MagneticMesh:
    def __init__(self, h=0.2, nx=8, Y0=2.0, N=143):
        self.Y0 = Y0; cx, faces, K, M, b = assemble(N, Y0, h, nx); self.cx, self.faces = cx, faces; G = K.shape[0]
        u, res, it = solve_liouville(K, M, b, Y0); self.u = u; self.Me = np.asarray(M @ np.exp(2 * u)).ravel(); self.Ac = float(self.Me.sum()); self.G = G
        self.TRI = []
        for fi, (nodes, tri, gid) in enumerate(faces):
            z = np.asarray(nodes)
            for T in tri: self.TRI.append((np.array([gid[T[0]], gid[T[1]], gid[T[2]]]), z[list(T)]))
        self.nT = len(self.TRI); self.edge_tris = {}
        for ti, (g3, z) in enumerate(self.TRI):
            for p, q in ((0, 1), (1, 2), (2, 0)): self.edge_tris.setdefault(self.edge_key(g3[p], g3[q]), []).append(ti)
        self.area_c = self._tri_area_c(); self.ev = Evaluator(cx, faces)
        width = [len(o) for o in cx.vertices]; self.cusp_nodes = {c: set() for c in range(4)}
        for fi, o in enumerate(cx.faces):
            nodes, tri, gid = faces[fi]; c_inf, c_0, c_1 = cx.vertex_of[o[0]], cx.vertex_of[cx.S(o[0])], cx.vertex_of[cx.S(o[1])]; Yi, Yz, Yo = Y0 * width[c_inf], Y0 * width[c_0], Y0 * width[c_1]
            for j in np.where(np.abs(nodes.imag - Yi) < 1e-7)[0]: self.cusp_nodes[c_inf].add(gid[j])
            for j in range(len(nodes)):
                if abs((1 - 1 / nodes[j]).imag - Yz) < 1e-6: self.cusp_nodes[c_0].add(gid[j])
                if abs((1 / (1 - nodes[j])).imag - Yo) < 1e-6: self.cusp_nodes[c_1].add(gid[j])
        self.bw = {width[c]: c for c in range(4)}
    @staticmethod
    def edge_key(a, b): return (a, b) if a < b else (b, a)
    def _tri_area_c(self):
        A = np.zeros(self.nT)
        for ti, (g3, z) in enumerate(self.TRI):
            x, y = z.real, z.imag; a2 = abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])); A[ti] = a2 / 2 / np.mean(y) ** 2 * np.exp(2 * np.mean(self.u[g3]))
        return A
    def tri_flux(self, theta):
        F = np.zeros(self.nT)
        for ti, (g3, z) in enumerate(self.TRI):
            for p, q in ((0, 1), (1, 2), (2, 0)): a, bb = g3[p], g3[q]; F[ti] += theta[self.edge_key(a, bb)] * (1.0 if a < bb else -1.0)
        return F
    def connection(self, target):
        """Coexact edge phases with prescribed triangle fluxes (dual-graph Laplacian solve)."""
        theta = {k: 0.0 for k in self.edge_tris}; rhs = target - self.tri_flux(theta); rows = []; cols = []; vals = []
        for key, ts in self.edge_tris.items():
            if len(ts) == 2: a, bb = ts; rows += [a, bb, a, bb]; cols += [a, bb, bb, a]; vals += [1, 1, -1, -1]
        Ld = sp.csr_matrix((vals, (rows, cols)), shape=(self.nT, self.nT)).tocsc(); psi = spl.spsolve(Ld + 1e-9 * sp.identity(self.nT), rhs - rhs.mean())
        for key, ts in self.edge_tris.items():
            if len(ts) == 2:
                a, bb = ts; g3a = self.TRI[a][0]; loc = [(p, q) for p, q in ((0, 1), (1, 2), (2, 0)) if self.edge_key(g3a[p], g3a[q]) == key][0]; sa = 1.0 if g3a[loc[0]] < g3a[loc[1]] else -1.0
                theta[key] = sa * (psi[a] - psi[bb])
        return theta, float(np.abs(self.tri_flux(theta) - target).max())
    def hamiltonian(self, theta):
        rows = []; cols = []; vals = []
        for (g3, z) in self.TRI:
            x, y = z.real, z.imag; P = np.stack([x, y], 1); Kloc = np.zeros((3, 3))
            for (p, q, r) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
                e1 = P[q] - P[p]; e2 = P[r] - P[p]; cot = np.dot(e1, e2) / abs(e1[0] * e2[1] - e1[1] * e2[0]); Kloc[q, r] -= cot / 2; Kloc[r, q] -= cot / 2; Kloc[q, q] += cot / 2; Kloc[r, r] += cot / 2
            for p in range(3):
                for q in range(3):
                    ph = 1.0 if p == q else np.exp(1j * theta[self.edge_key(g3[p], g3[q])] * (1.0 if g3[p] < g3[q] else -1.0)); rows.append(g3[p]); cols.append(g3[q]); vals.append(Kloc[p, q] * ph)
        KA = sp.csr_matrix((vals, (rows, cols)), shape=(self.G, self.G)); KA = (KA + KA.conj().T) / 2; Mi = sp.diags(1 / np.sqrt(self.Me)); H = (Mi @ KA @ Mi).tocsc()
        return (H + H.conj().T) / 2, Mi
    def eigenmodes(self, theta, nev):
        """Lowest eigenpairs; sections returned in the unitary frame, orthonormal for the compact area measure Me."""
        H, Mi = self.hamiltonian(theta); ev, V = spl.eigsh(H, k=nev, sigma=-1e-3, which="LM"); o = np.argsort(ev.real)
        return ev.real[o], np.asarray(Mi @ V[:, o])
    def divisor_modes(self, degree, lumps, nev):
        target = -2 * np.pi * degree / self.Ac * self.area_c
        for ti, n in lumps: target[ti] += 2 * np.pi * n
        theta, resid = self.connection(target); ev, psi = self.eigenmodes(theta, nev); return {"eigenvalues": ev, "modes": psi, "theta": theta, "landau": 2 * np.pi * abs(degree) / self.Ac, "flux_residual": resid}
    def eigenvalues(self, theta, nev):
        rows = []; cols = []; vals = []
        for (g3, z) in self.TRI:
            x, y = z.real, z.imag; P = np.stack([x, y], 1); Kloc = np.zeros((3, 3))
            for (p, q, r) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
                e1 = P[q] - P[p]; e2 = P[r] - P[p]; cot = np.dot(e1, e2) / abs(e1[0] * e2[1] - e1[1] * e2[0]); Kloc[q, r] -= cot / 2; Kloc[r, q] -= cot / 2; Kloc[q, q] += cot / 2; Kloc[r, r] += cot / 2
            for p in range(3):
                for q in range(3):
                    ph = 1.0 if p == q else np.exp(1j * theta[self.edge_key(g3[p], g3[q])] * (1.0 if g3[p] < g3[q] else -1.0)); rows.append(g3[p]); cols.append(g3[q]); vals.append(Kloc[p, q] * ph)
        KA = sp.csr_matrix((vals, (rows, cols)), shape=(self.G, self.G)); KA = (KA + KA.conj().T) / 2; Mi = sp.diags(1 / np.sqrt(self.Me)); H = (Mi @ KA @ Mi).tocsc(); H = (H + H.conj().T) / 2
        return np.sort(spl.eigsh(H, k=nev, sigma=-1e-3, which="LM", return_eigenvectors=False).real)
    def point_triangle(self, p):
        f_, lam, ids = locate(self.cx, self.faces, self.ev, p); return [ti for ti, (g3, z) in enumerate(self.TRI) if set(ids) <= set(g3)][0]
    def cusp_triangle(self, width):
        nodes = self.cusp_nodes[self.bw[width]]; return [ti for ti, (g3, z) in enumerate(self.TRI) if len(set(g3) & nodes) >= 2][0]
    def divisor_spectrum(self, degree, lumps, nev):
        """Bochner spectrum of O(D): lumps = [(triangle, multiplicity)], sum of multiplicities = degree.  Returns eigenvalues and
        the Landau value 2 pi degree/A_c (holomorphic sections sit there with multiplicity h^0)."""
        target = -2 * np.pi * degree / self.Ac * self.area_c
        for ti, n in lumps: target[ti] += 2 * np.pi * n
        theta, resid = self.connection(target); return {"eigenvalues": self.eigenvalues(theta, nev), "landau": 2 * np.pi * abs(degree) / self.Ac, "flux_residual": resid, "target_total_over_2pi": float(target.sum() / (2 * np.pi))}

def m1_towers(h=0.2, nx=8, nev_family=16, nev_higgs=28):
    """KK towers of M1: family bundle S0(sum P) (degree 15, h^0 = 3), spin structure S0 (degree 12, h^0 = 2), Higgs bundle
    K(2 sum P) (degree 30, h^0 = 18).  Cusp points 0 and 1/11 carry the theta-characteristic multiplicities."""
    mm = MagneticMesh(h, nx); cm = RR.cm_classes(); P = cm["P"]; T0 = {k: mm.point_triangle(P[k]) for k in (1, 2, 3)}; Tc0 = mm.cusp_triangle(143); Tc1 = mm.cusp_triangle(13)
    return {"A_c": mm.Ac,
            "family": mm.divisor_spectrum(15, [(Tc0, 6), (Tc1, 6), (T0[1], 1), (T0[2], 1), (T0[3], 1)], nev_family),
            "S0": mm.divisor_spectrum(12, [(Tc0, 6), (Tc1, 6)], 12),
            "higgs": mm.divisor_spectrum(30, [(Tc0, 12), (Tc1, 12), (T0[1], 2), (T0[2], 2), (T0[3], 2)], nev_higgs)}


def m1_fem_yukawa(h=0.2, nx=8, n_kk=3):
    """KK-TOWER-02: the M1 up-sector Yukawa tensor from kinetic-orthonormal FEM eigenmodes — Y_ijk = int psi_i psi_j conj(psi^H_k) dA_c
    with the three family zero modes of S0(sum P) and the 18 Higgs zero modes of K(2 sum P) = S0(sum P)^2 (the Higgs connection
    is exactly twice the family connection, so the integrand is gauge-invariant nodewise).  No Gram matrix or Cholesky enters:
    an independent check of the CC-26 normalisation.  Also returns the overlaps of the first n_kk family KK modes with
    (zero mode, Higgs mode), which control the KK decay widths KK -> f + H."""
    mm = MagneticMesh(h, nx); cm = RR.cm_classes(); P = cm["P"]; T0 = {k: mm.point_triangle(P[k]) for k in (1, 2, 3)}; Tc0 = mm.cusp_triangle(143); Tc1 = mm.cusp_triangle(13)
    F = mm.divisor_modes(15, [(Tc0, 6), (Tc1, 6), (T0[1], 1), (T0[2], 1), (T0[3], 1)], 3 + n_kk); Hg = mm.divisor_modes(30, [(Tc0, 12), (Tc1, 12), (T0[1], 2), (T0[2], 2), (T0[3], 2)], 20)
    psiF = F["modes"]; psiH = Hg["modes"][:, :18]
    Y = np.einsum("n,ni,nj,nk->ijk", mm.Me, psiF[:, :3], psiF[:, :3], np.conj(psiH))
    g = [np.einsum("n,n,ni,nk->ik", mm.Me, psiF[:, 3 + a], psiF[:, :3], np.conj(psiH)) for a in range(n_kk)]
    return {"Y": Y, "kk_overlaps": g, "family_eigenvalues": F["eigenvalues"], "higgs_eigenvalues": Hg["eigenvalues"], "A_c": mm.Ac,
            "family_landau": F["landau"], "higgs_landau": Hg["landau"]}

def ratio_distribution(Y, n=3000, seed=0):
    rng = np.random.default_rng(seed); R = []
    for _ in range(n):
        v = rng.standard_normal(Y.shape[2]) + 1j * rng.standard_normal(Y.shape[2]); s = np.sort(np.linalg.svd(np.einsum("ijk,k->ij", Y, v), compute_uv=False)); R.append(s / s[-1])
    return np.array(R)
