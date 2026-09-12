"""mtft.surface.spectral — hyperbolic Laplace spectrum of X0(N) by cusp-truncated FEM (SPEC-01, v0.28.2).

Construction: one truncated ideal triangle per Manin face, meshed in the upper-half-plane chart.
In two dimensions the Dirichlet energy is conformally invariant, so the P1 stiffness matrix of the
hyperbolic Laplacian is the Euclidean one in the chart; only the mass matrix carries dA = dx dy / y².
Cusps are truncated at horocycles with Y_c = Y0 * width(c) (a cusp form's tail decays like
e^{-2 pi Y / width}, so a common height fails at wide cusps — recorded design error of the first
prototype).  Faces are glued along the Manin edges with the reversal y -> 1/y, keyed by the
log-height parameter; Neumann conditions on the horocycles.

Reading the spectrum: with Neumann truncation the Eisenstein continuum becomes pseudo-modes that
MOVE with Y0 (toward 1/4 as Y0 grows); genuine Maass cusp forms are Y-stable.  ``maass_candidates``
returns the stable set from two truncations.  Gates: lambda_0 = 0 simple (connected surface);
smallest stable eigenvalue >= 975/4096 (Kim–Sarnak, proven) and compared with 1/4 (Selberg).
Results (h = 0.2, nx = 8): N = 11 stable {4.40, 6.47, 9.15, 9.54, 11.13, ...};
N = 143 stable {0.392, 0.451, 0.563, 0.972, 1.007, 1.011, 1.069, 1.202, 1.204, ...}, lambda_1 = 0.39.
Status: DIAGNOSTIC throughout (v0.28.4, KK-A07): mesh convergence (N = 11: 0.3%; N = 143: lambda_1
0.3923 -> 0.3910 between h = 0.2 and 0.15, 0.02% between truncation heights on the fine mesh) is a
diagnostic, not a continuum enclosure.  Scope (KK-A06): this is the untwisted scalar Laplacian on
the CUSPED quotient Y0(N) with Neumann horocycle truncation; its eigenvalues are m^2 (radius units)
for a minimally coupled scalar.  Fermion KK masses need the compact flux-twisted Dirac operator on
X0(N) (KK-04), a different operator on a different domain; these numbers must not be substituted.
"""
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spla
from scipy.spatial import Delaunay
from functools import lru_cache
import mtft.surface as S

R = lambda z: 1/(1-z)
R2 = lambda z: 1 - 1/z

def in_domain(z, Yi, Y0, Y1):
    x, y = z.real, z.imag
    return (x >= -1e-12) & (x <= 1+1e-12) & (abs(z-0.5) >= 0.5-1e-12) & (y <= Yi+1e-12) & \
           (x*x+y*y >= y/Y0-1e-12) & ((1-x)**2+y*y >= y/Y1-1e-12) & (y > 0)

@lru_cache(maxsize=None)
def face_mesh(Yi, Y0, Y1, h, nx):
    """Nodes/triangles of the ideal triangle 0,1,oo truncated at horocycles with parameters (Y_oo, Y_0, Y_1).
    Returns nodes, tri, and for the three edges the (node index, s) lists with s = log-height in the edge chart."""
    def edge_nodes(Ytop, Ybot):
        ks = np.arange(np.ceil(-np.log(Ybot)/h), np.floor(np.log(Ytop)/h)+1)
        s = np.unique(np.concatenate([[-np.log(Ybot)], ks*h, [np.log(Ytop)]]))
        return s
    s0 = edge_nodes(Yi, Y0); s1 = edge_nodes(Y0, Y1); s2 = edge_nodes(Y1, Yi)
    e0 = 1j*np.exp(s0); e1 = R(1j*np.exp(s1)); e2 = R2(1j*np.exp(s2))
    def horo(Y, f):  # interior horocycle nodes for a cusp at local height Y, mapped by f
        return f(np.linspace(0, 1, nx+2)[1:-1] + 1j*Y)
    env = lambda x: np.maximum(np.sqrt(np.maximum(1-x*x, 0)), np.sqrt(np.maximum(1-(1-x)**2, 0)))
    def third(Y, f):
        pts = []
        for x in np.linspace(0, 1, nx+2)[1:-1]:
            ny = max(3, int(np.log(Y/env(x))/h))
            for t in np.linspace(0, 1, ny+2)[1:-1]:
                pts.append(x + 1j*env(x)*(Y/env(x))**t)
        return f(np.array(pts))
    nodes = np.concatenate([e0, e1, e2, horo(Yi, lambda z: z), horo(Y0, R), horo(Y1, R2),
                            third(Yi, lambda z: z), third(Y0, R), third(Y1, R2)])
    ne = [len(e0), len(e1), len(e2)]
    P = np.column_stack([nodes.real, nodes.imag])
    tri = Delaunay(P).simplices
    cen = nodes[tri].mean(axis=1); tri = tri[in_domain(cen, Yi, Y0, Y1)]
    z = nodes[tri]; A2 = (z[:,1].real-z[:,0].real)*(z[:,2].imag-z[:,0].imag)-(z[:,2].real-z[:,0].real)*(z[:,1].imag-z[:,0].imag)
    tri = tri[np.abs(A2) > 1e-10]
    used = np.unique(tri); assert np.all(np.isin(np.arange(sum(ne)), used)), "edge node lost"
    remap = -np.ones(len(nodes), dtype=int); remap[used] = np.arange(len(used))
    edges = [(remap[np.arange(0, ne[0])], s0), (remap[np.arange(ne[0], ne[0]+ne[1])], s1), (remap[np.arange(ne[0]+ne[1], sum(ne))], s2)]
    return nodes[used], remap[tri], edges

def surface_spectrum(N, Y0=1.0, h=0.12, nx=14, k=20):
    cx = S.cell_complex(N)
    width = [len(o) for o in cx.vertices]
    Kblocks, Mblocks, rows, cols = [], [], [], []
    key_index = {}; nxt = 0
    for f, o in enumerate(cx.faces):
        c_inf, c_0, c_1 = cx.vertex_of[o[0]], cx.vertex_of[cx.S(o[0])], cx.vertex_of[cx.S(o[1])]
        Yi, Yz, Yo = Y0*width[c_inf], Y0*width[c_0], Y0*width[c_1]
        nodes, tri, edges = face_mesh(Yi, Yz, Yo, h, nx)
        Kl, Ml = assemble_local(nodes, tri)
        gid = -np.ones(len(nodes), dtype=int)
        for i, d in enumerate(o):
            e, sgn = cx.edge_of[d]
            idx, s = edges[i]
            for j, sv in zip(idx, s):
                key = (e, round(float(sgn*sv), 9))
                if key not in key_index: key_index[key] = nxt; nxt += 1
                gid[j] = key_index[key]
        for j in np.where(gid < 0)[0]: gid[j] = nxt; nxt += 1
        Ks, Ms = sps.coo_matrix(Kl), sps.coo_matrix(Ml)          # v0.28.4: scatter sparse triplets only
        rows.append(gid[Ks.row]); cols.append(gid[Ks.col]); Kblocks.append(Ks.data)
        rows.append(gid[Ms.row]); cols.append(gid[Ms.col]); Kblocks.append(np.zeros_like(Ms.data)); Mblocks.append(np.zeros_like(Ks.data)); Mblocks.append(Ms.data)
    G = nxt
    K = sps.coo_matrix((np.concatenate(Kblocks), (np.concatenate(rows), np.concatenate(cols))), shape=(G, G)).tocsc()
    M = sps.coo_matrix((np.concatenate(Mblocks), (np.concatenate(rows), np.concatenate(cols))), shape=(G, G)).tocsc()
    vals = spla.eigsh(K, k=k, M=M, sigma=-1e-2, which="LM", return_eigenvectors=False)
    return np.sort(vals), G, cx.inv


def assemble_local(nodes, tri):
    n = len(nodes); K = np.zeros((n, n)); M = np.zeros((n, n))
    # 6-point Gauss quadrature on triangle for the hyperbolic mass 1/y^2
    qw = np.array([0.109951743655322]*3 + [0.223381589678011]*3)/2
    qa = np.array([[0.816847572980459,0.091576213509771,0.091576213509771],[0.091576213509771,0.816847572980459,0.091576213509771],[0.091576213509771,0.091576213509771,0.816847572980459],
                   [0.108103018168070,0.445948490915965,0.445948490915965],[0.445948490915965,0.108103018168070,0.445948490915965],[0.445948490915965,0.445948490915965,0.108103018168070]])
    for t in tri:
        z = nodes[t]; x, y = z.real, z.imag
        A2 = (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0])
        if A2 < 0: t = t[[0,2,1]]; z = nodes[t]; x, y = z.real, z.imag; A2 = -A2
        b = np.array([y[1]-y[2], y[2]-y[0], y[0]-y[1]]); c = np.array([x[2]-x[1], x[0]-x[2], x[1]-x[0]])
        K[np.ix_(t,t)] += (np.outer(b,b)+np.outer(c,c))/(2*A2)          # Euclidean stiffness = hyperbolic (conformal invariance)
        for w, lam in zip(qw, qa):
            yq = lam @ y
            M[np.ix_(t,t)] += w*A2*np.outer(lam, lam)/yq**2
    return K, M



def maass_candidates(N, Y0a=1.0, Y0b=1.6, h=0.2, nx=8, k=18, tol=0.006):
    """Y-stable eigenvalues between two truncations, plus gates."""
    va, G, inv = surface_spectrum(N, Y0=Y0a, h=h, nx=nx, k=k)
    vb, _, _ = surface_spectrum(N, Y0=Y0b, h=h, nx=nx, k=k)
    stable = [float(x) for x in va[1:] if np.min(np.abs(vb - x)) / x < tol]
    moving = [float(x) for x in va[1:] if np.min(np.abs(vb - x)) / x >= tol]
    lam1 = min(stable) if stable else None
    gates = {"lambda0_zero": abs(float(va[0])) < 1e-6,
             "lambda0_simple": float(va[1]) > 1e-3,
             "kim_sarnak_975_4096": (lam1 is not None) and lam1 >= 975/4096,
             "selberg_quarter": (lam1 is not None) and lam1 >= 0.25}
    return {"N": N, "stable": stable, "pseudo_modes_Y0a": moving, "lambda_1": lam1, "dofs": G,
            "area_over_4pi": inv.index / 12, "gates": gates,
            "class": "DIAGNOSTIC",   # v0.28.4 (KK-A07): mesh convergence is not a continuum enclosure; no CERTIFIED tag without one
            "operator": "untwisted scalar hyperbolic Laplacian on the cusped Y0(N) (Neumann horocycle truncation); "
                        "NOT the compact flux-twisted Dirac operator of KK-04; a scalar eigenvalue is m^2 in radius units for a "
                        "minimally coupled scalar (KK-A06)"}


# ------------------------------------------------ SPEC-02: Hecke / Atkin–Lehner census of the modes
from . import transport as _TR
import mtft.surface as _S


def build_with_modes(N, Y0=1.0, h=0.2, nx=8, k=16):
    """Same assembly as surface_spectrum but returns faces (nodes, tri, gid), eigenpairs and M."""
    cx = _S.cell_complex(N); width = [len(o) for o in cx.vertices]
    faces = []; key_index = {}; nxt = 0; rows = []; cols = []; Kb = []; Mb = []
    for f, o in enumerate(cx.faces):
        c_inf, c_0, c_1 = cx.vertex_of[o[0]], cx.vertex_of[cx.S(o[0])], cx.vertex_of[cx.S(o[1])]
        nodes, tri, edges = face_mesh(Y0*width[c_inf], Y0*width[c_0], Y0*width[c_1], h, nx)
        Kl, Ml = assemble_local(nodes, tri); gid = -np.ones(len(nodes), dtype=int)
        for i, d in enumerate(o):
            e, sgn = cx.edge_of[d]; idx, s = edges[i]
            for j, sv in zip(idx, s):
                key = (e, round(float(sgn*sv), 9))
                if key not in key_index: key_index[key] = nxt; nxt += 1
                gid[j] = key_index[key]
        for j in np.where(gid < 0)[0]: gid[j] = nxt; nxt += 1
        Ks, Ms = sps.coo_matrix(Kl), sps.coo_matrix(Ml)
        rows.append(gid[Ks.row]); cols.append(gid[Ks.col]); Kb.append(Ks.data); Mb.append(np.zeros_like(Ks.data))
        rows.append(gid[Ms.row]); cols.append(gid[Ms.col]); Kb.append(np.zeros_like(Ms.data)); Mb.append(Ms.data)
        faces.append((nodes, tri, gid))
    G = nxt
    K = sps.coo_matrix((np.concatenate(Kb), (np.concatenate(rows), np.concatenate(cols))), shape=(G, G)).tocsc()
    M = sps.coo_matrix((np.concatenate(Mb), (np.concatenate(rows), np.concatenate(cols))), shape=(G, G)).tocsc()
    vals, vecs = spla.eigsh(K, k=k, M=M, sigma=-1e-2, which="LM"); order = np.argsort(vals)
    return cx, faces, vals[order], vecs[:, order], M

class Evaluator:
    def __init__(self, cx, faces):
        self.cx=cx; self.faces=faces; self.N=cx.N
        self.canon=_S.manin.canonical_factory(cx.N)
        self.dart_index={}
        for f,o in enumerate(cx.faces):
            for i,d in enumerate(o): self.dart_index[d]=(f,i)
        self.lifts=[np.array(_TR.lift_dart(o[0][0],o[0][1],cx.N)).reshape(2,2) for o in cx.faces]  # (a,b,c,d) -> [[a,b],[c,d]]
    def reduce(self, z):
        """z in H -> (w in T0, dart coset (c,d) mod N) with z ~ g_dart(w) mod Gamma_0(N)."""
        s=np.eye(2,dtype=np.int64)              # sigma with z_F = sigma z
        for _ in range(200):
            n=int(np.round(z.real)); z=z-n; s=np.array([[1,-n],[0,1]])@s
            if abs(z)<1-1e-14: z=-1/z; s=np.array([[0,-1],[1,0]])@s
            else: break
        if z.real<0: z=z+1; s=np.array([[1,1],[0,1]])@s
        # z_now = s z_orig, w = z_now in T0; coset of s^-1
        si=np.array([[s[1,1],-s[0,1]],[-s[1,0],s[0,0]]])
        return z, self.canon(int(si[1,0])%self.N, int(si[1,1])%self.N)
    def locate(self, z):
        w,dart=self.reduce(z); f,i=self.dart_index[dart]
        for _ in range(i): w=R(w)               # face-chart coordinate
        return f,w
    def value(self, u, z):
        f,w=self.locate(z); nodes,tri,gid=self.faces[f]
        P=nodes[tri]; x0,y0=P[:,0].real,P[:,0].imag; x1,y1=P[:,1].real,P[:,1].imag; x2,y2=P[:,2].real,P[:,2].imag
        det=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)
        l1=((w.real-x0)*(y2-y0)-(x2-x0)*(w.imag-y0))/det; l2=((x1-x0)*(w.imag-y0)-(w.real-x0)*(y1-y0))/det; l0=1-l1-l2
        lam=np.column_stack([l0,l1,l2]); ok=np.all(lam>-2e-3,axis=1)
        if not ok.any():                          # nearest triangle fallback (truncated cusp region)
            self.fallback_count = getattr(self, "fallback_count", 0) + 1
            j=int(np.argmin(np.min(lam,axis=1)*-1)); 
        else: j=int(np.where(ok)[0][0])
        lamj=np.clip(lam[j],0,1); lamj/=lamj.sum()
        return float(lamj@u[gid[tri[j]]])

def operators_on_mode(ev, u, samples, weights, N):
    vals=np.array([ev.value(u,z) for z in samples])
    norm=np.sum(weights*vals*vals)
    out={}
    for Q,(a,b,c,d) in (("W11",(66,5,143,11)),("W13",(78,7,143,13))):
        wz=[(a*z+b)/(c*z+d) for z in samples]; uw=np.array([ev.value(u,z) for z in wz]); out[Q]=float(np.sum(weights*uw*vals)/norm)
    t2=[]
    for z in samples:
        t2.append((ev.value(u,2*z)+ev.value(u,z/2)+ev.value(u,(z+1)/2))/np.sqrt(2))
    out["a2"]=float(np.sum(weights*np.array(t2)*vals)/norm)
    return out



def al_matrices(N):
    """Explicit Atkin–Lehner matrices [[Qa, b],[N c, Qd]] for each prime power exactly dividing N."""
    from .manin import factorize
    out = {}
    for p, a in factorize(N).items():
        Q = p ** a
        A, B, C, D = _TR.atkin_lehner_matrix(Q, N)
        out[f"W{Q}"] = (A, B, C, D)
    return out


def arithmetic_census(N, Y0=1.0, h=0.2, nx=8, k=16, samples_per_face=6, seed=1, stable_from=None):
    """For each Laplace eigenmode: AL parities, T2 eigenvalue, and the Eisenstein test
    a2 ≈ 2 cos(t log 2), t = sqrt(lambda − 1/4).  Cusp forms fail the Eisenstein test; pseudo-modes pass it."""
    cx, faces, vals, vecs, M = build_with_modes(N, Y0, h, nx, k)
    ev = Evaluator(cx, faces); rng = np.random.default_rng(seed); Md = M.diagonal()
    samples, weights = [], []
    for f, (nodes, tri, gid) in enumerate(faces):
        for j in rng.choice(len(nodes), size=min(samples_per_face, len(nodes)), replace=False):
            z = nodes[j]
            if z.imag < 0.3 or z.imag > 6: continue
            g = ev.lifts[f]; samples.append((g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])); weights.append(Md[gid[j]])
    weights = np.array(weights); W = al_matrices(N)
    rows = []
    for m in range(1, len(vals)):
        u = vecs[:, m]; v = np.array([ev.value(u, z) for z in samples]); norm = np.sum(weights*v*v)
        row = {"lambda": float(vals[m])}
        for name, (a, b, c, d) in W.items():
            uw = np.array([ev.value(u, (a*z+b)/(c*z+d)) for z in samples]); row[name] = float(np.sum(weights*uw*v)/norm)
        t2 = np.array([(ev.value(u, 2*z)+ev.value(u, z/2)+ev.value(u, (z+1)/2))/np.sqrt(2) for z in samples])
        row["a2"] = float(np.sum(weights*t2*v)/norm)
        t = np.sqrt(max(vals[m]-0.25, 0.0)); eis = 2*np.cos(t*np.log(2))
        row["eisenstein_a2"] = float(eis); row["eisenstein_like"] = abs(row["a2"]-eis) < 0.15   # coarse-mesh tolerance
        if stable_from is not None: row["Y_stable"] = any(abs(vals[m]-x) < 0.006*x for x in stable_from)
        rows.append(row)
    return rows


# ------------------------------------------------ WL-01: Wilson-line twisted spectrum (v0.28.5)
def wilson_spectrum(N, theta, q=1, Y0=1.0, h=0.2, nx=8, k=4):
    """Scalar Laplacian on cusped Y0(N) twisted by the flat U(1) connection with cycle-coordinates theta
    (holonomy phases e^{i q zeta(e)}, zeta = B theta on Manin edges, B = tree/cotree cycle matrix)."""
    cx = S.cell_complex(N); cb = S.tree_cotree(cx)
    zeta = q * (cb.basis_matrix.astype(float) @ np.asarray(theta, float))         # per Manin edge
    width = [len(o) for o in cx.vertices]
    rows, cols, Kv, Mv = [], [], [], []
    key_index = {}; nxt = 0
    for f, o in enumerate(cx.faces):
        c_inf, c_0, c_1 = cx.vertex_of[o[0]], cx.vertex_of[cx.S(o[0])], cx.vertex_of[cx.S(o[1])]
        nodes, tri, edges = face_mesh(Y0*width[c_inf], Y0*width[c_0], Y0*width[c_1], h, nx)
        Kl, Ml = assemble_local(nodes, tri)
        gid = -np.ones(len(nodes), dtype=int); phase = np.ones(len(nodes), dtype=complex)
        for i, d in enumerate(o):
            e, sgn = cx.edge_of[d]; idx, s = edges[i]
            for j, sv in zip(idx, s):
                key = (e, round(float(sgn*sv), 9))
                if key not in key_index: key_index[key] = nxt; nxt += 1
                gid[j] = key_index[key]
                if sgn == -1: phase[j] = np.exp(1j*zeta[e])                          # transition phase across the glued edge
        for j in np.where(gid < 0)[0]: gid[j] = nxt; nxt += 1
        Ks, Ms = sps.coo_matrix(Kl), sps.coo_matrix(Ml)
        pf = phase
        rows.append(gid[Ks.row]); cols.append(gid[Ks.col]); Kv.append(np.conj(pf[Ks.row])*Ks.data*pf[Ks.col])
        rows.append(gid[Ms.row]); cols.append(gid[Ms.col]); Kv.append(np.zeros(len(Ms.data)))
        Mv.append(np.zeros(len(Ks.data))); Mv.append(np.conj(pf[Ms.row])*Ms.data*pf[Ms.col])
    G = nxt
    K = sps.coo_matrix((np.concatenate(Kv), (np.concatenate(rows), np.concatenate(cols))), shape=(G, G)).tocsc()
    M = sps.coo_matrix((np.concatenate(Mv), (np.concatenate(rows), np.concatenate(cols))), shape=(G, G)).tocsc()
    vals = spla.eigsh(K, k=k, M=M, sigma=-1e-2, which="LM", return_eigenvectors=False)
    return np.sort(vals.real)



def wilson_line_mass_check(theta, q=1, Y0=1.0, h=0.2, nx=8):
    """Twisted ground state on X0(143) vs the Hodge prediction q^2 theta^T G theta / Area_truncated,
    with G = Jint J_true (frozen) in cycle coordinates and Area_truncated = 56 pi - 4/Y0 (four horocycle caps).
    A parameter-free comparison of the FEM instrument with the period-derived Hodge structure."""
    from .frozen import x0143
    d = x0143(); Jint = d["intersection_cycles"].astype(float); G = Jint @ d["J_true"]; G = (G + G.T) / 2
    G = -G if np.linalg.eigvalsh(G)[0] < 0 else G
    theta = np.asarray(theta, float)
    lam = wilson_spectrum(143, theta, q=q, Y0=Y0, h=h, nx=nx, k=3)
    area_trunc = 56 * np.pi - 4 / Y0
    pred = q * q * float(theta @ G @ theta) / area_trunc
    return {"lambda_0": float(lam[0]), "prediction": pred, "ratio": float(lam[0] / pred), "lambda_1": float(lam[1]),
            "class": "DIAGNOSTIC", "note": "vector-like pair lifted by any nontrivial Wilson line; mass^2 = q^2 <theta,theta>_Hodge / Area to O(theta^4)"}
