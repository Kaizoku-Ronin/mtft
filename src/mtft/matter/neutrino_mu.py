"""NU-02 (2026-10-02): the size of the lepton-number-violating mass mu — the arithmetic geometry of the localisation points and the
suppression it can supply.

1. The modular parametrisation, built from scratch (EXACT / NUMERICAL to 1e-12).  143a is the (+,+) newform: a_p = -sum_x chi(4x^3 - 4x^2 - 4x - 7)
   (the 2-division polynomial) for p not dividing 143, a_11 = a_13 = -1.  phi(tau) = sum a_n q^n / n; its periods along hyperbolic elements of
   Gamma_0(143) are integral in the basis (Omega+, Omega+ tau_E) with Omega+ = 2 int_{e1}^inf dx / sqrt(4x^3 - 4x^2 - 4x - 7) = 1.96992 and
   tau_E = 1/2 + 1.02327 i — Manin constant 1; phi o W_143 = phi to machine precision.
2. Where the fermionic partners live (NUMERICAL, certified by j-invariant pairs).  For c = 1 the fermionic gamma-modes sit on phi^-1(P_eta).  For each
   of the three nonzero 2-torsion points the fibre has exactly four points modulo Gamma_0(143), and the Atkin–Lehner involutions act on each fibre
   without fixed points: like the cusps, every fibre carries the regular representation (one mode per sector, (-,-) included).  The points lie deep
   in the bulk: hyperbolic distance rho = 3.589 (real spin structure) and 3.397 (the complex pair) below the unit horocycles of the cusps; the stretch
   of phi there (unit-area 143a1 length per hyperbolic length) is s = 0.1965 and 0.2049, identical on each fibre.
3. Suppression (ESTIMATE, two bulk methods with the same form).  A fermion zero mode decays away from its point like exp(-int m dl) with
   m = 2 pi |c| dist_E(phi(x), P_eta) / A_E; the Majorana mass needs it twice at the condensate:
       mu / M ~ exp(-beta sqrt(r)),   r = A_C / A_E,
       beta = 2 pi d^2 / (s sqrt(56 pi))  (constant stretch, d^2 = dist(O, P_eta)^2 / A_E)  = 0.588 (real), 0.732 (pair),
       beta = sigma_eta = pi d^2                      (uniform conformal stretch)            = 0.768 (real), 0.996 (pair).
   The arithmetic separation alone is about one width (sigma ~ 0.8 - 1.0): the large suppression comes from the size hierarchy r, which stretches it.
   The cusp region adds further, metric-dependent suppression (in the complete hyperbolic metric the cusps are infinitely far).  Targets (NU-01):
   second family 5.7e-5 needs r ~ 160 - 280; third family 3.8e-9 needs r ~ 640 - 1090 from the bulk alone.  At the r ~ 56 pi preferred by the QCD
   window (Matter-05), the bulk gives 1e-4 - 1e-5: enough for the second family; the third needs the cusp region or a larger r.
4. Flux (EXACT).  The exponent scales as 1/|c| (the fermion points spread to [c]^-1(P_eta)).  The cheap Pati–Salam flux of NU-01 puts c = 2 on the nu^c
   entry and halves it; (c, A1, L, B, A2) = (1, 1, -2, 0, 0) is also D-flat with the (4,1,2) condensate alone and puts c = 1 there, at 384 localised
   extras instead of 256.
OPEN: the smooth metric near the cusps (Arakelov is the canonical choice) and the exact overlap integrals; the Dirac pairing M; family structure — with an
Atkin–Lehner-symmetric background (all four cusps blown up), couplings conserve the Atkin–Lehner character and mu becomes block-diagonal in the four sectors."""
import math
import numpy as np
import mpmath as mp

A_INV = (0, -1, 1, -1, -2)
W = {11: ((66, 5), (143, 11)), 13: ((13, 7), (143, 78)), 143: ((0, -1), (143, 0))}

def _primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]

def newform_143a(nmax=3000):
    """a_1..a_nmax of the newform attached to 143a1 (a_p from the 2-division polynomial; a_11 = a_13 = -1; Hecke multiplicativity)."""
    ap = {}
    for p in _primes(nmax):
        if p in (11, 13): ap[p] = -1; continue
        if p == 2:
            ap[p] = 2 + 1 - (1 + sum(1 for x in range(2) for y in range(2) if (y * y + y - (x ** 3 - x * x - x - 2)) % 2 == 0)); continue
        tot = 0
        for x in range(p):
            v = (4 * x ** 3 - 4 * x * x - 4 * x - 7) % p
            if v: tot += 1 if pow(v, (p - 1) // 2, p) == 1 else -1
        ap[p] = -tot
    a = np.zeros(nmax + 1); a[1] = 1
    for n in range(2, nmax + 1):
        m, val, d = n, 1, 2
        while d * d <= m:
            k = 0
            while m % d == 0: m //= d; k += 1
            if k: val *= _apk(ap, d, k)
            d += 1
        if m > 1: val *= _apk(ap, m, 1)
        a[n] = val
    return a

def _apk(ap, p, k):
    if p in (11, 13): return ap[p] ** k
    A = [1, ap[p]]
    for _ in range(2, k + 1): A.append(ap[p] * A[-1] - p * A[-2])
    return A[k]

class ModularParametrisation:
    def __init__(self, nmax=3000):
        self.a = newform_143a(nmax); self.n = np.arange(1, nmax + 1); self.an = self.a[1:]
        mp.mp.dps = 30
        e1 = mp.findroot(lambda x: 4 * x ** 3 - 4 * x ** 2 - 4 * x - 7, 1.96)
        self.Omega = float(2 * mp.quad(lambda x: 1 / mp.sqrt(4 * x ** 3 - 4 * x ** 2 - 4 * x - 7), [e1, e1 + 1, mp.inf]))
        from mtft.research.theta_torus import tau_143a1
        self.tauE = complex(tau_143a1()); self.w1, self.w2 = self.Omega, self.Omega * self.tauE
        self.area = self.w1 * self.w2.imag
    def f(self, t): q = np.exp(2j * np.pi * t); return complex(np.sum(self.an * q ** self.n))
    def phi(self, t): q = np.exp(2j * np.pi * t); return complex(np.sum(self.an / self.n * q ** self.n))
    def coords(self, z):
        M = np.array([[self.w1.real, self.w2.real], [self.w1.imag, self.w2.imag]]); return np.linalg.solve(M, [z.real, z.imag])
    def reduce(self, z):
        c = self.coords(z); c = c - np.round(c); return c[0] * self.w1 + c[1] * self.w2
    def two_torsion(self):
        return {"real 1/2": self.w1 / 2, "tau/2": self.w2 / 2, "(1+tau)/2": (self.w1 + self.w2) / 2}
    def solve(self, t, P, it=80):
        for _ in range(it):
            F = self.reduce(self.phi(t) - P); dF = 2j * math.pi * self.f(t)
            t = t - F / dF
            if not 0.003 < t.imag < 3: return None
            if abs(F) < 1e-12: return t
        return None

def act(M, t): (a, b), (c, d) = M; return (a * t + b) / (c * t + d)

def _sl2_reduce(t):
    for _ in range(300):
        t = complex(t.real - round(t.real), t.imag)
        if abs(t) < 1 - 1e-12: t = -1 / t
        else: break
    return t

def jpair(t):
    """(j(tau), j(143 tau)): a complete invariant of the point of X0(143)."""
    mp.mp.dps = 30
    return complex(mp.kleinj(mp.mpc(_sl2_reduce(t)))), complex(mp.kleinj(mp.mpc(_sl2_reduce(143 * t))))

def _same(p, q): return all(abs(x - y) < 1e-6 * (1 + abs(x)) for x, y in zip(p, q))

# approximate fibre points found by Newton from random starts (session 2026-10-02); the tests re-converge them
KNOWN_FIBRES = {"real 1/2": [0.07893 + 0.027624j, 0.462946 + 0.019371j, 0.92107 + 0.027624j, 0.537054 + 0.019371j],
                "tau/2": [0.726247 + 0.018942j, 0.18363 + 0.033464j, -0.036858 + 0.006717j, 0.761248 + 0.008969j],
                "(1+tau)/2": [0.81637 + 0.033464j, 1.273753 + 0.018942j, 0.036858 + 0.006717j, 0.146259 + 0.005386j]}

def fibres(MP=None):
    """Re-converge the four points over each nonzero 2-torsion point, check distinctness and the Atkin–Lehner action."""
    MP = MP or ModularParametrisation(); out = {}
    for name, P in MP.two_torsion().items():
        pts = [MP.solve(t, P) for t in KNOWN_FIBRES[name]]
        jps = [jpair(t) for t in pts]
        distinct = all(not _same(jps[i], jps[k]) for i in range(4) for k in range(i + 1, 4))
        def idx(jp):
            for i, s in enumerate(jps):
                if _same(jp, s): return i
            return None
        perms = {Q: [idx(jpair(act(W[Q], t))) for t in pts] for Q in (11, 13)}
        perms[143] = [idx((jp[1], jp[0])) for jp in jps]
        free = all(all(p[i] is not None and p[i] != i for i in range(4)) for p in perms.values())
        heights = []
        for x in pts:
            h = x.imag
            for k in (1, 2, 3):
                for d in range(-600, 601):
                    if math.gcd(143 * k, d) == 1: h = max(h, x.imag / abs(143 * k * x + d) ** 2)
            heights.append(h)
        depth = max(heights)
        s = [2 * math.pi * abs(MP.f(x)) * x.imag / math.sqrt(MP.area) for x in pts]
        d2 = min(abs(P - (m * MP.w1 + n * MP.w2)) ** 2 for m in (-1, 0, 1) for n in (-1, 0, 1)) / MP.area
        out[name] = {"points": pts, "distinct": distinct, "AL_perms": perms, "AL_free": free, "heights": heights, "rho": math.log(1 / depth),
                     "stretch": s, "d2": d2, "sigma": math.pi * d2,
                     "beta_constant_stretch": 2 * math.pi * d2 / (s[0] * math.sqrt(56 * math.pi)), "beta_uniform": math.pi * d2}
    return out

def lattice_check(MP=None):
    """Periods of phi along hyperbolic elements [[a, b], [143, d]] (evaluated at height 1/143): integral lattice coordinates."""
    MP = MP or ModularParametrisation(); rows = []
    for G in (((7, 2), (143, 41)), ((5, 3), (143, 86))):
        t0 = -G[1][1] / 143 + 1j / 143
        rows.append(MP.coords(MP.phi(act(G, t0)) - MP.phi(t0)))
    w = MP.phi(act(W[143], 0.31 + 0.11j)) - MP.phi(0.31 + 0.11j)
    return {"periods": rows, "W143_shift": MP.coords(w)}

TARGETS = {"second family": 5.7e-5, "third family": 3.8e-9}

def required_r(beta, target): return (math.log(1 / target) / beta) ** 2

def requirements(fib=None):
    fib = fib or fibres(); out = {}
    for name in ("real 1/2", "tau/2"):
        b1, b2 = fib[name]["beta_constant_stretch"], fib[name]["beta_uniform"]
        out[name] = {t: (required_r(b2, v), required_r(b1, v)) for t, v in TARGETS.items()}
        out[name]["at_r_56pi"] = (math.exp(-b2 * math.sqrt(56 * math.pi)), math.exp(-b1 * math.sqrt(56 * math.pi)))
    return out

def flux_options():
    from mtft.matter import neutrino_gate as NG
    opts = {"cheap (1,1,-1,-1,-1)": NG.GAMMA_CHEAP, "PS-01 (1,1,0,-2,-2)": NG.GAMMA_PS01, "nu-optimal (1,1,-2,0,0)": {"c": 1, "A1": 1, "L": -2, "B": 0, "A2": 0}}
    return {k: {"d_flat_reps": [v["rep"] for v in NG.fi_check(g)["vacua"]], "nu_c_gamma_charge": NG.entry("A2", "A1", g)["gamma_charge"],
                "extras": sum(NG.localised_extras(g).values())} for k, g in opts.items()}
