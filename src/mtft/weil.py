"""mtft.weil -- Gabor-compressed Weil explicit-formula form (CANDIDATE, W1).

STATUS: staged for v0.15.0; does NOT enter the corpus until Kimi independently
(i) re-derives the rank-trace lemma, (ii) reproduces the E2 identity by an
implementation sharing no code with this one, (iii) confirms CC-02.

Source of the mathematics: "More than two thirds of the zeros of the Riemann
zeta function lie on the critical line" (Claude/Anthropic, Aug 2026; Lean repo
github.com/anthropics/zeta-23-lean). The asymptotic Theorems A-E of that paper
are [Ext] here -- external claims, not corpus-certified. Everything this module
computes is a FINITE identity or inequality certified independently of them:

  [Cert] prime-side/zero-side identity of the compressed form G (eq. 2.20 of
         the source): on I=[100,200], lam=1, eta=0.2, 371 zeros, max rel
         discrepancy 3.386e-6 (W1 study, 2026-08-10).
  [Cert] Poisson/Gabor frame identity (Lemma 2.2): rel err <= 6.3e-5 at
         lattice truncation K=4000.
  [Pr]   rank-trace inequality (Lemma 3.2 of the source; von Neumann trace
         inequality + x^2 >= cx - c^2/4). Numeric audit: 1e5 random
         instances, 0 violations; equality cases exact to 1.5e-14.
  [Pr]   Sum_{n>=1} w_n n^{-s} = F(s+1) = -zeta(s) * zeta'(s+1)   (CC-02).
         Three routes: Dirichlet algebra (w = Lambda_1 * 1, Lambda_1 = log/id),
         the shift w_n = f(n)/n with the settled F(s) = -zeta(s-1)zeta'(s),
         and numeric CERTIFIED(5.2e-12) at s=3, N=3e5.

WI-N1 (honest negative, filed): for an individual GL(2) L-function -- in
particular each X0(143) newform L-function -- this method certifies NOTHING.
Mechanism: the Montgomery-Vaughan mean-value step caps the arithmetic length
at X <= T^{1-eps}, so the normalised band-width cannot exceed Lambda* = 1/m
for degree m; at m = 2 the best constant is c = (1/2)/(1 + 1/12) = 6/13 and
the certified proportion 2 - 1/c = -1/6 < 0, for every window (source,
Remark 7.2(ii)). Do not attempt a direct port to the f1/f2/f3 L-functions.
The degree-1 route (Dirichlet characters mod 11, 13, 143 via Theorem E) is
open: see study proposal W2.

Conventions (source paper Sec. 1.8, 2.1): fhat(tau) = int f(u) e^{i tau u} du;
window I = [T, 2T]; l = log(T/2pi); L = lam*l; X = e^L; h = 2pi/L;
d = floor(T/h); tau_k = T + k h. Hard deps: numpy, mpmath. scipy optional
(complex digamma fast path).
"""
from __future__ import annotations
import numpy as np
import mpmath as mp

__all__ = [
    "Window", "gabor", "prime_powers", "nu_parts", "G_prime", "G_zero",
    "certificates", "rank_trace_gap", "audit_rank_trace", "mt_constant",
    "w_series_check", "selftest",
]

TWO_PI = 2.0 * np.pi


# ----------------------------------------------------------------------
# window / taper
# ----------------------------------------------------------------------
class Window:
    """C^2 flat-top taper phi(u) = ramp((L/2-|u|)/w), ramp(x)=x-sin(2pi x)/2pi,
    ramp width w = eta*L/2 (proportional taper of source Sec. 8). Provides the
    cosine transform phihat, the transform Phi of phi^2, taper constants a, b,
    the autocorrelation g = phi^2 * phi^2, and complex-argument phihat for
    synthetic off-line configurations."""

    def __init__(self, L: float, eta: float = 0.2, n_u: int = 6001, r_max: float = 801.0):
        if n_u % 2 == 0:
            n_u += 1
        self.L, self.eta, self.w = float(L), float(eta), eta * L / 2.0
        self.u = np.linspace(0.0, L / 2.0, n_u)
        du = self.u[1] - self.u[0]
        wS = np.ones(n_u); wS[1:-1:2] = 4.0; wS[2:-1:2] = 2.0
        self._wS = wS * du / 3.0                      # Simpson weights on [0, L/2]
        self.pu = self.phi(self.u)
        self.p2u = self.pu * self.pu
        self._rg = np.arange(0.0, r_max, 0.005)
        self._ph = self._cosT(self.pu, self._rg)
        self._Ph = self._cosT(self.p2u, self._rg)
        full_u = np.concatenate([-self.u[::-1], self.u[1:]])
        p2f = np.concatenate([self.p2u[::-1], self.p2u[1:]])
        self.a = float(np.trapezoid(p2f, full_u)) / L
        self.b = float(np.trapezoid(p2f * p2f, full_u)) / L
        gc = np.convolve(p2f, p2f) * (full_u[1] - full_u[0])
        yy = np.linspace(-L, L, gc.size)   # autocorrelation support of phi^2
        self._gy, self._gv = yy[yy >= 0], gc[yy >= 0]

    def phi(self, uu):
        x = np.clip((self.L / 2.0 - np.abs(uu)) / self.w, 0.0, 1.0)
        return x - np.sin(TWO_PI * x) / TWO_PI

    def _cosT(self, vals, r, step=4000):
        out = np.empty(np.shape(r), dtype=float)
        rf = np.ravel(r)
        of = out.ravel()
        for i in range(0, rf.size, step):
            rr = rf[i:i + step]
            of[i:i + step] = 2.0 * (np.cos(np.outer(rr, self.u)) @ (vals * self._wS))
        return out

    def phihat(self, r, direct: bool = False):
        """phihat on real arguments; grid+interp by default, direct Simpson if asked."""
        if direct:
            return self._cosT(self.pu, np.asarray(r, dtype=float))
        return np.interp(np.abs(r), self._rg, self._ph)

    def Phi(self, r):
        return np.interp(np.abs(r), self._rg, self._Ph)

    def g(self, y):
        return np.interp(np.abs(y), self._gy, self._gv)

    def phihat_complex(self, z):
        """phihat(z) = int phi(u) e^{i z u} du for complex z = x + i y (vectorised).
        Even phi: = int_0^{L/2} phi [e^{izu} + e^{-izu}] du."""
        z = np.asarray(z, dtype=complex)
        x, y = np.real(z), np.imag(z)
        ep = np.exp(-np.outer(y.ravel(), self.u))     # e^{-y u}
        em = np.exp(+np.outer(y.ravel(), self.u))     # e^{+y u}
        c = np.cos(np.outer(x.ravel(), self.u))
        s = np.sin(np.outer(x.ravel(), self.u))
        # e^{izu}+e^{-izu} = (e^{-yu}+e^{yu})cos(xu) + i (e^{-yu}-e^{yu}) sin(xu)
        re = (c * (ep + em)) @ self._wS
        im = (s * (ep - em)) @ self._wS
        return (re + 1j * im).reshape(z.shape)


# ----------------------------------------------------------------------
# grid, primes, density
# ----------------------------------------------------------------------
def gabor(T: float, lam: float = 1.0) -> dict:
    l = np.log(T / TWO_PI)
    L = lam * l
    h = TWO_PI / L
    d = int(np.floor(T / h))
    ell1 = l + 2 * np.log(2) - 1
    return dict(T=T, lam=lam, l=l, L=L, X=np.exp(L), h=h, d=d,
                tau_k=T + h * np.arange(d), ell1=ell1, lam1=L / ell1)


def prime_powers(x: float):
    """[(n, log p)] for prime powers n <= x. Sympy-free (release-gate rule)."""
    out, n = [], 2
    while n <= x:
        p = None
        for q in range(2, int(np.sqrt(n)) + 1):
            if n % q == 0:
                p = q
                break
        if p is None:
            out.append((n, float(np.log(n))))
        else:
            m = n
            while m % p == 0:
                m //= p
            if m == 1:
                out.append((n, float(np.log(p))))
        n += 1
    return out


def _mu_density(tau: np.ndarray) -> np.ndarray:
    """mu(tau) = Re psi(1/4 + i tau/2)/2pi - log(pi)/2pi. scipy fast path,
    mpmath coarse-grid + cubic refinement fallback (mu is smooth, mu'' ~ tau^-2)."""
    try:
        from scipy.special import digamma
        return np.real(digamma(0.25 + 0.5j * tau)) / TWO_PI - np.log(np.pi) / TWO_PI
    except Exception:
        lo, hi = float(np.min(tau)), float(np.max(tau))
        coarse = np.unique(np.concatenate([
            np.arange(lo - 1, hi + 1.25, 0.25),
            np.arange(-6, 6.05, 0.05),
        ]))
        vals = np.array([float(mp.re(mp.digamma(mp.mpc(0.25, 0.5 * t)))) for t in coarse])
        return np.interp(tau, coarse, vals) / TWO_PI - np.log(np.pi) / TWO_PI


def nu_parts(tau: np.ndarray, T: float, lam: float = 1.0):
    """(mu, Pi_X, P_X): the three parts of the unconditional prime-side density
    nu_X of source eqs. (2.3)-(2.6). Their sum integrates the compressed form."""
    g = gabor(T, lam)
    X, L = g["X"], g["L"]
    mu = _mu_density(tau)
    s = 0.5 + 1j * tau
    Pi = 1.0 / (TWO_PI * (0.25 + tau ** 2)) + np.real((np.sqrt(X) * np.exp(1j * tau * L) - 1) / s) / np.pi
    P = np.zeros_like(tau)
    for nq, lg in prime_powers(X):
        P -= (lg / np.sqrt(nq)) * np.cos(tau * np.log(nq)) / np.pi
    return mu, Pi, P


# ----------------------------------------------------------------------
# the compressed form, two independent routes
# ----------------------------------------------------------------------
def G_prime(T: float, lam: float = 1.0, eta: float = 0.2,
            pad: float = 250.0, dtau: float = 0.005, win: Window | None = None):
    """Prime-side route: G_kl = int phihat(tau-tau_k) phihat(tau-tau_l) nu_X dtau."""
    g = gabor(T, lam)
    win = win or Window(g["L"], eta)
    n = int(round((2 * T + pad - (T - pad)) / dtau))
    if n % 2 == 1:
        n += 1
    tau = np.linspace(T - pad, 2 * T + pad, n + 1)
    wS = np.ones(tau.size); wS[1:-1:2] = 4.0; wS[2:-1:2] = 2.0
    wS *= (tau[1] - tau[0]) / 3.0
    mu, Pi, P = nu_parts(tau, T, lam)
    nu = mu + Pi + P
    PhiM = win.phihat(tau[None, :] - g["tau_k"][:, None])
    G = (PhiM * (nu * wS)) @ PhiM.T
    return G, g, win


def G_zero(gammas, T: float, lam: float = 1.0, eta: float = 0.2,
           depths=None, mults=None, mirror: bool = True, win: Window | None = None):
    """Zero-side route: G = sum_rho m_rho phihat(gamma_rho - tau_k) conj(...).
    gammas: positive ordinates. depths[j] > 0 turns ordinate j into an off-line
    pair {rho, 1-rhobar} at |beta-1/2| = depths[j] (both members included).
    Returns (G, counts) with counts = dict(s1, s2, p, N)."""
    g = gabor(T, lam)
    win = win or Window(g["L"], eta)
    gammas = np.asarray(gammas, dtype=float)
    depths = np.zeros_like(gammas) if depths is None else np.asarray(depths, float)
    mults = np.ones_like(gammas) if mults is None else np.asarray(mults, float)
    tau_k = g["tau_k"]
    G = np.zeros((g["d"], g["d"]))
    s1 = s2 = p = 0
    N = 0.0
    on = depths == 0
    if np.any(on):
        gl = np.concatenate([gammas[on], -gammas[on]]) if mirror else gammas[on]
        ml = np.concatenate([mults[on], mults[on]]) if mirror else mults[on]
        U = win.phihat(gl[:, None] - tau_k[None, :], direct=True)
        G += (U * ml[:, None]).T @ U
        s1 = int(np.sum(mults[on] == 1)); s2 = int(np.sum(mults[on] > 1))
        N += float(np.sum(mults[on])) * (2 if mirror else 1)
        if mirror:
            s1 *= 2; s2 *= 2   # counts refer to the full gamma<->-gamma multiset
    off = ~on
    for gj, dj, mj in zip(gammas[off], depths[off], mults[off]):
        sides = (gj, -gj) if mirror else (gj,)
        for gs in sides:
            row = win.phihat_complex((gs - tau_k) - 1j * dj)
            # pair {rho, 1-rhobar}: summands A_k A_l + conj = 2m Re(a a^T),
            # the signature-(1,1) hyperbolic block of source Prop. 4.1 --
            # NOT Re(a conj(a)^T), which is PSD. (W1 correction record: the
            # first driver used the conjugated outer; exposed by the p < d
            # test case, where the inertia bound n_+ <= s1+s2+p is not
            # masked by the dimension cap n_+ <= d.)
            G += 2 * mj * np.real(np.outer(row, row))
            p += 1
            N += 2 * mj
    return G, dict(s1=s1, s2=s2, p=p, N=N)


def certificates(G: np.ndarray, g: dict, win: Window, N_I: int, N_Ip: int,
                 A_Ip: np.ndarray | None = None) -> dict:
    """Traces, the ratio C, the finite-T law of source eq. (7.2)/Sec. 8, the
    spectrum, and (if A_Ip, the form restricted to zeros in I', is given) the
    rank-trace and Cauchy-Schwarz certificates in units (4.4).
    Classes: every entry DIAGNOSTIC at finite T except the inequalities, which
    are exact statements about the supplied matrices."""
    L, a, b = g["L"], win.a, win.b
    Gt = G / L
    tr, tr2 = float(np.trace(Gt)), float(np.sum(Gt * Gt))
    C = tr * tr / tr2
    JT = 2 / L ** 3 * sum((lg ** 2 / nq) * float(win.g(np.log(nq)))
                          for nq, lg in prime_powers(g["X"]))
    law = g["lam1"] * a ** 2 / (b + g["lam1"] ** 2 * JT)
    out = dict(tr=tr, tr2=tr2, C=C, C_over_N=C / N_I, law=law, JT=JT,
               eigs=np.linalg.eigvalsh(Gt))
    if A_Ip is not None:
        Ah = A_Ip / (a * L ** 2)
        trA, frA = float(np.trace(Ah)), float(np.sum(Ah * Ah))
        out["cert_rank_trace"] = 4 * trA - 2 * N_Ip - frA
        out["cert_cs"] = 2 * trA ** 2 / frA - N_Ip
        out["trA_hat"] = trA
        out["frob_over_tr"] = frA / trA
    return out


# ----------------------------------------------------------------------
# the linear-algebra core (source Lemma 3.2) and its audit
# ----------------------------------------------------------------------
def rank_trace_gap(P, Q, r: int, b: int, c: float = 2.0) -> float:
    """||P+Q||_F^2 - [c tr P - c^2/4 r + 2c tr Q - c^2 b]; the lemma asserts
    this is >= 0 for Hermitian P >= 0 with rank <= r and Q with n_+(Q) <= b."""
    lhs = float(np.linalg.norm(P + Q, "fro") ** 2)
    return lhs - (c * float(np.trace(P)) - c * c / 4 * r
                  + 2 * c * float(np.trace(Q)) - c * c * b)


def audit_rank_trace(trials: int = 20000, seed: int = 143) -> dict:
    rng = np.random.default_rng(seed)
    viol, min_gap = 0, np.inf
    for t in range(trials):
        d = int(rng.integers(2, 9))
        r = int(rng.integers(0, d + 1))
        b = int(rng.integers(0, d + 1))
        P = np.zeros((d, d))
        if r:
            V = rng.standard_normal((d, r)) * rng.uniform(0.1, 3)
            P = V @ V.T
        Qw = rng.standard_normal((d, d)); Qw = (Qw + Qw.T) / 2
        ev, U = np.linalg.eigh(Qw)
        npos = int((ev > 0).sum())
        if npos > b:
            idx = np.argsort(ev)
            kill = idx[d - (npos - b):]
            ev[kill] = -np.abs(ev[kill]) * rng.uniform(0, 1, kill.size)
        Q = (U * ev) @ U.T
        c = 2.0 if t % 2 == 0 else float(rng.uniform(0.05, 4.0))
        gp = rank_trace_gap(P, Q, r, b, c)
        min_gap = min(min_gap, gp)
        if gp < -1e-9:
            viol += 1
    return dict(trials=trials, violations=viol, min_gap=min_gap)


def mt_constant(lam: float = 1.0) -> float:
    """Montgomery-Taylor c*_lam = sqrt2 tan(th)/(1+th tan th), th = lam/sqrt2.
    Certified proportion of the optimised certificate is 2 - 1/c*_lam. [Ext for
    the asymptotic meaning; the constant itself is a closed form.]"""
    th = lam / mp.sqrt(2)
    return float(mp.sqrt(2) * mp.tan(th) / (1 + th * mp.tan(th)))


# ----------------------------------------------------------------------
# CC-02 adjudicator
# ----------------------------------------------------------------------
def w_series_check(s: int = 3, N: int = 200_000) -> dict:
    """Sum w_n n^{-s} against three closed forms. Correct one [Pr]:
    -zeta(s)*zeta'(s+1) = F(s+1). The two printed corpus values (Paper 1
    Prop 1.5: -zeta'(s+1); AG Pr 4.1.4: -zeta'(s)/zeta(s-1)) are the series of
    the SUMMAND log(n)/n and a mistranscription respectively -- CC-02."""
    w = np.zeros(N + 1)
    for dd in range(2, N + 1):
        w[dd::dd] += np.log(dd) / dd
    n = np.arange(1, N + 1, dtype=float)
    S = float(np.sum(w[1:] / n ** s))
    mp.mp.dps = 30
    correct = float(-mp.zeta(s) * mp.zeta(s + 1, derivative=1))
    p1 = float(-mp.zeta(s + 1, derivative=1))
    ag = float(-mp.zeta(s, derivative=1) / mp.zeta(s - 1)) if s > 2 else float("nan")
    return dict(s=s, N=N, partial=S,
                correct=correct, diff_correct=abs(S - correct),
                paper1=p1, diff_paper1=abs(S - p1),
                ag=ag, diff_ag=abs(S - ag))


def selftest(verbose: bool = True) -> bool:
    ok = True
    r = w_series_check(3, 100_000)
    ok &= r["diff_correct"] < 1e-9 and r["diff_paper1"] > 1e-3
    a = audit_rank_trace(5000)
    ok &= a["violations"] == 0
    ok &= abs(mt_constant(1.0) - 0.753296067856) < 1e-9
    if verbose:
        print(f"weil.selftest: w-series diff {r['diff_correct']:.1e} (wrong forms "
              f"{r['diff_paper1']:.1e}/{r['diff_ag']:.1e}); rank-trace "
              f"{a['trials']} trials {a['violations']} viol; MT ok -> {ok}")
    return bool(ok)


if __name__ == "__main__":
    selftest()
