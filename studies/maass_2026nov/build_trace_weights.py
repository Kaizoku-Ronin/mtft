"""Recompute the geometric weights of the trace formulas of mtft.maass.trace, for lengths up to 17 (2026-11-15).

Needs cypari2 (PARI's quadclassunit gives the class number and regulator of a real quadratic field; in general that is
conditional on GRH, and mtft.maass.trace checks a sample against Dirichlet's finite sum, which is not).

For a level N (1, 11, 143), a divisor Q of N and a length bound LMAX, the Gamma_0(N)-classes of the coset W_Q Gamma_0(N)
(and of R W_Q Gamma_0(N), R the mirror) are counted through binary quadratic forms, as in the lab's record MP16:
  hyperbolic, det +Q, trace Q s, 2 cosh(l/2) = s sqrt(Q):  weight 2 sqrt(Q) W(Q^2 s^2 - 4 Q)
  glide,      det -Q, trace Q s, 2 sinh(l/2) = s sqrt(Q):  weight 2 sqrt(Q) W(Q^2 s^2 + 4 Q)
  W(Delta) = sum_{f | c0} f^-1 L(1, chi_D) prod_{p | c0/f} (1 - chi_D(p)/p) prod_{p | N/Q} rho_p,   Delta = D c0^2,
  rho_p = p + 1 if p | f, else 1 + chi_D(p) (taken as 1 when p | c0/f).
For the eta-twist at level 143 (record MP18) the term of index f changes sign when the trace is even and f is odd.

Usage, from the repository root:
    python studies/maass_2026nov/build_trace_weights.py [folder holding the lab's weights_*.npz, to compare with]
Output: src/mtft/maass/_data/trace_weights.npz, and its line in PROVENANCE.txt.  About a minute.
"""
import hashlib
import math
import sys
import time
from pathlib import Path

import cypari2
import numpy as np

pari = cypari2.Pari()
HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "src" / "mtft" / "maass" / "_data"
LMAX = 17.0
_L1 = {}


def primes_of(n):
    return [int(p) for p in pari.factor(n)[0]] if n > 1 else []


def divisors_from(ps_exp):
    divs = [1]
    for p, e in ps_exp:
        divs = [d * p ** j for d in divs for j in range(e + 1)]
    return divs


def kron(a, b):
    return int(pari.kronecker(a, b))


def L1(D):
    """L(1, chi_D) for a fundamental discriminant D > 0: 2 h log(eps) / sqrt(D)."""
    if D not in _L1:
        q = pari.quadclassunit(D)
        _L1[D] = 2 * int(q[0]) * float(q[3]) / math.sqrt(D)
    return _L1[D]


def W_of(Delta, N, Q, twist_trace=None):
    """W(Delta) for the coset W_Q at level N.  With twist_trace = t (level 143, Q = 1) the eta-twisted weight."""
    cd = pari.coredisc(Delta, 1)
    D, c0 = int(cd[0]), int(cd[1])
    L = L1(D)
    fac = [(int(p), int(e)) for p, e in zip(*pari.factor(c0))] if c0 > 1 else []
    tot = 0.0
    for f in divisors_from(fac):
        c = c0 // f
        term = L / f
        for p, _ in fac:
            if c % p == 0:
                term *= 1 - kron(D, p) / p
        for p in primes_of(N // Q):
            term *= (p + 1) if f % p == 0 else (1 + (0 if c % p == 0 else kron(D, p)))
        if twist_trace is not None and twist_trace % 2 == 0 and f % 2 == 1:
            term = -term
        tot += term
    return tot


def h_neg(D):
    """Number of reduced primitive positive definite forms of discriminant D < 0."""
    cnt, a = 0, 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (b < 0 and a == c) or math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            cnt += 1
        a += 1
    return cnt


def elliptic_count(N, Q):
    """Fixed points of W_Q on X0(N), each weighted by one over half the order of its stabiliser's centre-free part."""
    tot = 0.0
    for f in (1, 2):
        if (4 * Q) % (f * f):
            continue
        Dp = -4 * Q // (f * f)
        if Dp % 4 not in (0, 1):
            continue
        w = {-4: 2, -3: 3}.get(Dp, 1)
        rho = 1
        for p in primes_of(N // Q):
            rho *= (p + 1) if f % p == 0 else 1 + kron(Dp, p)
        tot += h_neg(Dp) * rho / w
    return tot


def classes(N, Q, lmax, twist=False):
    hyp, gli = [], []
    s = 1
    while True:
        x = s * math.sqrt(Q) / 2
        l_h = 2 * math.acosh(x) if x > 1 else None
        l_g = 2 * math.asinh(x)
        if l_g > lmax and (l_h is None or l_h > lmax):
            break
        if l_h is not None and l_h <= lmax:
            hyp.append((s, 2 * math.sqrt(Q) * W_of(Q * Q * s * s - 4 * Q, N, Q, s if twist else None)))
        if l_g <= lmax and not twist:
            gli.append((s, 2 * math.sqrt(Q) * W_of(Q * Q * s * s + 4 * Q, N, Q)))
        s += 1
    return hyp, gli


def main():
    t0 = time.time()
    data = {"lmax": LMAX}
    for N in (1, 11, 143):
        ps = primes_of(N)
        Qs = sorted({math.prod(ps[i] for i in range(len(ps)) if k >> i & 1) for k in range(2 ** len(ps))})
        for Q in Qs:
            hyp, gli = classes(N, Q, LMAX)
            key = f"N{N}_Q{Q}"
            data[key + "_hyp_s"] = np.array([r[0] for r in hyp], dtype=np.int64)
            data[key + "_hyp_w"] = np.array([r[1] for r in hyp])
            data[key + "_gli_s"] = np.array([r[0] for r in gli], dtype=np.int64)
            data[key + "_gli_w"] = np.array([r[1] for r in gli])
            if Q > 1:
                data[key + "_E"] = elliptic_count(N, Q)
                data[key + "_L"] = 2 * math.sqrt(Q) * W_of(4 * Q, N, Q)
            print(f"{key}: {len(hyp)} hyperbolic, {len(gli)} glide  [{time.time() - t0:.0f} s]", flush=True)
    tw, _ = classes(143, 1, LMAX, twist=True)
    data["N143_tw_hyp_s"] = np.array([r[0] for r in tw], dtype=np.int64)
    data["N143_tw_hyp_w"] = np.array([r[1] for r in tw])
    print(f"twisted: {len(tw)} hyperbolic; {len(_L1)} fundamental discriminants in all  [{time.time() - t0:.0f} s]", flush=True)

    if len(sys.argv) > 1:                                             # compare with the lab's stored tables
        lab = Path(sys.argv[1])
        worst = 0.0
        where = {(143, 1): "mp17/weights_143_1.npz", (143, 11): "mp16/weights_143_11.npz", (143, 13): "mp16/weights_143_13.npz",
                 (143, 143): "mp16/weights_143_143.npz", (11, 1): "mp17/weights_11_1.npz", (11, 11): "mp16/weights_11_11.npz", (1, 1): "mp19/weights_1_1.npz"}
        for (N, Q), name in where.items():
            z = np.load(lab / name)
            for kind in ("hyp", "gli"):
                if kind not in z.files:
                    continue
                a = z[kind]
                a = a[a[:, 1] <= LMAX]
                s, w = data[f"N{N}_Q{Q}_{kind}_s"], data[f"N{N}_Q{Q}_{kind}_w"]
                assert a.shape[0] == s.size and np.all(a[:, 0].astype(np.int64) == s), (N, Q, kind, a.shape, s.size)
                worst = max(worst, float(np.max(np.abs(a[:, 2] - w) / np.maximum(1.0, np.abs(w)))))
            if "E" in z.files:
                worst = max(worst, abs(float(z["E"]) - data[f"N{N}_Q{Q}_E"]), abs(float(z["L"]) - data[f"N{N}_Q{Q}_L"]))
        z = np.load(lab / "mp18/weights_tw.npz")["hyp"]
        z = z[z[:, 1] <= LMAX]
        assert z.shape[0] == data["N143_tw_hyp_s"].size
        worst = max(worst, float(np.max(np.abs(z[:, 2] - data["N143_tw_hyp_w"]) / np.maximum(1.0, np.abs(z[:, 2])))))
        print(f"against the lab's stored tables (lengths up to {LMAX}): largest relative difference {worst:.2e}")

    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "trace_weights.npz", **data)
    sha = hashlib.sha1((OUT / "trace_weights.npz").read_bytes()).hexdigest()
    prov = OUT / "PROVENANCE.txt"
    lines = [ln for ln in prov.read_text(encoding="utf-8").splitlines() if "trace_weights" not in ln] if prov.exists() else []
    while lines and not lines[-1].strip():
        lines.pop()
    lines += ["", f"{sha}  trace_weights.npz", "built by studies/maass_2026nov/build_trace_weights.py (cypari2; lengths up to 17)", ""]
    prov.write_text("\n".join(lines), encoding="utf-8")
    print(f"written {OUT / 'trace_weights.npz'} ({(OUT / 'trace_weights.npz').stat().st_size} bytes), sha1 {sha}; {time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
