"""Monte Carlo for the repaired SU(2) action on a small 4D lattice: where is the bulk transition in the (beta, kappa) plane?

s(U_p) = beta (1 - Tr U_p / 2) + (kappa / 2) sum_{n>=2} a_n (n + 1 - chi_{n/2}(U_p)),  chi_{n/2} = U_n(Tr U_p / 2) (Chebyshev, second kind).
Links are SU(2) quaternions; Metropolis with a near-identity proposal, checkerboard over (direction, site parity).  Hysteresis scan:
beta up from a hot start and down from a cold start at fixed kappa; a first-order bulk transition shows as a jump with hysteresis.
The Sym^2 coupling in Bhanot–Creutz's convention is beta_A = 3 kappa a_2 / 2 = 0.0530 kappa.
usage: python -I bulk_mc_su2.py L kappa out.json [sweeps] [betas...]
"""
import json
import math
import sys
import time

import numpy as np

L = int(sys.argv[1]); KAPPA = float(sys.argv[2]); OUT = sys.argv[3]
SWEEPS = int(sys.argv[4]) if len(sys.argv) > 4 else 120
BETAS = [float(b) for b in sys.argv[5:]] or list(np.round(np.arange(0.0, 3.01, 0.125), 3))
Y, NMAX = 0.18174, 20
rng = np.random.default_rng(143)


def weights(nmax):
    return np.array([0.0] + [sum(math.log(d) / d for d in range(1, n + 1) if n % d == 0) * math.exp(-2 * math.pi * Y * n) for n in range(1, nmax + 1)])
A = weights(NMAX)


def qmul(p, q):
    """Quaternion product, arrays (..., 4): U = a0 + i a.sigma."""
    a0, a1, a2, a3 = p[..., 0], p[..., 1], p[..., 2], p[..., 3]
    b0, b1, b2, b3 = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 - a2 * b3 + a3 * b2,
                     a0 * b2 + a2 * b0 - a3 * b1 + a1 * b3, a0 * b3 + a3 * b0 - a1 * b2 + a2 * b1], axis=-1)


def qconj(p):
    return p * np.array([1.0, -1.0, -1.0, -1.0])


def plaquette_action(c):
    """s as a function of c = Tr U_p / 2 (array); the tower via Chebyshev U_n(c)."""
    s = (1.0 - c)
    u0, u1 = np.ones_like(c), 2 * c
    tot = np.zeros_like(c)
    for n in range(2, NMAX + 1):
        u2 = 2 * c * u1 - u0
        tot += A[n] * ((n + 1) - u2)
        u0, u1 = u1, u2
    return s, tot                      # beta * s + (kappa / 2) * tot


def shift(arr, mu, sgn):
    return np.roll(arr, -sgn, axis=mu)


def staples(U, mu):
    """The six staples S_i for links in direction mu, shape (6, L, L, L, L, 4): U_mu(x) S_i^dagger... we return the matrices such that the
    plaquette is U_mu(x) @ S_i."""
    out = []
    for nu in range(4):
        if nu == mu:
            continue
        # forward: U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
        a = shift(U[nu], mu, +1); b = qconj(shift(U[mu], nu, +1)); c = qconj(U[nu])
        out.append(qmul(qmul(a, b), c))
        # backward: U_nu(x+mu-nu)^dag U_mu(x-nu)^dag U_nu(x-nu)
        a = qconj(shift(shift(U[nu], mu, +1), nu, -1)); b = qconj(shift(U[mu], nu, -1)); c = shift(U[nu], nu, -1)
        out.append(qmul(qmul(a, b), c))
    return np.stack(out)


def random_near_identity(shape, eps):
    v = rng.standard_normal(shape + (3,)) * eps
    ang = np.linalg.norm(v, axis=-1, keepdims=True)
    axis = v / np.maximum(ang, 1e-300)
    return np.concatenate([np.cos(ang), np.sin(ang) * axis], axis=-1)


PAR = (np.indices((L, L, L, L)).sum(0) % 2)


def sweep(U, beta, kappa, eps, hits=2):
    acc = 0; tot = 0
    for mu in range(4):
        S = staples(U, mu)                                                   # (6, ..., 4)
        for parity in (0, 1):
            m = PAR == parity
            Um = U[mu][m]                                                    # (n, 4)
            Sm = S[:, m]                                                     # (6, n, 4)
            for _ in range(hits):
                R = random_near_identity(Um.shape[:-1], eps)
                Un = qmul(R, Um)
                c_old = qmul(Um[None], Sm)[..., 0]; c_new = qmul(Un[None], Sm)[..., 0]
                so, to = plaquette_action(c_old); sn, tn = plaquette_action(c_new)
                dS = (beta * (sn - so) + 0.5 * kappa * (tn - to)).sum(0)
                ok = rng.random(dS.shape) < np.exp(-np.maximum(dS, -50))
                Um = np.where(ok[:, None], Un, Um); acc += ok.sum(); tot += ok.size
            U[mu][m] = Um
    return acc / tot


def avg_plaquette(U):
    tot = 0.0; n = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            a = U[mu]; b = shift(U[nu], mu, +1); c = qconj(shift(U[mu], nu, +1)); d = qconj(U[nu])
            tot += qmul(qmul(a, b), qmul(c, d))[..., 0].mean(); n += 1
    return tot / n


def run_scan(betas, kappa, start):
    U = np.zeros((4, L, L, L, L, 4))
    if start == "cold":
        U[..., 0] = 1.0
    else:
        v = rng.standard_normal((4, L, L, L, L, 4)); U = v / np.linalg.norm(v, axis=-1, keepdims=True)
    eps = 0.5; rows = []
    for beta in betas:
        meas = []
        for it in range(SWEEPS):
            a = sweep(U, beta, kappa, eps)
            eps = min(1.5, max(0.02, eps * (1.1 if a > 0.55 else 0.9 if a < 0.4 else 1.0)))
            if it >= SWEEPS // 3:
                meas.append(avg_plaquette(U))
        meas = np.array(meas)
        rows.append({"beta": float(beta), "P": float(meas.mean()), "err": float(meas.std() / math.sqrt(len(meas))), "eps": eps})
        print(f"kappa={kappa:6.1f} {start:4s} beta={beta:5.3f} P={meas.mean():.4f} +- {meas.std() / math.sqrt(len(meas)):.4f} eps={eps:.3f}", flush=True)
    return rows


t0 = time.time()
up = run_scan(BETAS, KAPPA, "hot"); down = run_scan(BETAS[::-1], KAPPA, "cold")
json.dump({"L": L, "kappa": KAPPA, "sweeps": SWEEPS, "beta_A_equivalent": 1.5 * KAPPA * A[2], "hot_up": up, "cold_down": down,
           "seconds": time.time() - t0}, open(OUT, "w"), indent=1)
print("done", time.time() - t0, "s")
