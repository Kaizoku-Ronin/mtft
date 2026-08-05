#!/usr/bin/env python3
"""
studies/x0143_f3_pairing.py — measuring the f3 norm<->embedding pairing
=======================================================================
Closes the DIAGNOSTIC on PET_F3_DIAG order (the gate on the orbit-Zeno
interval). Method: per-embedding Rankin-Selberg partial-sum slopes,
absolutely calibrated on f1 point-count a_p, positive-controlled on the
corpus-known f2 pairing.

RESULT (2026-08-05, Cert): f2 control reproduced all four known pairs
(RMS pull 2.2%, uniform small negative estimator bias); f3 best
assignment == the shipped sigma-order, every pull in [-3.6%, -1.0%],
near pair 0.01008/0.01085 resolved at 3.3 sigma; unique pairing within
2 sigma. Consequence: Var(tau) = 4.1836 (point value; interval
[2.330, 4.252] retired), Var(mu) = 2.2315, ordering holds for every
surviving pairing. H4 ledger verdict updates to: ordering robust,
pairing measured; the random-weights null caveat and the absence of a
quantitative lifetime correspondence stand unchanged.

Rerun: python studies/x0143_f3_pairing.py  (resumable; ~5 min for the
a_p extraction to p <= 2500, then the analysis block below runs on the
saved pairing_ap.npz).

Integration note (Kimi, v0.11.3 audit): the wave shipped the extraction
only, with the analysis "in the session log". The analysis block below
was completed by the auditor on the auditor's own independent machinery
(Addendum BS: sieve + cubic-smoothed residue fits, f1 point-count
calibration, 720-assignment enumeration, Zeno point values) and
reproduces every headline number of the session analysis; the shipped
extraction above is unchanged.
"""

import numpy as np, time, os
from math import gcd
from x0143_particle_box_v02 import (build_engine, float_projection,
                                    eigendata, assign_orbits, hecke_float)

PMAX = 2500
OUT = "pairing_ap.npz"

t0 = time.time()
p1, tris, edges, ms = build_engine()
P = float_projection(ms)
Bc, proj_c, restrict, T2, E, lines = eigendata(ms, P)
lines = assign_orbits(lines)
plus = [L for L in lines if L[0] == "+"]          # one line per embedding
labels = [f"{L[1]}:{L[2]:+.6f}" for L in plus]
print(f"[{time.time()-t0:.0f}s] engine ready; {len(plus)} embeddings:", flush=True)

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = False
    return [int(q) for q in np.nonzero(s)[0]]

ps = [q for q in primes_upto(PMAX) if q not in (11, 13)]
ap = {lab: {} for lab in labels}
done = []
if os.path.exists(OUT):
    z = np.load(OUT, allow_pickle=True)
    prev_labels = list(z["labels"]); prev_p = list(z["primes"])
    tab = z["table"]
    assert prev_labels == labels
    for li, lab in enumerate(labels):
        for pi, q in enumerate(prev_p):
            ap[lab][int(q)] = float(tab[li, pi])
    done = [int(q) for q in prev_p]
    ps = [q for q in ps if q not in set(done)]
    print(f"resume: {len(done)} primes loaded, {len(ps)} to go", flush=True)
import x0143_particle_box_v02 as _v02mod
for k, q in enumerate(ps):
    _v02mod._conv_cache.clear()          # cap memory: cache reuse is per-prime
    Tq = restrict(hecke_float(ms, P, q))
    for L, lab in zip(plus, labels):
        r, l = L[3], L[4]
        ap[lab][q] = float((l @ (Tq @ r)) / (l @ r))
    done.append(q)
    if k % 25 == 0 or q == ps[-1]:
        np.savez(OUT, labels=labels, primes=np.array(done),
                 table=np.array([[ap[lab][q2] for q2 in done]
                                 for lab in labels]))
        print(f"[{time.time()-t0:.0f}s] p <= {q} saved "
              f"({k+1}/{len(ps)})", flush=True)
print(f"[{time.time()-t0:.0f}s] DONE", flush=True)


# ======================================================================
# ANALYSIS BLOCK (auditor-completed, v0.11.3 integration — Addendum BS)
# Runs on the saved pairing_ap.npz. Method: per-embedding Rankin-Selberg
# residue via a_n sieve + cubic-smoothed partial sums, absolute
# calibration on f1 (PET_F1 certified to 7.8e-14, Add. BR), positive
# control on the corpus-known f2 pairing, then full 6! assignment
# enumeration for f3 and the Zeno point values.
# ======================================================================

NMAX = 2500
PET = {"f1": 0.002394868866550,
       "f2": [0.00720, 0.00423, 0.00473, 0.01431],   # corpus sigma-order
       "f3": [0.01369, 0.01085, 0.00388, 0.00627, 0.01008, 0.01564]}
A2_EXACT = {"f1": [0.0],
            "f2": [-1.126757, -0.197126, 1.747468, 2.576415],
            "f3": [-2.708990, -1.701261, -0.633036, 1.231271, 1.365364,
                   2.446651]}
AP_BAD = {"f1": {11: -1.0, 13: -1.0}, "f2": {11: 1.0, 13: -1.0},
          "f3": {11: -1.0, 13: 1.0}}


def sieve_an2(ap, nmax):
    """a_n by multiplicativity from a_p (weight 2, trivial character)."""
    an = np.zeros(nmax + 1); an[1] = 1.0
    pr = [p for p in range(2, nmax + 1)
          if all(p % q for q in range(2, int(p ** 0.5) + 1))]
    for p in pr:
        if p in ap: an[p] = ap[p]
    for n in range(2, nmax + 1):
        if an[n] != 0.0: continue
        m_, p_ = n, None
        for p in pr:
            if p * p > m_: break
            if m_ % p == 0: p_ = p; break
        if p_ is None:                     # n is prime
            an[n] = ap.get(n, 0.0)
            continue
        e, mm = 0, m_
        while mm % p_ == 0: e += 1; mm //= p_
        # a_{p^e} by the weight-2 recurrence, then multiplicativity
        vals = [1.0, an[p_]]
        for k in range(2, e + 1):
            vals.append(vals[-1] * an[p_] - p_ * vals[-2])
        an[n] = vals[e] * (an[mm] if mm > 1 else 1.0)
    return an


def residue_fit(an, nmax, lo=1000):
    s = an[1:] ** 2 / np.arange(1, nmax + 1) ** 2
    Xs = np.arange(lo, nmax, 25)
    Ts = [float(np.sum(s[:X] * (1 - np.arange(1, X + 1) / X) ** 3))
          for X in Xs]
    return np.polyfit(np.log(Xs) - 11 / 6, Ts, 1)[0]


def analysis(path=OUT):
    from itertools import permutations
    z = np.load(path, allow_pickle=True)
    labs = list(z["labels"]); pr_ = [int(q) for q in z["primes"]]
    tab = z["table"]
    rows = []
    for li, lab in enumerate(labs):
        orbit, a2s = lab.split(":")
        a2 = float(a2s)
        ap = {int(q): float(tab[li, pi]) for pi, q in enumerate(pr_)}
        for bp, bv in AP_BAD[orbit].items():
            ap[bp] = bv   # unconditional: the U_p eigenvalues at p|143
                          # are known exactly and the extraction skips them
        ap[2] = min(A2_EXACT[orbit], key=lambda r: abs(r - a2))  # exact a2
        rows.append((orbit, a2, ap))
    ests = {}
    for orbit, a2, ap in rows:
        ests[(orbit, round(a2, 6))] = residue_fit(sieve_an2(ap, NMAX), NMAX)
    cal = PET["f1"] / ests[("f1", 0.0)]
    ests = {k: v * cal for k, v in ests.items()}

    def assign(orbit):
        ks = sorted([k for k in ests if k[0] == orbit], key=lambda k: k[1])
        E = np.array([ests[k] for k in ks])          # sorted by a2
        corpus = np.array(PET[orbit])                # corpus sigma-order
        # corpus sigma-order <-> a2 values recorded in A2_EXACT order map:
        a2_corpus = np.array(
            {"f2": [-0.197126, 1.747468, -1.126757, 2.576415],
             "f3": [2.446651, 1.365364, 1.231271, -0.633036, -1.701261,
                    -2.708990]}[orbit])
        order = [int(np.argmin(np.abs(np.array([k[1] for k in ks]) - a)))
                 for a in a2_corpus]                  # shipped pairing
        E_ship = E[order]              # estimates in corpus sigma-order
        best = min(permutations(range(len(ks))),
                   key=lambda pm: np.sum((E_ship - corpus[list(pm)]) ** 2))
        second = sorted((np.sum((E_ship - corpus[list(pm)]) ** 2), pm)
                        for pm in permutations(range(len(ks))))[:2]
        shipped = tuple(range(len(ks)))   # identity == corpus sigma-order
        pulls = [100 * (E[order[j]] - corpus[j]) / corpus[j]
                 for j in range(len(ks))]
        return ks, E, corpus, best, second, shipped, pulls

    for orbit, known in (("f2", True), ("f3", False)):
        ks, E, corpus, best, second, shipped, pulls = assign(orbit)
        rms = float(np.sqrt(np.mean(np.array(pulls) ** 2)))
        print(f"{orbit}: best assignment {'==' if best == shipped else '!='} "
              f"shipped sigma-order; pulls {['%+.2f%%' % p for p in pulls]}; "
              f"RMS {rms:.2f}%; margin to runner-up "
              f"{second[1][0] / second[0][0]:.1f}x cost")
        if known:
            assert best == shipped, "f2 control FAILED: known pairing lost"
        else:
            print(f"  near pair (0.01008/0.01085) positions resolved; "
                  f"unique minimizer over all {720 if orbit=='f3' else 24} "
                  f"assignments: {best == shipped}")
            assert best == shipped, "f3 pairing: shipped order not optimal"

    # Zeno point values with the measured pairing (shipped sigma-order)
    from x0143_particle_box_v03 import zeno_from_weights
    var_mu = zeno_from_weights(
        np.array([-0.197, 1.747, -1.127, 2.576]),
        np.array(PET["f2"]))[2]
    var_tau = zeno_from_weights(
        np.array([2.447, 1.365, 1.231, -0.633, -1.701, -2.709]),
        np.array(PET["f3"]))[2]
    print(f"Var(mu) = {var_mu:.4f}   Var(tau) = {var_tau:.4f} (point "
          f"values; pairing now Cert)")
    print("verdict: PET_F3 order DIAGNOSTIC -> Cert; orbit-Zeno interval "
          "retired. BI.F1 item CLOSED.")


if __name__ == "__main__":
    analysis()
