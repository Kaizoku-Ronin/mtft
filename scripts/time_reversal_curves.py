"""Provenance script for TR-01: the thirteen star quotients X0(N)* of SQ-01 as curves, from PARI/GP's elliptic-curve database.

Requires cypari2 with the elldata package.  For each torus level N = pq the (+,+) isogeny class is the one whose curves have a_p = a_q = -1
(w_p = -a_p for p || N); in it the curve of modular degree 4 = |<W_p, W_q>| is the quotient.  Also recorded: torsion, root number, rank, a_2,
and Lambda(1 + i) (purely imaginary when the root number is -1).  Frozen as STAR_QUOTIENT_CURVES in mtft.matter.time_reversal.

    python scripts/time_reversal_curves.py            # JSON on stdout
    python scripts/time_reversal_curves.py --check    # compare with the frozen table; exit status 1 on any difference"""
import argparse
import json
import sys

import cypari2

TORI = (57, 58, 65, 77, 82, 91, 118, 123, 141, 142, 143, 145, 155)
SS = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
    out = {}
    for N in TORI:
        pr = [x for x in SS if N % x == 0]
        names = [str(r[0]) for r in pari.ellsearch(N)
                 if all(int(pari.ellap(pari.ellinit([int(x) for x in r[1]]), q)) == -1 for q in pr)]
        classes = sorted({n.rstrip("0123456789") for n in names})
        assert len(classes) == 1, (N, names)
        moddeg = {n: int(pari.ellmoddegree(pari(f'ellinit("{n}")'))) for n in names}
        quotient = [n for n, d in moddeg.items() if d == 4]
        assert len(quotient) == 1, (N, moddeg)
        E = pari(f'ellinit("{quotient[0]}")')
        L = pari.lfuncreate(E); lam = pari.lfunlambda(L, pari("1+I"))
        out[N] = {"curve": quotient[0], "a": [int(E[i]) for i in range(5)], "class_moddeg": moddeg,
                  "torsion": int(pari.elltors(E)[0]), "root_number": int(pari.ellrootno(E)), "rank": int(pari.ellrank(E)[0]),
                  "a2": int(pari.ellap(E, 2)) if N % 2 else None, "disc": int(E[11]),
                  "Lambda_1_plus_i": [float(pari.real(lam)), float(pari.imag(lam))]}
    print(json.dumps({"pari": ".".join(map(str, pari.version())), "curves": out}))
    if a.check:
        from mtft.matter import time_reversal as T
        bad = []
        for N, c in out.items():
            f = T.STAR_QUOTIENT_CURVES[N]
            if (c["curve"], tuple(c["a"]), c["class_moddeg"], c["torsion"], c["root_number"], c["rank"], c["a2"]) != \
               (f["curve"], f["a"], f["class_moddeg"], f["torsion"], f["root_number"], f["rank"], f["a2"]) or c["disc"] != T.discriminant(f["a"]):
                bad.append(N)
            if abs(c["Lambda_1_plus_i"][0]) > 1e-12:
                bad.append((N, "Lambda not imaginary"))
        print("differences from the frozen table:", bad or "none", file=sys.stderr)
        return 1 if bad else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
