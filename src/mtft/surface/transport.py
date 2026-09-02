"""mtft.surface.transport — EXACT Hecke / Atkin-Lehner transport onto the cycle lattice.

Route A: each representative dart (c:d) is lifted to g in SL2(Z); its
oriented edge is the modular symbol {g.oo, g.0}.  The cycle basis is pushed
through PARI ``mspathlog`` (weight 2, Gamma acts trivially).  T_p acts on
paths by the p+1 coset matrices; W_Q by an explicit det-Q matrix.  Solving
back in the cycle basis gives integer matrices on H_1(X0(N), Z).

Route B: PARI ``mshecke`` / ``msatkinlehner`` on the cuspidal symbol space
(dual), sharing no computational step with route A.  Gates: integrality,
charpoly(T_p) divides charpoly(mshecke(p)) with Eisenstein cofactor of
degree cusps-1, AL involutions and commutation, AL traces equal route B,
intersection form preserved by W_Q and T_p self-adjoint for it.

The same GP job evaluates ``mfsymbol`` of the rational cuspidal basis on
the same edge paths, giving the g x 2g period matrix over the cycles used
by :mod:`hodge_structure`.  Requires PARI/GP (``mtft.gprun.find_gp``).
"""
from __future__ import annotations

import re
import subprocess
import tempfile
from dataclasses import dataclass
from math import gcd
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp

from ..gprun import find_gp
from .cycles import CycleBasis
from .manin import ManinComplex, factorize


def egcd(a: int, b: int) -> Tuple[int, int, int]:
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def lift_dart(c: int, d: int, N: int) -> Tuple[int, int, int, int]:
    """(a, b, c', d') in SL2(Z) with (c', d') == (c, d) mod N."""
    if c % N == 0:
        return 1, 0, 0, 1
    if d % N == 0:
        return 0, -1, 1, 0
    cc = c % N
    for k in range(4 * N):
        dd = d % N + k * N
        if gcd(cc, dd) == 1:
            _, x, y = egcd(cc, dd)      # x cc + y dd = 1
            return y, -x, cc, dd
    raise ValueError(f"no coprime lift of ({c},{d}) mod {N}")


def cusp_class(den: int, N: int) -> int:
    return N if den == 0 else gcd(den, N)


def atkin_lehner_matrix(Q: int, N: int) -> Tuple[int, int, int, int]:
    """[[Q a, b], [N c, Q d]] with determinant Q, for Q exactly dividing N."""
    R = N // Q
    if gcd(Q, R) != 1:
        raise ValueError("Q must exactly divide N")
    _, a, b = egcd(Q, R)        # a Q + b R = 1  ->  Q a - R (-b) = 1
    return Q * a, -b, N, Q       # det = Q^2 a - N(-b)... check below
    

def _al_check(Q: int, N: int):
    A, B, C, D = atkin_lehner_matrix(Q, N)
    assert A * D - B * C == Q, (Q, N, A, B, C, D)
    return A, B, C, D


@dataclass
class Transport:
    N: int
    hecke: Dict[int, np.ndarray]           # p -> 2g x 2g integer matrix on cycles
    atkin_lehner: Dict[int, np.ndarray]    # Q -> matrix
    charpoly_route_A: Dict[int, sp.Poly]
    charpoly_route_B: Dict[int, sp.Poly]
    al_trace_route_B: Dict[int, int]
    period_Q: np.ndarray                   # 2g x 2g real [Re P; Im P]
    gates: List[dict]
    gp_stdout: str


def gp_script(cx: ManinComplex, cb: CycleBasis, primes: List[int], als: List[int], prec: int = 30) -> str:
    N = cx.N
    paths, bad = [], 0
    for e, d in enumerate(cx.edges):
        a, b, c, dd = lift_dart(d[0], d[1], N)
        if cusp_class(c, N) != cx.cusp_divisor_class[cx.vertex_of[d]] or \
           cusp_class(dd, N) != cx.cusp_divisor_class[cx.vertex_of[cx.S(d)]]:
            bad += 1
        paths.append((a, b, c, dd))
    if bad:
        raise ArithmeticError(f"cusp-class gate failed on {bad} edges")
    B = cb.basis_matrix
    lines = [f"default(parisizemax, 1500000000); default(realprecision, {prec});",
             f"N={N}; NE={len(paths)}; NC={B.shape[1]};",
             "PATHS=[" + ",".join(f"[{a},{b},{c},{d}]" for a, b, c, d in paths) + "];",
             "CYC=[" + ";".join(",".join(str(int(v)) for v in B[:, j]) for j in range(B.shape[1])) + "];",
             "M=msinit(N,2);",
             "DIM=msdim(M);",
             "cusp(n,d)=if(d==0,oo,n/d);",
             "pathof(P)=[cusp(P[1],P[3]),cusp(P[2],P[4])];",
             "pathvec(p)={my(t=mspathlog(M,p)); vector(DIM,i,if(t[i]==0,0,sum(j=1,matsize(t[i])[1],t[i][j,2])));}",
             "act(g,x)={my(n,d); if(x==oo, n=g[1,1]; d=g[2,1], n=g[1,1]*numerator(x)+g[1,2]*denominator(x); d=g[2,1]*numerator(x)+g[2,2]*denominator(x)); if(d==0,oo,n/d);}",
             "actpath(g,p)=[act(g,p[1]),act(g,p[2])];",
             "heckepath(p,path)={my(v=pathvec(actpath([p,0;0,1],path))); for(j=0,p-1,v+=pathvec(actpath([1,j;0,p],path))); v;}",
             "EV=vector(NE,e,pathvec(pathof(PATHS[e])));",
             "C=matrix(DIM,NC,i,j,sum(e=1,NE,CYC[j,e]*EV[e][i]));",
             "print(\"RANK=\",matrank(C));",
             "opmatrix(f)={my(TC=matrix(DIM,NC)); for(j=1,NC,my(v=vector(DIM)); for(e=1,NE,if(CYC[j,e]!=0,v+=CYC[j,e]*f(pathof(PATHS[e])))); TC[,j]=v~); matinverseimage(C,TC);}",
             "isint(A)={my(ok=1); for(i=1,matsize(A)[1],for(j=1,matsize(A)[2],if(denominator(A[i,j])!=1,ok=0))); ok;}",
             "H=mscuspidal(M)[1];"]
    for p in primes:
        lines += [f"T{p}=opmatrix(p->heckepath({p},p));",
                  f"print(\"INT_T{p}=\",isint(T{p}));",
                  f"print(\"CPA_{p}=\",charpoly(T{p}));",
                  f"print(\"CPB_{p}=\",charpoly(mshecke(M,{p})));",
                  f"print(\"MAT_T{p}=\",T{p});"]
    for Q in als:
        A, Bq, Cq, D = _al_check(Q, N)
        lines += [f"W{Q}=opmatrix(p->pathvec(actpath([{A},{Bq};{Cq},{D}],p)));",
                  f"print(\"INT_W{Q}=\",isint(W{Q}));",
                  f"print(\"TRB_W{Q}=\",trace(msatkinlehner(M,{Q},H)));",
                  f"print(\"MAT_W{Q}=\",W{Q});"]
    lines += ["mf=mfinit([N,2],1); BB=mfbasis(mf); g=#BB; P=matrix(g,NC);",
              "for(i=1,g,my(FS=mfsymbol(mf,BB[i]),pe=vector(NE)); for(e=1,NE,pe[e]=mfsymboleval(FS,pathof(PATHS[e]))); for(j=1,NC,P[i,j]=sum(e=1,NE,CYC[j,e]*pe[e])));",
              "Q=matrix(2*g,NC,i,j,if(i<=g,real(P[i,j]),imag(P[i-g,j])));",
              "print(\"MAT_Q=\",Q);",
              "print(\"DONE\");"]
    return "\n".join(lines) + "\n"


def _parse_matrix(text: str) -> np.ndarray:
    rows = text.strip().strip("[]").split(";")
    return np.array([[float(sp.Rational(x)) if "/" in x else float(x) for x in r.replace(" ", "").split(",")] for r in rows])


def _parse_int_matrix(text: str) -> np.ndarray:
    rows = text.strip().strip("[]").split(";")
    return np.array([[int(x) for x in r.replace(" ", "").split(",")] for r in rows], dtype=np.int64)


def run(cx: ManinComplex, cb: CycleBasis, Jint: np.ndarray, primes: Tuple[int, ...] = (2, 3),
        timeout: int = 600, prec: int = 30) -> Transport:
    gp = find_gp()
    if gp is None:
        raise RuntimeError("PARI/GP not found (set MTFT_GP)")
    N = cx.N
    primes = tuple(p for p in primes if N % p)
    als = [p ** a for p, a in factorize(N).items()]
    if len(als) > 1:
        als.append(N)
    script = gp_script(cx, cb, list(primes), als, prec)
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / f"surface_{N}.gp"
        path.write_text(script)
        proc = subprocess.run([gp, "-q", str(path)], capture_output=True, text=True, timeout=timeout)
    out = "\n".join(l for l in proc.stdout.splitlines() if "Warning" not in l)
    if "DONE" not in out:
        raise RuntimeError(f"GP job failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")

    def grab(key):
        m = re.search(rf"^{key}=(.*)$", out, re.M)
        if not m:
            raise RuntimeError(f"{key} missing from GP output")
        return m.group(1)

    x = sp.Symbol("x")
    hecke, cpa, cpb = {}, {}, {}
    gates = [{"name": "cycle_image_rank_2g", "status": "PASS" if int(grab("RANK")) == cb.rank else "FAIL"}]
    for p in primes:
        hecke[p] = _parse_int_matrix(grab(f"MAT_T{p}"))
        cpa[p] = sp.Poly(sp.sympify(grab(f"CPA_{p}")), x)
        cpb[p] = sp.Poly(sp.sympify(grab(f"CPB_{p}")), x)
        q, r = sp.div(cpb[p], cpa[p])
        gates.append({"name": f"T{p}_integral", "status": "PASS" if grab(f"INT_T{p}") == "1" else "FAIL"})
        gates.append({"name": f"T{p}_charpoly_divides_route_B_with_eisenstein_cofactor",
                      "status": "PASS" if (r.is_zero and q.degree() == cx.inv.cusps - 1) else "FAIL"})
    al, trb = {}, {}
    for Q in als:
        al[Q] = _parse_int_matrix(grab(f"MAT_W{Q}"))
        trb[Q] = int(grab(f"TRB_W{Q}"))
        gates.append({"name": f"W{Q}_integral", "status": "PASS" if grab(f"INT_W{Q}") == "1" else "FAIL"})
        gates.append({"name": f"W{Q}_involution", "status": "PASS" if np.array_equal(al[Q] @ al[Q], np.eye(cb.rank, dtype=np.int64)) else "FAIL"})
        gates.append({"name": f"W{Q}_trace_equals_route_B", "status": "PASS" if int(np.trace(al[Q])) == trb[Q] else "FAIL"})
        gates.append({"name": f"W{Q}_preserves_intersection", "status": "PASS" if np.array_equal(al[Q].T @ Jint @ al[Q], Jint) else "FAIL"})
    for p in primes:
        gates.append({"name": f"T{p}_selfadjoint_for_intersection", "status": "PASS" if np.array_equal(hecke[p].T @ Jint, Jint @ hecke[p]) else "FAIL"})
        for Q in als:
            gates.append({"name": f"T{p}_commutes_W{Q}", "status": "PASS" if np.array_equal(hecke[p] @ al[Q], al[Q] @ hecke[p]) else "FAIL"})
    if len(als) > 1:
        prod = np.eye(cb.rank, dtype=np.int64)
        for Q in als[:-1]:
            prod = prod @ al[Q]
        gates.append({"name": "product_of_prime_power_AL_is_W_N", "status": "PASS" if np.array_equal(prod, al[N]) else "FAIL"})
    Qm = _parse_matrix(grab("MAT_Q"))
    return Transport(N, hecke, al, cpa, cpb, trb, Qm, gates, out)


def sector_census(al: Dict[int, np.ndarray], Q1: int, Q2: int) -> Dict[str, int]:
    """Joint (w_Q1, w_Q2) eigen-multiplicities on H_1; halve for S_2."""
    n = al[Q1].shape[0]
    I = np.eye(n, dtype=np.int64)
    out = {}
    for a in (1, -1):
        for b in (1, -1):
            P = (I + a * al[Q1]) @ (I + b * al[Q2])
            out[f"({'+' if a > 0 else '-'},{'+' if b > 0 else '-'})"] = int(np.trace(P)) // 4
    return out
