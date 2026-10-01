"""Matter-02 exact cohomology and a holomorphic Yukawa rank witness.

Fixed representative: S=Bl_(infinity,0)(X0(143) x 143a1),
L=pi*pr_X*O_X(infinity) tensor O(-3E_exc), with no flat twists.
These are algebraic vector spaces and cubic couplings, not particle masses.
The split background has Higgs expectation value zero. Input provenance is
checked against the archived full weight-two cusp-form basis.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import argparse
import csv
import hashlib
import json

BACKGROUND = "Bl_(infinity,0)(X0(143) x 143a1); L=pi*pr_X*O_X(infinity) tensor O(-3E_exc); no flat twists"
STACK_POWERS = {"c":0, "L":1, "A1":2, "A2":-3, "B":-1}
_INPUTS = {"X0_143_S2_qexpansions.txt":"479316a0952d72dbd4731f50389003263cbfc4e1",
           "X0_143_AL_adapted_qexpansions.txt":"e71eaabe5cd23f3a3a35c204373ad6dfd5115358",
           "X0_143_AL_adapted_basis.txt":"5e4537c805a48e76027d4d12bdf0cea300cd1472"}


def _integer(q):
    if not isinstance(q,int) or isinstance(q,bool):
        raise TypeError("q must be an integer; booleans are not exponents")
    return q


def _h2(q):
    return 13 if q == 0 else max(13-q,0) if q > 0 else max(14+2*q,0)


@dataclass(frozen=True)
class Cohomology:
    q: int
    h0: int
    h1: int
    h2: int
    chi: int
    phi_dim: int
    holo_chiral_dim: int | None
    archived_net: int
    bhv_net: int


def cohomology(q: int) -> Cohomology:
    """Exact dimensions at the certified ordinary cusp; q0 gauge modes separate.

    phi_dim=h0(K_S tensor L^q)=h2(L^-q). For nonzero q, BHV bulk
    chiral dimension is h1(L^q)+phi_dim. Its difference from the conjugate
    is +3q, opposite the archived Euler-index orientation -3q.
    """
    q=_integer(q)
    chi=(-9*q*q-3*q)//2
    h0=1 if q == 0 else 0
    h2=_h2(q)
    h1=h0+h2-chi
    phi_dim=_h2(-q)
    return Cohomology(q,h0,h1,h2,chi,phi_dim,None if q == 0 else h1+phi_dim,-3*q,3*q)


def _rref(matrix):
    a=[[Fraction(x) for x in row] for row in matrix]
    pivots=[]
    if not a:
        return a,pivots
    rr=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(rr,len(a)) if a[i][col]),None)
        if pivot is None:
            continue
        a[rr],a[pivot]=a[pivot],a[rr]
        divisor=a[rr][col]
        a[rr]=[value/divisor for value in a[rr]]
        for i in range(len(a)):
            if i != rr and a[i][col]:
                scalar=a[i][col]
                a[i]=[u-scalar*v for u,v in zip(a[i],a[rr])]
        pivots.append(col)
        rr+=1
        if rr == len(a):
            break
    return a,pivots


def _det(matrix):
    a=[[Fraction(x) for x in row] for row in matrix]
    result=Fraction(1)
    for col in range(len(a)):
        pivot=next((i for i in range(col,len(a)) if a[i][col]),None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col],a[pivot]=a[pivot],a[col]
            result=-result
        result*=a[col][col]
        for i in range(col+1,len(a)):
            scale=a[i][col]/a[col][col]
            for j in range(col,len(a)):
                a[i][j]-=scale*a[col][j]
    return result


def _load_inputs():
    directory=Path(__file__).resolve().parents[1]/"canonical"/"_data"
    matrices={}; hashes={}
    for name,expected in _INPUTS.items():
        raw=(directory/name).read_bytes()
        sha1=hashlib.sha1(raw).hexdigest()
        if sha1 != expected:
            raise ValueError(f"Archived input provenance mismatch: {name}")
        rows=list(csv.reader(line for line in raw.decode().splitlines()
                             if line.strip() and not line.lstrip().startswith("#")))
        matrix=[[int(value) for value in row[1:]] for row in rows[1:]]
        expected_rows=13 if "_basis." in name else 141
        if len(matrix)!=expected_rows or any(len(row)!=13 for row in matrix):
            raise ValueError(f"Archived matrix shape mismatch: {name}")
        matrices[name]=matrix
        hashes[name]={"sha1":sha1,"sha256":hashlib.sha256(raw).hexdigest()}
    return matrices,hashes


def cusp_jet_certificate() -> dict:
    """Exact original/adapted jet ranks and basis-change certificate.

    f(z) dz/z has order >=m iff coefficients a1,...,am vanish. Rank13
    proves canonical orders0..12 and ordinary-cusp gaps1..13. Correctness
    of the archived modular-form input remains an explicit assumption.
    """
    inputs,hashes=_load_inputs()
    f=inputs["X0_143_S2_qexpansions.txt"]
    g=inputs["X0_143_AL_adapted_qexpansions.txt"]
    b=inputs["X0_143_AL_adapted_basis.txt"]
    transformed=[[sum(row[k]*b[k][j] for k in range(13)) for j in range(13)] for row in f]
    rows=[]
    for m in range(1,14):
        rank=len(_rref(f[1:m+1])[1]); rank_adapted=len(_rref(g[1:m+1])[1])
        rows.append({"vanishing_order":m,"rank_original":rank,"rank_adapted":rank_adapted,"h0_K_minus_mP":13-rank})
    passed=(f[0]==g[0]==[0]*13 and transformed==g
            and all(row["rank_original"]==row["rank_adapted"]==row["vanishing_order"] for row in rows))
    if not passed:
        raise ValueError("Archived basis failed the exact cusp-jet certificate")
    return {"status":"EXACT CONDITIONAL GEOMETRY","passed":passed,"P":"rational width1 cusp infinity",
            "jet_rows":rows,"full13_jet_determinant":str(_det(f[1:14])),
            "basis_change_determinant":str(_det(b)),"gap_sequence":list(range(1,14)),
            "input_hashes":hashes,"assumption":"Archived columns form the full S2(Gamma0(143)) basis."}


@dataclass(frozen=True)
class YukawaWitness:
    canonical_basis_coefficients: tuple[Fraction,...]
    canonical_q_coefficients: tuple[Fraction,...]
    phi_coefficients: tuple[tuple[int,Fraction],...]
    positive_jets: tuple[tuple[int,int],...]
    negative_principal_parts: tuple[tuple[int,int],...]
    residue_matrix: tuple[tuple[Fraction,...],...]
    determinant: Fraction
    root_bracket_trace: int
    full_matrix_shape: tuple[int,int] = (18,44)
    higgs_space_dimension: int = 10
    rank_lower_bound: int = 3


def _root_bracket_trace(a,b,c):
    """Tr([E_a,E_b] E_c) for matrix units, without physical charge assignment."""
    i,j=a; k,l=b; m,n=c
    return int(j==k and l==m and i==n)-int(l==i and j==m and k==n)


@lru_cache(maxsize=1)
def yukawa_witness() -> YukawaWitness:
    """Exact cubic rank witness; algebraic coordinates, no Higgs vev or metric.

    Solve a1=a2=a3=a5=0, a4=1 in the original cusp basis. The selected
    phi2 coefficient relative to K_X tensor O(2P) begins z^5. Residues
    pair three fat jets with three globally unobstructed exceptional
    principal parts. This is a lower-rank slice of an18x44 tensor with
    ten Higgs directions, not a normalized physical electron matrix.
    """
    cusp_jet_certificate()
    inputs,_=_load_inputs()
    f=inputs["X0_143_S2_qexpansions.txt"]
    equations=[row+[int(n==4)] for n,row in enumerate(f[1:6],1)]
    reduced,pivots=_rref(equations)
    if 13 in pivots or len(pivots)!=5:
        raise ValueError("Canonical order3 coefficient system is inconsistent")
    vector=[Fraction(0)]*13
    for row,pivot in enumerate(pivots):
        vector[pivot]=reduced[row][-1]
    coefficients=tuple(sum(Fraction(row[i])*vector[i] for i in range(13)) for row in f)
    if coefficients[1:6]!=(0,0,0,1,0):
        raise ValueError("Canonical witness has the wrong leading coefficients")
    phi={n+1:value for n,value in enumerate(coefficients) if value}
    positive=((0,1),(0,2),(1,1))
    negative=((6,2),(7,2),(6,3))
    residue=tuple(tuple(phi.get(ii-i-1,Fraction(0)) if jj-j==1 else Fraction(0)
                        for ii,jj in negative) for i,j in positive)
    det=_det(residue)
    # Weak index3, A1 index5, B index7 in the U8 parent. Their exponents
    # are1,-3,2, respectively, and their root bracket is nonzero.
    trace=_root_bracket_trace((5,3),(7,5),(3,7))
    if det != -1 or trace != -1:
        raise ValueError("Yukawa residue or root-bracket certificate failed")
    return YukawaWitness(tuple(vector),coefficients,tuple(phi.items()),positive,negative,residue,det,trace)


def _json_safe(value):
    if isinstance(value,Fraction):
        return str(value)
    if isinstance(value,dict):
        return {str(k):_json_safe(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [_json_safe(v) for v in value]
    return value


def report() -> dict:
    return _json_safe({"milestone":"Matter-02","status":"EXACT CONDITIONAL GEOMETRY",
        "background":BACKGROUND,"stack_powers":STACK_POWERS,
        "cohomology":[asdict(cohomology(q)) for q in range(-5,6)],
        "cusp_certificate":cusp_jet_certificate(),"yukawa_witness":asdict(yukawa_witness()),
        "chirality_note":"BHV N(q)-N(-q)=+3q; archived chi(q)-chi(-q)=-3q. Physical representation orientation OPEN.",
        "neutral_note":"q0 H0 gauge/vector modes are separate from H1 and H0(K) adjoint chiral modes.",
        "higgs_vev_selected":False,"physical_mass_matrix":False,
        "open":["physical chirality assignment","Higgs vacuum and paired-mode lifting",
                "L2 normalization and physical Yukawas","gauge kinetic and photon kernel matching"]})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=report()
    text=json.dumps(result,indent=2)+"\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text,end="")


if __name__=="__main__":
    main()
