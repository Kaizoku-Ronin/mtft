"""Physics-facing, but epistemically separated, tools built from X_0(143).

Everything through the Hodge metric, complex structure, and graph-channel
operator is mathematics on the certified curve/model.  Calling the resulting
quadratic form a *physical* mass or energy requires an external field theory
and scale; this module deliberately does not supply one.

The quantitative M8 tool is new: it transports the true period-derived Hodge
complex structure into the promoted ``mtft.hecke`` cuspidal frame, then splits
the canonical M7 graph coupling into complex-linear and J-antilinear pieces.
"""
from __future__ import annotations

from fractions import Fraction as Fr
import json, math

import mpmath as mp
import numpy as np

from mtft import hecke as H
from .core import hodge_complex_structure, hodge_metric, symplectic_form, data_path
from .bridge import hecke_to_symplectic_change


def _np(A):
    return np.array([[float(A[i,j]) for j in range(A.cols)] for i in range(A.rows)],float)


def hodge_structure_hecke(dps: int = 60) -> np.ndarray:
    """True period-derived J on the promoted 26-dim Hecke cuspidal basis."""
    P=np.array(hecke_to_symplectic_change(),float)  # s = P h
    J=_np(hodge_complex_structure(dps))
    return np.linalg.solve(P,J@P)


def hodge_metric_hecke(dps: int = 60) -> np.ndarray:
    """Positive Hodge metric in the promoted Hecke cuspidal basis."""
    P=np.array(hecke_to_symplectic_change(),float)
    G=_np(hodge_metric(dps))
    return P.T@G@P


def _harmonic_basis():
    rec=json.loads(data_path("X0_143_m7_harmonic_basis.json").read_text())
    return [[Fr(int(n),int(d)) for n,d in row] for row in rec["basis_26x84"]]


def _graph_potential(name: str):
    m=H.model(); E=m["E"]
    tris,tri_of,sS=m["tris"],m["tri_of"],m["sS"]
    erep,cusp_of,fans=m["erep"],m["cusp_of"],m["fans"]
    if name=="width":
        width={k:len(o) for k,o in enumerate(fans)}
        g=[sum(width[cusp_of[f]] for f in tris[t]) for t in range(56)]
    elif name=="degree":
        g=[0]*56
        for k in range(E):
            a,b=tri_of[erep[k]],tri_of[sS[erep[k]]];g[a]+=1;g[b]+=1
    elif name=="distance":
        adj=[[] for _ in range(56)]
        for k in range(E):
            a,b=tri_of[erep[k]],tri_of[sS[erep[k]]]
            if a!=b:adj[a].append(b);adj[b].append(a)
        g=[-1]*56;g[0]=0;front=[0]
        while front:
            nxt=[]
            for u in front:
                for v in adj[u]:
                    if g[v]<0:g[v]=g[u]+1;nxt.append(v)
            front=nxt
    else:
        raise ValueError("potential must be 'width', 'degree', or 'distance'")
    return g


def graph_coupling(potential: str = "width"):
    """M7 coupling V=G_graph^{-1}M in the promoted Hecke basis.

    The harmonic embedding is frozen as exact rational study data; G and M are
    recomputed from it at call time.  Returns ``(V, G_graph, M)`` as float64.
    The exact M7 study remains the certificate for the rational construction.
    """
    W=np.array([[float(x) for x in r] for r in _harmonic_basis()],float)
    m=H.model();g=_graph_potential(potential)
    gav=np.array([(g[m["tri_of"][m["erep"][k]]]
                  +g[m["tri_of"][m["sS"][m["erep"][k]]]])/2
                  for k in range(m["E"])],float)
    G=W@W.T;M=(W*gav[None,:])@W.T
    return np.linalg.solve(G,M),G,M


def complex_linear_decomposition(V, J):
    """Return (V_linear,V_antilinear) relative to J.

    V+ commutes with J and V- anticommutes with J:
        V+ = (V - J V J)/2,  V- = (V + J V J)/2.
    """
    V=np.asarray(V,float);J=np.asarray(J,float)
    return .5*(V-J@V@J), .5*(V+J@V@J)


def metric_hs_norm(A,G):
    r"""Basis-invariant Hilbert-Schmidt norm from a positive metric G.

      ||A||_G^2 = tr(G^{-1} A^T G A).
    """
    A=np.asarray(A,float);G=np.asarray(G,float)
    x=np.trace(np.linalg.solve(G,A.T@G@A))
    return math.sqrt(max(0.0,float(x)))


def cp_channel_report(potential: str = "width", dps: int = 60):
    """Quantitative M8 report using the *true period-derived* Hodge J.

    The reported ``antilinear_fraction`` is a new dimensionless diagnostic:
    ||V_-||_H/||V||_H under the Hodge metric.  It is NOT the historical M8b
    amplitude and is not a Standard-Model CP observable without an explicit
    physical coupling/identification.
    """
    J=hodge_structure_hecke(dps);G=hodge_metric_hecke(dps)
    V,Ggraph,M=graph_coupling(potential)
    Vp,Vm=complex_linear_decomposition(V,J)
    nV,nP,nM=(metric_hs_norm(A,G) for A in (V,Vp,Vm))
    I=np.array([[float(x) for x in r] for r in H.star_involution()])
    C=V@J-J@V
    return {
        "potential":potential,
        "norm_V":nV,
        "norm_complex_linear":nP,
        "norm_J_antilinear":nM,
        "antilinear_fraction":nM/nV,
        "antilinear_power_fraction":(nM/nV)**2,
        "commutator_norm_H":metric_hs_norm(C,G),
        "J_square_residual":float(np.linalg.norm(J@J+np.eye(26))),
        "V_star_commutator_rel":float(np.linalg.norm(V@I-I@V)/np.linalg.norm(V)),
        "commutator_star_odd_rel":float(np.linalg.norm(I@C+C@I)/np.linalg.norm(C)),
        "epistemic":"CERTIFIED(tol) mathematical diagnostic; physical CP interpretation is PHENO",
    }


def finite_charge_partition(charges, beta=1.0, dps: int = 60):
    """Finite diagnostic theta sum over explicitly supplied integer charges.

    Each charge is ``(n,m)`` with two length-13 sequences.  This intentionally
    refuses to pretend a truncated charge set is the full genus-13 theta
    function.
    """
    from .core import charge_energy
    with mp.workdps(dps):
        vals=[charge_energy(n,m,dps) for n,m in charges]
        return +sum(mp.e**(-mp.pi*mp.mpf(beta)*e) for e in vals)


__all__=["hodge_structure_hecke","hodge_metric_hecke","graph_coupling",
         "complex_linear_decomposition","metric_hs_norm","cp_channel_report",
         "finite_charge_partition"]
