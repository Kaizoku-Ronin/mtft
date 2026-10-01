#!/usr/bin/env python3
"""Matter-03: field-coordinate kinetic normalization, not mass prediction.

All Gram matrices, the two-mode tensor extension, and the physical Higgs
scale in this audit are explicitly supplied examples. The rank-three slice
comes from Matter-02 but is not an identified light-family block.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp


def canonical_frame(gram):
    gram = np.asarray(gram, dtype=complex)
    if gram.ndim != 2 or gram.shape[0] != gram.shape[1] or not gram.shape[0]:
        raise ValueError("Gram must be a nonempty square matrix")
    if not np.isfinite(gram).all():
        raise ValueError("Gram must be finite")
    if not np.allclose(gram, gram.conj().T, atol=1e-12, rtol=1e-12):
        raise ValueError("Gram must be Hermitian")
    try:
        chol = np.linalg.cholesky((gram + gram.conj().T) / 2)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Gram must be positive definite") from exc
    return np.linalg.solve(chol.conj().T, np.eye(len(gram), dtype=complex))


def normalize_tensor(tensor, grams):
    """Lower-index W=L^T Y_k E H_k, raw fields = C canonical fields."""
    frames = tuple(canonical_frame(gram) for gram in grams)
    normalized = np.einsum("ai,bj,ck,abc->ijk", *frames, tensor)
    return normalized, frames


def _array(matrix):
    return np.array(matrix.evalf().tolist(), dtype=complex)


def run_audit():
    checks = []

    def check(name, passed, value=None):
        checks.append({"name": name, "passed": bool(passed), "value": value})
        if not passed:
            raise AssertionError(name)

    I = sp.I
    Q = sp.Rational
    c_l = sp.Matrix([[1, (1 + I)/3, 0], [0, Q(1, 2), I/4], [0, 0, Q(2, 3)]])
    c_e = sp.Matrix([[Q(1, 2), -I/3, Q(1, 5)], [0, Q(2, 3), (1-I)/4], [0, 0, Q(3, 2)]])
    c_h = sp.Matrix([[Q(2, 3), I/5], [0, Q(3, 4)]])
    exact_frames = (c_l, c_e, c_h)
    exact_grams = tuple(sp.simplify(c.inv().conjugate().T * c.inv()) for c in exact_frames)
    for label, c, g in zip(("L", "E", "H"), exact_frames, exact_grams):
        check(f"exact {label} canonical kinetic identity", sp.simplify(c.conjugate().T*g*c)==sp.eye(c.rows))
        check(f"exact {label} positive Gram by invertible congruence", c.det()!=0)

    witness = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    # H canonical direction (1,0): raw Higgs direction is (2/3,0).
    exact_y = sp.simplify(c_l.T * witness * c_e * Q(2, 3))
    expected_det = c_l.det()*c_e.det()*Q(2, 3)**3*witness.det()
    check("exact ordinary-transpose holomorphic determinant", sp.simplify(exact_y.det()-expected_det)==0)
    check("exact rank-three slice survives positive normalization", exact_y.rank()==3)

    x1, y1, x2, y2 = sp.symbols("x1 y1 x2 y2", real=True)
    h = sp.Matrix([x1+I*y1, x2+I*y2])
    pauli = (sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-I],[I,0]]), sp.diag(1,-1))
    rho = sp.expand((h.conjugate().T*h)[0])
    pauli_sum = sum(sp.expand((h.conjugate().T*s*h)[0]**2) for s in pauli)
    check("exact Pauli D-term identity", sp.simplify(pauli_sum-rho**2)==0)
    g2, gy, u, d = sp.symbols("g2 gy u d", positive=True, real=True)
    neutral_d = (g2*g2+gy*gy)/8*(u*u-d*d)**2
    check("neutral paired-Higgs D-flat ray", sp.simplify(neutral_d.subs(d,u))==0)
    check("unbroken zero-FI D potential has no negative quadratic mass", sp.hessian(neutral_d,(u,d)).subs({u:0,d:0})==sp.zeros(2))

    grams = tuple(_array(g) for g in exact_grams)
    tensor = np.empty((3,3,2), dtype=complex)
    tensor[:,:,0] = _array(witness)
    tensor[:,:,1] = np.array([[0.2+0.1j,0,0.3j],[0.4,0.5j,0],[0,-0.2j,0.1]])
    yn, frames = normalize_tensor(tensor, grams)
    for label, c, g in zip(("L", "E", "H"), frames, grams):
        residual = float(np.linalg.norm(c.conj().T@g@c-np.eye(len(g)),2))
        check(f"numeric {label} canonical kinetic residual", residual<2e-13, residual)
    check("numeric/exact frame agreement", all(np.allclose(c,_array(exact),atol=2e-14) for c,exact in zip(frames,exact_frames)))
    check("numeric/exact rank witness normalization", np.allclose(yn[:,:,0],_array(exact_y),atol=2e-14))

    rng = np.random.default_rng(14303)
    n = np.array([1+0.2j, -0.3+0.5j]); n /= np.linalg.norm(n)
    # Scale one, in declared arbitrary physical units: no electroweak prediction.
    physical_v = 1.0
    matrix = physical_v/np.sqrt(2)*np.einsum("ijk,k->ij",yn,n)
    singular = np.linalg.svd(matrix, compute_uv=False)
    worst_mass = 0.0; worst_unitary = 0.0; worst_covariance = 0.0
    for trial in range(24):
        changes = tuple(2*np.eye(size)+0.4*(rng.normal(size=(size,size))+1j*rng.normal(size=(size,size))) for size in (3,3,2))
        changed_grams = tuple(b.conj().T@g@b for b,g in zip(changes,grams))
        changed_tensor = np.einsum("ai,bj,ck,abc->ijk",*changes,tensor)
        yn2, frames2 = normalize_tensor(changed_tensor,changed_grams)
        rotations = tuple(np.linalg.solve(c,b@cp) for c,b,cp in zip(frames,changes,frames2))
        unitary_res = max(float(np.linalg.norm(rot.conj().T@rot-np.eye(len(rot)),2)) for rot in rotations)
        n2 = rotations[2].conj().T@n
        matrix2 = physical_v/np.sqrt(2)*np.einsum("ijk,k->ij",yn2,n2)
        singular2 = np.linalg.svd(matrix2,compute_uv=False)
        mass_res = float(np.max(np.abs(singular2-singular))/singular[0])
        covariance_res = float(np.linalg.norm(matrix2-rotations[0].T@matrix@rotations[1],2)/np.linalg.norm(matrix,2))
        check(f"basis {trial}: canonical frames differ unitarily",unitary_res<4e-12,unitary_res)
        check(f"basis {trial}: same physical Higgs direction",abs(float(np.vdot(n2,n2).real)-1)<4e-12)
        check(f"basis {trial}: holomorphic mass covariance",covariance_res<4e-12,covariance_res)
        check(f"basis {trial}: physical singular values invariant",mass_res<4e-12,mass_res)
        worst_mass=max(worst_mass,mass_res); worst_unitary=max(worst_unitary,unitary_res); worst_covariance=max(worst_covariance,covariance_res)

    # A conjugate transpose on a holomorphic slot is not this convention.
    raw_h = frames[2]@n
    yh = np.einsum("abc,c->ab",tensor,raw_h)
    wrong = physical_v/np.sqrt(2)*frames[0].conj().T@yh@frames[1]
    wrong_s = np.linalg.svd(wrong,compute_uv=False)
    conjugation_error = float(np.max(np.abs(wrong_s-singular))/singular[0])
    check("incorrect conjugation detected on complex nonorthogonal example", conjugation_error>1e-3,conjugation_error)

    # Supplied metric data changes masses without changing holomorphic rank.
    y0 = _array(witness)
    diag_grams = (np.diag([1.,4.,25.]),np.diag([1.,9.,16.]),np.array([[1.]]))
    normalized_diag,_ = normalize_tensor(y0[:,:,None],diag_grams)
    diagonal_s = np.linalg.svd(normalized_diag[:,:,0],compute_uv=False)
    check("positive metric can change mass ratios while preserving rank", np.allclose(diagonal_s,[1.,1/8,1/15],atol=1e-14))
    for label,bad in (("indefinite",np.diag([1.,-1.])),("singular",np.diag([1.,0.])),("non-Hermitian",[[1,1j],[0,1]])):
        try:
            canonical_frame(bad)
        except ValueError:
            check(f"invalid {label} kinetic matrix rejected",True)
        else:
            check(f"invalid {label} kinetic matrix rejected",False)
    return {
        "status":"DIAGNOSTIC examples and EXACT conditional normalization identities; no mass prediction",
        "convention":"W=L^T Y_k E H_k; raw fields=C canonical fields; C†GC=I",
        "inputs":{"Grams":"supplied complex positive-definite examples, not B3 metrics", "tensor":"Matter-02 rank-three slice plus supplied second-Higgs example", "physical_v":physical_v,"scale_units":"arbitrary; INPUT", "light_family_selection":"OPEN"},
        "exact_witness_normalized_determinant":str(sp.factor(exact_y.det())),
        "diagnostic_singular_values":singular.tolist(),
        "diagonal_metric_singular_values":diagonal_s.tolist(),
        "basis_change_count":24,
        "max_mass_covariance_relative_residual":worst_covariance,
        "max_singular_value_relative_residual":worst_mass,
        "max_unitarity_residual":worst_unitary,
        "wrong_conjugation_relative_mass_error":conjugation_error,
        "check_count":len(checks),"passed":all(item["passed"] for item in checks),"checks":checks,
    }


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("kinetic_checks.json"))
    args=parser.parse_args()
    result=run_audit()
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(f"Matter-03 kinetics: {result['check_count']} checks passed; mass prediction OPEN")
