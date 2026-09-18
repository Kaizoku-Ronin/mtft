#!/usr/bin/env python3
"""Reproduce the v0.31.4 spin-lift and anomaly audit without editing MTFT.

Usage: python audit_v0314.py --source /path/to/mtft-0.31.4 --output results.json
Exact computations use SymPy; frozen Yukawa reconstruction is diagnostic.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import sympy as sp


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def simplify_matrix(m):
    return m.applyfunc(sp.simplify)


def encode_matrix(m):
    return [[str(x) for x in row] for row in m.tolist()]


def symmetric_square(m):
    a,b,c,d = list(m)
    return simplify_matrix(sp.Matrix([[a*a,a*b,b*b],
                                     [2*a*c,a*d+b*c,2*b*d],
                                     [c*c,c*d,d*d]]))


def group_stats(generators):
    identity = sp.eye(2)
    seen = {tuple(identity): identity}
    pending = [identity]
    while pending:
        a = pending.pop()
        for b in generators:
            c = simplify_matrix(a * b)
            if tuple(c) not in seen:
                seen[tuple(c)] = c
                pending.append(c)
                assert len(seen) <= 8
    orders = []
    for a in seen.values():
        orders.append(next(n for n in range(1, 9)
                           if simplify_matrix(a ** n) == identity))
    return {"order": len(seen), "element_order_counts": dict(Counter(orders))}


def dedekind_sum(d, c):
    # Equivalent to the sawtooth definition for coprime d,c, c>0.
    return sum((sp.Rational(r, c) * (sp.Rational(d * r, c)
                - sp.floor(sp.Rational(d * r, c)) - sp.Rational(1, 2))
                for r in range(1, c)), sp.Integer(0))


def eta_transform(exponents, Q, matrix):
    """Exact eta multiplier with sqrt(c*z+d) convention, c>0.

    eta(gamma*z) = exp(pi*i*((a+d)/(12*c)-s(d,c)-1/4))
                   * sqrt(c*z+d) * eta(z).
    Decompose delta*W_Q(tau) = gamma_delta(delta_prime*tau).
    Return C in F(W_Q tau)=C*(c*tau+d)^weight*F_permuted(tau).
    """
    a, b, c, d = matrix
    divisors = [1, 11, 13, 143]
    phase, scale = sp.Integer(0), sp.Integer(1)
    transformed = dict.fromkeys(divisors, 0)
    factors = []
    for delta, r in zip(divisors, exponents):
        delta_prime = Q * delta // int(sp.gcd(delta, Q)) ** 2
        lam = Q if delta % Q == 0 else 1
        ga = sp.Rational(delta * a, delta_prime * lam)
        gb = sp.Rational(delta * b, lam)
        gc = sp.Rational(c, delta_prime * lam)
        gd = sp.Rational(d, lam)
        assert all(x.is_integer for x in (ga, gb, gc, gd))
        assert ga * gd - gb * gc == 1 and gc > 0
        theta = (ga + gd) / (12 * gc) - dedekind_sum(int(gd), int(gc)) - sp.Rational(1, 4)
        phase += r * theta
        scale *= sp.Integer(lam) ** (-sp.Rational(r, 2))
        transformed[delta_prime] += r
        factors.append({"delta": delta, "delta_prime": delta_prime,
                        "gamma": [int(ga), int(gb), int(gc), int(gd)],
                        "phase_over_pi": str(theta)})
    constant = sp.simplify(scale * sp.exp(sp.I * sp.pi * phase))
    return constant, [transformed[x] for x in divisors], factors


def eta_series(exponents, precision):
    ds = [1, 11, 13, 143]
    valuation = sum(d*r for d, r in zip(ds, exponents)) // 24
    assert sum(d*r for d, r in zip(ds, exponents)) % 24 == 0
    coeff = [1] + [0] * precision
    for d, r in zip(ds, exponents):
        for _ in range(abs(r)):
            for n in range(1, precision // d + 1):
                step = d*n
                if r > 0:
                    for j in range(precision, step-1, -1):
                        coeff[j] -= coeff[j-step]
                else:
                    for j in range(step, precision+1):
                        coeff[j] += coeff[j-step]
    return [0] * valuation + coeff[:precision+1-valuation]


def convolve(a, b):
    return [sum(a[j]*b[k-j] for j in range(k+1))
            for k in range(min(len(a), len(b)))]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("results.json"))
    args = parser.parse_args()
    source = args.source.resolve()
    sys.path.insert(0, str(source / "src"))
    import mtft
    from mtft.surface import arithspin as AS, smflux as SF, hym as HY
    assert Path(mtft.__file__).resolve().is_relative_to(source)
    result = {"source": str(source), "python": platform.python_version(),
              "sympy": sp.__version__, "numpy": np.__version__,
              "scope": "targeted audit, not the full release suite"}
    tracked = ["src/mtft/surface/arithspin.py", "src/mtft/surface/smflux.py",
               "src/mtft/surface/hym.py", "src/mtft/surface/_data/x0143_weight2_basis.json",
               "src/mtft/surface/_data/x0143_m2_tensors_h02.npz"]
    result["input_sha256"] = {p: digest(source/p) for p in tracked}

    # Exact eta multipliers: independent of frozen modular-form matrices.
    unit = [-1, -1, 1, 1]
    form = [2, 2, 0, 0]
    eta_checks = {}
    for Q in (11, 13):
        cu, ru, detail = eta_transform(unit, Q, AS.AL_MATRICES[f"W{Q}"])
        cf, rf, _ = eta_transform(form, Q, AS.AL_MATRICES[f"W{Q}"])
        assert cu == (-1 if Q == 11 else -sp.Rational(1, 13))
        assert ru == (unit if Q == 11 else [-x for x in unit])
        assert Q * cf == (-1 if Q == 11 else 13)
        assert rf == (form if Q == 11 else [0, 0, 2, 2])
        eta_checks[f"W{Q}"] = {"unit_multiplier": str(cu),
                "unit_transformed_exponents": ru,
                "differential_multiplier": str(Q*cf),
                "form_transformed_exponents": rf,
                "SL2Z_decompositions": detail}
    ds = [1, 11, 13, 143]
    order_matrix = sp.Matrix([[sp.Rational(143*int(sp.gcd(d,c))**2,24*c*d)
                               for d in ds] for c in ds])
    omega_divisor = order_matrix * sp.Matrix(form) - sp.ones(4,1)
    assert list(omega_divisor) == [12,12,0,0]
    assert list(order_matrix*sp.Matrix(unit)) == [-6,-6,6,6]
    result["exact_eta_checks"] = eta_checks
    result["omega_divisor"] = list(map(int,omega_divisor))

    # Exact action on H^0(S0), and on its symmetric square mapped into H^0(K).
    A = sp.diag(1,-1)
    B = sp.Matrix([[0,-1/sp.sqrt(13)],[sp.sqrt(13),0]])
    A_spin = sp.I*A
    identity = sp.eye(2)
    assert A*A == identity and B*B == -identity
    assert simplify_matrix(A*B*A.inv()*B.inv()) == -identity
    assert (A*B)**2 == identity
    assert A_spin*A_spin == B*B == (A_spin*B)**2 == -identity
    assert simplify_matrix(A_spin*B*A_spin.inv()*B.inv()) == -identity
    D8, Q8 = group_stats([A,B]), group_stats([A_spin,B])
    assert D8["element_order_counts"] == {1:1,2:5,4:2}
    assert Q8["element_order_counts"] == {1:1,2:1,4:6}
    result["section_lift"] = {"A":encode_matrix(A),"B":encode_matrix(B),
                               "signs":[1,-1,-1,1],"group":D8}
    result["theta_compatible_lift"] = {"A":encode_matrix(A_spin),"B":encode_matrix(B),
                               "signs":[-1,-1,-1,-1],"group":Q8}

    # Separate check against the shipped exact q-expansion and AL dataset.
    basis_path = source / tracked[3]
    data = json.loads(basis_path.read_text())
    precision = 130
    F = sp.Matrix([x[:precision+1] for x in data["coefficients"]]).T
    WA = sp.Matrix(data["W11"]).applyfunc(sp.Rational)
    WB = sp.Matrix(data["W13"]).applyfunc(sp.Rational)
    assert WA*WA == WB*WB == sp.eye(13) and WA*WB == WB*WA
    conditions = sp.Matrix.vstack((F*WA*WB)[1:13,:],(F*WB)[1:13,:])
    ns = conditions.nullspace()
    assert len(ns) == 1
    v = ns[0] / (F*ns[0])[1]
    f = sp.Matrix(eta_series(form,precision))
    u = eta_series(unit,precision)
    assert F*v == f and WA*v == -v
    u2 = convolve(u,u)
    assert F*WB*v == 13*sp.Matrix(convolve(list(f),u2))
    product_forms = sp.Matrix.hstack(f,sp.Matrix(convolve(list(f),u)),
                                    sp.Matrix(convolve(list(f),u2)))
    product_coords = F.gauss_jordan_solve(product_forms)[0]
    natural_A = sp.diag(-1,1,-1)
    natural_B = sp.Matrix([[0,0,sp.Rational(1,13)],[0,-1,0],[13,0,0]])
    assert WA*product_coords == product_coords*natural_A
    assert WB*product_coords == product_coords*natural_B
    assert symmetric_square(A) != natural_A
    assert symmetric_square(A_spin) == natural_A
    assert symmetric_square(B) == natural_B
    result["exact_frozen_basis_checks"] = {
        "q_precision":precision,"weight2_Sturm_bound":28,
        "omega_basis_vector":list(map(str,v)),"square_map_rank":product_forms.rank(),
        "symmetric_square_W11":encode_matrix(natural_A),
        "symmetric_square_W13":encode_matrix(natural_B),
        "section_A_square_map_compatible":False,
        "spin_A_square_map_compatible":True,"B_square_map_compatible":True}
    fixed = {}
    for name,W in [("W11",WA),("W13",WB),("W143",WA*WB)]:
        trace = int(W.trace())
        fixed[name] = {"trace_on_holomorphic_differentials":trace,
                       "quotient_genus":(13+trace)//2,"fixed_points":2-2*trace}
    assert [x["fixed_points"] for x in fixed.values()] == [0,4,20]
    result["AL_geometry"] = fixed

    # Symbolic identity, stronger than evaluations at finitely many q vectors.
    S = SF.M1_STACKS
    names,N,m = list(S["N"]),S["N"],S["m"]
    alpha,beta,gamma = sp.symbols("alpha beta gamma")
    q = {x:alpha+beta*sp.Rational(S["Y"][x])+gamma*sp.Rational(S["B-L"][x]) for x in names}
    cubic = sp.expand(sum((m[a]-m[b])*N[a]*N[b]*(q[a]-q[b])**3
                         for i,a in enumerate(names) for b in names[i+1:]))
    grav = sp.expand(sum((m[a]-m[b])*N[a]*N[b]*(q[a]-q[b])
                         for i,a in enumerate(names) for b in names[i+1:]))
    mixed = {a:sp.expand(sum(N[b]*(m[a]-m[b])*(q[a]-q[b])
                         for b in names if b != a)) for a in names if N[a]>=2}
    assert cubic == grav == 0 and all(v == 0 for v in mixed.values())
    maj = SF.majorana_obstruction()
    assert maj["majorana_bilinear_charge"] == 2 and not maj["bare_majorana_allowed"]
    result["exact_anomalies"] = {"cubic_restricted":str(cubic),"gravitational_restricted":str(grav),
                                "mixed_nonabelian_restricted":{k:str(v) for k,v in mixed.items()},
                                "nu_c_BL":1,"majorana_bilinear_BL":2,
                                "phase_acts_trivially_on_bifundamentals":True}

    # Frozen numerical data: provenance and reproducibility, not an error bound.
    d = HY.load_m2_tensors()
    rng = np.random.default_rng(20260916)
    diagnostics = {}
    for label,key,raw,na,nb,nh in [
        ("M2_up","Yu","Y_up_raw","N_Q","N_u","N_Hu"),
        ("M2_down","Yd","Y_down_raw","N_L","N_L","N_Hd")]:
        yn = HY.normalise_yukawa(d[raw],d[na],d[nb],d[nh])
        err = float(np.linalg.norm(yn-d[key])/np.linalg.norm(d[key]))
        assert err < 1e-10
        basis_ok = HY.normalisation_is_basis_invariant(d[raw],d[na],d[nb],d[nh])
        assert basis_ok
        directions = rng.standard_normal((10000,yn.shape[2]))+1j*rng.standard_normal((10000,yn.shape[2]))
        matrices = np.einsum("ijk,nk->nij",yn,directions)
        singular = np.linalg.svd(matrices,compute_uv=False)
        ratios = singular[:,::-1]/singular[:,0,None]
        diagnostics[label] = {"relative_reconstruction_error":err,"basis_gate_passed":basis_ok,
                             "median_m1_m3":float(np.median(ratios[:,0])),
                             "median_m2_m3":float(np.median(ratios[:,1]))}
    result["frozen_data_diagnostic"] = {"sample_count_per_sector":10000,"seed":20260916,
                  "method":"isotropic complex Gaussian directions; ratios independent of vector norm",
                  "results":diagnostics,"keys":list(d),
                  "fit_vectors_present":any(k.startswith("v_") for k in d)}
    # Conditional minimal 6D gravity model: Einstein frame, compact negative
    # curvature, ordinary positive flux energy, nonnegative bulk vacuum energy.
    rho,a,b,c = sp.symbols("rho a b c", positive=True)
    V = a/rho**4+b/rho**6+c/rho**2
    dV = sp.diff(V,rho)
    assert dV == -4*a/rho**5-6*b/rho**7-2*c/rho**3
    result["conditional_radius_no_go"] = {
        "Einstein_frame_potential":str(V),"radius_derivative":str(dV),
        "conditions":"a>0, b>=0, c>=0, rho>0; unwarped product; two-derivative Einstein gravity; positive gauge kinetic terms; fixed quantized flux; no additional sources or corrections",
        "conclusion":"no stationary radius in this restricted added parent"}
    result["status"] = "PASS: exact identities and declared numerical gates verified"
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"groups":{"section":D8,"theta":Q8},
                      "fixed_points":fixed,"diagnostic":diagnostics,
                      "output":str(args.output)},indent=2))


if __name__ == "__main__":
    main()
