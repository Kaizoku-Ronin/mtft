#!/usr/bin/env python3
"""Exact non-involutive Schur and amplitude-feedback controls.

This synthetic rational 3x3 example has two active amplitudes and one
complementary amplitude. No MTFT matrix or future primary output is read.
An exact rational orthogonal diagonalizer makes all unitary checks symbolic.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def run() -> dict:
    z=sp.symbols("z")
    t,s=sp.symbols("t s", real=True)
    I=sp.I
    eye3=sp.eye(3)
    v=sp.Matrix([1,2,3])
    O=eye3-v*v.T/sp.Integer(7)
    spectrum=[-2,1,3]
    H=O*sp.diag(*spectrum)*O.T
    A=H[:2,:2]
    B=H[:2,2:3]
    D=H[2:3,2:3]
    d=D[0,0]
    eye2=sp.eye(2)
    checks={}

    def simp(x):
        return sp.simplify(sp.expand(x))

    def zero(name,expr):
        vals=list(expr) if isinstance(expr,sp.MatrixBase) else [expr]
        checks[name]=all(simp(x)==0 for x in vals)
        if not checks[name]:
            raise AssertionError(name)

    def truth(name,condition):
        checks[name]=bool(condition)
        if not checks[name]:
            raise AssertionError(name)

    zero("rational_orthogonal_frame",O.T*O-eye3)
    zero("Hermitian_generator",H-H.conjugate().T)
    truth("generator_is_not_an_involution",H*H!=eye3)
    zero("specified_characteristic_polynomial",(z*eye3-H).det()-(z+2)*(z-1)*(z-3))
    truth("complement_coupling_full_column_rank",B.rank()==1)
    truth("compression_pole_outside_full_spectrum",(d*eye3-H).det()!=0)
    sigma=B*B.T/(z-d)
    S=z*eye2-A-sigma
    Gfull=O*sp.diag(*[1/(z-lam) for lam in spectrum])*O.T
    G=Gfull[:2,:2]
    zero("full_resolvent",(z*eye3-H)*Gfull-eye3)
    zero("finite_Schur_identity_left",S*G-eye2)
    zero("finite_Schur_identity_right",G*S-eye2)
    zero("self_energy_pole_residue",sigma.applyfunc(lambda x:sp.cancel((z-d)*x).subs(z,d))-B*B.T)
    zero("projected_resolvent_null_direction_at_compression_pole",G.subs(z,d)*B)
    truth("projected_resolvent_nullity_matches_complement_pole",2-G.subs(z,d).rank()==1)
    zero("high_frequency_self_energy_moment",sigma.applyfunc(lambda x:sp.limit(z*x,z,sp.oo))-B*B.T)

    # Eigenvalue phases have integer frequencies. These checks are exact.
    def unitary(time):
        return O*sp.diag(*[sp.exp(-I*lam*time) for lam in spectrum])*O.T
    U=unitary(t)
    Us=unitary(s)
    Uts=unitary(t+s)
    zero("unitary_evolution",U.conjugate().T*U-eye3)
    zero("unitary_evolution_equation",I*U.diff(t)-H*U)
    zero("full_composition_law",U*Us-Uts)
    Upp=U[:2,:2]
    Upq=U[:2,2:3]
    Uqp=U[2:3,:2]
    zero("projected_composition_defect",Upp*Us[:2,:2]-Uts[:2,:2]+Upq*Us[2:3,:2])
    defect=Upp*Us[:2,:2]-Uts[:2,:2]
    zero("composition_defect_leading_ts",defect.diff(t).diff(s).subs({t:0,s:0})-B*B.T)
    zero("general_projected_norm_balance",Upp.conjugate().T*Upp+Uqp.conjugate().T*Uqp-eye2)

    # This is a normalized, initially dark active vector: B^dagger*x=0.
    x0=sp.Matrix([1,-11])/sp.sqrt(122)
    zero("dark_input_normalized",(x0.T*x0)[0,0]-1)
    zero("dark_input_in_coupling_kernel",B.T*x0)
    v2=B.T*A*x0
    truth("dark_input_has_nonzero_second_order_transfer",v2!=sp.zeros(1,1))
    xt=Upp*x0
    yt=Uqp*x0
    zero("dark_complement_initial_value",yt.subs(t,0))
    zero("dark_complement_first_derivative",yt.diff(t).subs(t,0))
    zero("dark_complement_second_derivative",yt.diff(t,2).subs(t,0)+v2)
    norm_y=(yt.conjugate().T*yt)[0,0]
    for order in range(4):
        zero(f"dark_norm_derivative_order_{order}",sp.diff(norm_y,t,order).subs(t,0))
    quartic=(v2.T*v2)[0,0]/4
    zero("dark_norm_quartic_coefficient",sp.diff(norm_y,t,4).subs(t,0)/sp.factorial(4)-quartic)
    zero("real_control_norm_even",norm_y-norm_y.subs(t,-t))

    # Duhamel memory convolution is evaluated term by term in the full spectral
    # resolution. The primitive below is valid here since d != each lambda.
    projected_residues=[O[:2,j:j+1]*O[:2,j:j+1].T for j in range(3)]
    memory=sp.zeros(2,1)
    for lam,Pj in zip(spectrum,projected_residues):
        integral=(sp.exp(-I*lam*t)-sp.exp(-I*d*t))/(I*(d-lam))
        memory += B*B.T*Pj*x0*integral
    zero("memory_equation_for_initially_active_input",I*xt.diff(t)-A*xt+I*memory)
    # Repeat with nonzero initial complementary amplitude to retain its forcing.
    y0=sp.Rational(2,3)
    x_full=Upp*x0+Upq*sp.Matrix([y0])
    memory_full=sp.zeros(2,1)
    psi0=sp.Matrix([x0[0],x0[1],y0])
    for j,lam in enumerate(spectrum):
        active_component=O[:2,j:j+1]*(O[:,j:j+1].T*psi0)[0,0]
        integral=(sp.exp(-I*lam*t)-sp.exp(-I*d*t))/(I*(d-lam))
        memory_full += B*B.T*active_component*integral
    zero("memory_equation_with_initial_complement_forcing",I*x_full.diff(t)-A*x_full-B*sp.exp(-I*d*t)*y0+I*memory_full)
    truth("minimal_auxiliary_Krylov_dimension",sp.Matrix.hstack(B.T,D*B.T).rank()==1)
    truth("no_static_active_generator_composition_defect_nonzero",B*B.T!=sp.zeros(2,2))

    # A secondary Hermitian example demonstrates that O(t^5), rather than
    # O(t^6), is the general remainder in the dark-state squared norm.
    H_complex=sp.Matrix([[0,1,I,1],[1,0,1,0],[-I,1,0,0],[1,0,0,0]])
    dark_complex=sp.Matrix([0,1,0,0])
    zero("quintic_control_Hermitian",H_complex-H_complex.conjugate().T)
    zero("quintic_control_initial_coupling_zero",(H_complex*dark_complex)[3,0])
    y_short=sum(((-I*t)**j/sp.factorial(j)*(H_complex**j*dark_complex)[3,0] for j in range(1,4)),sp.Integer(0))
    norm_short=sp.expand(sp.conjugate(y_short)*y_short)
    zero("quintic_control_amplitude",y_short+t*t/2+t*t*t/6)
    zero("quintic_control_norm_quartic",norm_short.coeff(t,4)-sp.Rational(1,4))
    zero("quintic_control_nonzero_fifth_order_term",norm_short.coeff(t,5)-sp.Rational(1,6))

    return {
        "status":"EXACT_SYMBOLIC_SYNTHETIC_CONTROL",
        "scope":"Rational 3x3 Hermitian non-involution; two active plus one complement; no MTFT data",
        "H":[[str(x) for x in row] for row in H.tolist()],
        "full_eigenvalues":spectrum,
        "compressed_D_eigenvalue":str(d),
        "dark_input":[str(x) for x in x0],
        "Bstar_A_dark_input":[str(x) for x in v2],
        "dark_norm_quartic_coefficient":str(quartic),
        "dark_norm_general_remainder":"O(t^5); this real example has stronger O(t^6) by evenness",
        "secondary_complex_Hermitian_control":{
            "H":[[str(x) for x in row] for row in H_complex.tolist()],
            "initial_active_state":"(0,1,0), initial complement 0",
            "complement_amplitude":"-t^2/2-t^3/6+O(t^4)",
            "complement_norm_squared":"t^4/4+t^5/6+O(t^6)",
            "purpose":"Shows that a nonzero fifth-order term is possible for Hermitian generators"
        },
        "checks":checks,
        "checks_passed":sum(checks.values()),
        "checks_total":len(checks),
        "domains":[
            "Schur expression requires z outside spec(D); projected full resolvent additionally requires z outside spec(H).",
            "A pole of Sigma at d is not a new full-system eigenvalue; the tested d is outside spec(H).",
            "Evolution time is dimensionless. This is a direct-sum amplitude projection, not a density-matrix partial trace."
        ],
        "sympy_version":sp.__version__,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("general_feedback_identity_results.json"))
    args=parser.parse_args()
    result=run()
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"checks_passed":result["checks_passed"],"checks_total":result["checks_total"],"output":str(args.output)},indent=2))
