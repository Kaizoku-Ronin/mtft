#!/usr/bin/env python3
"""Exact, synthetic two-channel control for the W143 feedback experiment.

No MTFT package, fitted numbers, or measured matrices are used. The positive
rational parameter r gives d=(1-r**2)/(1+r**2) in (-1,1) and b=2*r/(1+r**2),
so every coupled real Hermitian involution channel is represented up to basis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def run() -> dict:
    r = sp.symbols("r", positive=True)
    z = sp.symbols("z")
    t, s = sp.symbols("t s", real=True)
    d = (1-r**2)/(1+r**2)
    b = 2*r/(1+r**2)
    a = -d
    W = sp.Matrix([[a, b], [b, d]])
    eye = sp.eye(2)
    checks: dict[str, bool] = {}

    def zero(name: str, expression) -> None:
        values = list(expression) if isinstance(expression, sp.MatrixBase) else [expression]
        checks[name] = all(sp.simplify(sp.trigsimp(x)) == 0 for x in values)
        if not checks[name]:
            raise AssertionError(name)

    zero("involution", W*W-eye)
    zero("hermitian", W-W.conjugate().T)
    zero("block_a_defect", a*a+b*b-1)
    zero("block_d_defect", d*d+b*b-1)
    zero("block_intertwining", a*b+b*d)

    full_resolvent = (z*eye-W).inv()
    closed_resolvent = (z*eye+W)/(z*z-1)
    sigma = b*b/(z-d)
    heff = a+sigma
    projected = (z-d)/(z*z-1)
    zero("full_resolvent_closed_form", full_resolvent-closed_resolvent)
    zero("schur_projected_resolvent", 1/(z-heff)-projected)
    zero("effective_operator_closed_form", heff-(1-d*z)/(z-d))
    zero("self_energy_pushthrough", (z+a)*sigma-(1-a*a))
    zero("schur_polynomial_identity", (z-heff)*(z+a)-(z*z-1))
    zero("full_characteristic_polynomial", (z*eye-W).det()-(z*z-1))
    zero("compression_pole_projected_zero", projected.subs(z,d))
    zero("compression_pole_effective_residue", sp.cancel((z-d)*heff).subs(z,d)-b*b)
    zero("projected_resolvent_positive_pole_residue", sp.cancel((z-1)*projected).subs(z,1)-(1-d)/2)
    zero("projected_resolvent_negative_pole_residue", sp.cancel((z+1)*projected).subs(z,-1)-(1+d)/2)

    U = sp.cos(t)*eye-sp.I*sp.sin(t)*W
    zero("unitary_evolution", U.conjugate().T*U-eye)
    zero("unitary_evolution_equation", sp.I*U.diff(t)-W*U)
    u = U[0,0]
    v = U[1,0]
    zero("active_amplitude", u-(sp.cos(t)+sp.I*d*sp.sin(t)))
    zero("complement_amplitude", v+sp.I*b*sp.sin(t))
    zero("active_norm_loss", sp.conjugate(u)*u-(1-b*b*sp.sin(t)**2))
    zero("complement_norm_gain", sp.conjugate(v)*v-b*b*sp.sin(t)**2)
    zero("norm_conservation", sp.conjugate(u)*u+sp.conjugate(v)*v-1)
    zero("compressed_semigroup_defect", u*u.subs(t,s)-u.subs(t,t+s)-b*b*sp.sin(t)*sp.sin(s))
    # This antiderivative proves the eliminated-sector memory integral exactly:
    # integral_0^t exp(-i*d*(t-s))*u(s) ds = sin(t).
    zero("memory_integral_antiderivative", sp.diff(sp.exp(sp.I*d*s)*sp.sin(s),s)-sp.exp(sp.I*d*s)*u.subs(t,s))
    zero("memory_equation", sp.I*u.diff(t)-(a*u-sp.I*b*b*sp.sin(t)))
    zero("complete_revival_at_pi", U.subs(t,sp.pi)+eye)
    zero("no_transfer_at_zero", v.subs(t,0))
    return {
        "status": "EXACT_SYMBOLIC_CONTROL",
        "scope": "Synthetic canonical 2x2 Hermitian involution; no MTFT measurements",
        "parameter_domain": "r > 0, d=(1-r^2)/(1+r^2), b=2r/(1+r^2), a=-d",
        "domain_caveats": [
            "Full resolvent identities require z != -1,+1.",
            "Literal Schur-complement evaluation additionally requires z != d.",
            "At z=d the projected full resolvent is regular and zero in this channel; the Schur expression only has a meromorphic continuation there.",
            "Unitary evolution uses W as a dimensionless toy generator; no physical time, mass, or QFT interpretation is inferred."
        ],
        "checks": checks,
        "checks_passed": sum(checks.values()),
        "checks_total": len(checks),
        "formulas": {
            "H_eff": "(1-d*z)/(z-d)",
            "G_PP": "(z-d)/(z^2-1)",
            "self_energy": "(1-d^2)/(z-d)",
            "U_PP": "cos(t)+i*d*sin(t)",
            "norm_loss": "(1-d^2)*sin(t)^2",
            "semigroup_defect": "(1-d^2)*sin(t)*sin(s)",
            "memory_kernel": "(1-d^2)*exp(-i*d*t)"
        },
        "sympy_version": sp.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("feedback_identity_results.json"))
    args=parser.parse_args()
    results=run()
    args.output.write_text(json.dumps(results, indent=2)+"\n")
    print(json.dumps({"checks_passed":results["checks_passed"],"checks_total":results["checks_total"],"output":str(args.output)}, indent=2))
