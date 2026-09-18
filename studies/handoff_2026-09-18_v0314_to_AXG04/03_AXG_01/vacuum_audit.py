#!/usr/bin/env python3
"""AXG-01: exact, deliberately restricted vacuum checks.

Requires SymPy. Run with Python 3; writes vacuum_results.json beside this file.
No MTFT package files are modified. These are conditional calculations, not a
derivation of an axion sector, supersymmetric completion, or vacuum of MTFT.
"""
from pathlib import Path
import json
import sympy as s


checks = []


def check(name, condition, detail=""):
    passed = bool(condition)
    checks.append({"name": name, "passed": passed, "detail": detail})
    if not passed:
        raise AssertionError(name)


def zero(expr):
    return s.simplify(expr) == 0


rho, a = s.symbols("rho a", positive=True)
b, c = s.symbols("b c", nonnegative=True)
x = s.symbols("x", positive=True)
V = a / rho**4 + b / rho**6 + c / rho**2
derivative = s.diff(V, rho)
check("inherited_derivative", zero(derivative + 4*a/rho**5 + 6*b/rho**7 + 2*c/rho**3))
check("strict_monotonicity_for_assumed_signs", (-derivative).is_positive is True)
check("inverse_square_coordinate", zero(V.subs(rho, x**s.Rational(-1,2)) - (c*x+a*x**2+b*x**3)))

# The displayed integer matrix is a candidate charge map, not parent-derived.
K = s.Matrix([[0, -2, 1, 1, 0], [3, -4, 0, 0, 1]])
check("candidate_K_rank_two", K.rank() == 2)
check("surjective_real_shift_submatrix", K[:, [2,4]].det() == 1)
dv1, dv2 = s.symbols("partial_a1_V partial_a2_V")
invariance = K.T*s.Matrix([dv1,dv2])
check("pure_axion_potential_gradient_zero", s.solve(list(invariance), (dv1,dv2)) == {dv1:0,dv2:0})
f1,f2,f12 = s.symbols("G11 G22 G12", real=True)
G = s.Matrix([[f1,f12],[f12,f2]])
A = s.Matrix(s.symbols("A0:5", real=True))
da = s.Matrix(s.symbols("da0:2", real=True))
D = da + K*A
kinetic = (D.T*G*D)[0]/2
vacuum_subs = {symbol:0 for symbol in list(A)+list(da)}
check("stueckelberg_kinetic_vanishes_on_background", zero(kinetic.subs(vacuum_subs)))
check("radius_force_from_G_vanishes", all(zero(s.diff(kinetic,g).subs(vacuum_subs)) for g in [f1,f2,f12]))
mass = K.T*G*K
check("mass_matrix_is_kinetic_Hessian", s.hessian(kinetic,list(A)) == mass)
check("positive_metric_example_mass_rank_two", mass.subs({f1:1,f2:1,f12:0}).rank() == 2)
check("common_phase_massless", mass*s.ones(5,1) == s.zeros(5,1))

# Flavor-blind extra-U(1) D-term under the explicitly stated Kähler assumption.
h = s.Matrix(s.symbols("h1:4", real=True))
q, xi, gg, v = s.symbols("q xi g v", real=True)
norm = (h.T*h)[0]
VD = gg**2*(xi+q*norm)**2/2
R = s.Matrix([[s.Rational(3,5),-s.Rational(4,5),0],
              [s.Rational(4,5),s.Rational(3,5),0],[0,0,1]])
check("test_family_rotation_orthogonal", R.T*R == s.eye(3))
check("family_rotation_preserves_norm", zero(((R*h).T*(R*h))[0]-norm))
check("D_term_preserves_family_rotation", zero(VD.subs(dict(zip(list(h),list(R*h))), simultaneous=True)-VD))
stationary = {h[0]:v,h[1]:0,h[2]:0,xi:-q*v**2}
check("nonzero_D_flat_sample_stationary", all(zero(s.diff(VD,t).subs(stationary)) for t in h))
H = s.hessian(VD,list(h)).subs(stationary)
check("D_flat_hessian_radial_only", H == s.diag(4*gg**2*q**2*v**2,0,0))
check("two_real_tangent_directions_flat", H*s.Matrix([0,1,0]) == s.zeros(3,1) and H*s.Matrix([0,0,1]) == s.zeros(3,1))

# Positive control: change one sign explicitly; this is not an MTFT parameter.
C = s.symbols("C", positive=True)  # c = -C < 0
bp = s.symbols("b_positive", positive=True)
Vx = -C*x+a*x**2+bp*x**3
xstar = C/(a+s.sqrt(a*a+3*bp*C))
check("negative_c_control_stationary", zero(s.diff(Vx,x).subs(x,xstar)))
check("negative_c_control_root_positive", xstar.is_positive is True)
check("negative_c_control_other_root_negative", ((-a-s.sqrt(a*a+3*bp*C))/(3*bp)).is_negative is True)
check("negative_c_b_zero_root", zero((-C+2*a*x).subs(x,C/(2*a))))
check("negative_c_vacuum_value_identity", zero(Vx.subs(C,2*a*x+3*bp*x**2) - (-a*x**2-2*bp*x**3)))
check("negative_c_vacuum_value_strict_negative", (-a*x**2-2*bp*x**3).is_negative is True)
Vcontrol = a/rho**4+bp/rho**6-C/rho**2
stationary_C = 2*a/rho**2+3*bp/rho**4
control_hessian = s.diff(Vcontrol,rho,2).subs(C,stationary_C)
check("negative_c_radius_hessian", zero(control_hessian-(8*a/rho**6+24*bp/rho**8)))
check("negative_c_radius_hessian_positive", control_hessian.is_positive is True)
check("negative_c_b_zero_vacuum_value", zero((-C*x+a*x**2).subs(x,C/(2*a)) + C**2/(4*a)))
check("negative_c_b_zero_radius", zero(s.diff(a/rho**4-C/rho**2,rho).subs(rho,s.sqrt(2*a/C))))
check("negative_c_global_limit_small_radius", s.limit(Vcontrol,rho,0,dir="+") == s.oo)
check("negative_c_global_limit_large_radius", s.limit(Vcontrol,rho,s.oo) == 0)

# Optional two-power control illustrates why sign and exponent both matter.
n,p,AA,dd = s.symbols("n p A d", positive=True)
W = AA*rho**(-n)-dd*rho**(-p)
dstationary = (n*AA/p)*rho**(p-n)
check("two_power_stationarity_relation", zero(s.diff(W,rho).subs(dd,dstationary)))
check("two_power_vacuum_value", zero(W.subs(dd,dstationary) - AA*(p-n)/p*rho**(-n)))
check("two_power_second_derivative", zero(s.diff(W,rho,2).subs(dd,dstationary) - n*(n-p)*AA*rho**(-n-2)))

result = {
  "investigation": "AXG-01 vacuum subaudit",
  "status": "EXACT conditional algebra; no parent or vacuum derived",
  "assumptions": [
    "rho > 0; inherited Einstein-frame potential V=a*rho^-4+b*rho^-6+c*rho^-2, a>0, b,c>=0",
    "minimal extension contains only gauge kinetic terms and -G_IJ(rho) Da^I Da^J/2; G positive definite",
    "Lorentz-invariant classical vacuum uses A_mu=0, constant axions, and no gauge/axion field-strength background",
    "two axions shift by a rank-two real map K from five gauge parameters",
    "D-term flavor statement assumes a flavor-symmetric quadratic Kahler metric after canonical normalization, identical charge vectors within each family sector, and no additional family-sensitive interactions"
  ],
  "results": {
    "minimal_extension": {
      "classical_added_radius_potential": "0",
      "reason": "all covariant derivatives and field strengths vanish on the specified classical background; rho-dependent masses are not by themselves a potential",
      "pure_axion_local_potential": "constant in both axions: K^T grad_a V = 0 and rank(K)=2",
      "charged_dressing_exception": "gauge-invariant axion exponentials multiplied by appropriately charged matter operators are allowed; an example checked separately in AXG-01 is an axion-dressed Hu.Hd bilinear. Its flavor coefficients require additional parent input and are not fixed by the rank argument",
      "K": [list(K.row(i)) for i in range(2)],
      "local_axion_count": "both axions are eaten; global discrete gauge identifications require separate integral-lattice treatment",
      "remaining_runaway": "dV/drho=-4a/rho^5-6b/rho^7-2c/rho^3<0; no finite-radius critical point"
    },
    "D_term": {
      "general_formula": "V_D=1/2 D_a (Re f)^(-1)^{ab} D_b; D_a=xi_a(rho)+sum_s q_{a,s} ||H_s||^2 under the stated metric assumption",
      "selection": "can constrain radial norms, charge-sector competition, and rho if xi/f are supplied; does not select orientations among identical-charge family copies",
      "flat_sample": "V_D=g^2(xi+q(h1^2+h2^2+h3^2))^2/2; at (v,0,0), xi=-qv^2: Hessian=diag(4g^2q^2v^2,0,0)",
      "SUSY_caveat": "a supersymmetric completion can supply a saxion/Kahler moment-map potential even when A=0 and axions are constant; that potential is additional input absent from the minimal bosonic kinetic ansatz",
      "flavor_caveat": "higher-order non-universal Kahler terms, F-terms, Yukawas, soft terms, or nonperturbative interactions can break the flavor symmetry; they must be specified rather than assumed absent or effective"
    },
    "negative_c_control": {
      "status": "hypothetical sign change only; not a physical input established in MTFT",
      "definition": "c=-C<0, a>0, b>=0, x=rho^-2",
      "stationary_b_positive": "x*=C/(a+sqrt(a^2+3bC)); rho*=1/sqrt(x*)",
      "stationary_b_zero": "x*=C/(2a); rho*=sqrt(2a/C)",
      "vacuum_energy": "V*=-a*x*^2-2b*x*^3<0",
      "radius_hessian": "V''(rho*)=8a/rho*^6+24b/rho*^8>0",
      "interpretation": "unique global minimum along rho; with positive kinetic normalization it is stable in this direction; with 4D Einstein gravity and other fields stationary its negative potential gives AdS, not dS or Minkowski",
      "limitations": "does not establish stability of any other scalar direction, a controlled compactification, quantized flux compatibility, or a consistent negative-c source"
    },
    "two_power_control": {
      "definition": "W=A*rho^-n-d*rho^-p; A,d,n,p>0, p!=n",
      "stationarity": "rho*^(n-p)=n*A/(p*d)",
      "value": "W*=A*(p-n)/p*rho*^-n",
      "hessian": "W''*=n*(n-p)*A*rho*^(-n-2)",
      "classification": "n>p gives stable negative minimum; n<p gives unstable positive maximum; equal exponents give a single power or identically flat tuned potential",
      "scope": "two terms only; do not apply this classification to an unrestricted many-term potential"
    }
  },
  "necessary_next_data": [
    "parent-derived axion origin, periods, integral charge lattice, and anomaly-canceling couplings",
    "full gauge kinetic matrix f_ab(rho,other moduli), axion metric G_IJ, and scalar kinetic normalization",
    "if supersymmetry is assumed: actual Kahler potential, gauged Killing vectors/moment maps, superpotential and supersymmetry-breaking inputs",
    "charged Higgs/other scalar content, vacuum ansatz, and all family-sensitive interactions",
    "origin and signs of every classical potential contribution; any loop calculation also needs spectrum, regulator/counterterms and renormalization conditions",
    "all-moduli Hessian, scale separation, and higher-dimensional equations before claiming a gravitational vacuum"
  ],
  "checks": checks,
  "summary": {"passed":len(checks),"failed":0}
}

out = Path(__file__).with_name("vacuum_results.json")
out.write_text(json.dumps(result,indent=2,default=str)+"\n")
print(json.dumps(result["summary"]))
print(str(out))
