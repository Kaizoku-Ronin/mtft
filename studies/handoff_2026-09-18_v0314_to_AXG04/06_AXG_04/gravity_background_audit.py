"""AXG-04 exact unwarped Einstein--Yang--Mills--Higgs background audit.

Conditional action:
  L6 = M6^4 R6/2 - U - sum_i Tr(F_i^2)/(4 g6_i^2)
       - |D H|^2,   U = Lambda6 + V(H0).
Assume a smooth closed genus-13 curve, constant-curvature internal metric,
covariantly constant scalar vacuum, central internal magnetic flux, and an
unwarped maximally symmetric four-dimensional factor. No extra sources.

The proof is restricted to this action/ansatz, not all MTFT completions.
Run with Python and SymPy. No package import or network is required.
"""

from pathlib import Path
import json
import sympy as s


CHECKS = []


def check(name, condition):
    ok = bool(condition)
    CHECKS.append({"name": name, "passed": ok})
    if not ok:
        raise AssertionError(name)


def eq(name, lhs, rhs=0):
    check(name, s.simplify(lhs-rhs) == 0)


def main():
    M, rF, R, u, beta = s.symbols("M rho_F R u beta", positive=True)
    U, lam, k = s.symbols("U lambda4 k", real=True)
    # M denotes M6^4, not M6; rho_F is positive local magnetic energy.
    # Stress tensor follows from a 2D antisymmetric flux matrix.
    f, g6 = s.symbols("f g6", positive=True)
    Fmat = s.Matrix([[0, f], [-f, 0]])
    Fsq = sum(v*v for v in Fmat)
    flux_energy = Fsq/(4*g6**2)
    Tint = Fmat*Fmat.T/g6**2 - s.eye(2)*flux_energy
    check("2D magnetic stress is plus its energy density",
          s.simplify(Tint - s.eye(2)*flux_energy) == s.zeros(2))
    T4 = -(U+rF)
    T2 = rF-U
    traceT = 4*T4+2*T2
    eq("six-dimensional stress trace", traceT, -6*U-2*rF)
    eq("trace-reversed external Einstein equation",
       (T4-traceT/4)/M, (U-rF)/(2*M))
    eq("trace-reversed internal Einstein equation",
       (T2-traceT/4)/M, (U+3*rF)/(2*M))
    sol = s.solve([M*(-lam-k)-T4, M*(-2*lam)-T2], [lam, k])
    eq("direct Einstein external solution", sol[lam], (U-rF)/(2*M))
    eq("direct Einstein internal solution", sol[k], (U+3*rF)/(2*M))
    Ustar = -2*M/R**2-3*rF
    lamstar = -1/R**2-2*rF/M
    eq("negative internal curvature fixes U", sol[k].subs(U, Ustar), -1/R**2)
    eq("negative internal curvature forces AdS", sol[lam].subs(U, Ustar), lamstar)
    check("AdS expression is strictly negative", lamstar.is_negative is True)
    eq("Minkowski would require positive internal curvature",
       sol[k].subs(U, rF), 2*rF/M)
    L2 = 3/(1/R**2+2*rF/M)
    check("AdS curvature radius lacks parametric separation",
          s.simplify(3*R**2-L2).is_positive is True)

    # Flux quantization and action normalization.
    genus = 13
    area = 4*s.pi*(genus-1)*R**2
    eq("closed genus-13 hyperbolic area", area, 48*s.pi*R**2)
    eq("Gauss-Bonnet scalar curvature integral", area*(-2/R**2), -96*s.pi)
    N = (3, 2, 1, 1, 1)
    m = (0, -3, 3, 3, 0)
    weighted_m2 = sum(n*x*x for n, x in zip(N, m))
    check("M1 trace-weighted flux norm", weighted_m2 == 36)
    rho_common = 2*s.pi**2*weighted_m2/(g6**2*area**2)
    eq("common-coupling magnetic energy in stated trace convention",
       rho_common, 1/(32*g6**2*R**4))
    beta_common = s.simplify(rho_common*R**4)
    eq("beta coefficient", beta_common, 1/(32*g6**2))
    # C3X uses a distinct direct-product action with only one U1_X flux,
    # not the trace norm36 of either five-stack control.
    gX=s.symbols("gX6",positive=True)
    rho_C3X=2*s.pi**2/(gX*gX*area*area)
    eq("C3X unit-flux magnetic energy",rho_C3X,1/(1152*gX*gX*R**4))
    eq("C3X beta",rho_C3X*R**4,1/(1152*gX*gX))

    # Fixed integer flux: rho_F=beta/R^4. Negative U=-u gives one R^2>0.
    t = s.symbols("t", positive=True)
    radius_poly = u*t**2-2*M*t-3*beta
    tstar = (M+s.sqrt(M**2+3*u*beta))/u
    tminus = (M-s.sqrt(M**2+3*u*beta))/u
    eq("radius quadratic follows from internal Einstein equation",
       s.expand(2*M*t**2*(sol[k].subs({U:-u,rF:beta/t**2})+1/t)),
       -radius_poly)
    eq("positive radius root", radius_poly.subs(t,tstar))
    check("positive radius root has positive sign", tstar.is_positive is True)
    negative_form = -3*beta/(M+s.sqrt(M**2+3*u*beta))
    eq("negative root rationalization",tminus,negative_form)
    check("other root is negative",negative_form.is_negative is True)
    eq("zero-flux radius limit", tstar.subs(beta,0), 2*M/u)

    # Einstein-frame reduction: ds6^2=e^-2sigma ds4E^2+e^2sigma dsX0^2.
    sigma = s.symbols("sigma", real=True)
    A0, R0, MP2 = s.symbols("A0 R0 MP2", positive=True)
    a, b, c = s.symbols("a b c", real=True)
    V = a*s.exp(-4*sigma)+b*s.exp(-6*sigma)+c*s.exp(-2*sigma)
    cstar = -2*a-3*b
    eq("radion stationarity condition", s.diff(V,sigma).subs({sigma:0,c:cstar}))
    eq("stationary Einstein-frame potential", V.subs({sigma:0,c:cstar}), -a-2*b)
    eq("radion coordinate Hessian", s.diff(V,sigma,2).subs({sigma:0,c:cstar}),8*a+24*b)
    # For n=2: Lkin=-2 MP^2 (d sigma)^2 and phi=2 MP sigma.
    eq("canonical radion normalization n=2", s.Rational(2*(2+2),2),4)
    mrad2 = s.simplify((8*a+24*b)/(4*MP2))
    physical_subs = {a:M*A0/R**2,b:A0*rF,MP2:M*A0}
    physical_mrad2 = s.simplify(mrad2.subs(physical_subs))
    eq("canonical radion squared mass",physical_mrad2,2/R**2+6*rF/M)
    check("radion squared mass positive",physical_mrad2.is_positive is True)
    eq("4D Einstein equation equals direct 6D result",
       ((-a-2*b)/MP2).subs(physical_subs),lamstar)
    eq("4D stationary condition equals internal 6D equation",
       cstar.subs(physical_subs)/A0,Ustar)
    z = s.symbols("z", positive=True)
    mL2 = s.simplify((physical_mrad2*L2).subs(rF,z*M/R**2))
    eq("dimensionless radion mass in AdS",mL2,6*(1+3*z)/(1+2*z))
    check("radion mass above lower endpoint",s.simplify(mL2-6).is_positive is True)
    check("radion mass below upper endpoint",s.simplify(9-mL2).is_positive is True)

    # Potential monotonicity without a negative constant potential.
    ap,bp,cp = s.symbols("ap bp cp",positive=True)
    positive_derivative = s.diff(V,sigma).subs({a:ap,b:bp,c:cp})
    check("positive potential terms cannot stabilize radius",
          positive_derivative.is_negative is True)

    # Explicit arbitrary-units control: M=beta=1,U=-5 gives R=1.
    control = {M:1,beta:1,u:5}
    eq("exact control radius",tstar.subs(control),1)
    eq("exact control external curvature",sol[lam].subs({M:1,U:-5,rF:1}),-3)
    eq("exact control internal curvature",sol[k].subs({M:1,U:-5,rF:1}),-1)
    eq("exact control radion mass",physical_mrad2.subs({M:1,R:1,rF:1}),8)
    eq("exact control AdS radius",L2.subs({M:1,R:1,rF:1}),1)

    # Covariant constancy is stronger than a coordinate-constant Higgs.
    qHu = (0,1,-1,0,0)
    qHd = (0,1,0,-1,0)
    hu_flux = sum(x*y for x,y in zip(qHu,m))
    hd_flux = sum(x*y for x,y in zip(qHd,m))
    check("Hu cannot have nonzero covariantly constant section",hu_flux == -6)
    check("Hd cannot have nonzero covariantly constant section",hd_flux == -6)
    # The new AXG-04 candidate reverses the a,b flux signs. Its magnetic
    # stress is unchanged, but a constant Higgs is no longer excluded by
    # degree. It still requires identification of the actual connections.
    mprime = (0,-3,-3,-3,0)
    weighted_mprime2 = sum(n*x*x for n,x in zip(N,mprime))
    hu_flux_prime = sum(x*y for x,y in zip(qHu,mprime))
    hd_flux_prime = sum(x*y for x,y in zip(qHd,mprime))
    check("candidate flux has identical magnetic stress",weighted_mprime2 == weighted_m2)
    check("candidate Hu passes degree-zero necessity",hu_flux_prime == 0)
    check("candidate Hd passes degree-zero necessity",hd_flux_prime == 0)
    # At fixed area, curvature, top-form magnetic energy and constant U
    # are shape-independent. Thus the pure background sector cannot
    # select a point in the 6g-6 real-dimensional Teichmuller space.
    check("real hyperbolic shape dimension at genus 13",6*genus-6 == 72)

    result = {
      "study":"AXG-04 unwarped gravity background gate",
      "classification":"Exact conditional action-and-ansatz result; no complete vacuum construction",
      "action":"M6^4 R6/2 - U - sum Tr(F_i^2)/(4 g6_i^2) - |DH|^2",
      "trace_convention":"Tr(I_N^2)=N; differs by factor two in gauge coefficient from AXG02 kappa_i=2N_i/g6^2",
      "assumptions":["smooth closed genus-13 curve","Gaussian curvature -1/R^2","unwarped maximally symmetric M4 times X2","purely internal central constant magnetic curvature","covariantly constant scalar vacuum at a critical point of V","positive Einstein, gauge and scalar kinetic coefficients","no added sources or higher-derivative terms"],
      "symbols":{"M":"M6^4","lambda4":"R_mu_nu=lambda4*g_mu_nu; in AdS lambda4=-3/L4^2","k":"R_mn=k*g_mn=-g_mn/R^2","rho_F":"positive magnetic energy density in six dimensions","U":"Lambda6+V(H0), constant on the vacuum"},
      "einstein_equations":{"external":str(sol[lam]),"internal":str(sol[k]),"required_U":str(Ustar),"external_at_negative_internal_curvature":str(lamstar)},
      "flux":{"N":N,"original_m":m,"candidate_m":mprime,"weighted_sum_Nm2_both":weighted_m2,"area":str(area),"general_beta":"sum_i N_i*m_i^2/g6_i^2 / 1152","common_coupling_rho_F":str(rho_common),"common_coupling_beta":str(beta_common)},
      "radius":{"condition":"U=-u<0","R_squared":str(tstar),"positive_roots":1,"magnetic_energy":"beta/R^4","AdS_radius_squared":str(L2),"scale_separation":"L4^2<3 R^2 for nonzero magnetic energy"},
      "C3X_specialization":{"group":"SU3 x SU2 x U1_h x U1_X","only_flux":"one unit of U1_X, Higgs X charge zero","Maxwell_normalization":"-F_X^2/(4 gX6^2)","rho_F":str(rho_C3X),"beta":str(s.simplify(rho_C3X*R**4)),"scope":"Same Einstein/radion formulas with this beta. On the maximally symmetric product and zero other gauge backgrounds the GS four-form sources vanish as differential forms; a zero tensor background is compatible with its local equations."},
      "radion":{"metric_ansatz":"ds6^2=exp(-2 sigma)ds4E^2+exp(2 sigma)dsX0^2","R":"R0 exp(sigma)","MP_squared":"M6^4 A0","canonical_field":"phi=2 MP sigma","potential":str(V),"coefficients":{"a":"M6^4 A0/R0^2","b":"2*pi^2 sum_i(N_i*m_i^2/g6_i^2)/A0","c":"A0 U"},"at_reference_radius_equal_vacuum":{"c":"-2a-3b","V":"-a-2b","mass_squared":str(physical_mrad2),"mass_squared_times_AdS_radius_squared":str(mL2)}},
      "control_arbitrary_units":{"M":1,"beta":1,"U":-5,"R_squared":1,"lambda4":-3,"AdS_radius_squared":1,"canonical_radion_mass_squared":8,"status":"chosen action parameters, not derived MTFT values"},
      "higgs":{"original_Hu_flux_degree":hu_flux,"original_Hd_flux_degree":hd_flux,"candidate_Hu_flux_degree":hu_flux_prime,"candidate_Hd_flux_degree":hd_flux_prime,"result":"Original M1 excludes nonzero covariantly constant Hu,Hd. Candidate sign-reversed a,b flux passes the degree-zero necessity, but matching line bundles and flat holonomies/connections are still required. Nonconstant profiles require gradient stress and their full field equations."},
      "shape_moduli":{"real_dimension":72,"result":"The stated curvature+central-top-form-flux+constant-potential sector is independent of shape at fixed area; arithmetic complex structure is not selected."},
      "limits":["Radion positivity is not full scalar/KK stability","No Einstein action or Planck/gauge scale is derived from arithmetic","Changing warp factors, localized sources, matter gradients or corrections requires a new field-equation analysis","No discrete-anomaly or six-dimensional anomaly cancellation follows from the background equations"],
      "sources":[{"url":"https://arxiv.org/pdf/hep-th/0405173","location":"sections 2 and 3","use":"Primary flux compactification framework and radion/AdS cross-check"},{"url":"https://people.math.harvard.edu/~ctm/papers/home/text/class/notes/rs/course.pdf","location":"section 3.2","use":"Fenchel-Nielsen coordinates and dimension of Teichmuller space"}],
      "checks":CHECKS,"checks_passed":sum(c["passed"] for c in CHECKS),"checks_failed":sum(not c["passed"] for c in CHECKS),
    }
    target=Path(__file__).with_name("gravity_background_results.json")
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"checks_passed":result["checks_passed"],"checks_failed":result["checks_failed"],"output":str(target)}))


if __name__ == "__main__":
    main()
