"""Independent exact exponent ledger and numerical mode-basis checks."""
from pathlib import Path
import json
import math
import numpy as np
import sympy as sp
from mtft.matter.normalization import (
    canonical_u1_vertex, breathing_gauge_kinetics, electroweak_matching,
    electromagnetic_invariants,
)

ROOT=Path(__file__).resolve().parents[1]


def main():
    d,sigma,x,q,g=sp.symbols('d sigma x q g',positive=True)
    exact=[]
    def gate(name, expression):
        residual=sp.simplify(expression)
        assert residual==0,(name,residual)
        exact.append({'name':name,'class':'EXACT conditional algebra','residual':str(residual)})
    # Four external dimensions: external sqrt(det) e^-2dσ, internal e^dσ.
    gate('Einstein-frame volume times two external inverse metrics',-2*d+d+2*d-d)
    gate('canonical breathing exponent',d/sp.sqrt(d*(d+2)/2)-sp.sqrt(2*d/(d+2)))
    gate('C3X d=2 gauge coefficient exponent',sp.sqrt(2*2/sp.Integer(2+2))-1)
    gate('B3 d=4 uniform breathing coefficient exponent',sp.sqrt(2*4/sp.Integer(4+2))-sp.sqrt(sp.Rational(4,3)))
    # Explicit rational nonorthogonal mode matrix: exact K=B†B and C=B^-1.
    b=sp.Matrix([[1,sp.I],[2,3]])
    k=b.conjugate().T*b;c=b.inv()
    matrix=c.conjugate().T*k*c-sp.eye(2)
    for i in range(2):
        for j in range(2):gate(f'exact kinetic normalization entry {i}{j}',matrix[i,j])
    matrix=c.conjugate().T*(q*g*k)*c-q*g*sp.eye(2)
    for i in range(2):
        for j in range(2):gate(f'exact U1 charge vertex entry {i}{j}',matrix[i,j])
    g2,gy=sp.symbols('g2 gy',positive=True)
    e=g2*gy/sp.sqrt(g2*g2+gy*gy)
    gate('electroweak inverse coupling matching',1/e**2-1/g2**2-1/gy**2)
    rng=np.random.default_rng(14301)
    residuals=[]
    for n in (2,3,6):
        b=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))+4*np.eye(n)
        gram=b.conj().T@b
        for _ in range(10):
            a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))+4*np.eye(n)
            norm=canonical_u1_vertex(a.conj().T@gram@a,-1,0.3)
            residuals.append(max(norm['kinetic_residual'],norm['vertex_residual']))
    assert max(residuals)<2e-11
    null=electromagnetic_invariants([1,0,0],[0,1,0])
    assert null['F_squared']==0 and null['energy_density']==1
    weak=electroweak_matching(1,math.sqrt(3/5))
    assert abs(weak['sin2_theta_w']-3/8)<1e-15
    breathing_samples=[]
    for dimension in (2,4):
        exponent=math.sqrt(2*dimension/(dimension+2))
        for value in (-1.0,0.0,0.1,1.0):
            result=breathing_gauge_kinetics(value,dimension)
            expected_z=math.exp(exponent*value)
            expected_g=math.exp(-exponent*value/2)
            assert math.isclose(result['Z'],expected_z,rel_tol=2e-15)
            assert math.isclose(result['g_eff_over_g0'],expected_g,rel_tol=2e-15)
            assert math.isclose(result['mp_dZ_dchi'],exponent*expected_z,rel_tol=2e-15)
            assert result['implemented_evolution'] is False
            breathing_samples.append({'internal_dimension':dimension,'chi_over_mp':value,
                                      'Z':result['Z'],'g_eff_over_g0':result['g_eff_over_g0'],
                                      'mp_dZ_dchi':result['mp_dZ_dchi'],'status':'PASS'})
    out={'study':'Matter-01 U1 conditional reduction','exact_gates':exact,
         'numerical_class':'DIAGNOSTIC, finite basis-change sample',
         'basis_change_samples':len(residuals),'max_basis_residual':max(residuals),
         'C3X_breathing':breathing_gauge_kinetics(0.1,2),
         'B3_uniform_breathing':breathing_gauge_kinetics(0.1,4),
         'breathing_implementation_samples':breathing_samples,
         'null_plane_wave':null,'unified_input_matching':weak,
         'open':['unbroken massless photon in selected MTFT vacuum',
                 'parent coupling and physical scale','electron identification and mass',
                 'coupled gauge-fermion-breathing evolution','Einstein metric feedback']}
    path=ROOT/'results/u1_reduction_checks.json';path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'exact_gates':len(exact),'basis_samples':len(residuals),'max_basis_residual':max(residuals)}))


if __name__=='__main__':main()
