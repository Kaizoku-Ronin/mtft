#!/usr/bin/env python3
"""Compare inherited/new monodromy about the primitive period-three parameter."""
from pathlib import Path
from fractions import Fraction
import json
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'geometry'/'transport_results.json'


def coeffs(raw):
    return np.array([float(Fraction(v)) for v in raw])


def load_connection():
    data=json.loads(SOURCE.read_text())
    pairs=[[(coeffs(v['numerator_ascending']),coeffs(v['denominator_ascending']))
            for v in row] for row in data['connection']]
    def evaluate(c):
        return np.array([[np.polynomial.polynomial.polyval(c,n)/
                          np.polynomial.polynomial.polyval(c,d) for n,d in row]
                         for row in pairs],dtype=complex)
    return evaluate


def inherited(c):
    return np.array([[-1/(4*c),1/(4*c*(c+1))],[-.25,1/(4*c)]],dtype=complex)


def line(a,b):
    return lambda t:(a+(b-a)*t,b-a)


def circle(a,r):
    def evaluate(t):
        e=np.exp(2j*np.pi*t)
        return a-r*e,-r*2j*np.pi*e
    return evaluate


def integrate(connection,dimension,path,tolerance):
    mat=np.eye(dimension,dtype=complex)
    nfev=0
    for segment in path:
        def rhs(t,y):
            c,dc=segment(t)
            return (dc*connection(c)@y.reshape(dimension,dimension)).ravel()
        sol=solve_ivp(rhs,(0,1),mat.ravel(),method='DOP853',rtol=tolerance,
                      atol=tolerance*.1)
        if not sol.success:raise RuntimeError(sol.message)
        mat=sol.y[:,-1].reshape(dimension,dimension)
        nfev+=sol.nfev
    return mat,nfev


def serialize(a):
    return [[[float(z.real),float(z.imag)] for z in row] for row in a]


def main():
    a=brentq(lambda c:c**3+2*c*c+c+1,-1.8,-1.7,xtol=5e-15)
    radius=.05;base=-2;start=a-radius
    path=[line(base,start),circle(a,radius),line(start,base)]
    new_connection=load_connection()
    result={'curve_new':'D2: v^2=(u-c)*((u^2+c)^2+c)',
            'curve_inherited':'C2: y^2=(x^2+c)^2+c',
            'source_connection':str(SOURCE.relative_to(HERE.parent)),
            'parameter_root_polynomial_ascending':[1,1,2,1],
            'parameter_root':a,'basepoint':base,'radius':radius,
            'orientation':'counterclockwise, circle starts at its leftmost point',
            'precision':'IEEE complex128 DOP853, not interval certified',
            'rank_warning':'Singular values are a numerical rank diagnostic, not an exact rank proof.',
            'runs':[]}
    previous=None
    for tol in [1e-8,1e-10,1e-13]:
        md,nd=integrate(new_connection,4,path,tol)
        me,ne=integrate(inherited,2,path,tol)
        d=md-np.eye(4)
        sv=np.linalg.svd(d,compute_uv=False)
        run={'rtol':tol,'evaluations_new':nd,'evaluations_inherited':ne,
             'new_monodromy_real_imag_pairs':serialize(md),
             'inherited_monodromy_real_imag_pairs':serialize(me),
             'new_trace_minus_four_abs':float(abs(np.trace(md)-4)),
             'new_determinant_minus_one_abs':float(abs(np.linalg.det(md)-1)),
             'new_minus_identity_norm':float(np.linalg.norm(d)),
             'new_minus_identity_singular_values':[float(v) for v in sv],
             'second_to_first_singular_value_ratio':float(sv[1]/sv[0]),
             'new_nilpotent_square_norm':float(np.linalg.norm(d@d)),
             'new_relative_nilpotent_square_norm':float(np.linalg.norm(d@d)/np.linalg.norm(d)**2),
             'inherited_minus_identity_norm':float(np.linalg.norm(me-np.eye(2)))}
        if previous is not None:
            run['new_matrix_change_from_previous_norm']=float(np.linalg.norm(md-previous))
        result['runs'].append(run);previous=md
    (HERE/'new_component_monodromy_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result['runs'][-1].items() if 'pairs' not in k},indent=2))


if __name__=='__main__':main()
