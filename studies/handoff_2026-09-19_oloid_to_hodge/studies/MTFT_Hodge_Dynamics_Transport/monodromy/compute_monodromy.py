#!/usr/bin/env python3
"""Based complex period transport, with explicit paths and gauge check.

Requires numpy/scipy. No integral homology basis is asserted. Numerical
monodromy is integrated in the displayed algebraic de Rham frame.
"""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
I = np.eye(2, dtype=complex)
J = np.array([[0, 1], [-1, 0]], dtype=complex)


def connection(c):
    return np.array([[-1/(4*c), 1/(4*c*(c+1))], [-.25, 1/(4*c)]], dtype=complex)


def gauge(c):
    return np.diag([np.exp(.31j*c), np.exp(-.31j*c)])


def gauged_connection(c):
    g = gauge(c)
    return np.diag([.31j, -.31j]) + g @ connection(c) @ np.linalg.inv(g)


def line(a, b):
    return lambda t: (a+(b-a)*t, b-a)


def circle(center, radius=.2):
    def evaluate(t):
        e = np.exp(2j*np.pi*t)
        return center+radius*e, 2j*np.pi*radius*e
    return evaluate


def based_loop(center):
    if center == 0:
        vertices = [1+0j, .2+0j]
    elif center == -1:
        # The upper-half-plane stem specifies the based homotopy class.
        vertices = [1+0j, 1+.75j, -.8+.75j, -.8+0j]
    else:
        raise ValueError(center)
    stem = [line(a,b) for a,b in zip(vertices, vertices[1:])]
    retrace = [line(b,a) for a,b in reversed(list(zip(vertices, vertices[1:])))]
    return stem + [circle(center)] + retrace


def transport(path, rtol, transformed=False):
    mat = I.copy()
    evaluations = 0
    conn = gauged_connection if transformed else connection
    for segment in path:
        def rhs(t, flat):
            c, dc = segment(t)
            return (dc * conn(c) @ flat.reshape(2, 2)).ravel()
        sol = solve_ivp(rhs, (0.,1.), mat.ravel(), method='DOP853',
                        rtol=rtol, atol=rtol*.1)
        if not sol.success:
            raise RuntimeError(sol.message)
        mat = sol.y[:,-1].reshape(2,2)
        evaluations += sol.nfev
    return mat, evaluations


def serialize(z):
    if isinstance(z, np.ndarray):
        return [[serialize(w) for w in row] for row in z] if z.ndim == 2 else [serialize(w) for w in z]
    z = complex(z)
    return {'real':z.real, 'imag':z.imag}


def metrics(m):
    return {'matrix':serialize(m), 'trace':serialize(np.trace(m)),
            'eigenvalues':serialize(np.linalg.eigvals(m)),
            'determinant_error':float(abs(np.linalg.det(m)-1)),
            'symplectic_error':float(np.linalg.norm(m.T@J@m-J)),
            'operator_norm':float(np.linalg.norm(m,2)),
            'spectral_radius':float(max(abs(np.linalg.eigvals(m))))}


def main():
    result = {'basepoint':1, 'loop_radii':.2,
              'minus_one_stem':[[1,0],[1,.75],[-.8,.75],[-.8,0]],
              'orientation':'counterclockwise',
              'multiplication':'first loop 0 then loop -1 gives M_minus_one @ M_zero',
              'precision':'IEEE complex128; adaptive DOP853, not interval certified',
              'runs':[]}
    previous = None
    finest = None
    for tolerance in [1e-8, 1e-10, 1e-13]:
        m0, e0 = transport(based_loop(0), tolerance)
        m1, e1 = transport(based_loop(-1), tolerance)
        product = m1@m0
        hyperbolic = m1@m1@m0
        run = {'rtol':tolerance, 'evaluations':e0+e1,
               'around_zero':metrics(m0), 'around_minus_one':metrics(m1),
               'zero_then_minus_one':metrics(product),
               'zero_then_minus_one_twice':metrics(hyperbolic),
               'zero_square_plus_identity_error':float(np.linalg.norm(m0@m0+I)),
               'minus_one_nilpotent_square_error':float(np.linalg.norm((m1-I)@(m1-I))),
               'product_plus_identity_nilpotent_square_error':float(np.linalg.norm((product+I)@(product+I)))}
        if previous is not None:
            run['matrix_change_from_previous'] = max(float(np.linalg.norm(m0-previous[0])),
                                                      float(np.linalg.norm(m1-previous[1])))
        result['runs'].append(run)
        previous = (m0,m1)
        finest = (m0,m1,product,hyperbolic)
    m0,m1,product,hyperbolic = finest
    g = gauge(1)
    transformed, _ = transport(based_loop(-1), 1e-13, transformed=True)
    result['nonconstant_gauge_check'] = {
        'gauge':'G(c)=diag(exp(0.31 i c),exp(-0.31 i c))',
        'law':'A_tilde=G_prime G_inverse+G A G_inverse; M_tilde=G(1) M G(1)^(-1)',
        'monodromy_conjugacy_error':float(np.linalg.norm(transformed-g@m1@np.linalg.inv(g))),
        'trace_change':float(abs(np.trace(transformed)-np.trace(m1)))}
    result['expected_exact_character_data'] = {
        'trace_zero_loop':0, 'trace_minus_one_loop':2, 'trace_product':-2,
        'trace_hyperbolic_word':-4,
        'hyperbolic_word_eigenvalues':[-2-np.sqrt(3),-2+np.sqrt(3)],
        'hyperbolic_word_growth_per_repetition':float(np.log(2+np.sqrt(3))),
        'qualification':'Character identities follow from scalar Fuchsian equation plus specified based loop topology; numerical matrices do not identify an integral homology frame.'}
    rows=[]
    for count in [1,2,4,8,16,24,32]:
        for name,m in [('elliptic_zero',m0),('parabolic_minus_one',m1),
                       ('parabolic_alternation',product),('hyperbolic_word',hyperbolic)]:
            norm=float(np.linalg.norm(np.linalg.matrix_power(m,count),2))
            rows.append({'word':name,'repetitions':count,'operator_norm':norm,
                         'log_norm_per_repetition':float(np.log(norm)/count)})
    result['growth_caution']='Finite complex128 powers are illustrative; exact Cayley-Hamilton formulas give bounded/linear/exponential asymptotic classifications.'
    (HERE/'monodromy_results.json').write_text(json.dumps(result,indent=2)+'\n')
    with (HERE/'monodromy_growth.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps({'finest':result['runs'][-1], 'gauge':result['nonconstant_gauge_check']},indent=2))


if __name__ == '__main__':
    main()
