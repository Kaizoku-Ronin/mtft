"""Rational CRT projectors and exact symmetry-selection certificates."""
import sympy as sp
from common import X, POLYS, BLOCKS, DIMS, SIGNS, frozen, dump, matrix_eval, rational_matrix

def main():
    d = frozen()
    T, E, I, zero = sp.Matrix(d['T2']), sp.Matrix(d['intersection_cycles']), sp.eye(26), sp.zeros(26)
    W11, W13, T3 = (sp.Matrix(d[k]) for k in ['W11','W13','T3'])
    factors = [sp.Poly(X,X), sp.Poly(X+2,X), POLYS['q4'], POLYS['q6']]
    m = sp.Poly(sp.prod(f.as_expr() for f in factors), X)
    assert matrix_eval(m,T) == zero and sp.gcd(m,m.diff()).degree() == 0
    assert sp.Poly(T.charpoly(X).as_expr(),X) == sp.Poly(X**2*(X+2)**4*POLYS['q4'].as_expr()**2*POLYS['q6'].as_expr()**2,X)
    assert E.T == -E and abs(E.det()) == 1
    for W in [W11,W13]:
        assert W*W == I and W.T*E*W == E
        assert W*T == T*W and W*T3 == T3*W
    assert W11*W13 == W13*W11 == sp.Matrix(d['W143'])
    assert T*T3 == T3*T and T.T*E == E*T and T3.T*E == E*T3
    Ps, crt = {}, {}
    for name, f, dim in zip(BLOCKS, factors, DIMS):
        q = m.exquo(f)
        polynomial = (q*sp.invert(q,f)).rem(m)
        P = matrix_eval(polynomial,T)
        assert P*P == P and P.trace() == int(dim) and P.T*E == E*P
        assert P*T3 == T3*P and P*W11 == W11*P and P*W13 == W13*P
        Ps[name], crt[name] = P, str(polynomial.as_expr())
    assert sum(Ps.values(),zero) == I
    for a in BLOCKS:
        for b in BLOCKS:
            assert Ps[a]*Ps[b] == (Ps[a] if a==b else zero)
    # Independent projector construction from the packaged integral block bases.
    bases = [sp.Matrix(d['block_ghost' if k=='old' else f'block_{k}']) for k in BLOCKS]
    B = sp.Matrix.hstack(*bases); Bi = B.inv(); offset = 0
    for name,dim in zip(BLOCKS,DIMS):
        diagonal = sp.zeros(26)
        for j in range(offset,offset+int(dim)):
            diagonal[j,j] = 1
        assert B*diagonal*Bi == Ps[name]
        offset += int(dim)
    Qs = [(I+a*W11)*(I+b*W13)/4 for a,b in SIGNS]
    ranks = []
    for name in BLOCKS:
        ranks.append([int((Ps[name]*Q).trace()) for Q in Qs])
        for Q in Qs:
            assert (Ps[name]*Q)**2 == Ps[name]*Q
    allowed = [[any(ranks[i][s]>0 and ranks[j][s]>0 for s in range(4)) for j in range(4)] for i in range(4)]
    assert ranks == [[2,0,0,0],[0,0,2,2],[0,0,8,0],[0,12,0,0]]
    # An exact counterexample: block-trace balance does not enforce invariance.
    p = sp.diag(1,0); v = sp.Matrix([[0,1],[1,0]])
    assert p*v*p == sp.zeros(2) and (sp.eye(2)-p)*v*(sp.eye(2)-p) == sp.zeros(2)
    assert p*v-v*p != sp.zeros(2)
    dump('sector_certificate.json', {
        'block_order':BLOCKS, 'real_dimensions':[int(v) for v in DIMS],
        'minimal_polynomial':str(m.as_expr()), 'CRT_polynomials':crt,
        'projectors':{k:rational_matrix(v) for k,v in Ps.items()},
        'AL_sign_order':[list(s) for s in SIGNS], 'AL_intersection_real_dimensions':ranks,
        'AL_allowed_direct_coupling':allowed,
        'checks':{'minimal_polynomial_annihilates':True,'characteristic_polynomial':True,
                  'idempotent_complete_orthogonal_projectors':True,'exact_ranks':True,
                  'independent_integral_basis_projectors':True,'Hecke_and_AL_commutation':True,
                  'exact_intersection_form_and_adjointness':True,'AL_character_intersections':True,
                  'balanced_but_mixing_counterexample':True},
        'matrix_recursion_lemma':{
            'assumptions':'For each sector projector P, [P,Z0]=[P,C]=0.',
            'identity':'[P,Z^2+C] = [P,Z] Z + Z [P,Z] + [P,C]',
            'conclusion':'Induction gives [P,Zk]=0 for every nonnegative integer k. If Z0,C belong to Q[T2], all iterates belong to Q[T2]. This matrix lift is a proposed model of the scalar image, not an identification with Mandelbrot dynamics.',
            'counterexample':{'P':[[1,0],[0,0]],'V':[[0,1],[1,0]],
                              'property':'Both diagonal blocks vanish, yet [P,V] is nonzero. In Schrödinger evolution exp(-it V), the transfer probability is sin(t)^2.'}}})
    print('Exact AL rank table:',ranks,flush=True)

if __name__ == '__main__':
    main()
