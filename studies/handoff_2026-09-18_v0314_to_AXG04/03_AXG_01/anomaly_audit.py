#!/usr/bin/env python3
"""AXG-01: exact 4D anomaly/Stueckelberg algebra for the declared M1 spectrum.

Run with Python + sympy (+ numpy if checking the copied package source):
    python anomaly_audit.py

No model parameter is fitted, and no 6D completion is asserted. Anomaly
normalization: left Weyl fermions, T(fundamental)=1/2, f=F/(2*pi),
ch_2(SU fundamental)=-c_2 and Ahat=1-p_1/24+....
"""
from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'anomaly_results.json'
SOURCE = HERE / 'input' / 'smflux_v0314.py'
STACKS = ('c', 'L', 'a', 'b', 'd')
N = (3, 2, 1, 1, 1)
m = (0, -3, 3, 3, 0)
q = s.Matrix(s.symbols('q_c q_L q_a q_b q_d'))
c, ell, a, b, d = q
Y = s.Matrix([s.Rational(1, 6), 0, -s.Rational(1, 2), s.Rational(1, 2), -s.Rational(1, 2)])
BL = s.Matrix([s.Rational(1, 3), 0, 0, 0, -1])
phase = s.ones(5, 1)
K = s.Matrix([[0, -2, 1, 1, 0], [3, -4, 0, 0, 1]])
checks = []


def check(label, condition):
    assert bool(condition), label
    checks.append(label)


def equal(label, lhs, rhs):
    if isinstance(lhs, s.MatrixBase) or isinstance(rhs, s.MatrixBase):
        diff = (lhs - rhs).applyfunc(s.simplify)
        check(label, diff == s.zeros(*diff.shape))
    else:
        check(label, s.simplify(lhs-rhs) == 0)


def matrix(A):
    return [[str(x) for x in row] for row in A.tolist()]


def expression(p):
    return str(s.expand(p))


def directional(p, v):
    return s.expand(sum(v[i]*s.diff(p, q[i]) for i in range(5)))


source_record = {
    'source': 'mtft-0.31.4/src/mtft/surface/smflux.py:M1_STACKS',
    'package_sha256': '2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5',
    'copied_source_present': SOURCE.exists(),
    'fallback': 'M1 constants are recorded literally in this script; no network or package installation needed.',
}
module = None
if SOURCE.exists():
    source_record['source_sha256'] = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    spec = importlib.util.spec_from_file_location('axg_smflux', SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source_record['M1_STACKS_checked'] = True
    check('source stack dimensions', tuple(module.M1_STACKS['N'][k] for k in STACKS) == N)
    check('source flux degrees', tuple(module.M1_STACKS['m'][k] for k in STACKS) == m)
    equal('source hypercharge', s.Matrix([s.Rational(module.M1_STACKS['Y'][k]) for k in STACKS]), Y)
    equal('source B-L', s.Matrix([s.Rational(module.M1_STACKS['B-L'][k]) for k in STACKS]), BL)

# Oriented bifundamental index n_ij=m_i-m_j accounts for charge conjugation.
P = s.expand(sum((m[i]-m[j])*N[i]*N[j]*(q[i]-q[j])**3
                 for i in range(5) for j in range(i+1, 5)))
Agr = s.expand(sum((m[i]-m[j])*N[i]*N[j]*(q[i]-q[j])
                   for i in range(5) for j in range(i+1, 5)))
A3 = s.expand(sum(s.Rational(1, 2)*N[j]*(m[0]-m[j])*(q[0]-q[j])
                  for j in range(1, 5)))
A2 = s.expand(sum(s.Rational(1, 2)*N[j]*(m[1]-m[j])*(q[1]-q[j])
                  for j in range(5) if j != 1))
mixed = s.Matrix([[s.diff(p, z) for z in q] for p in (A3, A2, Agr)])

x, y, z, w = c-ell, a-ell, b-ell, d-ell
r1, r2 = K*q
Q1 = 24*(y*y-y*z+z*z)+27*x*x+9*w*w
Q2 = -9*(y*y+z*z)
Q = [s.expand(Q1), s.expand(Q2)]
equal('cubic exact ideal factorization', P, r1*Q1+r2*Q2)
equal('mixed SU3 factorization', A3, s.Rational(3, 2)*r1)
equal('mixed SU2 factorization', A2, 3*r1+s.Rational(3, 2)*r2)
equal('mixed gravity factorization', Agr, 24*r1)
check('mixed anomaly rank two', mixed.rank() == 2)
equal('massless named generators lie in K kernel', K*s.Matrix.hstack(phase, Y, BL), s.zeros(2, 3))
check('named kernel generators independent', s.Matrix.hstack(phase, Y, BL).rank() == 3)
alpha, beta, gamma = s.symbols('alpha beta gamma')
kernel_sub = dict(zip(q, alpha*phase+beta*Y+gamma*BL))
equal('entire restricted cubic polynomial vanishes', P.subs(kernel_sub, simultaneous=True), 0)
equal('common phase decouples from full cubic', directional(P, phase), 0)

# A(u,v,w)=tr_chiral(Q_u Q_v Q_w), not the polynomial coefficient with multinomial factors.
AYqq = directional(P, Y)/3
ABLqq = directional(P, BL)/3
AYYq = directional(AYqq, Y)/2
AYBLq = directional(AYqq, BL)/2
ABLBLq = directional(ABLqq, BL)/2
equal('Yqq anomaly', AYqq, r1*(12*(z-y)+3*(x-w))+3*r2*(y-z))
equal('BLqq anomaly', ABLqq, 6*r1*(x-w))
equal('YYq anomaly', AYYq, 7*r1-s.Rational(3, 2)*r2)
equal('YBLq anomaly', AYBLq, 2*r1)
equal('BLBLq anomaly', ABLBLq, 4*r1)
check('massless Y has mixed anomalies with massive directions', AYYq != 0)
check('massless BL has mixed anomalies with massive directions', ABLBLq != 0)

tensor = []
for i, j, k in itertools.combinations_with_replacement(range(5), 3):
    val = s.diff(P, q[i], q[j], q[k])/6
    if val:
        tensor.append({'indices': [STACKS[i], STACKS[j], STACKS[k]], 'value': str(val)})

# Verify the local Chern-Simons descent change explicitly with the graded Leibniz rule.
# P stands for the unnormalized Abelian anomaly six-form. f_i are commuting 2-forms.
# C4=(1/3) sum_I (K A)_I wedge [sum_j A_j partial_j Q_I(f)].
# d(A_i wedge A_j)=f_i wedge A_j-A_i wedge f_j.
C4 = {}
for i, j in itertools.combinations(range(5), 2):
    coeff = s.expand(sum((K[I, i]*s.diff(Q[I], q[j])-K[I, j]*s.diff(Q[I], q[i]))/3 for I in range(2)))
    if coeff:
        C4[(i, j)] = coeff
dC4 = s.zeros(5, 1)
for (i, j), coeff in C4.items():
    dC4[j] += q[i]*coeff
    dC4[i] -= q[j]*coeff
I5sym = s.Matrix([s.diff(P, qi)/3 for qi in q])
I5fac = K.T*s.Matrix(Q)
equal('graded descent difference equals dC4', I5sym-I5fac, dC4)
equal('Euler verifies dI5 symmetric', (q.T*I5sym)[0], P)
equal('Euler verifies dI5 factorized', (q.T*I5fac)[0], P)
equal('factorized descent preserves all massless gauge parameters',
      s.Matrix.hstack(phase, Y, BL).T*I5fac, s.zeros(3, 1))

c2_3, c2_2, p1 = s.symbols('c2_3 c2_2 p1')
I6 = P/6-2*A3*c2_3-2*A2*c2_2-Agr*p1/24
X41 = Q1/6-3*c2_3-6*c2_2-p1
X42 = Q2/6-3*c2_2
equal('full normalized anomaly six-form factorization', I6, r1*X41+r2*X42)

# Limited period check: on a CLOSED SPIN four-manifold, any integral degree-two
# class has even square. Integral mixed coefficients + half-integral diagonal
# coefficients suffice here. Assumes a split SU3 x SU2 x U1^5 product cover.
period_checks = []
for I, Qi in enumerate(Q):
    for exps, coeff in s.Poly(Qi/6, *q).terms():
        is_square = max(exps) == 2
        integral = (2*coeff if is_square else coeff).is_Integer
        check(f'X4_{I+1} spin-period coefficient {exps}', integral)
        period_checks.append({'axion': I+1, 'monomial_exponents': list(exps),
                              'coefficient': str(coeff), 'uses_even_square': is_square})

SU3_cubic = sum(N[j]*(m[0]-m[j]) for j in range(1, 5))
SU2_doublets = sum(N[j]*abs(m[1]-m[j]) for j in range(5) if j != 1)
equal('SU3 cubic gauge anomaly zero', SU3_cubic, 0)
check('SU2 Witten parity even in stated bifundamental spectrum', SU2_doublets % 2 == 0)

snf = smith_normal_form(K, domain=ZZ)
check('K Smith invariant factors one one', [abs(snf[i, i]) for i in range(2)] == [1, 1])
minors = {','.join(STACKS[i] for i in inds): int(K[:, inds].det())
          for inds in itertools.combinations(range(5), 2)}
check('gcd maximal K minors one', math.gcd(*minors.values()) == 1)
kernel_lattice = s.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 2, -1], [-3, 4, 0]])
equal('primitive integer kernel basis', K*kernel_lattice, s.zeros(2, 3))
check('primitive kernel free parameters have unit coordinate minor', kernel_lattice[:3, :].det() == 1)
Kdet = s.Matrix([[0, -1, 1, 1, 0], [1, -2, 0, 0, 1]])
equal('K descends to determinant characters at charge level', Kdet*s.diag(*N), K)

KKt = K*K.T
identity_eigs = [16-2*s.sqrt(41), 16+2*s.sqrt(41)]
for eigen in identity_eigs:
    equal('identity metric characteristic eigenvalue '+str(eigen), (KKt-eigen*s.eye(2)).det(), 0)
check('identity metric positive eigenvalues', all(v.is_positive for v in identity_eigs))
H0 = 2*s.diag(*N)
stack_gram = K*H0.inv()*K.T
equal('stack metric reduced mass matrix', stack_gram, s.Matrix([[2, 2], [2, 6]]))
stack_eigs = [4-2*s.sqrt(2), 4+2*s.sqrt(2)]
for eigen in stack_eigs:
    equal('stack metric characteristic eigenvalue '+str(eigen), (stack_gram-eigen*s.eye(2)).det(), 0)
check('stack metric positive eigenvalues', all(v.is_positive for v in stack_eigs))
check('mass matrix rank two', (K.T*K).rank() == 2)

if module:
    model = {'degrees': dict(zip(STACKS, m)), 'y': {name: module.Fr(str(Y[i])) for i, name in enumerate(STACKS)}}
    ledger = module.anomaly_ledger(model)
    for i, name in enumerate(STACKS):
        for j, label in enumerate(('SU3^2', 'SU2^2', 'grav')):
            equal(f'package mixed ledger {name} {label}', s.Rational(ledger['ledger'][name][label]), mixed[j, i])
    # Function mixed_nonabelian_anomalies uses twice the T(fund)=1/2 ledger convention.
    for t in range(12):
        point = [s.Rational(((i+2)*(t+3)) % 11-5, (t % 3)+1) for i in range(5)]
        qdict = {name: module.Fr(str(point[i])) for i, name in enumerate(STACKS)}
        nc = module.M1_STACKS['N']; mc = module.M1_STACKS['m']
        equal(f'package cubic value {t}', P.subs(dict(zip(q, point))), s.Rational(module.abelian_anomaly_polynomial(nc, mc, qdict)))
        mix = module.mixed_nonabelian_anomalies(nc, mc, qdict)
        for label, p in (('c', A3), ('L', A2)):
            equal(f'package doubled mixed convention {t} {label}', 2*p.subs(dict(zip(q, point))), s.Rational(mix[label]))

report = {
    'test_id': 'AXG-01-ANOMALY', 'status': 'EXACT algebra, conditional 4D EFT interpretation',
    'source': source_record,
    'normalization': {
        'P': 'tr_left_chiral Q(q)^3; fundamental non-Abelian Dynkin index 1/2',
        'I6': 'P(f)/6 - 2 A3(f) c2_3 - 2 A2(f) c2_2 - Agr(f) p1/24',
        'package_note': 'mixed_nonabelian_anomalies returns twice the SU(N)^2 coefficients in anomaly_ledger.',
        'charges': 'stack fundamental charge +1; stack order c,L,a,b,d',
    },
    'input': {'N': list(N), 'm': list(m), 'Y': matrix(Y), 'B-L': matrix(BL), 'phase': matrix(phase), 'K': matrix(K)},
    'cubic': {'expanded': expression(P), 'symmetric_tensor_nonzero_sorted_indices': tensor,
              'difference_coordinates': {'x': 'q_c-q_L', 'y': 'q_a-q_L', 'z': 'q_b-q_L', 'w': 'q_d-q_L'},
              'r1': 'y+z', 'r2': '3*x+w',
              'Q1': '24*(y^2-y*z+z^2)+27*x^2+9*w^2', 'Q2': '-9*(y^2+z^2)',
              'Q1_expanded': expression(Q1), 'Q2_expanded': expression(Q2),
              'Q_symmetric_matrices': [matrix(s.hessian(Qi, q)/2) for Qi in Q],
              'factorization': 'P=r1*Q1+r2*Q2',
              'restriction_to_phase_Y_BL': 0},
    'mixed': {'matrix_rows_SU3_SU2_gravity': matrix(mixed),
              'SU3': '3*r1/2', 'SU2': '3*r1+3*r2/2', 'gravity': '24*r1',
              'SU3_cubic': SU3_cubic, 'SU2_doublets': SU2_doublets, 'SU2_global_parity': SU2_doublets % 2},
    'massless_mixed_cubic': {'A(Y,q,q)': expression(AYqq), 'A(B-L,q,q)': expression(ABLqq),
                            'A(Y,Y,q)': '7*r1-3*r2/2', 'A(Y,B-L,q)': '2*r1', 'A(B-L,B-L,q)': '4*r1'},
    'descent': {
        'necessity': 'Symmetric consistent descent has nonzero massless-current anomalies with massive backgrounds; generalized Chern-Simons/Bardeen counterterms or an equivalent preserving scheme are required.',
        'local_identity': 'I5_symmetric-I5_factorized=dC4, with C4=(1/3)sum_I(KA)_I wedge sum_a A_a partial_a Q_I(f), for unnormalized P.',
        'normalized_abelian_counterterm': 'Divide the displayed C4 by 6 for the P/6 anomaly convention; overall action sign tracks anomaly convention.',
        'counterterm_coefficients_Ai_wedge_Aj': [{'i': STACKS[i], 'j': STACKS[j], 'coefficient_in_f': expression(value)} for (i, j), value in C4.items()],
        'X4_1': 'Q1/6-3*c2_3-6*c2_2-p1', 'X4_2': 'Q2/6-3*c2_2',
        'factorized_full_I6': 'r1*X4_1+r2*X4_2',
        'axion_variation': 'With D a=d a+K A, delta a=-K lambda. A Wess-Zumino term with appropriate sign cancels the factorized local anomaly.',
    },
    'limited_period_check': {
        'assumptions': 'closed spin 4-manifold; split SU3 x SU2 x U1^5 cover bundles; f_a integral; axions of period 2*pi and conventional exponentiated WZ normalization',
        'result': 'X4_1 and X4_2 have integral periods under these assumptions: every diagonal half-integral term multiplies an even self-intersection; all mixed coefficients and c2,p1 coefficients are integral.',
        'coefficients': period_checks,
        'limitation': 'Does not establish differential-cohomological/global consistency on nontrivial U(N) quotient bundles, the full parent anomaly structure, or a UV completion.',
    },
    'integer_lattice': {'SNF': matrix(snf), 'maximal_minors': minors,
                        'primitive_kernel_basis_columns': matrix(kernel_lattice),
                        'determinant_character_matrix': matrix(Kdet),
                        'determinant_note': 'K rows are characters det(U2)^(-1)*U1a*U1b and det(U3)*det(U2)^(-2)*U1d; this checks axion shifts descend at character level only.',
                        'remnant': 'For the standard primitive U1^5 lattice, kernel torus is connected U1^3; no additional finite Stueckelberg remnant.',
                        'limitation': 'Full gauge group global quotient, non-Abelian centers, and admissible bundles need independent treatment.'},
    'mass_spectrum': {'identity_metrics_KKt': matrix(KKt), 'identity_nonzero_m2_in_unit_scale': [str(v) for v in identity_eigs],
                       'identity_zero_multiplicity': 3, 'stack_H0': matrix(H0), 'stack_reduced_matrix': matrix(stack_gram),
                       'stack_nonzero_m2': ['g4^2*f_axion^2*('+str(v)+')' for v in stack_eigs],
                       'stack_zero_multiplicity': 3,
                       'assumptions': 'H=H0/g4^2, G_axion=f_axion^2 I2, no other Higgs/Stueckelberg masses or kinetic mixing; illustrative inputs, not a GeV prediction.'},
    'unresolved': ['Origin, number and periodicities of axions', 'Axion and gauge kinetic matrices and their moduli dependence',
                   'Choice of global gauge group and charge lattice', 'Full nonlinear counterterms and anomaly conventions in action',
                   '6D parent local/global anomalies and UV completion', 'Vacuum stabilization and B-L breaking'],
    'reference': {'authors': 'Anastasopoulos, Bianchi, Dudas, Kiritsis', 'title': 'Anomalies, Anomalous U(1)s and generalized Chern-Simons terms', 'url': 'https://arxiv.org/abs/hep-th/0605225'},
    'checks': checks, 'check_count': len(checks), 'all_passed': True,
}
OUTPUT.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({'output': str(OUTPUT), 'check_count': len(checks), 'all_passed': True,
                  'source_checked': SOURCE.exists(), 'K_rank': K.rank(), 'SU2_doublets': SU2_doublets}, indent=2))
