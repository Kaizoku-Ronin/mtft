#!/usr/bin/env python3
"""R2C-01: exact, finite 6D chirality audit of the frozen M1 sector list.

Requires Python >=3.10 and SymPy. No MTFT installation or floating-point
eigensolver is used. Run: python chirality_audit.py --output results
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import csv
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import platform
import runpy
import sys

import sympy as s

ROOT = Path(__file__).resolve().parent
names = ('c', 'L', 'a', 'b', 'd')
pairs = tuple(combinations(range(5), 2))
labels = tuple(names[i] + names[j] for i, j in pairs)
pair_by_label = dict(zip(labels, pairs))
positions = {name: k for k, name in enumerate(labels)}
N = (3, 2, 1, 1, 1)
m0 = (0, -3, 3, 3, 0)
Y = s.Matrix([s.Rational(1, 6), 0, -s.Rational(1, 2), s.Rational(1, 2), -s.Rational(1, 2)])
BL = s.Matrix([s.Rational(1, 3), 0, 0, 0, -1])
degrees = tuple(m0[i] - m0[j] for i, j in pairs)
dims = tuple(N[i] * N[j] for i, j in pairs)
families = ('cL', 'ca', 'cb', 'Ld', 'ad', 'bd')
targets = dict(zip(families, (3, -3, -3, -3, 3, 3)))
vertices = {'up': ('cL', 'ca'), 'down': ('cL', 'cb'),
            'neutrino': ('Ld', 'ad'), 'charged_lepton': ('Ld', 'bd')}
f = s.symbols('f_c f_L f_a f_b f_d')
C, D3, W, p1, p2 = s.symbols('C D3 W p1 p2')
alpha, beta, gamma = s.symbols('alpha beta gamma')
AHAT8 = (7*p1**2 - 4*p2)/s.Integer(5760)
checks: list[str] = []


def require(label, condition):
    # Deliberately not Python assert: checks remain active under python -O.
    if not bool(condition):
        raise AssertionError(label)
    checks.append(label)


def equal(label, a, b):
    if isinstance(a, s.MatrixBase) or isinstance(b, s.MatrixBase):
        delta = s.Matrix(a) - s.Matrix(b)
        require(label, all(s.expand(v) == 0 for v in delta))
    else:
        require(label, s.expand(a-b) == 0)


def matrix_json(M):
    return [[str(v) for v in row] for row in M.tolist()]


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, default=str) + '\n')


def pushforward(poly, flux):
    return s.expand(sum(flux[i]*s.diff(poly, f[i]) for i in range(5)))


def charge(i, j):
    return s.eye(5)[:, i] - s.eye(5)[:, j]


def validate_inputs():
    provenance = json.loads((ROOT/'PROVENANCE.json').read_text())
    for row in provenance['snapshots']:
        observed = hashlib.sha256((ROOT/row['bundle_path']).read_bytes()).hexdigest()
        require('frozen source SHA256: '+row['bundle_path'], observed == row['sha256'])
    parent = runpy.run_path(str(ROOT/'input/parents_v0320.py'))['M1']
    require('source parent stack/rank/degree record',
            parent['stacks'] == names and parent['ranks'] == N and parent['degrees'] == m0)
    require('parent explicitly leaves 6D chirality unassigned',
            parent['six_d_chiralities'] == 'not assigned in the specification')
    tree = ast.parse((ROOT/'input/smflux_v0320.py').read_text())
    stacks = next(ast.literal_eval(node.value) for node in tree.body
                  if isinstance(node, ast.Assign)
                  and any(isinstance(t, ast.Name) and t.id == 'M1_STACKS' for t in node.targets))
    require('independent source stack constants', tuple(stacks['m'][v] for v in names) == m0
            and tuple(stacks['N'][v] for v in names) == N)
    equal('source hypercharge', s.Matrix([s.Rational(stacks['Y'][v]) for v in names]), Y)
    equal('source B-L', s.Matrix([s.Rational(stacks['B-L'][v]) for v in names]), BL)
    return provenance


def clifford_audit():
    sx = s.Matrix([[0, 1], [1, 0]])
    sy = s.Matrix([[0, -s.I], [s.I, 0]])
    sz = s.diag(1, -1)
    I2, I8 = s.eye(2), s.eye(8)
    kron = s.kronecker_product
    gs = [kron(sx,I2,I2), kron(sy,I2,I2), kron(sz,sx,I2),
          kron(sz,sy,I2), kron(sz,sz,sx), kron(sz,sz,sy)]
    for a in range(6):
        for b in range(a, 6):
            equal(f'Clifford anticommutator {a+1},{b+1}', gs[a]*gs[b]+gs[b]*gs[a],
                  2*I8 if a == b else s.zeros(8))
    G = s.I*gs[0]*gs[1]*gs[2]*gs[3]*gs[4]*gs[5]
    CP = gs[1]*gs[3]*gs[5]
    equal('chirality squares to one', G*G, I8)
    projectors = {1:(I8+G)/2, -1:(I8-G)/2}
    for g in gs:
        equal('chirality anticommutes with generator', G*g+g*G, s.zeros(8))
        equal('charge-conjugation intertwiner', CP*g*CP.inv(), -g.T)
    contraction_ranks = []
    for a, b in product((-1, 1), repeat=2):
        scalar = projectors[-a]*projectors[b]
        transposed = projectors[a].T*CP*projectors[b]
        # An internal odd Clifford insertion changes the chirality selection.
        vector = projectors[-a]*gs[4]*projectors[b]
        require(f'Dirac scalar selection {a},{b}', scalar.rank() == (4 if a != b else 0))
        require(f'transpose scalar selection {a},{b}', transposed.rank() == (4 if a != b else 0))
        require(f'odd Clifford selection {a},{b}', vector.rank() == (4 if a == b else 0))
        contraction_ranks.append({'chi1':a, 'chi2':b, 'Dirac_scalar_rank':scalar.rank(),
                                  'transpose_scalar_rank':transposed.rank(), 'vector_rank':vector.rank()})
    return {'construction':'Complexified Euclidean Cl(6) tests algebraic Spin(1,5) invariant tensors; not a Wick-rotated dynamical model.',
            'gamma_matrices':[matrix_json(g) for g in gs], 'chirality_matrix':matrix_json(G),
            'charge_conjugation_matrix':matrix_json(CP), 'ranks':contraction_ranks}


def geometry_and_charges():
    # Arithmetic input H0(S0)=span(1,u), u=+/-i/sqrt(13) at the chosen CM points.
    v = s.I/s.sqrt(13)
    E = s.Matrix([[1,v], [1,v], [1,-v]])
    require('CM evaluation matrix has rank two', E.rank() == 2)
    equal('opposite-sign evaluation minor', E.extract([0,2],[0,1]).det(), -2*v)
    h = {3:(3,0), -3:(0,3), 6:(6,0), -6:(0,6), 0:(2,2)}
    ledger = []
    for name, (i,j), d, dimension in zip(labels, pairs, degrees, dims):
        h0, h1 = h[d]
        equal('Riemann-Roch spin index '+name, h0-h1, d)
        ledger.append({'sector':name, 'pair':[names[i],names[j]], 'dimension':dimension,
                       'degree':d, 'h0':h0, 'h1':h1, 'hypercharge':Y[i]-Y[j],
                       'B_minus_L':BL[i]-BL[j], 'charge_vector':list(charge(i,j)),
                       'family_target_index':targets.get(name)})
    Q, uc, dc = charge(0,1), charge(2,0), charge(3,0)
    L, ec, nc = charge(4,1), charge(3,4), charge(2,4)
    Hu, Hd = charge(1,2), charge(1,3)
    for name, expression in [('up',Q+uc+Hu),('down',Q+dc+Hd),
                             ('charged_lepton',L+ec+Hd),('neutrino',L+nc+Hu)]:
        equal('4D gauge triangle '+name, expression, s.zeros(5,1))
    for name, expression in [('up',-charge(0,1)+charge(0,2)-Hu),
                             ('down',-charge(0,1)+charge(0,3)-Hd),
                             ('neutrino',-charge(1,4)+charge(2,4)+Hu),
                             ('charged_lepton',-charge(1,4)+charge(3,4)+Hd)]:
        equal('6D vector Dirac triangle '+name, expression, s.zeros(5,1))
    equal('original Higgs degree Hu', (Hu.T*s.Matrix(m0))[0], -6)
    equal('original Higgs degree Hd', (Hd.T*s.Matrix(m0))[0], -6)
    return {'sector_ledger':ledger, 'evaluation_matrix':matrix_json(E),
            'evaluation_rank':E.rank(), 'spin_degree':12, 'genus':13,
            'purity_status':'Conditional exact deduction from frozen arithmetic input; not an independent reconstruction of S0 or the CM points.',
            'zero_degree_warning':'With these explicit untwisted stack bundles, cd and ab each have (h0,h1)=(2,2).',
            'Higgs_H1_dimension':18}


def anomaly_sectors():
    # Chern character components in degree-2 units: ch0, ch1, ..., ch4.
    ch = [[s.Integer(3),0,-C,D3/2,C**2/12],
          [s.Integer(2),0,-W,0,W**2/12],
          [s.Integer(1),0,0,0,0]]
    ch += [ch[2],ch[2]]
    roots = ((alpha,beta,-alpha-beta),(gamma,-gamma),(0,),(0,),(0,))
    cartan_sub = {C:-(alpha**2+alpha*beta+beta**2), D3:-alpha*beta*(alpha+beta), W:-gamma**2}
    sectors = {}
    for label, (i,j) in zip(labels, pairs):
        d = f[i]-f[j]
        coeff = []
        for k in range(5):
            coeff.append(s.expand(sum(ch[i][a]*(-1)**b*ch[j][b]*d**e/s.factorial(e)
                         for a in range(k+1) for b in range(k-a+1) for e in [k-a-b])))
        I8 = s.expand(coeff[4]-p1*coeff[2]/24+coeff[0]*AHAT8)
        I6 = s.expand(coeff[3]-p1*coeff[1]/24)
        weights = [d+a-b for a in roots[i] for b in roots[j]]
        direct8 = sum(t**4/24-p1*t**2/48+AHAT8 for t in weights)
        direct6 = sum(t**3/6-p1*t/24 for t in weights)
        equal('independent Cartan I8: '+label, I8.subs(cartan_sub), direct8)
        equal('independent Cartan I6: '+label, I6.subs(cartan_sub), direct6)
        # Stronger than one numerical flux: compare each formal flux derivative.
        for k in range(5):
            equal(f'index pushforward per stack: {label},{names[k]}',
                  s.diff(I8,f[k]), ((1 if k == i else 0)-(1 if k == j else 0))*I6)
        equal('irreducible gravitational coefficient '+label, s.diff(I8,p2), -s.Rational(N[i]*N[j],1440))
        sectors[label] = {'I8':I8, 'I6':I6, 'dimension':N[i]*N[j],
                          'c3_coefficient':s.diff(I8,D3)}
    return sectors


def enumerate_assignments(sectors, out):
    rows = []
    baseline8 = sum(sectors[label]['I8'] for label in labels)
    baseline6 = pushforward(baseline8,m0)
    aY, aBL, common = s.symbols('aY aBL common')
    safe_sub = dict(zip(f, Y*aY+BL*aBL+s.ones(5,1)*common))
    safe6 = {k:s.expand(v['I6'].subs(safe_sub)) for k,v in sectors.items()}
    for idx, signs in enumerate(product((-1,1), repeat=10)):
        signed = dict(zip(labels,signs))
        indices = tuple(a*b for a,b in zip(signs,degrees))
        family_ok = all(indices[positions[k]] == target for k,target in targets.items())
        full_ok = indices == degrees
        allowed = {name:signed[a] != signed[b] for name,(a,b) in vertices.items()}
        signed_dim = sum(a*b for a,b in zip(signs,dims))
        color = s.expand(sum(signed[label]*sectors[label]['c3_coefficient'] for label in labels))
        invariant = s.expand(sum(indices[positions[k]]*safe6[k] for k in labels)) == 0
        row = {'assignment_id':idx, **{'chi_'+k:signed[k] for k in labels},
               **{'index_'+k:indices[positions[k]] for k in labels},
               'preserves_six_family_indices':family_ok, 'preserves_all_M1_indices':full_ok,
               **{'scalar_'+k:v for k,v in allowed.items()},
               'all_four_scalar_bilinears':all(allowed.values()),
               'signed_6D_dimension':signed_dim, 'I8_p2_coefficient':str(-s.Rational(signed_dim,1440)),
               'I8_c3_coefficient':str(color), '4D_Y_BL_phase_anomaly_free':invariant}
        rows.append(row)
    with (out/'assignments.csv').open('w',newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    survivors = [r for r in rows if r['preserves_six_family_indices']]
    full = [r for r in rows if r['preserves_all_M1_indices']]
    scalar = [r for r in rows if r['all_four_scalar_bilinears']]
    both = [r for r in survivors if r['all_four_scalar_bilinears']]
    # Combinatorial expectations follow from six fixed signs, eight fixed signs,
    # or two connected three-vertex chirality graphs plus four isolated signs.
    require('scan exhaustive unique sign vectors', len(rows) == 2**10 and
            len({tuple(r['chi_'+k] for k in labels) for r in rows}) == 2**10)
    require('six fixed signs leave four free signs', len(survivors) == 2**4)
    require('eight fixed signs leave two free signs', len(full) == 2**2)
    require('two Yukawa components and four isolated sectors', len(scalar) == 2**6)
    require('no fixed-flux family/scalar intersection', len(both) == 0)
    require('each family-preserving assignment forbids all four scalars',
            all(not r['scalar_'+k] for r in survivors for k in vertices))
    require('all family-preserving signed dimensions strictly positive',
            min(r['signed_6D_dimension'] for r in survivors) > 0)
    require('all assignments have a nonzero color cubic coefficient',
            all(r['I8_c3_coefficient'] != '0' for r in rows))
    full_polys = []
    for row in full:
        poly = s.expand(sum(row['chi_'+label]*sectors[label]['I8'] for label in labels))
        equal('full-ledger survivor has unchanged 4D anomaly '+str(row['assignment_id']),
              pushforward(poly,m0),baseline6)
        full_polys.append({'assignment_id':row['assignment_id'], 'chi_cd':row['chi_cd'],
                           'chi_ab':row['chi_ab'], 'I8':str(poly), 'I6':str(baseline6),
                           'signed_dimension':row['signed_6D_dimension'],
                           'p2_coefficient':row['I8_p2_coefficient'], 'c3_coefficient':row['I8_c3_coefficient']})
    return {'counts':{'all':len(rows), 'six_family_indices':len(survivors),
                      'all_M1_indices':len(full), 'all_four_scalar_bilinears':len(scalar),
                      'families_and_all_four_scalar_bilinears':len(both),
                      'family_preserving_and_Y_BL_phase_anomaly_free':sum(r['4D_Y_BL_phase_anomaly_free'] for r in survivors)},
            'family_preserving_assignments':survivors,
            'family_signed_dimension_histogram':dict(sorted(Counter(r['signed_6D_dimension'] for r in survivors).items())),
            'full_ledger_survivors':full_polys, 'baseline_I6':str(baseline6),
            'all_assignment_scalar_count_histogram':dict(sorted(Counter(sum(r['scalar_'+k] for k in vertices) for r in rows).items())),
            'quark_only_witness':'Fixed Q and up-antiquark indices already force chi_cL=chi_ca=+1; a scalar needs them opposite.'}


def relaxed_flux_controls(sectors):
    unknown = s.symbols('m_c m_L m_a m_b m_d')
    controls = []
    cochar = s.Matrix([0,0,1,-1,0])
    k1 = s.Matrix([[0,-2,1,1,0]])
    t = s.symbols('t')
    for aq, al in product((-1,1),repeat=2):
        core = {'cL':aq,'ca':-aq,'cb':-aq,'Ld':al,'ad':-al,'bd':-al}
        equations = [unknown[0]]
        for label,target in targets.items():
            i,j = pair_by_label[label]
            equations.append(core[label]*(unknown[i]-unknown[j])-target)
        solutions = s.solve(equations,unknown,dict=True)
        require(f'flux equations have unique solution up to common shift: {aq},{al}',len(solutions)==1)
        flux = tuple(int(solutions[0][v]) for v in unknown)
        expected = (0,-3*aq,-3*aq,-3*aq,3*(al-aq))
        require(f'closed-form relaxed flux: {aq},{al}',flux==expected)
        require(f'flux control changes M1: {aq},{al}',flux!=m0)
        I8 = s.expand(sum(core[k]*sectors[k]['I8'] for k in families))
        equal(f'core pure gravitational anomaly cancels: {aq},{al}',s.diff(I8,p2),0)
        equal(f'core cubic color direction: {aq},{al}',s.diff(I8,D3),aq*(k1*s.Matrix(f))[0]/2)
        equal(f'color direction annihilates primitive cocharacter: {aq},{al}',(k1*cochar)[0],0)
        equal(f'color direction annihilates new flux: {aq},{al}',(k1*s.Matrix(flux))[0],0)
        absub = {C:0,D3:0,W:0,p1:0,p2:0,**dict(zip(f,cochar*t))}
        quartic = s.expand(I8.subs(absub)).coeff(t,4)
        equal(f'primitive quartic from direct charged-field sum: {aq},{al}',
              quartic, sum(core[k]*dims[positions[k]]*((cochar[pair_by_label[k][0]]-cochar[pair_by_label[k][1]])**4)
                           for k in families)/s.Integer(24))
        equal(f'closed-form primitive quartic: {aq},{al}',quartic,-s.Rational(3*aq+al,12))
        require(f'ordinary tensor lattice norm fails: {aq},{al}',not (8*quartic).is_Integer)
        n2 = sum(core[k]*(N[j] if i==1 else N[i]) for k in families
                 for i,j in [pair_by_label[k]] if 1 in (i,j))
        # Do not silently delete index-zero sectors: enumerate all sixteen
        # completions of this altered-flux core back to the ten-sector field list.
        other = ('cd','La','Lb','ab')
        completions = []
        for extra in product((-1,1),repeat=4):
            signs = {**core,**dict(zip(other,extra))}
            fullpoly = I8 + sum(signs[k]*sectors[k]['I8'] for k in other)
            totaldim = sum(signs[k]*dims[positions[k]] for k in labels)
            quartic_full = s.expand(fullpoly.subs(absub)).coeff(t,4)
            require(f'pure-gravity completion classification: {aq},{al},{extra}',
                    (totaldim==0) == (extra[0]==extra[3] and extra[1]==extra[2]==-extra[0]))
            if totaldim == 0:
                equal(f'closed-form full-completion norm: {aq},{al},{extra}',
                      8*quartic_full, -2*aq-s.Rational(2*al,3)+4*extra[0])
                require(f'pure-gravity completion fails this tensor norm: {aq},{al},{extra}',
                        not (8*quartic_full).is_Integer)
            completions.append({'extra_chiralities':dict(zip(other,extra)), 'signed_dimension':totaldim,
                                'a_minus_b_required_norm':str(8*quartic_full),
                                'passes_this_norm_test':bool((8*quartic_full).is_Integer)})
        require(f'two pure-gravity completions per changed flux: {aq},{al}',
                sum(r['signed_dimension']==0 for r in completions)==2)
        controls.append({'chi_Q':aq,'chi_L':al,'core_chiralities':core,'flux':flux,
                         'nonfamily_degrees':{k:flux[pair_by_label[k][0]]-flux[pair_by_label[k][1]] for k in other},
                         'Higgs_degrees':[flux[1]-flux[2],flux[1]-flux[3]],
                         'six_sector_I8':str(I8), 'net_complex_SU2_doublets':n2, 'SU2_mod6':n2%6,
                         'primitive_cocharacter':list(cochar),'primitive_quartic':str(quartic),
                         'required_tensor_norm':str(8*quartic),'ten_sector_completions':completions,
                         'ten_sector_gravity_and_this_norm_pass':sum(r['signed_dimension']==0 and r['passes_this_norm_test'] for r in completions)})
    # A genuinely independent check of the published C3X *local polynomial*:
    # start with the six-field charge sum above, then restrict and multiply.
    plus = next(c for c in controls if c['chi_Q']==1 and c['chi_L']==1)
    core8 = sum(plus['core_chiralities'][k]*sectors[k]['I8'] for k in families)
    h,x = s.symbols('h x')
    mu = s.Matrix([0,-1,-1,-1,0])
    embedded = s.expand(3*core8.subs(dict(zip(f,6*Y*h+mu*x))))
    factorized = (W+9*h**2)*(3*C+W+p1/2-27*h**2-6*x**2)
    equal('C3X factorization independently from spectrum',embedded,factorized)
    return {'branches':controls,'C3X_independent_I8':str(embedded),
            'C3X_factors':[str(W+9*h**2),str(3*C+W+p1/2-27*h**2-6*x**2)],
            'C3X_scope':'Local anomaly polynomial only. Three input copies times one mode; not a one-parent three-mode solution.',
            'norm_scope':'Ordinary integral tensor-charge lattice with I8=Omega(X4,X4)/2 and X4=(b_tt/2)t^2 along the primitive cocharacter; requires b_tt^2=8[t^4]I8 integral. Additional matter or other anomaly-cancellation structures change the problem.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results')
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    provenance = validate_inputs()
    clifford = clifford_audit()
    geometry = geometry_and_charges()
    sectors = anomaly_sectors()
    enumeration = enumerate_assignments(sectors,out)
    controls = relaxed_flux_controls(sectors)
    write_json(out/'clifford_matrices.json',clifford)
    write_json(out/'sector_anomaly_basis.json',{'stack_order':names,'sector_order':labels,
              'symbols':{'C':'c2(SU3)','D3':'c3(SU3)','W':'c2(SU2)','p1':'p1(T)','p2':'p2(T)'},
              'normalization':'positive complex Weyl: [Ahat(T) ch(R)]8',
              'sectors':sectors})
    result = {'study':'R2C-01','status':'EXACT within the explicitly fixed finite field/operator class',
              'runtime':{'python':platform.python_version(),'sympy':s.__version__},
              'baseline_archive_sha256':provenance['archive_sha256'],
              'sign_convention':'nL-nR in oriented R_ij = epsilon_ij*(m_i-m_j); positive internal chirality is H0.',
              'sector_order':labels,'geometry':geometry,'enumeration':enumeration,
              'relaxed_flux_controls':controls,
              'checks':checks,'summary':{'passed':len(checks),'failed':0}}
    write_json(out/'audit_results.json',result)
    write_json(out/'verification.json',{'passed':len(checks),'failed':0,'checks':checks,
              'independence':'Chern characters versus direct Cartan weights, exact Clifford matrices, representation-by-representation flux pushforward; scan counts also derived by free-sign counting.'})
    print(json.dumps({'counts':enumeration['counts'],
                     'family_signed_dimensions':enumeration['family_signed_dimension_histogram'],
                     'full_survivors':[{'cd':r['chi_cd'],'ab':r['chi_ab'],'signed_dimension':r['signed_dimension']}
                                       for r in enumeration['full_ledger_survivors']],
                     'relaxed_flux':[{'chi_Q':c['chi_Q'],'chi_L':c['chi_L'],'flux':c['flux'],
                                      'required_norm':c['required_tensor_norm'],
                                      'ten_sector_gravity_and_norm_pass':c['ten_sector_gravity_and_this_norm_pass']}
                                     for c in controls['branches']],
                     'checks_passed':len(checks),'output':str(out)},indent=2))


if __name__ == '__main__':
    main()
