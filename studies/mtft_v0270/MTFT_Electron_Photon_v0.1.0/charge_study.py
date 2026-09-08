"""Exact sector algebra and conditional electromagnetic charge constraints."""
import numpy as np
import sympy as sp
from scipy.linalg import expm
from common import ROOT, BLOCKS, COMPLEX_DIMS, dump, read, comm, herm, internal_stage, gamma_matrices, slash


def exact_certificate():
    data = read('inputs/x0143_certified.json')['arrays']
    T = sp.Matrix(data['T2'])
    cert = read('inputs/sector_certificate.json')
    Ps = [sp.Matrix(cert['projectors'][b]) for b in BLOCKS]
    checks = {
        'completeness': sum(Ps, sp.zeros(26)) == sp.eye(26),
        'orthogonal_idempotents': all(P*Q == (P if i == j else sp.zeros(26))
                                      for i,P in enumerate(Ps) for j,Q in enumerate(Ps)),
        'commute_with_T2': all(P*T == T*P for P in Ps),
        'ranks': [int(P.trace()) for P in Ps] == [2,4,8,12],
    }
    # Incidence constraints are exact for the stated nonzero mixing pattern.
    # P_i [Q,M] P_j = (q_i-q_j) P_i M P_j.
    edge_sets = {'unmixed': [], 'old_quartic': [(1,2)],
                 'fully_mixed': [(i,j) for i in range(4) for j in range(i+1,4)]}
    constraints = {}
    for label, edges in edge_sets.items():
        C = sp.zeros(len(edges), 4)
        for k,(i,j) in enumerate(edges):
            C[k,i],C[k,j] = 1,-1
        constraints[label] = {'edges': [[BLOCKS[i],BLOCKS[j]] for i,j in edges],
                              'incidence_matrix': [[int(v) for v in row] for row in C.tolist()],
                              'charge_space_dimension': 4-C.rank(),
                              'nullspace_basis': [[int(v) for v in x] for x in C.nullspace()]}
    assert all(checks.values())
    return {'status': 'EXACT finite rational algebra, conditional on the stated mixing pattern',
            'checks': checks, 'constraints': constraints,
            'scope': 'Four scalar sector charges only; not the full commutant or a charge quantization theorem.'}


def mass_charge_dimension(M, Ps):
    C = np.column_stack([comm(P,M).ravel() for P in Ps])
    s = np.linalg.svd(np.vstack([C.real, C.imag]), compute_uv=False)
    return int(4 - np.count_nonzero(s > 1e-9)), s.tolist()


def main():
    exact = exact_certificate()
    dump('charge_exact_certificate.json', exact)
    stage = internal_stage()
    Ps = np.array([herm(P) for P in stage['P']])
    T = herm(stage['T2'])
    I = np.eye(13)
    M0 = I + T / (4*np.linalg.norm(T, 2))
    W11,W13 = herm(stage['W11']),herm(stage['W13'])
    residuals = dict(stage['residuals'])
    residuals.update({'charge_commutator_norm_formula': 0., 'equal_charge_commutator': 0.,
                      'balanced_diagonal_blocks': 0., 'AL_commutator': 0., 'Ward_identity_mass_defect': 0.})
    scans, models = [], []
    saved = {}
    min_mass = 100.
    for seed in range(2026090700,2026090716):
        rng = np.random.default_rng(seed)
        X = herm(rng.normal(size=(13,13)) + 1j*rng.normal(size=(13,13)))
        W = W11@W13
        Xal = (X+W11@X@W11+W13@X@W13+W@X@W.conj().T)/4
        V = herm(Ps[1]@Xal@Ps[2] + Ps[2]@Xal@Ps[1])
        V /= np.linalg.norm(V)
        residuals['balanced_diagonal_blocks'] = max(residuals['balanced_diagonal_blocks'],
            float(max(np.linalg.norm(P@V@P) for P in Ps)))
        residuals['AL_commutator'] = max(residuals['AL_commutator'],
            float(np.linalg.norm(comm(V,W11))),float(np.linalg.norm(comm(V,W13))))
        for epsilon in [0., .05, .2, .4]:
            M = M0+epsilon*V
            min_mass = min(min_mass, float(np.linalg.eigvalsh(M)[0]))
            for delta in [0.,.5,1.,2.]:
                charges = np.array([0.,-1.,-1.+delta,2.])
                Q = np.einsum('b,bij->ij',charges,Ps)
                actual = float(np.linalg.norm(comm(Q,M)))
                expected = epsilon*abs(delta)  # ||V||_F = 1, only old↔q4 entries.
                residuals['charge_commutator_norm_formula'] = max(
                    residuals['charge_commutator_norm_formula'],abs(actual-expected))
                if delta == 0:
                    residuals['equal_charge_commutator'] = max(residuals['equal_charge_commutator'], actual)
                scans.append({'seed':seed,'epsilon':epsilon,'charge_difference':delta,
                              'commutator_norm':actual,'exact_pattern_prediction':expected})
        for label,M in [('unmixed', M0),('old_quartic', M0+.2*V),
                        ('fully_mixed',M0+.2*X/np.linalg.norm(X))]:
            dimension,singular_values = mass_charge_dimension(M,Ps)
            expected_dimension = exact['constraints'][label]['charge_space_dimension']
            assert dimension == expected_dimension
            edge_norms = [[float(np.linalg.norm(P@M@Q)) for Q in Ps] for P in Ps]
            models.append({'seed':seed,'mixing':label,'charge_space_dimension':dimension,
                           'singular_values':singular_values,'block_norms':edge_norms})
        if seed == 2026090700:
            saved = {'M0':M0,'balanced_V':V,'full_V':X/np.linalg.norm(X),'projectors':Ps}

    # Same minimal vertex and mass term: Ward residual detects forbidden mixing.
    gamma = gamma_matrices()
    p = np.array([1.2,.3,-.1,.4]); q = np.array([.4,-.2,.5,.1])
    M = saved['M0']+.2*saved['balanced_V']
    S1 = np.kron(slash(p+q,gamma),I)-np.kron(np.eye(4),M)
    S0 = np.kron(slash(p,gamma),I)-np.kron(np.eye(4),M)
    examples = []
    for label,charges in [('common_charge',[-1,-1,-1,-1]),
                          ('equal_mixed_charges',[0,-1,-1,2]),
                          ('different_mixed_charges',[0,-1,0,2])]:
        Q = np.einsum('b,bij->ij',charges,Ps)
        Q4 = np.kron(np.eye(4),Q)
        lhs = np.kron(slash(q,gamma),Q)
        rhs = S1@Q4-Q4@S0
        defect = lhs-rhs
        predicted = np.kron(np.eye(4),comm(M,Q))
        residuals['Ward_identity_mass_defect'] = max(residuals['Ward_identity_mass_defect'],
                                                   float(np.linalg.norm(defect-predicted)))
        U = expm(1j*.37*Q)
        examples.append({'case':label,'charges':charges,'mass_commutator_norm':float(np.linalg.norm(comm(Q,M))),
                         'stripped_Ward_residual':float(np.linalg.norm(defect)),
                         'finite_mass_covariance_defect':float(np.linalg.norm(U.conj().T@M@U-M))})
    np.savez_compressed(ROOT/'charge_operators.npz',**saved)
    assert max(residuals.values()) < 1e-10
    assert examples[0]['stripped_Ward_residual'] < 1e-10
    assert examples[1]['stripped_Ward_residual'] < 1e-10
    assert abs(examples[2]['stripped_Ward_residual']-.4) < 1e-10
    summary = {'status':'DIAGNOSTIC numerical realization of exact conditional identities',
               'seed_count':16,'scan_count':len(scans),'complex_dimensions':COMPLEX_DIMS,
               'minimum_trial_mass_eigenvalue':min_mass,
               'mass_units':'dimensionless; affine T2 baseline is a chosen positive trial operator',
               'gauge_identity_examples':examples,'residuals':residuals,
               'trial_charge_dimensions':{k:v['charge_space_dimension'] for k,v in exact['constraints'].items()},
               'scans':scans,'models':models}
    dump('charge_results.json',summary)
    print('Charge study:',len(scans),'scan points; max identity residual',max(residuals.values()))


if __name__ == '__main__':
    main()
