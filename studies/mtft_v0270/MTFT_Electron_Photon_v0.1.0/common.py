"""Frozen-input readers and explicitly recorded physical reference inputs."""
import json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
BLOCKS = ['ell', 'old', 'q4', 'q6']
REAL_DIMS = [2, 4, 8, 12]
COMPLEX_DIMS = [1, 2, 4, 6]
ALPHA_INV_REFERENCE = 137.035999177
ALPHA_INV_REFERENCE_UNCERTAINTY = 0.000000021
ELECTRON_REST_ENERGY_EV = 510998.95069
REFERENCE_URL = 'https://physics.nist.gov/cuu/Constants/Table/allascii.txt'


def dump(name, value):
    (ROOT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def read(name):
    return json.loads((ROOT / name).read_text())


def comm(a, b):
    return a @ b - b @ a


def herm(a):
    return (a + a.conj().T) / 2


def internal_stage():
    """Return operators on the +i eigenspace of the Hodge complex structure.

    All adjoints below use the frozen Hodge metric, whitened to the identity.
    A 13-dimensional internal space does not specify 13 physical particles.
    """
    data = read('inputs/x0143_certified.json')['arrays']
    E = np.array(data['intersection_cycles'], float)
    J = np.array(data['J_true'], float)
    G = herm(E @ J)
    if np.linalg.eigvalsh(G)[0] < 0:
        J, G = -J, -G
    eig, vec = np.linalg.eigh(G)
    assert eig[0] > 0
    R = (vec * np.sqrt(eig)) @ vec.T
    Ri = (vec / np.sqrt(eig)) @ vec.T
    Jw = R @ J @ Ri
    values, basis = np.linalg.eigh(1j * (Jw - Jw.T) / 2)
    B = basis[:, :13]
    assert np.max(np.abs(values[:13] + 1)) < 1e-10
    certificate = read('inputs/sector_certificate.json')
    Ps = [np.array(sp.Matrix(certificate['projectors'][b]), float) for b in BLOCKS]
    convert = lambda A: B.conj().T @ R @ A @ Ri @ B
    result = {name: convert(np.array(data[name], float)) for name in ['T2', 'T3', 'W11', 'W13']}
    result['P'] = np.array([convert(P) for P in Ps])
    result['residuals'] = {
        'J_squared_plus_identity': float(np.linalg.norm(Jw @ Jw + np.eye(26))),
        'J_skew': float(np.linalg.norm(Jw + Jw.T)),
        'complex_basis_isometry': float(np.linalg.norm(B.conj().T @ B - np.eye(13))),
        'J_eigenspace': float(np.linalg.norm(Jw @ B - 1j * B)),
        'complex_subspace_invariance': float(max(
            np.linalg.norm((np.eye(26) - B @ B.conj().T) @ R @ A @ Ri @ B)
            for A in Ps + [np.array(data['T2'], float)]
        )),
        'projector_hermiticity': float(max(np.linalg.norm(P-P.conj().T) for P in result['P'])),
        'projector_completeness': float(np.linalg.norm(result['P'].sum(axis=0)-np.eye(13))),
        'T2_hermiticity': float(np.linalg.norm(result['T2']-result['T2'].conj().T)),
    }
    return result


def gamma_matrices():
    """Supplied Dirac representation, metric signature (+---)."""
    pauli = [np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]), np.diag([1,-1])]
    zero = np.zeros((2,2))
    return np.array([np.diag([1,1,-1,-1])] +
                    [np.block([[zero, s], [-s, zero]]) for s in pauli], complex)


def slash(p, gamma):
    return np.einsum('m,mab->ab', np.array([1,-1,-1,-1])*p, gamma)
