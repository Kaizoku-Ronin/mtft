"""Shared constants and input readers; no MTFT installation is required to rerun."""
import json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
X = sp.Symbol('x')
POLYS = {'q4': sp.Poly(X**4-3*X**3-X**2+5*X+1, X),
         'q6': sp.Poly(X**6-10*X**4+2*X**3+24*X**2-7*X-12, X)}
BLOCKS = ['ell', 'old', 'q4', 'q6']
DIMS = np.array([2, 4, 8, 12])
SIGNS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
CLASSES = ['hecke_commuting', 'block_preserving', 'AL_preserving',
           'balanced_old_q4', 'unrestricted']
CLASS_LABELS = ['Hecke commuting', 'Within sectors', 'Atkin–Lehner preserving',
                'Balanced old–quartic', 'Unrestricted']

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def read(name):
    return json.loads((ROOT / name).read_text())

def frozen():
    return {k: np.array(v) for k, v in read('inputs/x0143_certified.json')['arrays'].items()}

def matrix_eval(poly, matrix):
    answer = sp.zeros(matrix.rows)
    for c in sp.Poly(poly, X).all_coeffs():
        answer = answer * matrix + c * sp.eye(matrix.rows)
    return answer

def rational_matrix(matrix):
    return [[str(v) for v in row] for row in matrix.tolist()]

def numeric_stage():
    data = frozen()
    E = data['intersection_cycles'].astype(float)
    J = data['J_true']
    G = (E @ J + (E @ J).T) / 2
    if np.linalg.eigvalsh(G)[0] < 0:
        J = -J
        G = -G
    values, vectors = np.linalg.eigh(G)
    assert values[0] > 0
    R = (vectors * np.sqrt(values)) @ vectors.T
    Ri = (vectors / np.sqrt(values)) @ vectors.T
    cert = read('sector_certificate.json')
    Ps = np.array([R @ np.array(sp.Matrix(cert['projectors'][k]), float) @ Ri for k in BLOCKS])
    return {'J': R @ J @ Ri, 'E': Ri @ E @ Ri, 'G': G,
            'T2': R @ data['T2'] @ Ri, 'T3': R @ data['T3'] @ Ri,
            'W11': R @ data['W11'] @ Ri, 'W13': R @ data['W13'] @ Ri,
            'P': Ps, 'R': R, 'Ri': Ri}
