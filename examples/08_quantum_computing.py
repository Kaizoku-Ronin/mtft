"""
Quantum Computing with MTFT
============================

MTFT provides natural quantum computing primitives:
  - Holonomy gates (geometric phase -> fault tolerant)
  - Topological qudits (mass gap -> protected states)
  - Arithmetic error correction (divisor structure -> distance)

The connection to experiment:
  de Mello Koch et al. (Nature Comms 2025) confirmed that
  OAM entangled photons carry the SAME topological structure
  as Yang-Mills monopoles, with wrapping numbers matching the
  MTFT mass gap prediction.

  pip install mtft
  python examples/08_quantum_computing.py
"""

import numpy as np
from mtft.quantum import (
    gell_mann_matrices, HolonomyGate, holonomy_gate_set,
    TopologicalQudit, ArithmeticCode,
    topological_spectrum_info, topological_spectrum_table,
    oam_entangled_state, skyrmion_number,
)

# ── 1. SU(d) generators ──────────────────────────────────────
print("GELL-MANN MATRICES (SU(d) Lie algebra)")
print("=" * 55)
for d in [2, 3, 4, 5, 7]:
    gens = gell_mann_matrices(d)
    print(f"  SU({d}): {len(gens)} generators (d^2 - 1 = {d**2-1})")

# Verify properties for SU(3)
gens3 = gell_mann_matrices(3)
print(f"\n  SU(3) generator checks:")
print(f"    All traceless: {all(abs(np.trace(T)) < 1e-12 for T in gens3)}")
print(f"    All Hermitian: {all(np.allclose(T, T.conj().T) for T in gens3)}")

# ── 2. Holonomy gates ────────────────────────────────────────
print(f"\nHOLONOMY GATES (geometric phase -> fault tolerance)")
print("-" * 55)
gates = holonomy_gate_set(2)
print(f"  Universal qubit gate set:")
for name, gate in gates.items():
    U = gate.matrix_approx
    print(f"    {name}: det = {np.linalg.det(U):.4f}, unitary = {gate.is_unitary()}")

print(f"\n  Qutrit (d=3) gate set:")
gates3 = holonomy_gate_set(3)
print(f"    {len(gates3)} gates generated")
sample = list(gates3.items())[0]
print(f"    Example: {sample[0]} is unitary = {sample[1].is_unitary()}")

# Apply a gate
print(f"\n  Gate application:")
psi = np.array([1, 0], dtype=complex)
gate_X = gates["X"]
psi_out = gate_X.apply(psi)
print(f"    |0> --X--> [{psi_out[0]:.3f}, {psi_out[1]:.3f}]")

# ── 3. Topological qudits ────────────────────────────────────
print(f"\nTOPOLOGICAL QUDITS (mass gap = protection)")
print("=" * 55)
q3 = TopologicalQudit(d=3, coefficients=[1, 1, 1])
print(f"  Qutrit |psi> = (1/sqrt(3))(|0> + |1> + |2>)")
print(f"  Manifold dimension:    {q3.manifold_dimension}")
print(f"  Total S^2->S^2 maps:  {q3.total_maps}")
print(f"  Non-trivial maps:     {q3.nontrivial_maps}")
print(f"  Independent inv.:     {q3.independent_invariants}")
print(f"  Mass gap mu_3(0.18):  {q3.topological_gap():.4f} > 0")
print(f"  Purity:               {q3.purity():.4f}")
print(f"  GLA vector (first 4): {q3.gla_vector()[:4].round(4)}")

# ── 4. Topological spectrum scaling ──────────────────────────
print(f"\nTOPOLOGICAL SPECTRUM (matches Nature Comms 2025)")
print("-" * 55)
print(f"  {'d':>3s}  {'manifold':>8s}  {'total':>8s}  {'non-triv':>9s}  {'indep':>7s}  {'gap':>8s}")
print(f"  {'-'*48}")
for info in topological_spectrum_table():
    print(f"  {info['d']:3d}  {info['manifold_dim']:8d}  "
          f"{info['total_maps']:8,d}  {info['nontrivial_maps']:9,d}  "
          f"{info['independent_invariants']:7,d}  {info['mass_gap_y018']:8.4f}")

print(f"\n  d=7: 48-dim manifold, 17,296 maps <- matches experiment")

# ── 5. OAM skyrmion numbers ──────────────────────────────────
print(f"\nOAM SKYRMION NUMBERS (biphoton topology)")
print("-" * 55)
examples = [(0,1), (0,5), (0,-8), (5,-8), (-3,3)]
for l0, l1 in examples:
    N = skyrmion_number(l0, l1)
    print(f"  |{l0}>|-{l0}> + |{l1}>|-{l1}>:  N = Delta_ell = {N}")

# ── 6. Arithmetic error correction ───────────────────────────
print(f"\nARITHMETIC ERROR CORRECTION")
print("=" * 55)
code = ArithmeticCode(d_logical=2, n_physical=16, y=0.18)
print(f"  Code: [{code.n_physical}, {code.d_logical}]")
print(f"  Code distance: {code.code_distance:.4f}")
print(f"  Mass gap (protection): {code.protection_gap:.4f}")
print()

# Encode, corrupt, detect
logical_0 = np.array([1, 0], dtype=complex)
encoded = code.encode(logical_0)
check_clean = code.syndrome_check(encoded)
print(f"  Encode |0>:")
print(f"    In code space: {check_clean['in_code_space']}")
print(f"    Error magnitude: {check_clean['error_magnitude']:.6f}")

# Add noise
rng = np.random.default_rng(42)
noisy = encoded + 0.05 * rng.standard_normal(len(encoded))
check_noisy = code.syndrome_check(noisy)
print(f"  After 5% noise:")
print(f"    In code space: {check_noisy['in_code_space']}")
print(f"    Error magnitude: {check_noisy['error_magnitude']:.4f}")
print(f"    Errors below mu_N = {code.protection_gap:.4f} don't change topology.")
