"""
tests/test_quantum.py
======================
Regression tests for quantum computing module + original 5 projects.
"""
import math
import numpy as np
import pytest
from mtft.quantum import ArithmeticCode, holonomy_gate_set, skyrmion_number


class TestArithmeticCode:
    @pytest.mark.parametrize("n", [7, 12, 16, 24, 60])
    def test_encode_decode_roundtrip(self, n):
        code = ArithmeticCode(d_logical=2, n_physical=n, y=0.18)
        for state in [np.array([1,0], dtype=complex), np.array([0,1], dtype=complex)]:
            enc = code.encode(state)
            dec = code.decode(enc)
            fid = abs(np.vdot(state, dec/np.linalg.norm(dec)))**2
            assert fid > 0.99, f"n={n}, fid={fid}"

    def test_gap_positive(self):
        for n in [7, 13, 16, 60]:
            code = ArithmeticCode(d_logical=2, n_physical=n, y=0.18)
            assert code.protection_gap > 0

    def test_divisor_scaling_differs(self):
        c_prime = ArithmeticCode(d_logical=2, n_physical=13, y=0.18)
        c_composite = ArithmeticCode(d_logical=2, n_physical=12, y=0.18)
        assert c_prime.divisor_density != c_composite.divisor_density
        assert c_prime.code_distance != c_composite.code_distance

    def test_composite_stronger_gap(self):
        c13 = ArithmeticCode(d_logical=2, n_physical=13, y=0.18)
        c60 = ArithmeticCode(d_logical=2, n_physical=60, y=0.18)
        assert c60.protection_gap > c13.protection_gap

    def test_code_info(self):
        info = ArithmeticCode(d_logical=2, n_physical=16, y=0.18).code_info()
        assert 'sigma_0' in info
        assert 'divisor_density' in info
        assert info['sigma_0'] == 5

    def test_syndrome_below_gap(self):
        code = ArithmeticCode(d_logical=2, n_physical=16, y=0.18)
        psi = code.encode(np.array([1, 0], dtype=complex))
        rng = np.random.default_rng(42)
        nv = rng.standard_normal(len(psi)) + 1j * rng.standard_normal(len(psi))
        nv /= np.linalg.norm(nv)
        noisy = psi + 0.04 * nv
        noisy /= np.linalg.norm(noisy)
        check = code.syndrome_check(noisy)
        assert check['error_magnitude'] < code.protection_gap


class TestHolonomyGates:
    def test_gate_set_complete(self):
        gates = holonomy_gate_set(2)
        assert set(gates.keys()) == {"X", "Y", "Z", "H", "T"}

    def test_all_unitary(self):
        gates = holonomy_gate_set(2)
        for name, gate in gates.items():
            assert gate.is_unitary(), f"{name} not unitary"

    def test_X_gate_flips(self):
        gates = holonomy_gate_set(2)
        result = gates["X"].apply(np.array([1, 0], dtype=complex))
        assert abs(abs(result[1]) - 1.0) < 1e-6


class TestSkyrmion:
    def test_skyrmion_basic(self):
        assert skyrmion_number(0, 1) == -1
        assert skyrmion_number(0, 5) == -5
        assert skyrmion_number(3, 3) == 0


class TestQuantumProjects:
    """Regression tests for the 5 original quantum simulation projects."""

    def test_p1_ft_circuit(self):
        code = ArithmeticCode(d_logical=2, n_physical=16, y=0.18)
        gates = holonomy_gate_set(2)
        psi = code.encode(gates["X"].apply(np.array([1, 0], dtype=complex)))
        clean = psi.copy()
        rng = np.random.default_rng(42)
        for _ in range(10):
            nv = rng.standard_normal(len(psi)) + 1j * rng.standard_normal(len(psi))
            nv /= np.linalg.norm(nv)
            noisy = psi + 0.04 * rng.uniform(0.5, 2.0) * nv
            noisy /= np.linalg.norm(noisy)
            check = code.syndrome_check(noisy)
            assert check['error_magnitude'] < code.protection_gap
            psi = clean.copy()

    def test_p3_oam_qkd(self):
        rng = np.random.default_rng(143)
        alice = rng.integers(0, 2, size=16)
        enc = {0: (0, 1), 1: (0, 5)}
        bob = [0 if skyrmion_number(*enc[b]) == -1 else 1 for b in alice]
        assert list(alice) == bob

    def test_p5_crosstalk(self):
        cA = ArithmeticCode(d_logical=2, n_physical=16, y=0.18)
        cB = ArithmeticCode(d_logical=2, n_physical=16, y=0.18)
        pA = cA.encode(np.array([0, 1], dtype=complex))
        pB = cB.encode(np.array([1, 0], dtype=complex))
        cA_c, cB_c = pA.copy(), pB.copy()
        rng = np.random.default_rng(999)
        for _ in range(10):
            shared = rng.standard_normal(len(pA)) + 1j * rng.standard_normal(len(pA))
            shared /= np.linalg.norm(shared)
            burst = shared * rng.uniform(0.15, 0.25)
            nA = pA + burst; nA /= np.linalg.norm(nA)
            nB = pB + burst; nB /= np.linalg.norm(nB)
            assert cA.syndrome_check(nA)['error_magnitude'] < cA.protection_gap
            assert cB.syndrome_check(nB)['error_magnitude'] < cB.protection_gap
            pA, pB = cA_c.copy(), cB_c.copy()
