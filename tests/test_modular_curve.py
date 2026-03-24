"""
tests/test_modular_curve.py
============================
Pytest suite for mtft.modular_curve.

Validates genus, cusps, index, vortex energy, Hecke spectrum,
and homology against known values from LMFDB / Shimura theory.

Run:  pytest tests/test_modular_curve.py -v
"""

import math
import numpy as np
import pytest
from mtft.modular_curve import ModularCurve, X0, VortexConfig, HomologyData


# ── Known invariants from LMFDB ─────────────────────────────
# Format: (level, genus, num_cusps, index)

KNOWN_CURVES = [
    (1,    0,  1,   1),
    (2,    0,  2,   3),
    (3,    0,  2,   4),
    (5,    0,  2,   6),
    (7,    0,  2,   8),
    (11,   1,  2,  12),
    (13,   0,  2,  14),
    (14,   1,  4,  24),
    (15,   1,  4,  24),
    (17,   1,  2,  18),
    (23,   2,  2,  24),
    (26,   2,  4,  42),
    (37,   2,  2,  38),
    (43,   3,  2,  44),
    (67,   5,  2,  68),
    (143, 13,  4, 168),
    (169,   8, 14, 182),
]


class TestModularCurveInvariants:
    """Test genus, cusps, and index against known values."""

    @pytest.mark.parametrize("N,g,c,idx", KNOWN_CURVES)
    def test_genus(self, N, g, c, idx):
        X = ModularCurve(N)
        assert X.genus == g, f"X₀({N}): expected genus {g}, got {X.genus}"

    @pytest.mark.parametrize("N,g,c,idx", KNOWN_CURVES)
    def test_num_cusps(self, N, g, c, idx):
        X = ModularCurve(N)
        assert X.num_cusps == c, f"X₀({N}): expected {c} cusps, got {X.num_cusps}"

    @pytest.mark.parametrize("N,g,c,idx", KNOWN_CURVES)
    def test_index(self, N, g, c, idx):
        X = ModularCurve(N)
        assert X.index == idx, f"X₀({N}): expected index {idx}, got {X.index}"


class TestX0_143:
    """Detailed tests for the MTFT-critical level 143 = 11 × 13."""

    @pytest.fixture
    def X(self):
        return ModularCurve(143)

    def test_factorization(self, X):
        assert X._factors == {11: 1, 13: 1}

    def test_genus_13(self, X):
        assert X.genus == 13

    def test_four_cusps(self, X):
        assert X.num_cusps == 4

    def test_index_168(self, X):
        assert X.index == 168

    def test_no_elliptic_points(self, X):
        """143 = 11 × 13, both ≡ −1 mod 4, so ε₂ = 0.
        11 ≡ 2 mod 3 ⟹ (−3/11) = −1, so ε₃ = 0."""
        assert X.num_elliptic_2 == 0
        assert X.num_elliptic_3 == 0

    def test_cusp_labels(self, X):
        labels = {c.representative for c in X.cusps}
        # Should contain ∞, 0, and representatives for d=11, d=13
        assert "∞" in labels
        assert "0" in labels

    def test_cusp_widths_sum(self, X):
        """Sum of cusp widths must equal the index."""
        total = sum(c.width for c in X.cusps)
        assert total == X.index, (
            f"Sum of cusp widths {total} ≠ index {X.index}"
        )

    def test_euler_characteristic(self, X):
        assert 2 - 2 * X.genus == -24

    def test_homology_rank(self, X):
        H = X.homology()
        assert H.rank == 26
        assert H.monodromy_cycles == 26

    def test_homology_surplus(self, X):
        """26 cycles − 12 SM generators = 14 surplus."""
        H = X.homology()
        assert H.surplus_cycles == 14

    def test_summary_dict(self, X):
        s = X.summary()
        assert s["level"] == 143
        assert s["genus"] == 13
        assert s["index"] == 168
        assert s["num_cusps"] == 4

    def test_repr(self, X):
        r = repr(X)
        assert "143" in r
        assert "13" in r


class TestVortexConfig:
    """Test τ-vortex energy computations."""

    def test_bps_single_charge(self):
        v = VortexConfig(winding_numbers={"∞": 1, "0": 0})
        assert v.total_charge == 1
        assert abs(v.bps_bound - 2 * math.pi) < 1e-10
        assert abs(v.energy - 2 * math.pi) < 1e-10
        assert abs(v.bps_ratio - 1.0) < 1e-10

    def test_bps_cancel(self):
        """Opposite winding numbers: Q=0 but E>0."""
        v = VortexConfig(winding_numbers={"∞": 1, "0": -1})
        assert v.total_charge == 0
        assert v.bps_bound == 0.0
        assert v.energy > 0

    def test_bps_ratio_geq_1(self):
        """E/E_BPS ≥ 1 always (Bogomolny bound)."""
        for ns in [(1, 1), (2, -1), (1, 0, -1, 1), (3, 0, 0, 0)]:
            labels = [f"c{i}" for i in range(len(ns))]
            v = VortexConfig(winding_numbers=dict(zip(labels, ns)))
            if v.total_charge != 0:
                assert v.bps_ratio >= 1.0 - 1e-10

    def test_zero_config(self):
        v = VortexConfig(winding_numbers={"∞": 0, "0": 0})
        assert v.total_charge == 0
        assert v.energy == 0.0

    def test_default_vortex(self):
        X = ModularCurve(143)
        v = X.vortex_energy()
        assert v.total_charge == 1
        assert "∞" in v.winding_numbers


class TestHeckeSpectrum:
    """Test Hecke eigenvalue generation."""

    def test_ramanujan_bound(self):
        """All normalised eigenvalues must satisfy |a_p/(2√p)| ≤ 1."""
        X = ModularCurve(143)
        hs = X.hecke_spectrum(max_prime=200, seed=42)
        assert hs.ramanujan_violations == 0
        assert np.all(np.abs(hs.eigenvalues) <= 1.0 + 1e-10)

    def test_primes_exclude_bad(self):
        """Primes dividing N should be excluded."""
        X = ModularCurve(143)
        hs = X.hecke_spectrum(max_prime=200, seed=42)
        assert 11 not in hs.primes
        assert 13 not in hs.primes

    def test_newform_count_positive(self):
        X = ModularCurve(143)
        hs = X.hecke_spectrum(seed=42)
        assert hs.num_newforms >= 1

    def test_reproducible(self):
        X = ModularCurve(143)
        h1 = X.hecke_spectrum(seed=42)
        h2 = X.hecke_spectrum(seed=42)
        np.testing.assert_array_equal(h1.eigenvalues, h2.eigenvalues)


class TestHomology:
    """Test symplectic intersection form."""

    def test_symplectic_form_shape(self):
        X = ModularCurve(143)
        H = X.homology()
        assert H.intersection_matrix.shape == (26, 26)

    def test_symplectic_antisymmetric(self):
        X = ModularCurve(143)
        J = X.homology().intersection_matrix
        np.testing.assert_array_equal(J, -J.T)

    def test_symplectic_determinant(self):
        """det(J) = 1 for standard symplectic form."""
        X = ModularCurve(143)
        J = X.homology().intersection_matrix
        assert abs(np.linalg.det(J) - 1.0) < 1e-10

    def test_J_squared(self):
        """J² = -I for the symplectic form."""
        X = ModularCurve(143)
        J = X.homology().intersection_matrix
        np.testing.assert_array_equal(J @ J, -np.eye(26, dtype=int))


class TestFordCircles:
    """Test Ford circle generation."""

    def test_nonempty(self):
        X = ModularCurve(143)
        circles = X.ford_circles(max_denom=20)
        assert len(circles) > 0

    def test_positive_radius(self):
        X = ModularCurve(143)
        for cx, cy, r, cls in X.ford_circles(max_denom=20):
            assert r > 0


class TestConvenienceConstructor:
    def test_X0_shorthand(self):
        X = X0(143)
        assert X.genus == 13
        assert isinstance(X, ModularCurve)


class TestEdgeCases:
    def test_level_1(self):
        X = ModularCurve(1)
        assert X.genus == 0
        assert X.num_cusps == 1

    def test_prime_level(self):
        X = ModularCurve(11)
        assert X.genus == 1  # Elliptic curve!

    def test_invalid_level(self):
        with pytest.raises(ValueError):
            ModularCurve(0)
        with pytest.raises(ValueError):
            ModularCurve(-5)
