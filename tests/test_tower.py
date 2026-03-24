"""
tests/test_tower.py
====================
Regression tests for the Multi-N tower module.
"""
import math
import numpy as np
import pytest
from mtft.tower import (
    tower_stiffness, even_n_universality, phase_transition_scaling,
    boundary_tracking, arithmetic_genome, character_orthogonality,
    arithmetic_periodic_table,
)
from mtft.arithmetic import mass_gap_stiffness


class TestMassGapPositive:
    """μ_N(y) > 0 unconditionally for all N ≥ 2, y > 0."""

    @pytest.mark.parametrize("N", [2, 3, 5, 7, 10, 13])
    def test_gap_positive_moderate_y(self, N):
        for y in [0.05, 0.10, 0.18, 0.50, 1.0]:
            mu = mass_gap_stiffness(y, N=N, n_max=200)
            assert mu > 0, f"μ_{N}({y}) = {mu} ≤ 0"


class TestEvenNUniversality:
    def test_all_even_collapse(self):
        univ = even_n_universality(y=0.10, N_max=30)
        assert univ['universal'], f"Max deviation: {univ['max_deviation']}"

    def test_deviation_machine_precision(self):
        univ = even_n_universality(y=0.10, N_max=20)
        assert univ['max_deviation'] < 1e-10


class TestN2Scaling:
    def test_scaling_converges(self):
        pts = phase_transition_scaling(list(range(3, 16)))
        n2_vals = [p['N2_yc'] for p in pts if p['y_c'] is not None]
        assert len(n2_vals) > 5
        # Should converge toward ~1.6-2.2
        assert all(0.5 < v < 5.0 for v in n2_vals)


class TestBoundaryTracking:
    def test_all_confined(self):
        results = boundary_tracking()
        for r in results:
            assert r['status'] == 'CONFINED'
            assert r['mu_N'] > 0

    def test_gap_grows_with_N(self):
        results = boundary_tracking([2, 5, 10, 20])
        mus = [r['mu_N'] for r in results]
        # Gap should grow dramatically under tracking
        assert mus[-1] > mus[0] * 100


class TestArithmeticGenome:
    def test_su3_suppresses_p3(self):
        g = arithmetic_genome(3, y=0.10)
        p3_frac = g['prime_contributions'][3]['fraction']
        assert p3_frac < 0.01, f"SU(3) should suppress p=3, got {p3_frac}"

    def test_su5_suppresses_p5(self):
        g = arithmetic_genome(5, y=0.10)
        p5_frac = g['prime_contributions'][5]['fraction']
        assert p5_frac < 0.01

    def test_prime_power_dominant(self):
        g = arithmetic_genome(2, y=0.10)
        assert g['prime_power_total_fraction'] > 0.90


class TestCharacterOrthogonality:
    @pytest.mark.parametrize("N", [3, 5, 7, 13])
    def test_identity_exact(self, N):
        result = character_orthogonality(N, y=0.10)
        assert result['identity_verified']

    def test_convergence_to_mu_inf(self):
        result = character_orthogonality(13, y=0.10)
        assert result['convergence_ratio'] > 0.99


class TestPeriodicTable:
    def test_three_tiers(self):
        table = arithmetic_periodic_table(N_max=10, y=0.10)
        tiers = {r['tier'] for r in table}
        assert 'floor' in tiers
        assert 'peak' in tiers

    def test_even_are_floor(self):
        table = arithmetic_periodic_table(N_max=10, y=0.10)
        for r in table:
            if r['is_even']:
                assert r['tier'] == 'floor'
