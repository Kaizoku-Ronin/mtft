"""
tests/test_predictions.py
==========================
Regression tests for all MTFT zero-parameter predictions.
"""
import math
import pytest
from mtft.constants import (
    FEIGENBAUM_DELTA as D, T_INF, TORQUE_FULL, EULER_GAMMA as G,
    LAMBERT_OMEGA as O, XI, PI, GAUGE, SM, PDG, QUARKS, LEPTONS,
)


class TestGaugeSector:
    def test_alpha_inverse_leading(self):
        assert abs(2 * PI * D**2 - 137.0) < 0.1

    def test_alpha_inverse_3term(self):
        a = 2*PI*D**2 + 1/(4*D) - XI*T_INF/D**6
        assert abs(a - 137.035999) < 0.001

    def test_alpha_s(self):
        assert abs(T_INF/4 - 0.1172) < 0.002

    def test_alpha_s_13(self):
        assert abs(13**(-5/6) - 0.1180) < 0.002

    def test_weinberg_angle(self):
        assert abs(3/13 - 0.2308) < 0.001

    def test_W_Z_ratio(self):
        assert abs(1/(2*O) - PDG.m_W/PDG.m_Z) < 0.001

    def test_charge_gaussian(self):
        e_pred = math.sqrt(2) / D
        e_obs = math.sqrt(4 * PI / PDG.alpha_inv)
        assert abs(e_pred - e_obs) / e_obs < 0.001

    def test_coupling_shift_jacobian(self):
        from mtft.falsify import coupling_shift
        cs = coupling_shift(1e-6)
        assert cs.delta_alpha_s != 0
        assert cs.delta_sin2_tW != 0
        assert abs(cs.delta_sin2_tW) < abs(cs.delta_alpha_s)


class TestHiggsSector:
    def test_higgs_mass(self):
        m_H = 246.22 * G / (2*O)
        assert abs(m_H - 125.25) / 125.25 < 0.001

    def test_higgs_quartic(self):
        lam = G**2 / (8*O**2)
        lam_obs = PDG.m_H**2 / (2 * 246.22**2)
        assert abs(lam - lam_obs) / lam_obs < 0.002

    def test_higgs_W_ratio(self):
        ratio = math.sqrt(PI) * TORQUE_FULL**2
        obs = PDG.m_H / PDG.m_W
        assert abs(ratio - obs) / obs < 0.001

    def test_hosotani_mtft(self):
        from mtft.hosotani import HosotaniMTFT
        hm = HosotaniMTFT()
        m = hm.gauge_masses()
        assert abs(m['m_H'] - 125.3) < 0.5
        assert abs(m['sin2_theta_W'] - 3/13) < 1e-6


class TestLeptonSector:
    def test_koide_ratio(self):
        from mtft.koide import koide_ratio
        m_e, m_mu, m_tau = LEPTONS.e*1e3, LEPTONS.mu*1e3, LEPTONS.tau*1e3
        Q = koide_ratio(m_e, m_mu, m_tau)
        assert abs(Q - 2/3) < 1e-4

    def test_tau_from_koide(self):
        from mtft.koide import predict_tau_mass
        pred = predict_tau_mass()
        assert pred['error_percent'] < 0.01

    def test_generation_count(self):
        from mtft.x0_143 import generation_count
        assert generation_count() == 3


class TestQuarkSector:
    def test_md_mu_ratio(self):
        ratio = math.sqrt(D)
        obs = QUARKS.d / QUARKS.u
        assert abs(ratio - obs) / obs < 0.01

    def test_mc_mmu_ratio(self):
        obs = QUARKS.c / LEPTONS.mu
        assert abs(12 - obs) / obs < 0.02

    def test_mb_mtau_ratio(self):
        pred = 13**(1/3)
        obs = QUARKS.b / LEPTONS.tau
        assert abs(pred - obs) / obs < 0.01


class TestCosmology:
    def test_confinement_depth(self):
        y_pred = (G + O) / (2 * PI)
        assert abs(y_pred - 0.18174) < 0.001

    def test_desert_check(self):
        from mtft.falsify import desert_check
        result = desert_check()
        assert not result['falsified']

    def test_friedmann_today(self):
        from mtft.cosmology import FriedmannMTFT
        f = FriedmannMTFT()
        H2 = f.H_squared(1.0)
        assert abs(H2 - 1.0) < 0.1

    def test_expansion_acceleration(self):
        from mtft.cosmology import FriedmannMTFT
        hist = FriedmannMTFT().expansion_history()
        # q should cross zero (acceleration onset)
        assert any(hist['q'] < 0) and any(hist['q'] > 0)
