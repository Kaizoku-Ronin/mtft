"""
Verified regression tests for the Speiser–Hadamard lab in riemann.py
(audit Addendum I, July 2026 — every anchor independently recomputed).

Covers:
  * negative real zeros of ζ′ via the exact H2 functional-equation solver
  * the certified 19-zero nontrivial census (Speiser condition Re > 1/2)
  * the Hadamard identity ∂²log(−ζ′) = 2/(s−1)² − Σ(s−ρ′)⁻²
  * the decomposition lemma g_D(β) and its von Mangoldt cross-check
  * the exact shift identity μ_N(y) = (1/4π²)[T″(y) − Re T″(y − i/N)]
  * the modularity no-go: −1/4 · ln²(1/X) cusp coefficient
"""

import math

import pytest

from mtft.riemann import (
    ZETAPRIME_ZEROS,
    ZETAPRIME_CENSUS_HEIGHT,
    zetaprime_zero_count_berndt,
    zetaprime_negative_zero,
    zetaprime_refine,
    zetaprime_logcurvature,
    hadamard_zetaprime_check,
    dirichlet_curvature,
    von_mangoldt_curvature,
    divisor_log_weights,
    weighted_theta,
    filtered_moment_identity,
    weighted_theta_cusp_fit,
    _mp,
)


class TestNegativeZeros:
    """One zero of ζ′ per interval (−2n−2, −2n), solved exactly via H2."""

    ANCHORS = {
        1: -2.7172628292,
        2: -4.9367621086,
        3: -7.0745971450,
    }

    def test_first_three_anchors(self):
        for n, val in self.ANCHORS.items():
            assert abs(zetaprime_negative_zero(n) - val) < 1e-9

    def test_one_per_interval(self):
        for n in range(1, 21):
            r = zetaprime_negative_zero(n)
            assert -2 * n - 2 < r < -2 * n

    def test_actual_zetaprime_sign_change(self):
        mp = _mp()
        with mp.workdps(20):
            for n in range(1, 6):
                r = mp.mpf(zetaprime_negative_zero(n))
                f_lo = mp.zeta(r - mp.mpf("1e-8"), derivative=1)
                f_hi = mp.zeta(r + mp.mpf("1e-8"), derivative=1)
                assert mp.sign(f_lo) != mp.sign(f_hi)


class TestNontrivialCensus:
    """The certified 19-zero census to height 100."""

    def test_census_size(self):
        assert len(ZETAPRIME_ZEROS) == 19

    def test_all_within_height(self):
        assert all(0 < z.imag <= ZETAPRIME_CENSUS_HEIGHT
                   for z in ZETAPRIME_ZEROS)

    def test_speiser_condition(self):
        """All zeros have Re > 1/2 (RH-consistent; Speiser 1935)."""
        assert min(z.real for z in ZETAPRIME_ZEROS) > 0.5

    def test_first_zero_anchor(self):
        z0 = ZETAPRIME_ZEROS[0]
        assert abs(z0.real - 2.4631618694543213) < 1e-12
        assert abs(z0.imag - 23.298320492762858) < 1e-12

    def test_refine_recovers_census_zero(self):
        mp = _mp()
        r = zetaprime_refine(complex(2.46, 23.3))
        assert abs(r - ZETAPRIME_ZEROS[0]) < 1e-10
        with mp.workdps(20):
            assert abs(mp.zeta(mp.mpc(r.real, r.imag), derivative=1)) < 1e-12

    def test_berndt_count_consistent(self):
        """N′(100) = 17.1 + O(log T); certified count is 19."""
        assert abs(zetaprime_zero_count_berndt(100.0) - 19) < 4


class TestHadamardIdentity:
    """∂²log(−ζ′(s)) = 2/(s−1)² − Σ_{ρ′} (s−ρ′)⁻²."""

    def test_logcurvature_anchor(self):
        assert abs(zetaprime_logcurvature(8.0)
                   - 0.014052477562801085548) < 1e-14
        assert abs(zetaprime_logcurvature(30.0)
                   - 1.3597768711309190e-6) < 1e-16

    def test_pole_coefficient_is_two(self):
        """Direct from ζ derivatives, independent of the Hadamard product."""
        assert abs(zetaprime_logcurvature(1.001) * 1e-6 - 2.0) < 2e-5

    @pytest.mark.parametrize("s", [3.0, 5.0, 8.0])
    def test_identity_battery(self, s):
        """Residuals ~1e-5: 20x below the smooth tail's ±2e-4 uncertainty."""
        out = hadamard_zetaprime_check(s, n_neg=120)
        assert abs(out["residual"]) < 1e-4


class TestDecompositionLemma:
    """g_D(β) = ∂²log ζ(β) + ∂²log(−ζ′(β+1))."""

    def test_gD_anchor_beta3(self):
        out = dirichlet_curvature(3.0)
        assert abs(out["g_D"] - 0.33510387864414189) < 1e-14
        assert abs(out["zetaprime_share"] - 0.48588864) < 1e-6

    def test_gD_anchors_beta25_beta4(self):
        out25 = dirichlet_curvature(2.5)
        out4 = dirichlet_curvature(4.0)
        assert abs(out25["g_D"] - 0.60337094150219630) < 1e-14
        assert abs(out25["zetaprime_share"] - 0.41917573) < 1e-6
        assert abs(out4["g_D"] - 0.13359050881584953) < 1e-14
        assert abs(out4["zetaprime_share"] - 0.58039014) < 1e-6

    def test_von_mangoldt_crosscheck(self):
        zeta_piece = dirichlet_curvature(3.0)["zeta_piece"]
        vm = von_mangoldt_curvature(3.0, n_max=20000)
        assert abs(vm - zeta_piece) < 1e-6


class TestShiftIdentity:
    """μ_N(y) = (1/4π²)[T″(y) − Re T″(y − i/N)] — exact."""

    def test_identity_holds(self):
        out = filtered_moment_identity(0.05, N=3, n_max=2000)
        assert out["rel_diff"] < 1e-12

    def test_identity_other_N(self):
        for N in (4, 5):
            out = filtered_moment_identity(0.08, N=N, n_max=2000)
            assert out["rel_diff"] < 1e-12

    def test_N3_filter_is_center_projector(self):
        """(1 − cos(2πn/3)) = 0 on 3|n, 3/2 else — the SU(3) filter."""
        for n in range(1, 31):
            expect = 0.0 if n % 3 == 0 else 1.5
            assert abs((1 - math.cos(2 * math.pi * n / 3)) - expect) < 1e-12


class TestModularityNoGo:
    """Double pole at s=0 forces −1/4·ln²(1/X): Θ̃ cannot be modular."""

    def test_weights_match_dirichlet_series(self):
        """w_n = Σ_{d|n} (log d)/d  ⟺  W(s) = −ζ(s)ζ′(s+1)."""
        w = divisor_log_weights(12)
        # w_6 = log2/2 + log3/3 + log6/6
        assert abs(w[6] - (math.log(2)/2 + math.log(3)/3
                           + math.log(6)/6)) < 1e-15
        assert w[4] == pytest.approx(math.log(2)/2 + math.log(4)/4)

    def test_cusp_log_squared_coefficient(self):
        out = weighted_theta_cusp_fit(y_points=(1e-2, 1e-3, 1e-4),
                                      n_max=200000)
        assert abs(out["A_fit"] - (-0.25)) < 5e-3
        assert abs(out["C_fit"] - out["C_predicted_approx"]) < 0.05
