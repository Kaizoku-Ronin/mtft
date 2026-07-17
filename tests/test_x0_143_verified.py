"""
tests/test_x0_143_verified.py
==============================
Regression suite for the independently verified X₀(143) arithmetic spine
(audit + Correction Sessions 1–6, July 2026).  Every test here encodes a
TRUTH that was cross-checked at least two independent ways (PARI polmods,
curve point-counts, trace-form totals, Newton identities, Chebotarev).

These tests would have caught bugs B1 (wrong f₁ column), B11 (phantom
polynomial) and the curve-eigenvalue errors at 37/41/43/47/53/59/61.
"""
import math

import numpy as np
import pytest

from mtft.x0_143 import (
    CURVE_143A1, ORBIT_TRACES_VERIFIED, ORBIT_TRACE_F1, ORBIT_TRACE_F2,
    ORBIT_TRACE_F3, TRACE_TOTALS_50, FIELD_POLY_F2, FIELD_POLY_F3,
    FIELD_DISCRIMINANTS, GALOIS_GROUPS, ROOT_NUMBERS_LIST,
    hecke_polynomial_f2_T2, hecke_polynomial_f3_T2,
    hecke_polynomial_f2_T3, hecke_polynomial_f3_T3,
    rankin_selberg_epsilon, rankin_selberg_Q, RS_COUPLING, ORBIT_DIMENSIONS,
)


# ── curve point-counting (independent of the stored table) ───────────────

def _ap_143a1(p: int) -> int:
    """a_p = p + 1 − #E(F_p) for y² + y = x³ − x² − x − 2, brute force."""
    cnt = 1  # point at infinity
    for x in range(p):
        rhs = (x**3 - x**2 - x - 2) % p
        for y in range(p):
            if (y * y + y - rhs) % p == 0:
                cnt += 1
    return p + 1 - cnt


class TestCurveEigenvalues:
    def test_stored_table_matches_point_counts(self):
        for p, a_stored in CURVE_143A1.hecke_eigenvalues.items():
            assert a_stored == _ap_143a1(p), f"a_{p}: stored {a_stored} ≠ point count"

    def test_hasse_bounds(self):
        for p, a in CURVE_143A1.hecke_eigenvalues.items():
            assert abs(a) <= 2 * math.sqrt(p) + 1e-9


class TestOrbitTraceSumRule:
    """B1 regression: per-orbit traces must sum to the trace form."""

    def test_prime_sum_rule(self):
        tot = {n: TRACE_TOTALS_50[n - 1] for n in range(1, 51)}
        for p, (a1, t2, t3) in ORBIT_TRACES_VERIFIED.items():
            assert a1 + t2 + t3 == tot[p], f"sum rule fails at p={p}"

    def test_full_sum_rule_to_50(self):
        for n in range(50):
            assert (ORBIT_TRACE_F1[n] + ORBIT_TRACE_F2[n] + ORBIT_TRACE_F3[n]
                    == TRACE_TOTALS_50[n]), f"sum rule fails at n={n+1}"

    def test_trace_dims_at_1(self):
        assert (ORBIT_TRACE_F1[0], ORBIT_TRACE_F2[0], ORBIT_TRACE_F3[0]) == (1, 4, 6)

    def test_multiplicativity_f1(self):
        a = [0] + ORBIT_TRACE_F1  # a[n] = a_n(f1)
        for m in range(2, 51):
            for n in range(2, 50 // m + 1):
                if math.gcd(m, n) == 1:
                    assert a[m * n] == a[m] * a[n], f"a_{m*n} ≠ a_{m}·a_{n}"

    def test_prime_power_recursion_f1(self):
        a = [0] + ORBIT_TRACE_F1
        for p in [2, 3, 5, 7]:
            k = 1
            while p ** (k + 2) <= 50:
                n1, n2, n3 = p ** k, p ** (k + 1), p ** (k + 2)
                assert a[n3] == a[p] * a[n2] - p * a[n1]
                k += 1


class TestCharpolys:
    """B11 regression: charpolys from Newton-identity moments of the traces."""

    @staticmethod
    def _newton_moment(coeffs_desc, k):
        """Sum of k-th powers of roots via Newton's identities."""
        d = len(coeffs_desc) - 1
        e = [0.0] + [-(-1) ** i * coeffs_desc[i] for i in range(1, d + 1)]
        # e_i = (-1)^i a_i (monic poly x^d + a_1 x^{d-1} + ...)
        e = [1.0] + [(-1) ** i * coeffs_desc[i] for i in range(1, d + 1)]
        p = [0.0] * (k + 1)
        p[0] = float(d)
        for m in range(1, k + 1):
            p[m] = (-1) ** (m - 1) * m * e[m] + sum(
                (-1) ** (i - 1) * e[i] * p[m - i] for i in range(1, m)
            )
        return p[k]

    def test_f2_T2_moments(self):
        poly = list(hecke_polynomial_f2_T2())
        assert self._newton_moment(poly, 1) == pytest.approx(3.0)   # Tr_f2(a2) = 3
        assert self._newton_moment(poly, 2) == pytest.approx(11.0)  # Tr(a2²) = Tr(a4)+2·4

    def test_f3_T2_moments(self):
        poly = list(hecke_polynomial_f3_T2())
        assert self._newton_moment(poly, 1) == pytest.approx(0.0)   # Tr_f3(a2) = 0
        assert self._newton_moment(poly, 2) == pytest.approx(20.0)  # Tr(a2²) = Tr(a4)+2·6

    def test_all_roots_real_and_ramanujan(self):
        for poly, bound in [(hecke_polynomial_f2_T2(), 2 * math.sqrt(2)),
                            (hecke_polynomial_f3_T2(), 2 * math.sqrt(2)),
                            (hecke_polynomial_f2_T3(), 2 * math.sqrt(3)),
                            (hecke_polynomial_f3_T3(), 2 * math.sqrt(3))]:
            roots = np.roots(list(poly))
            assert max(abs(r.imag) for r in roots) < 1e-6
            assert max(abs(r) for r in roots) < bound + 1e-6

    def test_phantom_polynomial_retired(self):
        """The wrong Paper-26 polynomial is gone from the module."""
        roots = np.roots(list(hecke_polynomial_f3_T2()))
        phantom = 0.5732 + 0.3564j
        assert min(abs(r - phantom) for r in roots) > 0.5


class TestFieldsAndGalois:
    @staticmethod
    def _bareiss_det(mat):
        """Exact integer determinant (Bareiss algorithm)."""
        M = [row[:] for row in mat]
        n = len(M)
        sign = 1
        prev = 1
        for k in range(n - 1):
            if M[k][k] == 0:
                for r in range(k + 1, n):
                    if M[r][k] != 0:
                        M[k], M[r] = M[r], M[k]
                        sign = -sign
                        break
                else:
                    return 0
            for i in range(k + 1, n):
                for j in range(k + 1, n):
                    M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
            prev = M[k][k]
        return sign * M[n - 1][n - 1]

    @classmethod
    def _discriminant(cls, coeffs_desc):
        """disc(P) = (−1)^{d(d−1)/2} · Res(P, P′) for monic P — exact integers."""
        d = len(coeffs_desc) - 1
        p = list(coeffs_desc)                       # descending
        dp = [coeffs_desc[i] * (d - i) for i in range(d)]  # derivative, descending
        # Sylvester matrix of P (deg d) and P′ (deg d−1), size 2d−1
        size = 2 * d - 1
        syl = [[0] * size for _ in range(size)]
        for i in range(d - 1):                      # rows of shifted P
            for j, c in enumerate(p):
                syl[i][i + j] = c
        for i in range(d):                          # rows of shifted P′
            for j, c in enumerate(dp):
                syl[d - 1 + i][i + j] = c
        res = cls._bareiss_det(syl)
        return res * (-1) ** (d * (d - 1) // 2)

    def test_discriminants(self):
        assert self._discriminant(FIELD_POLY_F2) == FIELD_DISCRIMINANTS["f2"]
        assert self._discriminant(FIELD_POLY_F3) == FIELD_DISCRIMINANTS["f3"]

    def test_fields_unramified_at_level_primes(self):
        for disc in FIELD_DISCRIMINANTS.values():
            assert disc % 11 != 0 and disc % 13 != 0

    def test_fields_totally_real(self):
        for poly in (FIELD_POLY_F2, FIELD_POLY_F3):
            roots = np.roots(list(poly))
            assert max(abs(r.imag) for r in roots) < 1e-8

    def test_galois_labels(self):
        assert GALOIS_GROUPS == {"f2": "S4", "f3": "S6"}


class TestRootNumbers:
    def test_up_eigenvalues_rational(self):
        """For p ‖ N the U_p eigenvalue is ±1: Tr(a_p) = ±dim per orbit."""
        dims = ORBIT_DIMENSIONS
        for traces, d in ((ORBIT_TRACE_F1, dims[0]), (ORBIT_TRACE_F2, dims[1]),
                          (ORBIT_TRACE_F3, dims[2])):
            for p in (11, 13):
                assert abs(traces[p - 1]) == d

    def test_root_number_vector(self):
        assert ROOT_NUMBERS_LIST == (-1, +1, +1)
        # ε = −w₁₁·w₁₃ with w_p = −a_p
        for i, traces in enumerate((ORBIT_TRACE_F1, ORBIT_TRACE_F2, ORBIT_TRACE_F3)):
            d = ORBIT_DIMENSIONS[i]
            a11 = traces[10] // d
            a13 = traces[12] // d
            eps = -((-a11) * (-a13))
            assert eps == ROOT_NUMBERS_LIST[i]


class TestRankinSelberg:
    def test_epsilon_matrix_N50(self):
        eps = rankin_selberg_epsilon(50)
        stored = RS_COUPLING["eps_matrix_N50"]
        for i in range(3):
            for j in range(3):
                assert eps[i][j] == pytest.approx(stored[i][j], abs=1e-5)

    def test_Q_values(self):
        q = rankin_selberg_Q(50)
        assert q["Q"] == pytest.approx(RS_COUPLING["Q_N50"], abs=1e-5)
        assert q["Q_corr"] == pytest.approx(RS_COUPLING["Q_corr_N50"], abs=1e-5)
        assert q["strict"] is True

    def test_Q_positive(self):
        """The Session-4 universal positivity: Q(143) > 0."""
        assert rankin_selberg_Q(50)["Q"] > 0


class TestCesaroIdentity:
    """B3 lock: (1/N) Σ w_n → −ζ′(2) = TORQUE_FULL (proved limit)."""

    def test_cesaro_converges_to_minus_zeta_prime_2(self):
        from mtft.constants import TORQUE_FULL
        N = 50000
        w = [0.0] * (N + 1)
        for d in range(2, N + 1):
            ld = math.log(d) / d
            for k in range(d, N + 1, d):
                w[k] += ld
        total = sum(w)
        assert abs(total / N - TORQUE_FULL) < 2e-3
