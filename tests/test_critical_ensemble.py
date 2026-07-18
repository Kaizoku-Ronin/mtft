#!/usr/bin/env python3
"""
Test Suite for critical_ensemble.py
====================================

Anchors:
  - lambda_1 exact closed form 1 + gamma/2 - (1/2)ln(4 pi)
  - series vs Cauchy coefficient agreement (independent methods)
  - Cauchy radius independence
  - zero-sum partials: nonnegative terms, monotone lower bounds,
    gap within tail model
  - self-computed lambda table snapshot (certified by three-leg agreement)
  - Bombieri-Lagarias caveat present in every report object

Roger Tano — MTFT Research Program — July 2026
"""

import sys
import os
import unittest


from mpmath import mp, mpf, log, pi, euler

from mtft.critical_ensemble import (
    lambda_1_closed_form, li_lambda, li_lambda_batch,
    logxi_taylor, logxi_taylor_cauchy, li_lambda_cauchy,
    li_lambda_zero_sum, certify, li_criterion_report,
    BOMBIERI_LAGARIAS_CAVEAT, XI_ANALYTICITY_RADIUS,
    THREE_ENSEMBLE_TABLE, LEVEL, GENUS,
)


class TestClosedFormAnchor(unittest.TestCase):
    """LEG A — lambda_1 against the exact closed form."""

    def test_lambda1_closed_form_value(self):
        """1 + gamma/2 - ln(4 pi)/2 = 0.0230957089661210..."""
        mp.dps = 30
        cf = lambda_1_closed_form()
        self.assertAlmostEqual(float(cf), 0.023095708966121034, places=16)

    def test_lambda1_series_vs_closed_form(self):
        """Series algebra reproduces the closed form to ~dps precision."""
        mp.dps = 30
        err = abs(li_lambda(1, dps=30) - lambda_1_closed_form())
        self.assertLess(float(err), 1e-25)

    def test_lambda1_cauchy_vs_closed_form(self):
        """Cauchy method independently reproduces the closed form."""
        mp.dps = 25
        err = abs(li_lambda_cauchy(1, r=2.0, dps=25) - lambda_1_closed_form())
        self.assertLess(float(err), 1e-20)


class TestMethodIndependence(unittest.TestCase):
    """LEG B — series algebra vs Cauchy integrals, coefficient level."""

    @classmethod
    def setUpClass(cls):
        mp.dps = 30
        cls.n = 6
        cls.a_series = logxi_taylor(cls.n, dps=30)
        cls.a_r15 = logxi_taylor_cauchy(cls.n, r=1.5, dps=30)
        cls.a_r25 = logxi_taylor_cauchy(cls.n, r=2.5, dps=30)

    def test_series_vs_cauchy(self):
        for k in range(1, self.n + 1):
            err = float(abs(self.a_series[k] - self.a_r15[k]))
            self.assertLess(err, 1e-30,
                            f"coefficient a_{k} disagrees: {err}")

    def test_radius_independence(self):
        for k in range(1, self.n + 1):
            err = float(abs(self.a_r15[k] - self.a_r25[k]))
            self.assertLess(err, 1e-30,
                            f"radius dependence at a_{k}: {err}")

    def test_bad_radius_rejected(self):
        with self.assertRaises(ValueError):
            logxi_taylor_cauchy(3, r=1.0)     # s = 0 on contour
        with self.assertRaises(ValueError):
            logxi_taylor_cauchy(3, r=3.5)     # trivial zero inside


class TestLambdaTable(unittest.TestCase):
    """Snapshot anchors, self-computed and certified by three-leg agreement
    (series/Cauchy coefficient agreement < 1e-30 at two radii; lambda_1
    matches its closed form to 1e-25; lambda_3 independently reproduced
    via Cauchy at r = 2.2 to full precision). Consistent with Keiper
    (1992) to the ~7 digits comparable from literature:
    0.0230957, 0.0923457, 0.2076389, 0.3687904, 0.5755427."""

    SNAPSHOT = {
        1: 0.023095708966121034,
        2: 0.09234573522804667,
        3: 0.2076389205543248,
        4: 0.36879047949224164,
        5: 0.57554271446117745,
        6: 0.8275660122823793,
        7: 1.1244601175709595,
        8: 1.4657556771470606,
    }

    def test_snapshot(self):
        mp.dps = 30
        vals = li_lambda_batch(8, dps=30)
        for n, expected in self.SNAPSHOT.items():
            self.assertAlmostEqual(float(vals[n - 1]), expected, places=12,
                                   msg=f"lambda_{n} drifted")

    def test_batch_matches_single(self):
        mp.dps = 25
        batch = li_lambda_batch(5, dps=25)
        for n in range(1, 6):
            err = float(abs(batch[n - 1] - li_lambda(n, dps=25)))
            self.assertLess(err, 1e-20)

    def test_positivity_prefix(self):
        """Finite prefix positive — expected, and logically empty (caveat)."""
        mp.dps = 25
        vals = li_lambda_batch(20, dps=25)
        self.assertTrue(all(float(v) > 0 for v in vals))

    def test_invalid_n(self):
        with self.assertRaises(ValueError):
            li_lambda(0)


class TestZeroSum(unittest.TestCase):
    """LEG C — zero-sum diagnostics."""

    @classmethod
    def setUpClass(cls):
        mp.dps = 20
        cls.results = {n: li_lambda_zero_sum(n, num_pairs=60, dps=20)
                       for n in (1, 2, 3)}

    def test_partial_nonnegative(self):
        for n, r in self.results.items():
            self.assertGreaterEqual(r.partial, 0)

    def test_partial_is_lower_bound(self):
        """|1 - 1/rho| = 1 on the line => terms >= 0 => partial <= analytic."""
        for n, r in self.results.items():
            self.assertLessEqual(r.partial, r.analytic + 1e-12,
                                 f"partial exceeds analytic at n={n}")

    def test_gap_within_tail_model(self):
        for n, r in self.results.items():
            self.assertLess(r.gap, 3 * r.tail_estimate + 1e-6,
                            f"gap {r.gap} vs tail {r.tail_estimate} at n={n}")

    def test_monotone_in_pairs(self):
        mp.dps = 20
        r30 = li_lambda_zero_sum(2, num_pairs=30, compare_analytic=False, dps=20)
        r60 = self.results[2]
        self.assertLessEqual(r30.partial, r60.partial + 1e-12)

    def test_caveat_attached(self):
        for r in self.results.values():
            self.assertEqual(r.caveat, BOMBIERI_LAGARIAS_CAVEAT)


class TestCertification(unittest.TestCase):
    """Full three-leg certification."""

    def test_certify_passes(self):
        res = certify(n_test=6, dps=25, zero_pairs=60)
        self.assertLess(res.lambda1_closed_form_err, 1e-19)
        self.assertLess(res.cauchy_agreement_err, 1e-20)
        self.assertLess(res.cauchy_radius_independence_err, 1e-20)
        self.assertTrue(res.zero_sum_bracketing_ok)
        self.assertTrue(res.passed)


class TestReport(unittest.TestCase):
    """The on-the-fly report object."""

    def test_report_structure(self):
        rep = li_criterion_report(n_max=10, dps=25)
        self.assertEqual(rep.n_max, 10)
        self.assertEqual(len(rep.values), 10)
        self.assertTrue(rep.all_positive)
        self.assertEqual(rep.min_at, 1)   # lambda_1 is the smallest
        self.assertEqual(rep.caveat, BOMBIERI_LAGARIAS_CAVEAT)

    def test_caveat_wording(self):
        self.assertIn("ALL n", BOMBIERI_LAGARIAS_CAVEAT)
        self.assertIn("Bombieri-Lagarias", BOMBIERI_LAGARIAS_CAVEAT)

    def test_three_ensemble_table_mentions_all(self):
        for word in ("Laplace", "Dirichlet", "Critical",
                     "Speiser", "Li", "Th 1"):
            self.assertIn(word, THREE_ENSEMBLE_TABLE)


class TestConstants(unittest.TestCase):
    def test_mtft_constants(self):
        self.assertEqual(LEVEL, 143)
        self.assertEqual(GENUS, 13)

    def test_analyticity_radius(self):
        """First zero at 1/2 + 14.1347i; |1 - rho_1| = sqrt(1/4 + gamma^2)."""
        import math
        d = math.sqrt(0.25 + 14.134725 ** 2)
        self.assertAlmostEqual(XI_ANALYTICITY_RADIUS, d, delta=0.02)


if __name__ == "__main__":
    unittest.main(verbosity=2)
