"""Tier 11 (JC certificate, estimator standards) + Legend tests."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import unittest
from fractions import Fraction as Fr

from mtft.jc_counterexample import (
    verify_all, verify_jacobian, verify_tautological_identity,
    verify_missed_curve, apply_F, COLLISION_FIBER, COLLISION_TARGET)
from mtft.estimator_standards import (
    binned_log_slope, stride_resonance_check)
from mtft.legend import (REGISTRY, legend_map, card, trace, status, search)


class TestJCCertificate(unittest.TestCase):
    def test_full_certificate(self):
        cert = verify_all()
        self.assertTrue(cert.valid)
        self.assertEqual(cert.component_degrees, (7, 6, 4))
    def test_jacobian_identity(self):
        self.assertTrue(verify_jacobian())
    def test_collision_exact(self):
        for p in COLLISION_FIBER:
            self.assertEqual(apply_F(p), COLLISION_TARGET)
    def test_tautological_identity(self):
        self.assertTrue(verify_tautological_identity())
    def test_missed_curve_empty(self):
        self.assertTrue(verify_missed_curve())
    def test_ag_d5_flag(self):
        self.assertIn("AG-D5", verify_all().notes["flag_AG_D5"])


class TestEstimatorStandards(unittest.TestCase):
    def test_gamma3_resonance_flagged(self):
        c, f, r = stride_resonance_check(25.0109, 6.0)
        self.assertTrue(r); self.assertLess(f, 0.01)
    def test_nonresonant_clean(self):
        _, _, r = stride_resonance_check(25.0109, 5.5)
        self.assertFalse(r)
    def test_binned_slope_recovers_powerlaw(self):
        import random; random.seed(143)
        ys = [10 ** (-6 + 5 * i / 400) for i in range(401)]
        vals = [y ** (-0.25) * (1 + 0.05 * random.uniform(-1, 1)) for y in ys]
        s, used, dropped = binned_log_slope(ys, vals, min_bin=10)
        self.assertAlmostEqual(s, -0.25, delta=0.02)
        self.assertGreaterEqual(used, 3)
    def test_min_bin_guard(self):
        with self.assertRaises(ValueError):
            binned_log_slope([1.0, 2.0], [1.0, 1.0])


class TestLegend(unittest.TestCase):
    def test_registry_integrity(self):
        for name, e in REGISTRY.items():
            self.assertEqual(name, e.name)
            self.assertIn(e.tag, ("Df","Pp","Pr","Conj","Heur","Cert"))
            for up in e.upstream:
                self.assertIn(up, REGISTRY)
    def test_all_chains_reach_integers(self):
        def reaches(n, path=()):
            if n == "integers": return True
            if n in path: return False
            e = REGISTRY[n]
            return any(reaches(u, path+(n,)) for u in e.upstream) if e.upstream else False
        for name in REGISTRY:
            if name != "integers":
                self.assertTrue(reaches(name), f"{name} !-> N")
    def test_trace_terminates(self):
        self.assertIn("always the integers", trace("alpha_inverse"))
    def test_map_covers_tier11(self):
        m = legend_map()
        self.assertIn("jc_counterexample", m)
        self.assertIn("Certificates & Standards", m)
    def test_card_and_search(self):
        self.assertIn("Speiser", card("dirichlet_curvature"))
        self.assertIn("li_lambda", search("Li coefficients"))


if __name__ == "__main__":
    unittest.main(verbosity=1)
