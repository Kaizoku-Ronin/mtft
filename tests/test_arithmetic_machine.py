#!/usr/bin/env python3
"""
Test Suite for arithmetic_machine.py
=====================================

Verifies the five-primitive decomposition of computation,
configuration space geometry, computational stiffness,
primitive complexity classes, and halting surface topology.

Roger Tano — MTFT Research Program — April 2026
"""

import math
import unittest

from mtft.arithmetic_machine import (
    # Constants
    LEVEL, GENUS, INDEX, DIM_NEW, ORBIT_DIMS, CANONICAL_DEG,
    MONSTER_DIM, HECKE_TRACES,
    # Primitives
    Primitive, PrimitiveLevel, PrimitiveDecomposition,
    decompose_turing_machine, decompose_lambda_calculus,
    decompose_recursive_function,
    # Configuration space
    ConfigPoint, config_space_size, hamming_distance,
    modular_distance, hecke_weighted_norm,
    # Stiffness
    StiffnessResult, hecke_sign, hecke_constraint_index,
    search_space_compression, computational_stiffness,
    # Complexity
    PrimitiveClassification, PRIMITIVE_CLASSIFICATIONS,
    classify_computation, level_hierarchy,
    # Halting surface
    HaltingSurface, analyze_halting_surface,
    # Entropy
    ArithmeticEntropy, arithmetic_entropy,
    # Bridge
    ComputationPhysicsBridge, COMPUTATION_PHYSICS_DICTIONARY,
    bridge_lookup,
    # Verification
    verify_config_space_identity, verify_search_compression_bounds,
    verify_primitive_hierarchy, verify_stiffness_nonnegativity,
)


class TestMTFTConstants(unittest.TestCase):
    """Verify all MTFT structural constants."""

    def test_level_factorization(self):
        """N = 143 = 11 × 13."""
        self.assertEqual(LEVEL, 143)
        self.assertEqual(LEVEL, 11 * 13)

    def test_genus(self):
        """genus(X₀(143)) = 13."""
        self.assertEqual(GENUS, 13)

    def test_index(self):
        """[SL(2,ℤ) : Γ₀(143)] = 168 = |PSL(2,7)|."""
        self.assertEqual(INDEX, 168)

    def test_dim_new(self):
        """dim S₂ᶰᵉʷ(Γ₀(143)) = 11."""
        self.assertEqual(DIM_NEW, 11)

    def test_orbit_dimensions(self):
        """Galois orbits: [1, 4, 6] summing to 11."""
        self.assertEqual(ORBIT_DIMS, (1, 4, 6))
        self.assertEqual(sum(ORBIT_DIMS), DIM_NEW)

    def test_canonical_degree(self):
        """2g - 2 = 24."""
        self.assertEqual(CANONICAL_DEG, 2 * GENUS - 2)
        self.assertEqual(CANONICAL_DEG, 24)

    def test_monster_dimension(self):
        """dim V♮ = 196883."""
        self.assertEqual(MONSTER_DIM, 196_883)

    def test_hecke_traces_length(self):
        """200 Hecke traces precomputed."""
        self.assertEqual(len(HECKE_TRACES), 200)

    def test_hecke_trace_1(self):
        """a_1 = dim S₂ᶰᵉʷ = 11 (trace of identity)."""
        self.assertEqual(HECKE_TRACES[0], 11)


class TestPrimitiveDecomposition(unittest.TestCase):
    """Test the five-primitive decomposition of computation."""

    def test_turing_machine_decomposition(self):
        """TM decomposes into all five primitives."""
        decomp = decompose_turing_machine()
        self.assertIsInstance(decomp, PrimitiveDecomposition)
        self.assertEqual(decomp.level, PrimitiveLevel.GEOMETRIC)
        # Check that iterate component mentions iteration/repeated
        self.assertTrue(
            "iteration" in decomp.iterate_component.lower()
            or "repeatedly" in decomp.iterate_component.lower()
            or "δ" in decomp.iterate_component
        )

    def test_lambda_calculus_decomposition(self):
        """Lambda calculus decomposes into all five primitives."""
        decomp = decompose_lambda_calculus()
        self.assertEqual(decomp.level, PrimitiveLevel.GEOMETRIC)
        self.assertIn("EXTRACT", decomp.description.upper())

    def test_recursive_function_decomposition(self):
        """Recursive functions decompose into all five primitives."""
        decomp = decompose_recursive_function()
        self.assertEqual(decomp.level, PrimitiveLevel.GEOMETRIC)
        self.assertIn("DIVIDE", decomp.description.upper())

    def test_three_formulations_same_level(self):
        """Church-Turing: all three are Level 4 (GEOMETRIC)."""
        tm = decompose_turing_machine()
        lc = decompose_lambda_calculus()
        rf = decompose_recursive_function()
        self.assertEqual(tm.level, lc.level)
        self.assertEqual(lc.level, rf.level)

    def test_primitive_enum_values(self):
        """Primitives numbered 1-5 matching AG convention."""
        self.assertEqual(Primitive.ITERATE, 1)
        self.assertEqual(Primitive.DIVIDE, 2)
        self.assertEqual(Primitive.ASSEMBLE, 3)
        self.assertEqual(Primitive.EXTRACT, 4)
        self.assertEqual(Primitive.CURVE, 5)

    def test_level_enum_values(self):
        """Levels numbered 0-4."""
        self.assertEqual(PrimitiveLevel.PRIM_REC, 0)
        self.assertEqual(PrimitiveLevel.GENERAL_REC, 1)
        self.assertEqual(PrimitiveLevel.ANALYTIC, 2)
        self.assertEqual(PrimitiveLevel.INVERSIVE, 3)
        self.assertEqual(PrimitiveLevel.GEOMETRIC, 4)


class TestConfigurationSpace(unittest.TestCase):
    """Test configuration space geometry."""

    def test_config_space_size_identity(self):
        """|C_13(1)| = 2^13 × 13 = 106496."""
        self.assertEqual(config_space_size(1, 13), 106496)
        self.assertEqual(config_space_size(1, 13), (2**13) * 13)

    def test_config_space_linearity(self):
        """|C_g(n)| = n × |C_g(1)|."""
        base = config_space_size(1, GENUS)
        for n in range(1, 5):
            self.assertEqual(config_space_size(n, GENUS), n * base)

    def test_hamming_distance_self(self):
        """d_H(c, c) = 0."""
        c = ConfigPoint(0, 3, (0, 1, 0, 1, 0))
        self.assertEqual(hamming_distance(c, c), 0)

    def test_hamming_distance_symmetry(self):
        """d_H(c1, c2) = d_H(c2, c1)."""
        c1 = ConfigPoint(0, 3, (0, 1, 0, 1, 0))
        c2 = ConfigPoint(1, 3, (0, 0, 0, 1, 1))
        self.assertEqual(hamming_distance(c1, c2), hamming_distance(c2, c1))

    def test_hamming_distance_computation(self):
        """d_H counts differing components."""
        c1 = ConfigPoint(0, 0, (0, 0, 0))
        c2 = ConfigPoint(1, 1, (1, 0, 1))
        # state diff (1) + head diff (1) + tape diff (2) = 4
        self.assertEqual(hamming_distance(c1, c2), 4)

    def test_modular_distance_nonneg(self):
        """d_M ≥ 0 always."""
        c1 = ConfigPoint(0, 0, (0, 0, 0))
        c2 = ConfigPoint(1, 1, (1, 1, 1))
        self.assertGreaterEqual(modular_distance(c1, c2, 2, 3), 0)

    def test_modular_distance_self_zero(self):
        """d_M(c, c) = 0."""
        c = ConfigPoint(0, 2, (1, 0, 1, 0))
        self.assertEqual(modular_distance(c, c, 1, 4), 0)

    def test_hecke_weighted_norm(self):
        """Hecke norm weights tape by |a_n|."""
        tape = (1, 0, 0, 0, 0)  # only position 0 is 1
        # |a_1| = 11
        self.assertEqual(hecke_weighted_norm(tape), 11)

        tape2 = (1, 1, 0, 0, 0)  # positions 0 and 1
        # |a_1| + |a_2| = 11 + 3 = 14
        self.assertEqual(hecke_weighted_norm(tape2), 14)

    def test_config_point_index(self):
        """ConfigPoint indexing is injective."""
        indices = set()
        for s in range(2):
            for h in range(3):
                for bits in range(8):
                    tape = tuple((bits >> i) & 1 for i in range(3))
                    c = ConfigPoint(s, h, tape)
                    idx = c.to_index(2, 3)
                    self.assertNotIn(idx, indices)
                    indices.add(idx)


class TestHeckeConstraints(unittest.TestCase):
    """Test Hecke sign oracle and constraint functions."""

    def test_hecke_sign_positive(self):
        """a_1 = 11 > 0 → bosonic (+1)."""
        self.assertEqual(hecke_sign(1), 1)

    def test_hecke_sign_negative(self):
        """a_6 = -4 < 0 → fermionic (-1)."""
        self.assertEqual(hecke_sign(6), -1)

    def test_hecke_sign_zero(self):
        """a_5 = 0 → free (0)."""
        self.assertEqual(hecke_sign(5), 0)

    def test_hecke_sign_out_of_range(self):
        """Beyond table: return 0 (unconstrained)."""
        self.assertEqual(hecke_sign(201), 0)
        self.assertEqual(hecke_sign(0), 0)

    def test_constraint_index_mapping(self):
        """(state, symbol) → unique Hecke index."""
        self.assertEqual(hecke_constraint_index(0, 0), 1)
        self.assertEqual(hecke_constraint_index(0, 1), 2)
        self.assertEqual(hecke_constraint_index(1, 0), 3)
        self.assertEqual(hecke_constraint_index(1, 1), 4)

    def test_search_compression_bounds(self):
        """0 < compression ≤ 1 for all n."""
        for n in range(1, 6):
            c = search_space_compression(n)
            self.assertGreater(c, 0)
            self.assertLessEqual(c, 1)

    def test_search_compression_n1(self):
        """n=1: a_1=11>0, a_2=3>0, both bosonic → compression = 1/4."""
        c = search_space_compression(1)
        # Both slots constrained: each loses factor of 2 in write choices
        # 2 slots × half = (1/2)^2 = 1/4
        self.assertAlmostEqual(c, 0.25)


class TestComputationalStiffness(unittest.TestCase):
    """Test the computational stiffness function μ_C."""

    def test_stiffness_n1(self):
        """μ_C(1) is computable and non-negative."""
        result = computational_stiffness(1, max_steps=5000, tape_len=5)
        self.assertGreaterEqual(result.stiffness, 0)
        self.assertGreater(result.computation_time, 0)

    def test_stiffness_result_structure(self):
        """StiffnessResult has all required fields."""
        result = computational_stiffness(1, max_steps=1000, tape_len=5)
        self.assertEqual(result.n_states, 1)
        self.assertIsInstance(result.bb_unconstrained, int)
        self.assertIsInstance(result.bb_hecke, int)
        self.assertIsInstance(result.compression_ratio, float)
        self.assertIsInstance(result.stiffness, float)

    def test_bb_hecke_leq_bb_unc(self):
        """BB_Hecke ≤ BB_unconstrained always."""
        result = computational_stiffness(1, max_steps=5000, tape_len=5)
        if result.bb_unconstrained >= 0 and result.bb_hecke >= 0:
            self.assertLessEqual(result.bb_hecke, result.bb_unconstrained)

    def test_search_space_smaller(self):
        """Hecke search space ≤ unconstrained search space."""
        result = computational_stiffness(1, max_steps=5000, tape_len=5)
        self.assertLessEqual(result.search_compression, 1.0)


class TestPrimitiveComplexity(unittest.TestCase):
    """Test primitive complexity classification."""

    def test_all_levels_populated(self):
        """Each level has at least one classified computation."""
        hier = level_hierarchy()
        for level in PrimitiveLevel:
            self.assertIn(level, hier,
                          f"Level {level.name} has no classified computations")

    def test_level_ordering(self):
        """Lower-level computations use fewer primitives."""
        for name, cls in PRIMITIVE_CLASSIFICATIONS.items():
            n_prims = len(cls.primitives_used)
            self.assertLessEqual(n_prims, cls.level.value + 1)

    def test_classify_known_computation(self):
        """Can classify named computations."""
        c = classify_computation("addition")
        self.assertIsNotNone(c)
        self.assertEqual(c.level, PrimitiveLevel.PRIM_REC)

    def test_classify_level4(self):
        """Coupling constants require all five primitives."""
        c = classify_computation("coupling_constants")
        self.assertIsNotNone(c)
        self.assertEqual(c.level, PrimitiveLevel.GEOMETRIC)
        self.assertEqual(len(c.primitives_used), 5)

    def test_mass_gap_is_level4(self):
        """Mass gap is a Level 4 computation."""
        c = classify_computation("mass_gap")
        self.assertIsNotNone(c)
        self.assertEqual(c.level, PrimitiveLevel.GEOMETRIC)

    def test_monster_hash_is_level3(self):
        """MonsterHash is Level 3 (INVERSIVE)."""
        c = classify_computation("monster_hash")
        self.assertIsNotNone(c)
        self.assertEqual(c.level, PrimitiveLevel.INVERSIVE)

    def test_halting_problem_level1(self):
        """Halting problem appears at Level 1."""
        c = classify_computation("halting_problem")
        self.assertIsNotNone(c)
        self.assertEqual(c.level, PrimitiveLevel.GENERAL_REC)


class TestHaltingSurface(unittest.TestCase):
    """Test halting surface topology analysis."""

    def test_halting_surface_structure(self):
        """HaltingSurface has required fields."""
        hs = analyze_halting_surface(1, tape_len=5)
        self.assertEqual(hs.n_states, 1)
        self.assertEqual(hs.tape_len, 5)
        self.assertGreater(hs.total_configs, 0)
        self.assertGreaterEqual(hs.halting_configs, 0)

    def test_boundary_fraction_bounds(self):
        """0 ≤ boundary fraction ≤ 1."""
        hs = analyze_halting_surface(1, tape_len=5)
        self.assertGreaterEqual(hs.boundary_fraction, 0)
        self.assertLessEqual(hs.boundary_fraction, 1)

    def test_genus_bound(self):
        """Genus bound ≤ GENUS = 13."""
        hs = analyze_halting_surface(2, tape_len=5)
        self.assertLessEqual(hs.genus_bound, GENUS)

    def test_halting_configs_leq_total(self):
        """Halting configs ≤ total configs."""
        hs = analyze_halting_surface(1, tape_len=5)
        self.assertLessEqual(hs.halting_configs, hs.total_configs)


class TestArithmeticEntropy(unittest.TestCase):
    """Test arithmetic entropy computation."""

    def test_entropy_nonneg(self):
        """Shannon entropy ≥ 0."""
        # Build a simple halting machine
        transitions = {
            (0, 0): (1, 1, -1),  # write 1, move right, halt
            (0, 1): (0, 0, -1),  # write 0, move left, halt
        }
        ent = arithmetic_entropy(transitions, 1, tape_len=5)
        self.assertGreaterEqual(ent.shannon_entropy, 0)

    def test_hecke_entropy_nonneg(self):
        """Hecke entropy ≥ 0."""
        transitions = {
            (0, 0): (1, 1, -1),
            (0, 1): (0, 0, -1),
        }
        ent = arithmetic_entropy(transitions, 1, tape_len=5)
        self.assertGreaterEqual(ent.hecke_entropy, 0)

    def test_orbit_length_positive(self):
        """Orbit length > 0 for any machine."""
        transitions = {
            (0, 0): (1, 1, -1),
            (0, 1): (0, 0, -1),
        }
        ent = arithmetic_entropy(transitions, 1, tape_len=5)
        self.assertGreater(ent.orbit_length, 0)


class TestComputationPhysicsBridge(unittest.TestCase):
    """Test the computation-physics bridge dictionary."""

    def test_bridge_entries_exist(self):
        """All bridge entries are populated."""
        self.assertGreater(len(COMPUTATION_PHYSICS_DICTIONARY), 0)

    def test_bridge_lookup(self):
        """Can look up bridge entries."""
        b = bridge_lookup("halting")
        self.assertIsNotNone(b)
        self.assertIn("confinement", b.physics.lower())

    def test_bridge_structure(self):
        """All bridge entries have required fields."""
        for key, bridge in COMPUTATION_PHYSICS_DICTIONARY.items():
            self.assertIsInstance(bridge, ComputationPhysicsBridge)
            self.assertIsInstance(bridge.computation, str)
            self.assertIsInstance(bridge.physics, str)
            self.assertIsInstance(bridge.bridge_mechanism, str)


class TestVerification(unittest.TestCase):
    """Test all verification invariants."""

    def test_config_space_identity(self):
        """|C_13(1)| = 2^genus × genus."""
        self.assertTrue(verify_config_space_identity())

    def test_search_compression_bounds(self):
        """Search compression in (0, 1] for all tested n."""
        self.assertTrue(verify_search_compression_bounds())

    def test_primitive_hierarchy(self):
        """Primitive hierarchy is well-formed."""
        self.assertTrue(verify_primitive_hierarchy())

    def test_stiffness_nonnegativity(self):
        """μ_C(n) ≥ 0 for tested n."""
        self.assertTrue(verify_stiffness_nonnegativity(max_n=1))


if __name__ == "__main__":
    unittest.main(verbosity=2)
