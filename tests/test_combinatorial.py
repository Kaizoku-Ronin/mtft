"""
tests/test_combinatorial.py
============================
Gate for the combinatorial ancestry module (2024 lineage).  Every
headline claim is certified by two routes sharing no steps (E2).
"""
import math
from fractions import Fraction as Fr

import numpy as np
import pytest

from mtft import combinatorial as C


# ── 1. Faulhaber engine ─────────────────────────────────────────────

class TestFaulhaber:
    @pytest.mark.parametrize("p", [1, 2, 3, 4, 5, 6, 7, 9, 11])
    def test_two_routes_agree(self, p):
        """Route A direct summation vs route B Bernoulli closed form,
        at p+3 points — pins the degree-(p+1) polynomial.  EXACT."""
        coeffs = C.faulhaber_coeffs(p)
        for n in range(0, p + 3):
            assert C._peval(coeffs, Fr(n)) == C.power_sum(p, n)

    def test_known_closed_forms(self):
        assert C.faulhaber_coeffs(1) == [Fr(0), Fr(1, 2), Fr(1, 2)]
        # S_2 = n(n+1)(2n+1)/6 = n/6 + n^2/2 + n^3/3
        assert C.faulhaber_coeffs(2) == [Fr(0), Fr(1, 6), Fr(1, 2), Fr(1, 3)]

    def test_C_formula_and_vanishing_corrections(self):
        assert C.C_figurate(1) == 1 and C.C_figurate(3) == 1
        assert C.C_figurate(5) == Fr(4, 3) and C.C_figurate(7) == 2
        # R_1 = R_3 = 0: the 2024 anchor cases
        assert not C.figurate_decomposition(1).R_T
        assert not C.figurate_decomposition(3).R_T


# ── 2. Sigma involution ─────────────────────────────────────────────

class TestSigmaInvolution:
    @pytest.mark.parametrize("p", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    def test_parity_two_routes(self, p):
        """Odd p sigma-even, even p sigma-odd; recurrence route vs
        coefficient-substitution route.  EXACT."""
        assert C.sigma_parity_check(p)

    def test_invariant_ring_and_odd_sector(self):
        # T itself round-trips; (2n+1)^2 = 8T + 1
        assert C.to_T_basis(C._T_POLY) == {1: Fr(1)}
        sq = C._pmul([Fr(1), Fr(2)], [Fr(1), Fr(2)])  # (2n+1)^2
        assert C.to_T_basis(sq) == {0: Fr(1), 1: Fr(8)}
        # S_2 = (2n+1)·T/3 ;  S_4 = (2n+1)·T(6T−1)/15
        assert C.odd_sector_factor(C.faulhaber_coeffs(2)) == {1: Fr(1, 3)}
        assert C.odd_sector_factor(C.faulhaber_coeffs(4)) == {1: Fr(-1, 15), 2: Fr(6, 15)}

    def test_non_sigma_even_rejected(self):
        with pytest.raises(ValueError):
            C.to_T_basis([Fr(0), Fr(1)])  # f(n) = n is not sigma-even


# ── 3. Figurate decomposition + jump rule ───────────────────────────

class TestFigurate:
    @pytest.mark.parametrize("p", [5, 7, 9, 11, 13])
    def test_degree_reduction(self, p):
        d = C.figurate_decomposition(p)
        assert d.deg_T_R <= d.q - 1                       # T-degree bound
        Rn = C.R_polynomial_coeffs(p)
        assert len(Rn) - 1 <= p - 1                       # n-degree bound

    def test_known_R5(self):
        # R_5 = T^2/3 (verified by hand in the ancestry review)
        assert C.figurate_decomposition(5).R_T == {2: Fr(1, 3)}

    @pytest.mark.parametrize("p", [5, 7, 9])
    def test_jump_rule_telescopes_to_closed_form(self, p):
        """Route A closed-form polynomial vs route B exact recurrence,
        plus the archived initial value R(1) = C(p) − 1.  EXACT."""
        Rn = C.R_polynomial_coeffs(p)
        rec = C.jump_rule_R(p, 30)
        assert rec[1] == C.C_figurate(p) - 1
        for m in range(31):
            assert rec[m] == C._peval(Rn, Fr(m))


# ── 4. s-gonal absorption defect (the honest negative) ──────────────

class TestSgonalDefect:
    @pytest.mark.parametrize("s", [3, 4, 5, 6, 8, 12])
    @pytest.mark.parametrize("p", [3, 5, 7, 9])
    def test_defect_closed_form(self, s, p):
        assert C.sgonal_defect_measured(s, p) == C.sgonal_defect(s)

    def test_triangular_unique_sigma_invariant(self):
        for s in [3, 4, 5, 6, 8, 12]:
            P = C.sgonal_coeffs(s)
            diff = C._padd(C.sigma_reflect(P), C._pscale(P, Fr(-1)))
            # P_s(−1−n) − P_s(n) = (2n+1)(s−3)
            assert diff == C._ptrim([Fr(s - 3), Fr(2 * (s - 3))])


# ── 5. Graph uncertainty ────────────────────────────────────────────

class TestGraphUncertainty:
    def _graphs(self):
        yield C.graph_cycle(7), True
        yield C.graph_complete(6), True
        yield C.graph_hypercube(3), True
        yield C.graph_path(6), False
        yield C.graph_star(7), False
        yield C.graph_erdos_renyi(12, 0.4, seed=143), None

    def test_commutator_two_routes(self):
        for A, _ in self._graphs():
            assert np.array_equal(C.commutator_DL(A), C.degree_gradient_matrix(A))

    def test_regularity_iff_zero(self):
        for A, regular in self._graphs():
            zero = not C.commutator_DL(A).any()
            if regular is not None:
                assert C.is_regular(A) is regular
            assert zero == C.is_regular(A)

    def test_support_is_degree_gradient_edges(self):
        A = C.graph_star(7)
        d = C.degrees(A)
        comm = C.commutator_DL(A)
        for u in range(7):
            for v in range(7):
                expect = bool(A[u, v]) and d[u] != d[v]
                assert bool(comm[u, v]) == expect

    def test_robertson_holds_on_complex_states(self):
        rng = np.random.default_rng(11 * 13)
        for A, _ in self._graphs():
            n = len(A)
            D = np.diag(C.degrees(A)).astype(float)
            L = D - A
            for _ in range(5):
                psi = rng.standard_normal(n) + 1j * rng.standard_normal(n)
                r = C.robertson_margin(D, L, psi)
                assert r["margin"] >= -1e-9


# ── 6. Number-phase and entropic uncertainty ────────────────────────

class TestNumberPhase:
    def test_robertson_random_states(self):
        rng = np.random.default_rng(168)
        for d in [4, 6, 9]:
            K = np.diag(np.arange(d, dtype=float))
            Th = C.phase_operator(d)
            assert np.allclose(Th, Th.conj().T)
            for _ in range(6):
                psi = rng.standard_normal(d) + 1j * rng.standard_normal(d)
                assert C.robertson_margin(K, Th, psi)["margin"] >= -1e-9

    def test_archived_2024_regression(self):
        """Delta_K = sqrt(5)/2 reproduces (EXACT).  The archived
        Theta-side triple reproduces under the ANGULAR convention
        (spectrum 2*pi*j/d), pinned by the exact 2*pi/6 factor the
        index-units pass missed by.  Cert(5e-4); the index rows keep
        the discovery route and must NOT match."""
        r = C.number_phase_regression()
        assert abs(r["angular"]["dX"] - math.sqrt(5) / 2) < 1e-12
        assert abs(C.ARCHIVED_2024["dK"] - math.sqrt(5) / 2) < 5e-4
        assert r["angular"]["matches_archived"]
        assert not r["index_plus"]["matches_archived"]
        for conv in ("index_plus", "index_minus", "angular"):
            assert r[conv]["margin"] >= -1e-9        # the inequality itself

    def test_maassen_uffink(self):
        rng = np.random.default_rng(29)
        for d in [3, 6, 8, 13]:
            # saturating state: a basis vector (H_K = 0, H_Theta = log d)
            e0 = np.zeros(d, dtype=complex); e0[0] = 1.0
            m = C.entropic_margin(e0)
            assert abs(m["margin"]) < 1e-9
            for _ in range(6):
                psi = rng.standard_normal(d) + 1j * rng.standard_normal(d)
                assert C.entropic_margin(psi)["margin"] >= -1e-9
        b = C.entropic_margin(C.binomial_state(5, 0.5).astype(complex))
        assert b["margin"] >= -1e-9


# ── 7. q-combinatorics and thermodynamics ───────────────────────────

class TestQThermo:
    def test_gaussian_binomial_two_routes(self):
        for n in range(0, 9):
            for k in range(0, n + 1):
                assert C.gaussian_binomial(n, k) == C.gaussian_binomial_product(n, k)

    def test_galois_anchors_q2(self):
        for n, g in C.GALOIS_ANCHORS_Q2.items():
            assert C.galois_number(n, 2) == g

    def test_thermo_identities_two_routes(self):
        E = np.log(np.arange(1, 40, dtype=float))          # small primon box
        for beta in (1.5, 2.5, 4.0):
            ens = C.Ensemble(energies=E)
            assert abs(ens.U(beta) - ens.U_fd(beta)) < 1e-6
            assert abs(ens.var_E(beta) - ens.var_fd(beta)) < 1e-4
            assert abs(ens.entropy(beta) - ens.entropy_identity(beta)) < 1e-9
            assert abs(ens.fisher_beta(beta) - ens.var_fd(beta)) < 1e-4
            assert ens.heat_capacity(beta) >= 0

    def test_multiplicity_energy_specials(self):
        # binomial multiplicities: Z(1) = 2^n exactly
        om = [math.comb(10, k) for k in range(11)]
        sp = C.multiplicity_Z_exact_specials(om)
        assert sp == {"Z0": 11, "Z1": 2 ** 10}
        # Gaussian multiplicities at q=2: Z(1) = Galois number
        om2 = [sum(c * 2 ** d for d, c in enumerate(C.gaussian_binomial(5, k)))
               for k in range(6)]
        assert C.multiplicity_Z_exact_specials(om2)["Z1"] == C.GALOIS_ANCHORS_Q2[5]
        # float Z at beta = 1 agrees with the exact integer
        ens = C.multiplicity_ensemble(om)
        assert abs(math.exp(ens.log_Z(1.0)) - 2 ** 10) < 1e-6


# ── 8. The bridge ───────────────────────────────────────────────────

class TestBridge:
    def test_sieve_matches_package_weights(self):
        """Route B for the sieve itself: mtft.arithmetic.weight_array."""
        mtft = pytest.importorskip("mtft")
        w = C.tano_weight_sieve(300)
        pkg = mtft.weight_array(300)
        pkg = np.asarray(pkg, dtype=float).ravel()
        mine = w[1:301] if len(pkg) == 300 else w[: len(pkg)]
        assert np.allclose(mine, pkg[: len(mine)], atol=1e-12)

    @pytest.mark.parametrize("s,N,cap", [(2.0, 200_000, 2e-3), (3.0, 60_000, 1e-6)])
    def test_weight_dirichlet_identity(self, s, N, cap):
        r = C.weight_dirichlet_identity_check(s, N)
        assert r["tail_bound"] < cap                       # bound is meaningful
        assert r["gap"] <= r["tail_bound"]                 # Cert(gap <= tail)

    def test_mean_weight_closed_form_and_endpoint(self):
        import mpmath as mp
        # <w>_2 = −zeta'(3): truncated mean within combined tails
        tr = C.mean_tano_weight(2.0, 200_000)
        ex = C.mean_tano_weight_exact(2.0)
        assert abs(ex - float(-mp.zeta(3, derivative=1))) < 1e-12
        assert abs(tr - ex) < 5e-3
        # beta -> 1+ endpoint: −zeta'(2) = cold-gas alpha = 2·T_INF
        end = C.mean_tano_weight_exact(1.0)
        assert abs(end - 0.9375482543158438) < 1e-12
        try:
            from mtft.constants import T_INF
            # bounds the STORED constant's precision (13 digits), not math
            assert abs(end - 2.0 * T_INF) < 1e-9
        except ImportError:
            pass

    def test_primon_fisher_is_logzeta_curvature(self):
        r = C.primon_fisher_check(beta=3.0, n_max=50_000)
        assert r["gap"] < 1e-5
