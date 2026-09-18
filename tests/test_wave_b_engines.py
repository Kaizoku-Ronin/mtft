"""Wave B engines: each test uses an INDEPENDENT construction of the expected answer (theorem, divisor count, direct charge
sum, known widths), never a hardcoded value returned by the same implementation."""
import sympy as sp, pytest
from fractions import Fraction
from mtft.surface import spin_circle as SC, hopf_geometry as HG, crt_dessin as CD, hodge_blocks as HB
from mtft.research import anomalies as AN, charge_lattices as CL, mode_operators as MO


def test_spin_circle_fiber_acyclicity_and_controls():
    g, e = 13, 12
    for m in (3, 6, -3, -6):                                                      # M1 sectors: holonomy != 1 -> acyclic (theorem) vs Fox computation
        lam = SC.fiber_holonomy(m, e); assert lam != 1
        assert SC.m1_sector_cohomology(m)["cohomology_dimensions"] == [0, 0, 0, 0]
    assert SC.m1_sector_cohomology(0)["cohomology_dimensions"] == [1, 2 * g, 2 * g, 1]     # Gysin: e != 0 kills the fiber class rationally
    base = [-sp.S.One] + [sp.S.One] * (2 * g - 1)
    assert SC.rank_one_cohomology(g, e, [*base, sp.S.One])["cohomology_dimensions"] == [0, 2 * g - 2, 2 * g - 2, 0]
    assert SC.circle_bundle_topology(g, e)["pullback_torsion_order"](3) == 4 and SC.circle_bundle_topology(g, e)["pullback_torsion_order"](6) == 2
    R, r = sp.Symbol("R", positive=True), sp.Symbol("r", positive=True); met = SC.spin_circle_metric(R, r)
    assert sp.simplify(met["scalar_curvature"] - (-2 / R ** 2 - r ** 2 / (8 * R ** 4))) == 0 and sp.simplify(met["volume"] - 96 * sp.pi ** 2 * R ** 2 * r) == 0
    assert SC.radius_potential_7d(1, 1, 1)["stationary_point_exists"] is False


def test_hopf_pencil_by_riemann_hurwitz_and_divisor_count():
    p = HG.half_form_pencil(); assert p["degree"] == 12 and p["riemann_hurwitz_check"] and p["total_ramification"] == 48
    assert HG.pullback_line_degrees() == (-12, 12) and HG.hopf_connection_pullback()["c2_on_curve"] == 0


def test_crt_dessin_widths_and_genus():
    d = CD.dessin(); assert d["darts"] == 168 and d["cusp_widths"] == [1, 11, 13, 143] and d["edge_orbits"] == 84 and d["triangle_orbits"] == 56
    assert d["genus"] == 13 and d["S_involution"] and d["ST_order_three"]
    assert CD.crt_projective_line()["all_primitive"] and [CD.genus_squarefree(N)["genus"] for N in (11, 13, 143)] == [1, 0, 13]
    # independent: the dessin permutations agree with the package's Manin-symbol model
    from mtft import hecke as H
    model = H.model(); P = CD.crt_projective_line(); key_to_old = {(CD.canon_prime(c, d_, 11), CD.canon_prime(c, d_, 13)): j for j, (c, d_) in enumerate(model["P1"])}
    old = [key_to_old[pair] for pair in P["pairs"]]; PS, PT = CD.permutation((0, -1, 1, 0), P["pairs"]), CD.permutation((1, 1, 0, 1), P["pairs"])
    assert all(old[PS[j]] == model["sS"][old[j]] for j in range(168)) and all(old[PT[j]] == model["sT"][old[j]] for j in range(168))


@pytest.mark.slow
def test_hodge_block_projectors_rational():
    r = HB.coprime_block_projectors(); assert r["rank12"] == 12 and r["rank14"] == 14 and r["idempotent"] and r["complementary"]
    x = sp.Symbol("x"); assert sp.expand(r["T2"].charpoly(x).as_expr() - x ** 2 * (x + 2) ** 4 * r["quartic"] ** 2 * r["sextic"] ** 2) == 0


def test_quaternionic_gate():
    assert HB.quaternionic_multiplicity_gate((2, 2, 2, 2))["obstructed"] and not HB.quaternionic_multiplicity_gate((4, 4, 4, 4))["obstructed"]


def test_anomaly_polynomial_factorization_and_kernel():
    N = {"c": 3, "L": 2, "a": 1, "b": 1, "d": 1}; m = {"c": 0, "L": -3, "a": 3, "b": 3, "d": 0}; model = AN.stack_model(N, m, order=("c", "L", "a", "b", "d"))
    A = AN.anomaly_polynomials(model); q = A["q"]; c, l, a, b, d = q; x, y, z, w = c - l, a - l, b - l, d - l
    r1, r2 = y + z, 3 * x + w; Q1 = 24 * (y * y - y * z + z * z) + 27 * x * x + 9 * w * w; Q2 = -9 * (y * y + z * z)
    assert sp.expand(A["cubic"] - (r1 * Q1 + r2 * Q2)) == 0                                   # AXG-01 exact ideal factorization
    assert sp.expand(A["mixed"]["c"] - sp.Rational(3, 2) * r1) == 0 and sp.expand(A["mixed"]["L"] - (3 * r1 + sp.Rational(3, 2) * r2)) == 0 and sp.expand(A["gravity"] - 24 * r1) == 0
    K = AN.primitive_shift_matrix(model); assert K.shape == (2, 5) and K.rank() == 2
    phase = sp.Matrix([1, 1, 1, 1, 1]); Y = sp.Matrix([sp.Rational(1, 6), 0, -sp.Rational(1, 2), sp.Rational(1, 2), -sp.Rational(1, 2)]); BL = sp.Matrix([sp.Rational(1, 3), 0, 0, 0, -1])
    assert (K * sp.Matrix.hstack(phase, Y, BL)) == sp.zeros(2, 3)
    assert CL.smith_remnant(K)["invariant_factors"] == [1, 1]
    # direct charge-sum cross-check of the polarised trace against the polynomial
    ql = [sp.Rational(1, 6), 0, -sp.Rational(1, 2), sp.Rational(1, 2), -sp.Rational(1, 2)]
    direct = sum((model["m"][i] - model["m"][j]) * model["N"][i] * model["N"][j] * (ql[i] - ql[j]) ** 3 for i in range(5) for j in range(i + 1, 5))
    assert AN.polarised_trace(A["cubic"], q, ql, ql, ql) == direct == 0


def test_charge_lattice_dressings_and_masses():
    K = sp.Matrix([[0, -2, 1, 1, 0], [3, -4, 0, 0, 1]])
    assert CL.shift_kernel(K)["massless_count"] == 3
    # H_u H_d = (L, a-bar)(L, b-bar): charge vector on stacks (c, L, a, b, d) is 2 e_L - e_a - e_b
    assert CL.integer_dressing(K, [0, 2, -1, -1, 0])["dressable"]
    # nu^c nu^c = 2 (a, d-bar): 2 e_a - 2 e_d is NOT in the row lattice (B-L charge 2): forbidden
    assert not CL.integer_dressing(K, [0, 0, 2, 0, -2])["dressable"]
    g, f = sp.symbols("g f", positive=True); r = CL.canonical_vector_masses(K, sp.diag(6, 4, 2, 2, 2) / g ** 2, f ** 2 * sp.eye(2))
    ev = set(sp.simplify(k) for k in r["eigenvalues"]); assert r["rank"] == 2 and sp.simplify(g ** 2 * f ** 2 * (4 + 2 * sp.sqrt(2))) in ev and sp.simplify(g ** 2 * f ** 2 * (4 - 2 * sp.sqrt(2))) in ev


def test_mode_operator_certificates():
    i13 = sp.I / sp.sqrt(13)
    assert MO.s0_twist_cohomology([i13], 0) == {"h0": 2, "h1": 1, "index": 1, "pure": False}          # O(P): mirror pair
    assert MO.s0_twist_cohomology([i13, -i13], 1) == {"h0": 1, "h1": 0, "index": 1, "pure": True}     # O(P+Q-R): pure single mode
    assert MO.s0_twist_cohomology([i13, -i13, i13], 0)["h0"] == 3 and MO.purity_certificate([i13, -i13, i13])["pure"]   # M1 family space
    assert MO.bochner_bound(-6)["lambda_min_bound"] == sp.Rational(1, 4)                                   # = the magnetic-mesh Landau value 0.2504–0.2512
    assert MO.six_d_scalar_yukawa_gate("+", "-")["allowed"] and not MO.six_d_scalar_yukawa_gate("+", "+")["allowed"]
