"""Test suite for mtft v0.2 — arithmetic-first constants + lattice."""
import math, cmath
import numpy as np
import pytest

# ── MTFT Arithmetic Constants ────────────────────────────────

def test_alpha_inverse():
    from mtft.constants import GAUGE, PDG
    assert abs(GAUGE.alpha_inv - PDG.alpha_inv) / PDG.alpha_inv < 1e-7

def test_weinberg_angle():
    from mtft.constants import GAUGE
    assert GAUGE.sin2_theta_W == pytest.approx(3/13)
    assert GAUGE.cos2_theta_W == pytest.approx(10/13)

def test_strong_coupling():
    from mtft.constants import GAUGE, PDG
    assert abs(GAUGE.alpha_s - PDG.alpha_s_mZ) / PDG.alpha_s_mZ < 0.01

def test_W_Z_ratio():
    from mtft.constants import GAUGE, PDG
    obs = PDG.m_W / PDG.m_Z
    assert abs(GAUGE.W_Z_ratio - obs) / obs < 0.005

def test_higgs_mass():
    from mtft.constants import HIGGS, PDG
    assert abs(HIGGS.m_H - PDG.m_H) / PDG.m_H < 0.005

def test_sm_masses():
    from mtft.constants import SM, PDG
    assert abs(SM.m_W - PDG.m_W) / PDG.m_W < 0.04  # MTFT uses sin²θ=3/13
    assert abs(SM.m_Z - PDG.m_Z) / PDG.m_Z < 0.04

def test_critical_depth_from_arithmetic():
    from mtft.constants import CriticalDepths
    y_arith = CriticalDepths.y_conf_from_arithmetic()
    assert abs(y_arith - CriticalDepths.y_conf) / CriticalDepths.y_conf < 0.01

# ── Arithmetic Weights ───────────────────────────────────────

def test_weight_values():
    from mtft.arithmetic import weight
    assert weight(1) == 0.0
    assert abs(weight(2) - math.log(2)/2) < 1e-12
    assert abs(weight(3) - math.log(3)/3) < 1e-12
    w6 = math.log(2)/2 + math.log(3)/3 + math.log(6)/6
    assert abs(weight(6) - w6) < 1e-12

def test_torque_convergence():
    from mtft.arithmetic import torque_partial_sum
    from mtft.constants import T_INF
    T_approx = torque_partial_sum(10000)
    assert abs(T_approx - T_INF) / T_INF < 0.002

def test_S1_at_spectral_depth():
    from mtft.arithmetic import S1
    from mtft.constants import GAUGE, CriticalDepths
    s1 = S1(CriticalDepths.y_spec)
    alpha_inv = 2.0 / s1
    assert abs(alpha_inv - GAUGE.alpha_inv) / GAUGE.alpha_inv < 0.01

def test_mass_gap_positive():
    from mtft.arithmetic import mass_gap_stiffness
    for N in [2, 3, 4, 5]:
        mu = mass_gap_stiffness(0.15, N=N)
        assert mu > 0, f"μ_{N}(0.15) must be > 0"

def test_su3_hessian_isotropy():
    from mtft.arithmetic import su3_hessian_eigenvalues
    lam1, lam2 = su3_hessian_eigenvalues(0.18174)
    ratio = lam1 / lam2
    assert abs(ratio - 1.0) < 0.01  # eigenvalues nearly equal at y_conf

def test_confinement_depth():
    from mtft.arithmetic import find_confinement_depth
    y_c = find_confinement_depth()
    assert abs(y_c - 0.1817) < 0.002

# ── Modular Geometry ─────────────────────────────────────────

def test_sl2z_transform():
    from mtft.modular import sl2z_transform, S_MATRIX, T_MATRIX
    tau = 0.3 + 1.2j
    assert abs(sl2z_transform(tau, S_MATRIX) - (-1/tau)) < 1e-12
    assert abs(sl2z_transform(tau, T_MATRIX) - (tau + 1)) < 1e-12

def test_fundamental_domain():
    from mtft.modular import reduce_to_fundamental
    tau_red, _ = reduce_to_fundamental(3.7 + 0.4j)
    assert abs(tau_red) >= 1.0 - 1e-8
    assert abs(tau_red.real) <= 0.5 + 1e-8

# ── Modular Forms ────────────────────────────────────────────

def test_eta_modular_property():
    from mtft.forms import dedekind_eta
    tau = 0.2 + 1.0j
    ratio = dedekind_eta(tau + 1) / (cmath.exp(1j*math.pi/12) * dedekind_eta(tau))
    assert abs(ratio - 1.0) < 1e-8

def test_j_invariant_at_i():
    from mtft.forms import j_invariant
    j_i = j_invariant(1j)
    assert abs(j_i.real - 1728) < 1.0

# ── Hosotani ─────────────────────────────────────────────────

def test_hosotani_vacuum():
    from mtft.hosotani import HosotaniPotential
    hp = HosotaniPotential(fermion_fraction=0.4)
    theta0 = hp.find_vacuum()
    assert 0 <= theta0 <= math.pi
    assert hp.second_derivative(theta0) > -1e-6  # at or near minimum

# ── Particles ────────────────────────────────────────────────

def test_particle_catalog():
    from mtft.particles import StandardModel
    sm = StandardModel()
    assert len(sm.all) >= 17
    assert sm.by_name("Top").mass_GeV == pytest.approx(172.69, abs=0.1)

def test_kappa_hierarchy():
    from mtft.particles import StandardModel
    kappas = [k for _, k in StandardModel().kappa_hierarchy()]
    assert kappas == sorted(kappas)

# ── Dark Sector ──────────────────────────────────────────────

def test_flat_rotation():
    from mtft.dark_sector import TauVortexHalo
    halo = TauVortexHalo(A=1e20, r0=1e30)
    r = np.array([1e33, 5e33, 1e34, 5e34])
    v = halo.v_circular(r)
    assert np.all(np.abs(v - v.mean()) / v.mean() < 0.15)

def test_v_infinity():
    from mtft.dark_sector import TauVortexHalo
    from mtft.constants import PhysicalConstants as PC
    A = 1e20
    halo = TauVortexHalo(A=A)
    expected = math.sqrt(2*math.pi * PC.G_N * A**2)
    assert abs(halo.v_infinity - expected) / expected < 1e-10

# ── Information Geometry ─────────────────────────────────────

def test_lyapunov():
    from mtft.info_geometry import lyapunov_exponent
    assert lyapunov_exponent(4.0) > 0
    assert lyapunov_exponent(2.5) < 0

def test_fisher_metric_positive():
    from mtft.info_geometry import fisher_rao_metric
    assert fisher_rao_metric(3.8) > 0

# ── Cosmology ────────────────────────────────────────────────

def test_friedmann():
    from mtft.cosmology import FriedmannMTFT
    cosmo = FriedmannMTFT()
    assert abs(cosmo.H_squared(1.0) - 1.0) < 0.1

def test_expansion_history():
    from mtft.cosmology import FriedmannMTFT
    hist = FriedmannMTFT().expansion_history()
    assert len(hist["a"]) == 500

# ── Lattice ──────────────────────────────────────────────────

def test_lattice_cold_start():
    from mtft.lattice import LatticeConfig, avg_plaquette
    cfg = LatticeConfig(N=2, L=2)
    assert cfg.volume == 16
    plaq = avg_plaquette(cfg)
    assert abs(plaq - 1.0) < 1e-10  # identity links → plaquette = 1

def test_mtft_action_cold():
    from mtft.lattice import LatticeConfig, MTFTAction
    cfg = LatticeConfig(N=2, L=2)
    action = MTFTAction(beta=4.0, kappa=1.0, y=0.18)
    # Cold start → all plaquettes = I → Wilson = 0, arith = 0
    assert action.wilson_action(cfg) < 1e-10
    assert action.arithmetic_action(cfg) < 1e-10

def test_mass_gap_bound():
    from mtft.lattice import MTFTAction
    action = MTFTAction(y=0.15, n_max=100)
    mu3 = action.mass_gap_bound(N=3)
    assert mu3 > 0

def test_metropolis_runs():
    from mtft.lattice import LatticeConfig, MTFTAction, metropolis_sweep
    cfg = LatticeConfig(N=2, L=2)
    action = MTFTAction(beta=4.0, kappa=0.5, y=0.18)
    acc = metropolis_sweep(cfg, action, n_hits=3, epsilon=0.3)
    assert 0 <= acc <= 1

# ── X₀(143) ──────────────────────────────────────────────────

def test_generation_count():
    from mtft.x0_143 import generation_count
    assert generation_count() == 3

def test_tano_mass_formula():
    from mtft.x0_143 import tano_mass_predictions
    pred = tano_mass_predictions()
    assert pred["errors_percent"]["muon"] < 2.0
    assert pred["errors_percent"]["tau"] < 2.0

def test_hecke_polynomial_roots():
    from mtft.x0_143 import hecke_roots_f3, A2_COMPLEX
    roots = hecke_roots_f3()
    complex_roots = [r for r in roots if abs(r.imag) > 0.01]
    assert len(complex_roots) >= 2
    best = min(complex_roots, key=lambda r: abs(r - A2_COMPLEX))
    assert abs(best - A2_COMPLEX) < 0.01

def test_koide_angle_vs_experiment():
    from mtft.x0_143 import koide_angle_tano, koide_angle_experimental
    theta_tano = koide_angle_tano()
    theta_exp = koide_angle_experimental()
    assert abs(math.degrees(theta_tano) - math.degrees(theta_exp)) < 0.5

# ── Koide ────────────────────────────────────────────────────

def test_koide_leptons():
    from mtft.koide import koide_leptons
    K = koide_leptons()
    assert abs(K - 2/3) < 1e-4

def test_koide_manifold_closed():
    from mtft.koide import koide_manifold_point, koide_ratio
    # Only phases where all √m > 0: need 1+√2 cos(φ+2πi/3) > 0 for all i
    # Valid window: φ ∈ roughly (-0.3, 0.9) and shifted copies
    for phi in [0.22, 0.3, 0.6]:
        w = koide_manifold_point(phi, r=2.0)
        if all(w > 0):
            m = w ** 2
            K = koide_ratio(m[0], m[1], m[2])
            assert abs(K - 2/3) < 1e-10, f"K={K} at φ={phi}"

def test_tau_prediction():
    from mtft.koide import predict_tau_mass
    result = predict_tau_mass()
    assert result["error_percent"] < 0.01  # 0.006% claimed

def test_bisector_gives_two_thirds():
    from mtft.koide import bisector_stability_K
    assert abs(bisector_stability_K() - 2/3) < 1e-14

# ── Burning Ship ─────────────────────────────────────────────

def test_burning_ship_bounded():
    from mtft.burning_ship import burning_ship_iterate
    assert burning_ship_iterate(0.0 + 0.0j) == 200  # origin is bounded
    assert burning_ship_iterate(10.0 + 0.0j) < 10   # escapes fast

def test_fold_budget():
    from mtft.burning_ship import ANISOTROPIC
    assert abs(ANISOTROPIC.fold_budget - 0.5) < 0.01

def test_jacobian_det():
    from mtft.burning_ship import jacobian_determinant
    z = 1.0 + 2.0j
    assert abs(jacobian_determinant(z) - 4*abs(z)**2) < 1e-10

# ── Dimensional Bridge ───────────────────────────────────────

def test_charge_identity():
    from mtft.dimensional_bridge import charge_from_feigenbaum
    result = charge_from_feigenbaum()
    assert result["error_percent"] < 0.05

def test_product_lock():
    from mtft.dimensional_bridge import feigenbaum_product_lock
    result = feigenbaum_product_lock()
    assert result["error_percent"] < 0.05

def test_electron_mass_formula():
    from mtft.dimensional_bridge import electron_mass_from_eta
    result = electron_mass_from_eta()
    # The ln(m_P/m_e) should be within a few percent
    assert result["error_percent"] < 10  # rough agreement at this depth

# ── Decay ────────────────────────────────────────────────────

def test_modular_decay():
    from mtft.decay import ModularDecay
    d = ModularDecay(Gamma_0=1e-10, delta_eff=0.001)
    assert d.half_life() > 0
    assert d.half_life_correction >= 1.0  # correction makes it longer-lived

# ── Verify Scorecard ─────────────────────────────────────────

def test_full_report():
    from mtft.verify import full_report
    report = full_report()
    assert len(report) >= 20
    sub1 = sum(1 for p in report if p.error_percent < 1.0)
    assert sub1 >= 10  # at least 10 predictions better than 1%

# ── Viz (data export, no plot) ───────────────────────────────

def test_viz_stiffness_data():
    from mtft.viz import stiffness_landscape
    data = stiffness_landscape(n_points=20, N_values=(2, 3))
    assert "y" in data
    assert "mu_2" in data
    assert len(data["y"]) == 20

def test_viz_rotation_curve_data():
    from mtft.viz import rotation_curve_plot
    data = rotation_curve_plot(n_points=50)
    assert "v_total" in data
    assert len(data["v_total"]) == 50

def test_viz_export_plotly():
    from mtft.viz import stiffness_landscape, export_for_plotly
    data = stiffness_landscape(n_points=10, N_values=(3,))
    exported = export_for_plotly(data)
    assert isinstance(exported["y"], list)

# ── Crypto ───────────────────────────────────────────────────

def test_arithmetic_hash_deterministic():
    from mtft.crypto import ArithmeticHash
    h = ArithmeticHash()
    assert h.hexdigest(b"MTFT") == h.hexdigest(b"MTFT")

def test_arithmetic_hash_avalanche():
    from mtft.crypto import ArithmeticHash
    h = ArithmeticHash()
    av = h.verify_avalanche(b"test message for hash")
    assert av["bits_changed"] > 30  # significant diffusion

def test_hash_different_inputs():
    from mtft.crypto import ArithmeticHash
    h = ArithmeticHash()
    assert h.hexdigest(b"hello") != h.hexdigest(b"Hello")

def test_burning_ship_prng():
    from mtft.crypto import BurningShipPRNG
    prng = BurningShipPRNG(seed=42)
    vals = [prng.random_float() for _ in range(100)]
    assert all(0 <= v < 1 for v in vals)
    assert len(set(vals)) > 90  # should be nearly all unique

def test_sl2z_key_exchange():
    from mtft.crypto import ModularKeyExchange
    w = ModularKeyExchange.generate_secret_word(10, seed=7)
    pub = ModularKeyExchange.public_key(w)
    # Public key should still be in SL(2,Z): det = 1
    det = int(pub[0,0]*pub[1,1] - pub[0,1]*pub[1,0])
    assert det == 1

def test_arithmetic_lattice():
    from mtft.crypto import ArithmeticLattice
    lat = ArithmeticLattice(dimension=8, y=0.18)
    assert lat.shortest_vector_estimate() > 0
    lwe = lat.generate_lwe_instance(secret_dim=4)
    assert lwe["A"].shape == (4, 4)

# ── Quantum Computing ────────────────────────────────────────

def test_gell_mann_count():
    from mtft.quantum import gell_mann_matrices
    assert len(gell_mann_matrices(2)) == 3   # Pauli matrices
    assert len(gell_mann_matrices(3)) == 8   # Gell-Mann matrices
    assert len(gell_mann_matrices(4)) == 15

def test_gell_mann_traceless():
    from mtft.quantum import gell_mann_matrices
    for T in gell_mann_matrices(3):
        assert abs(np.trace(T)) < 1e-12

def test_holonomy_gate_unitary():
    from mtft.quantum import HolonomyGate
    gate = HolonomyGate(d=3, generator_index=0, angle=0.5)
    assert gate.is_unitary()

def test_topological_qudit():
    from mtft.quantum import TopologicalQudit
    q = TopologicalQudit(d=3, coefficients=[1, 1, 1])
    assert q.manifold_dimension == 8
    assert q.total_maps == 56
    assert q.nontrivial_maps == 18
    assert q.topological_gap() > 0

def test_topological_spectrum_d7():
    from mtft.quantum import topological_spectrum_info
    info = topological_spectrum_info(7)
    assert info["manifold_dim"] == 48
    assert info["total_maps"] == 17296

def test_arithmetic_code_roundtrip():
    from mtft.quantum import ArithmeticCode
    code = ArithmeticCode(d_logical=2, n_physical=16)
    logical = np.array([1, 0], dtype=complex)
    encoded = code.encode(logical)
    decoded = code.decode(encoded)
    # Should approximately recover the logical state
    assert abs(abs(decoded[0]) - 1.0) < 0.2

def test_skyrmion_number():
    from mtft.quantum import skyrmion_number
    assert skyrmion_number(0, 1) == -1
    assert skyrmion_number(0, -5) == 5
    assert skyrmion_number(5, -8) == 13
