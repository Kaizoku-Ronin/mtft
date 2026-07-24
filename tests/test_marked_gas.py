"""test_marked_gas.py — the marked primon gas module (mtft 0.9.0).

The 13-gate green suite (note v0.1.1 §7) as pytest, plus module-API
regression tests.  Classes per the Legend: EXACT / CERTIFIED(bound) /
DIAGNOSTIC.  Run times target < 60 s total (quick-mode N's).
"""
import math

import numpy as np
import pytest
import mpmath as mp

import mtft.marked_gas as MG
from mtft.marked_gas import (ALPHA_COLD, B_COLD, Certified,
                             bc_deformation, cold_gas_report, correlator,
                             edge_mass, flow_phase, gates, kms_check,
                             psi_coefficients, spectral_function, spectrum,
                             weights_sieve, z1, z2, zD_certified_interval)


# ── constants (pinned to mpmath) ────────────────────────────────

def test_alpha_cold_constant():
    with mp.workdps(30):
        assert abs(ALPHA_COLD - float(-mp.zeta(2, derivative=1))) < 1e-15


def test_b_cold_constant():
    with mp.workdps(30):
        assert abs(B_COLD - float(mp.e ** (-mp.zeta(2, derivative=2)))) < 1e-15


def test_amplitude_target_value():
    with mp.workdps(30):
        a = float(-mp.zeta(2, derivative=1))
        tgt = float(mp.e ** (-mp.zeta(2, derivative=2)) / mp.gamma(a + 1))
    assert abs(tgt - 0.14027492490189596) < 1e-15


# ── traces ──────────────────────────────────────────────────────

def test_z1_trace_resid():
    r = z1(2)
    assert r.detail["resid"] <= r.bound + 1e-27
    assert isinstance(r, Certified) and r.err_class.startswith("CERTIFIED")


def test_z1_trace_beta3():
    r = z1(3)
    assert abs(r.value - 1.202056903159594) < 1e-12


def test_z2_trace_resid():
    r = z2(2)
    assert r.detail["resid"] <= r.bound + 1e-27


def test_z2_target():
    r = z2(2)
    assert abs(r.value - r.detail["target"]) < 1e-25


# ── the convolution interval ────────────────────────────────────

def test_zD_interval_contains_target():
    r = zD_certified_interval(2, N=1_000_000)
    assert r.detail["inside"]
    assert r.bound < 1e-10


def test_zD_interval_width_scaling():
    """Width ~ (log N)²/(2N²) across N (audit U.1)."""
    w1 = zD_certified_interval(2, N=500_000).bound
    w2 = zD_certified_interval(2, N=2_000_000).bound
    pred = (math.log(500_000) / math.log(2_000_000)) ** 2 * (2_000_000 / 500_000) ** 2
    assert 0.6 < (w1 / w2) / pred < 1.7


def test_zD_reorder_identity():
    """G3b: the (d,k)-swapped head equals the direct sieve sum (exact)."""
    N, b = 50_000, 2.0
    w = weights_sieve(N)
    nn = np.arange(2, N + 1, dtype=np.float64)
    direct = float(np.sum(w[2:] * nn ** (-b)))
    d = np.arange(2, N + 1, dtype=np.int64)
    M = N // d
    H = MG._kahan_cumjumps(np.arange(1, N + 1, dtype=np.float64) ** (-b),
                           list(np.unique(M)))
    swap = float(np.sum(np.log(d) * d.astype(np.float64) ** (-b - 1.0)
                        * np.array([H[m] for m in M])))
    assert abs(direct - swap) < 5e-13


def test_zD_target_constant():
    """−ζ(2)ζ′(3) matches the census value (cross-engine anchor)."""
    with mp.workdps(30):
        t = float(-mp.zeta(2) * mp.zeta(3, derivative=1))
    assert abs(t - 0.32590460645923011912) < 1e-18


# ── kinematics / flow / KMS ─────────────────────────────────────

def test_spectrum_values():
    E = spectrum(beta=2, nmax=10)
    n = np.arange(2, 11, dtype=float)
    assert np.allclose(E, 3.0 * np.log(n) - np.log(np.log(n)), rtol=0, atol=1e-14)


def test_flow_phase_unit_modulus():
    for p in (2, 3, 5):
        for t in (0.7, 1.9, math.pi):
            assert abs(abs(flow_phase(p, 7, t)) - 1.0) < 1e-14


def test_flow_phase_vs_matrix():
    """G4's real leg (v0.9.1, audit R2): the closed form against the
    DEFINITION path — literal matrix conjugation U mu U† with
    K = −log ρ̂ from the Gibbs weights (no spectral-formula input).
    v0.9.0's version of this test re-typed the closed form and compared
    it to itself."""
    from mtft.marked_gas import flow_phase_matrix
    for p in (2, 3, 5):
        for t in (0.7, 1.9, math.pi):
            lev, ph = flow_phase_matrix(p, t, beta=2.0, Nb=400)
            worst = float(np.max(np.abs(ph - flow_phase(p, lev, t,
                                                         beta=2.0))))
            assert worst < 1e-12


def test_kms_termwise_identity():
    """All four KMS legs (matrix, spectral, restricted spectral, and
    the apples-to-apples cross-check) at the 1e-12 bound."""
    r = kms_check(2, 2, 0.7, nbasis=2000)
    assert r.value < r.bound
    for leg in ("matrix_residual", "spectral_residual",
                "spectral_residual_sub", "cross_residual"):
        assert r.detail[leg] < 1e-12, leg


def test_kms_wrong_sign_fails():
    """The t − i convention must NOT satisfy the identity (G5's story:
    the convention is pinned by the gate, not chosen)."""
    r = kms_check(2, 2, 0.7, nbasis=2000)
    assert r.detail["wrong_sign_control"] > 1e-2


def test_bc_deformation_uv_restoration():
    """UV restoration is asymptotic: |D_p(t) − 1| ~ t·log p/log n."""
    p, t = 2, 1.3
    n_small = abs(bc_deformation(p, 100, t) - 1.0)
    n_large = abs(bc_deformation(p, 1_000_000, t) - 1.0)
    assert n_large < n_small
    rate = n_large / (t * math.log(p) / math.log(1_000_000))
    assert 0.9 < rate < 1.1


# ── cold gas ────────────────────────────────────────────────────

def test_weights_sieve_spot():
    w = weights_sieve(20)
    assert abs(w[2] - math.log(2) / 2) < 1e-16
    assert abs(w[12] - (math.log(2)/2 + math.log(3)/3 + math.log(4)/4
                        + math.log(6)/6 + math.log(12)/12)) < 1e-14


def test_psi_recurrence_spot():
    w, a, A = psi_coefficients(200)
    # n a_n = Σ_{k≤n} w_k a_{n−k} re-checked directly at n = 137
    n = 137
    lhs = n * a[n]
    rhs = sum(w[k] * a[n - k] for k in range(1, n + 1))
    assert abs(lhs - rhs) < 1e-10 * max(1.0, abs(rhs))


def test_cold_gas_slope_quick():
    rep = cold_gas_report(N=8_000, lo=2_000)
    assert abs(rep["slope_two_point"]["value"] - rep["alpha"]) < 2e-2


def test_cold_gas_amplitude_quick():
    rep = cold_gas_report(N=8_000, lo=2_000)
    assert abs(rep["amplitude_endpoint"]["value"]
               / rep["amplitude_target"] - 1) < 0.05


def test_cold_gas_amplitude_trend():
    """The honest trend converges to B/Γ(α+1) from above (U.3 item 3)."""
    rep = cold_gas_report(N=100_000, lo=20_000)
    cps = rep["summatory_over_n_alpha"]
    tgt = rep["amplitude_target"]
    keys = sorted(cps)
    vals = [cps[k] for k in keys]
    assert all(v > tgt for v in vals)
    assert all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))
    assert abs(vals[-1] / tgt - 1) < 2e-3


def test_cold_gas_pointwise_diagnostic():
    rep = cold_gas_report(N=8_000, lo=2_000)
    assert abs(rep["pointwise_exponent"]["value"]
               - rep["pointwise_exponent"]["target"]) < 0.05


# ── correlator / spectral function / edge law ───────────────────

def test_correlator_bound():
    r = correlator(2, 1.0, N=2000)
    assert abs(r.value) <= 1.0 + r.bound
    assert r.bound < 1e-5


def test_correlator_kms_consistency():
    """F(t+i) equals the damped correlator (termwise; same content as
    G5) — evaluated through the module API at the KMS point (v0.9.1:
    correlator accepts complex t with Im t ≥ 0; audit R3)."""
    p, t, beta = 2, 0.9, 2.0
    rF = correlator(p, t + 1j, beta, N=2000)
    n = np.arange(2, 2001, dtype=np.float64)
    rho = np.log(n) * n ** (-(beta + 1.0))
    dE = (beta + 1.0) * math.log(p) - np.log(np.log(p * n) / np.log(n))
    Z2 = float(-mp.zeta(3, derivative=1))
    Gt = np.sum(rho / Z2 * np.exp(-dE) * np.exp(1j * t * dE))
    assert abs(rF.value - Gt) < 1e-12
    with pytest.raises(ValueError, match="Im t"):
        correlator(p, t - 1j, beta, N=100)


def test_spectral_function_monotone_to_edge():
    sf = spectral_function(p=2, beta=2, nmax=5_000)
    assert sf["increasing"] and sf["last_line_below_edge"]
    assert abs(sf["edge"] - 3.0 * math.log(2)) < 1e-15


def test_edge_mass_pinned_convention():
    """Per-level evaluation matches the asymptote at ε = 0.1 (U.4/V.3)."""
    r = edge_mass(p=2, beta=2, eps=0.1, nmax=2_000_000)
    assert 0.95 < r.detail["ratio"] < 1.05


def test_edge_mass_convention_is_per_level():
    """The pinned convention must NOT round M (the ε = 0.2 artifact).
    v0.9.1: the reference is now genuinely independent — the exact
    Hurwitz-zeta tail Σ_{n≥a} (log n) n^{-3} = −ζ′(3, a), not a
    transcription of the function's numpy internals (audit R2)."""
    r = edge_mass(p=2, beta=2, eps=0.2, nmax=2_000_000)
    m0 = 2
    while math.log1p(math.log(2) / math.log(m0)) >= 0.2:
        m0 += 1
    assert m0 == r.detail["m0"] == 23
    with mp.workdps(30):   # context form: newer mpmath guards module.dps
        ref = float((-mp.zeta(3, m0, derivative=1)
                     + mp.zeta(3, 2_000_001, derivative=1))
                    / (-mp.zeta(3, derivative=1)))
    assert abs(r.value - ref) < 1e-12


def test_edge_mass_em_predictor_subppm():
    """Audit R1: the m₀-anchored EM predictor is deterministic and
    sub-ppm where the bare law oscillates ±f(M)/2 (5.5% at ε = 0.3).
    Shipped deviations from 1 (mass_plus_tail/pred_em): −4.6e-9,
    −1.8e-12, +1.1e-15 at ε = 0.3, 0.2, 0.1 — sub-ppm throughout,
    deterministic in m₀ (Addendum X)."""
    for eps, tol in ((0.3, 1e-5), (0.2, 1e-6), (0.1, 1e-6)):
        d = edge_mass(p=2, beta=2, eps=eps, nmax=2_000_000).detail
        assert abs(d["ratio_em_corrected"] - 1.0) < tol, eps
    # and the bare law must show the old scatter (guard against a
    # future regression that re-inflates the predictor)
    d03 = edge_mass(p=2, beta=2, eps=0.3).detail
    assert abs(d03["ratio"] - 1.0) > 0.01


def test_edge_mass_tail_correction():
    """At ε = 0.07 the nmax = 2×10⁶ truncation is 7.5e-5 of the mass
    (audit R1); mass_plus_tail must cancel it to sub-ppm."""
    d = edge_mass(p=2, beta=2, eps=0.07, nmax=2_000_000).detail
    assert abs(d["ratio_em"] - 0.999925) < 2e-6          # truncation visible
    assert abs(d["ratio_em_corrected"] - 1.0) < 1e-6     # and corrected


# ── the gate suite ──────────────────────────────────────────────

def test_gates_quick_all_green():
    g = gates(quick=True)
    assert g["all_green"], [k for k, v in g.items() if not v.get("ok", True)]


def test_gates_classes_present():
    g = gates(quick=True)
    for name, v in g.items():
        if name == "all_green":
            continue
        assert v["err_class"].split("(")[0] in (
            "EXACT", "CERTIFIED", "DIAGNOSTIC", "PHENO", "GIVEN")


def test_certified_dataclass_shape():
    r = z1(2)
    assert hasattr(r, "value") and hasattr(r, "err_class")
    assert hasattr(r, "bound") and hasattr(r, "detail")


def test_legend_registration():
    from mtft.legend import REGISTRY
    assert "marked_gas" in REGISTRY
    assert "kms_flow" in REGISTRY
    assert "cold_gas_amplitude" in REGISTRY
    assert "spectral_edge_soft" in REGISTRY
