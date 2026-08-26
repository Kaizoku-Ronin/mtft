"""v0.21: native X0(143) periods, Hodge bridge, and physics-facing diagnostics."""
import math
import numpy as np
import mpmath as mp

import mtft
import mtft.periods as P
from mtft.periods import gates as G


def test_version_and_data_ship():
    assert mtft.__version__ == "0.21.0"
    for name in ("X0_143_period_data_v6.json","X0_143_period_basis_v6.json",
                 "X0_143_m7_harmonic_basis.json","PROVENANCE.txt"):
        assert P.data_path(name).exists()


def test_integral_symplectic_exact():
    r=G.gate_integral_symplectic();assert r["det_S"]==1


def test_period_reconstruction_and_serialization_guard():
    r=G.gate_period_reconstruction(50)
    assert mp.mpf(r["tau_frozen_residual"]) < mp.mpf("1e-45")
    assert mp.mpf(r["legacy_Omega_sym_discrepancy"]) > 1


def test_exact_basis_bridge():
    r=G.gate_basis_bridge(45)
    assert (r["det_R"],r["det_C"],r["det_P"]) == (1,1,1)


def test_true_hodge_structure_on_hecke_stage():
    r=G.gate_hodge_bridge(50)
    assert r["J_square"] < 1e-10
    assert max(r["Hecke_commutators"].values()) < 1e-10
    assert r["star_anticommutator_rel"] < 1e-10


def test_qexpansion_native_span_and_tail():
    r=G.gate_qexpansion_span(40)
    assert r["rank"]==13 and r["max_abs_span_residual"]<1e-8
    assert r["q140_tail_bound_at_y_1_sqrt143"] < 1e-27


def test_bergman_density_positive():
    y=1/math.sqrt(143)
    b=P.bergman_density(1j*y,nmax=140,dps=45)
    assert b>0


def test_charge_energy_positive_and_even():
    n=[0]*13;m=[0]*13;n[0]=1
    a=P.charge_energy(n,m,45)
    n[0]=-1;b=P.charge_energy(n,m,45)
    assert a>0 and abs(a-b)<mp.mpf("1e-35")


def test_quantitative_m8_true_period_J():
    r=G.gate_quantitative_m8(50)
    # First v0.21 measured anchor; not a particle observable.
    assert abs(r["antilinear_fraction"]-0.06796942853668321) < 1e-10
    assert r["commutator_star_odd_rel"] < 1e-10
