"""mtft.surface v0.26.0: Ising two-route gates, frozen certified data, dynamics controls (no GP needed)."""
import numpy as np
import pytest

from mtft.surface import ising as I, frozen as FR, dynamics as D


@pytest.mark.parametrize("N", [6, 11, 15, 35])
def test_ising_pfaffian_sum_equals_brute_force(N):
    r = I.gate_report(N)
    assert r["status"] == "PASS" and r["relative_difference"] < 1e-12
    assert all(g["status"] == "PASS" for g in r["structural_gates"])


def test_ising_143_setup_and_sample():
    S = I.IsingSurface(143)
    assert S.genus == 13 and S.F == 56 and S.E == 84
    assert all(g["status"] == "PASS" for g in S.K.gates)
    assert S.K.twists.shape == (26, 252)
    with pytest.raises(RuntimeError):
        S.dimer_sum(I.T_CRITICAL_HONEYCOMB)          # 4^13 refused without allow_long
    smp = S.sample(I.T_CRITICAL_HONEYCOMB, n_samples=8, seed=1)
    assert len(smp["even"]) + len(smp["odd"]) == 8


def test_frozen_x0143_gates_without_gp():
    d = FR.x0143()
    assert all(d["gates"].values())
    assert d["T2"].shape == (26, 26) and d["provenance"]["gp_script_sha1"]


def test_dynamics_controls_orthonormal_frame():
    st = D.stage_from_frozen()
    g = st.gates()
    assert g["J_squared"] < 1e-12 and g["G_min_eigenvalue"] > 0.05
    c = D.closure_controls(st)
    assert c["random_pair"]["dimension"] == 351             # generic: sp(26)
    assert c["hecke_commuting"]["dimension"] == 4            # abelian
    assert c["block_diagonal"]["dimension"] == 127 == c["block_diagonal_expected"]
    F = D.hamiltonian(st, np.diag(np.arange(1.0, 27.0)))
    assert D.flow_gates(st, F)["hamiltonian_residual"] < 1e-12
