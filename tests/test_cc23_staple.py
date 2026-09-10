"""CC-23: local Wilson action change must equal full recomputation (the test that was missing)."""
import numpy as np

from mtft import lattice as L


def test_local_delta_action_matches_full_recomputation():
    rng = np.random.default_rng(23)
    cfg = L.LatticeConfig(N=2, L=3, L_t=3)
    for t in range(3):
        for x in range(3):
            for y in range(3):
                for z in range(3):
                    for mu in range(4):
                        cfg.set_link((t, x, y, z), mu, L.random_su_n(2, rng))
    act = L.MTFTAction(beta=2.0)
    for _ in range(12):
        site = tuple(int(v) for v in rng.integers(0, 3, 4)); mu = int(rng.integers(4))
        S0 = act.wilson_action(cfg)
        staple = L._staple_sum(cfg, site, mu, act.beta)
        Uold = cfg.get_link(site, mu).copy(); Unew = L.su_n_near_identity(2, 0.3, rng) @ Uold
        d_local = -(act.beta / 2) * np.real(np.trace((Unew - Uold) @ staple))
        cfg.set_link(site, mu, Unew); d_full = act.wilson_action(cfg) - S0; cfg.set_link(site, mu, Uold)
        assert abs(d_local - d_full) < 1e-10
