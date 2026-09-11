"""CC-24 readout coefficient; CC-25 unit-extent guard; exact-input hardening (from Astra's v0.28.0 audit)."""
import numpy as np
import pytest
import sympy as sp

from mtft import lattice as L
from mtft.surface import marked as MK, oldsector as OS


def test_cc24_readout_is_sector_independent():
    a = sp.Symbol("alpha")
    for p, q, ap, aq in ((13, 17, 4, -2), (17, 13, -2, 4), (13, 19, 4, 0)):
        S = OS.two_prime_sector(p, q, ap, aq)
        assert MK.readout_coefficient(S) == (q + 1) ** 2 - aq ** 2
        assert MK.two_objectives(S)["readout_alpha"] == a


def test_exact_input_guards():
    with pytest.raises(TypeError):
        OS.local_block(13.0, 4)
    with pytest.raises(TypeError):
        OS.local_block(13, True)
    S = OS.two_prime_sector(13, 17, 4, -2)
    with pytest.raises(TypeError):
        MK.decomposition(sp.Matrix(4, 4, lambda i, j: 0.5 if i == j else 0), S["local_p"]["G"], S["local_q"]["G"])


def test_cc25_unit_extent_rejected_before_mutation():
    rng = np.random.default_rng(0)
    cfg = L.LatticeConfig(N=2, L=1, L_t=2)
    before = cfg.links.copy()
    with pytest.raises(ValueError):
        L.metropolis_sweep(cfg, L.MTFTAction(beta=2.0), rng=rng) if "rng" in L.metropolis_sweep.__code__.co_varnames else L.metropolis_sweep(cfg, L.MTFTAction(beta=2.0))
    assert np.array_equal(cfg.links, before)
