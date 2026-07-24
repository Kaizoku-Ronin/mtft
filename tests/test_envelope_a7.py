"""test_envelope_a7.py — regression tests for the envelope_slope A.7
patch (mtft 0.9.0, audit S.4 / T-E4).

Three facts under test:
  1. the min_bin guard alone (legacy 91-point grid) discards all eleven
     bins — guard and density had to ship together;
  2. the upgraded default density (A.7 recommended_samples_per_decade)
     is healthy, densification-stable, and keeps the corrected Th 1
     anchors (on-line slope 0, off-line quadruplet 1/2 − β₀);
  3. terminal-bin leverage is guarded and the bin width does not alias
     the γ₃ stride resonance.
"""
import math

import numpy as np
import pytest

from mtft.estimator_standards import (recommended_samples_per_decade,
                                      stride_resonance_check)
from mtft.riemann import (RIEMANN_ZEROS, envelope_slope, offline_quadruplet,
                          on_line_zeros)

GAMMA1 = float(RIEMANN_ZEROS[0])   # RIEMANN_ZEROS is a list of float
GAMMA3 = float(RIEMANN_ZEROS[2])   # ordinates (float.imag == 0.0 — audit catch)


# ── 1. the naive-rewire trap ────────────────────────────────────

def test_legacy_density_guard_discards_all_bins():
    """91 points over 5.2 decades ≈ 8.8/half-decade bin < min_bin=10:
    the bare guard discards every bin (ChatGPT's observation, R)."""
    on = on_line_zeros(4)
    with pytest.raises(ValueError, match="usable bins"):
        envelope_slope(on, n_points=91, min_bin=10)


# ── 2. the upgraded default ─────────────────────────────────────

def test_default_density_is_healthy():
    on = on_line_zeros(4)
    r = envelope_slope(on)
    spd = recommended_samples_per_decade(GAMMA1)
    assert r["samples_per_decade"] >= spd - 1
    assert r["n_bins"] >= 10
    assert r["n_bins_dropped"] <= 2


def test_online_slope_anchor():
    on = on_line_zeros(4)
    r = envelope_slope(on)
    assert abs(r["slope"]) < 0.05


def test_offline_quadruplet_anchor():
    on = on_line_zeros(4)
    for beta0 in (0.6, 0.75, 0.9):
        zs = offline_quadruplet(beta0, RIEMANN_ZEROS[0]) + on[1:]
        r = envelope_slope(zs)
        assert abs(r["slope"] - (0.5 - beta0)) < 0.05, \
            f"β₀={beta0}: slope {r['slope']:.4f}"


def test_density_stability():
    """Slope is invariant under 2× densification (anti-aliasing)."""
    on = on_line_zeros(4)
    r1 = envelope_slope(on, samples_per_period=10)
    r2 = envelope_slope(on, samples_per_period=20)
    assert abs(r1["slope"] - r2["slope"]) < 5e-3


# ── 3. terminal-bin leverage and γ₃ aliasing ────────────────────

def test_terminal_bin_removal():
    """Same range, same grid: the terminal bin (12 samples over
    [-2.0, -1.7]) is kept at min_bin=10 and dropped at min_bin=15 —
    and the guarded slope moves by less than the bin's leverage."""
    on = on_line_zeros(4)
    r_keep = envelope_slope(on, -7.0, -1.7, n_points=200, min_bin=10)
    assert r_keep["n_bins_dropped"] == 0
    r_drop = envelope_slope(on, -7.0, -1.7, n_points=200, min_bin=15)
    assert r_drop["n_bins_dropped"] == 1      # the sparse terminal bin
    assert abs(r_drop["slope"] - r_keep["slope"]) < 0.02


def test_bin_width_not_gamma3_resonant():
    """The half-decade envelope bin width must not alias γ₃'s period
    (the 54.994-cycles/6-decades trap lives at window stride 6)."""
    _, _, resonant = stride_resonance_check(GAMMA3, 0.5)
    assert not resonant
    # and the trap itself still flags where it should (regression):
    assert stride_resonance_check(GAMMA3, 6.0)[2]
