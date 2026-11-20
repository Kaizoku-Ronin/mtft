"""WL-01: walls between twists — the exact wall ball, the wall distance, the systole band and one-way bottleneck costs
(records FF-14, FF-17, FF-20)."""
import math
import random

import pytest

from mtft.lines import twists as T
from mtft.lines import walls as WL

# FF-20's exact ball sizes |B(r)| on its grid to r = 12 (record ff20_numbers.json, T2).
FF20_BALL = {0.0: 1, 2.5: 1, 2.75: 5, 4.25: 5, 4.5: 7, 5.0: 7, 5.25: 11, 5.5: 33, 8.0: 289, 10.0: 1020, 12.0: 5559}


def _words(tmax):
    """Second route to the level-1 classes: every positive word in R, L with trace ≤ tmax, by cyclic minimal rotation."""
    R, L = ((1, 1), (0, 1)), ((1, 0), (1, 1))

    def mul(A, B):
        return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)) for i in range(2))

    out, stack = set(), [("", ((1, 0), (0, 1)))]
    while stack:
        w, P = stack.pop()
        if w and "R" in w and "L" in w:
            rot = min(w[i:] + w[:i] for i in range(len(w)))
            if all(rot != (rot[:d] * (len(rot) // d)) for d in range(1, len(rot)) if len(rot) % d == 0):
                out.add(rot)
        if len(w) >= tmax - 1:                               # a word with both letters and length n has trace ≥ n + 1
            continue
        for c, M in (("R", R), ("L", L)):
            Q = mul(P, M)
            if Q[0][0] + Q[1][1] <= tmax:
                stack.append((w + c, Q))
    return out


def test_level1_classes_two_routes():
    for tmax in (3, 4, 7, 12, 20):
        cls = WL.level1_classes(tmax)
        words = _words(tmax)
        assert len(cls) == len(words)
        assert sorted(A[0][0] + A[1][1] for _, A in cls) == sorted(_trace(w) for w in words)
    assert len(WL.level1_classes(3)) == 1 and len(WL.level1_classes(4)) == 3


def _trace(w):
    a, b, c, d = 1, 0, 0, 1
    for ch in w:
        a, b, c, d = (a, a + b, c, c + d) if ch == "R" else (a + b, b, c + d, d)
    return a + d


def test_systoles():
    P = WL.primitive_classes(3.0)
    assert len(P) == 8 and all(abs(L - 2 * math.acosh(2)) < 1e-12 for L, _, _ in P)
    assert WL.SYSTOLE == pytest.approx(2.633915793849633, abs=1e-15)
    chk = WL.check_systoles()
    assert chk["systole classes"] == 4 and chk["rank mod 2"] == 4
    assert chk["pairwise intersections"] == [[0] * 4 for _ in range(4)]
    assert list(chk["band sizes over all twists"].values()) == [2 ** 22 * k for k in (1, 4, 6, 4, 1)]
    bands = {k: WL.systole_band(v) for k, v in T.NAMED.items()}
    assert bands["minimum"] == bands["its W13 image"] == 4
    assert bands["φ*η"] == bands["φ*(0,1)"] == bands["φ*(1,1)"] == bands["FF-14 weakest"] == bands["O"] == 0
    for y in WL.systoles():                                  # a wall along a systole never changes the band
        for z in (T.NAMED["minimum"], T.NAMED["lowest symmetric"], 12345678):
            assert WL.systole_band(z ^ T.wall_twist(y)) == WL.systole_band(z)


def test_the_wall_ball_against_ff20():
    sizes = WL.ball_sizes(sorted(FF20_BALL))
    assert sizes == FF20_BALL
    c = WL.compare_routes(10.0)
    assert c["classes in the ball"] == 1020 and c["walls of 1-3 geodesics"] == 1019
    assert c["sums that beat the ball"] == 0 and c["walls of 1-3 geodesics where the routes disagree"] == 0


def test_ff17s_walls():
    eta, worst = T.NAMED["φ*η"], T.NAMED["FF-14 weakest"]
    assert eta ^ worst == 6533120 and T.dual_class(6533120) == 12664832          # FF-17: δ and its dual class
    assert WL.wall_distance(eta, worst, 10.0) == pytest.approx(9.897419174603995, abs=1e-12)
    assert WL.wall(eta, worst, 10.0) == (pytest.approx(9.897419174603995, abs=1e-12),)
    assert T.dual_class(T.NAMED["minimum"]) == 40693698                          # FF-17: O to the minimum, no wall ≤ 10
    assert WL.wall_distance(0, T.NAMED["minimum"], 10.0) is None
    assert WL.wall_distance(T.NAMED["minimum"], T.NAMED["its W13 image"], 12.0) is None    # FF-20: more than 16


def test_the_wall_distance_is_a_metric():
    ball = WL.wall_ball(10.0)
    pts = [T.wall_twist(y) for y in ball]
    rng = random.Random(143)
    for _ in range(300):
        a, b, c = rng.sample(pts, 3)
        t = rng.getrandbits(26)
        dab, dbc, dac = (WL.wall_distance(a, b, 10.0), WL.wall_distance(b, c, 10.0), WL.wall_distance(a, c, 10.0))
        assert dab == WL.wall_distance(b, a, 10.0) == WL.wall_distance(a ^ t, b ^ t, 10.0)
        if dab is not None and dbc is not None and dab + dbc <= 10.0:
            assert dac is not None and dac <= dab + dbc + 1e-9
    assert WL.wall_distance(12345, 12345, 1.0) == 0.0


def test_one_way_costs():
    assert WL.thin_wall_exponent(1.0, 0.0) == math.inf and WL.thin_wall_exponent(1.0, -0.5) == math.inf
    assert WL.thin_wall_exponent(1.0, 1.0) == pytest.approx(27 * math.pi ** 2 / 2)
    inf = math.inf
    cost = [[0, 5, inf, inf], [inf, 0, 2, inf], [9, inf, 0, 3], [1, inf, inf, 0]]
    c = WL.bottleneck(cost)
    assert c[0][3] == 5 and c[3][2] == 5 and c[1][0] == 3 and c[2][1] == 5
    n = len(c)
    assert all(c[i][k] <= max(c[i][j], c[j][k]) for i in range(n) for j in range(n) for k in range(n))
