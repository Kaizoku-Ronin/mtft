"""Exception-Spacing Curvature Law: fast structural gates + slow instances."""
import pytest
from mpmath import mp, mpf, zeta, ln

from mtft import exception_spectrum as S


def test_rigidity_quarter():
    mp.dps = 40
    K = S.K_atoms([mpf(0), mpf(1), mpf(2)], [0, 0, 1],
                  [mpf(1), mpf("0.4"), mpf("0.2")])
    assert abs(K - mpf(1) / 4) < mpf(10) ** -30


def test_line_atom_lemma_matches_numeric():
    mp.dps = 60
    line = [(mpf(0), mpf(1)), (ln(2), mpf("0.31")), (ln(3), mpf("0.08"))]
    pred = S.line_atom_variation(line, ln(7), mpf("0.17"), ln(15))
    xs = [x for x, _ in line] + [ln(7), ln(15)]
    mks = [0, 0, 0, 1, 1]
    h = mpf(10) ** -25
    K0 = S.K_atoms(xs, mks, [w for _, w in line] + [mpf("0.17"), mpf(0)])
    K1 = S.K_atoms(xs, mks, [w for _, w in line] + [mpf("0.17"), h])
    assert abs((K1 - K0) / h - pred) / abs(pred) < mpf(10) ** -15


def test_flat_base_lemma_and_ratio():
    mp.dps = 60
    cloud = [(mpf("-0.7"), mpf("0.9")), (mpf("0.4"), mpf("0.3")),
             (mpf("1.3"), mpf("1.1"))]
    V, t, r, lam = mpf("0.37"), mpf("1.21"), mpf("2.05"), mpf("0.3")
    xs = [x for x, _ in cloud] + [x + t for x, _ in cloud]
    mks = [0, 0, 0, 1, 1, 1]
    w0 = [q for _, q in cloud] + [V * q for _, q in cloud]
    h = mpf(10) ** -25
    K0 = S.K_atoms(xs, mks, w0, lam=lam)
    assert abs(K0) < mpf(10) ** -50
    dU = (S.K_atoms(xs + [r], mks + [0], w0 + [h], lam=lam) - K0) / h
    dM = (S.K_atoms(xs + [t + r], mks + [1], w0 + [h], lam=lam) - K0) / h
    predU = S.flat_base_variation(cloud, V, r, lam=lam)
    assert abs(dU - predU) / abs(predU) < mpf(10) ** -15
    assert abs(dM / dU + 1 / V) * V < mpf(10) ** -15


def test_phase_classifier_and_channels():
    assert S.classify_phase(7, 15).startswith("spherical")
    assert S.classify_phase(3, 6).startswith("critical")
    assert S.classify_phase(3, 5).startswith("hyperbolic")
    chans = S.defect_spectrum([3, 6, 10, 15, 21, 28, 36], 3, 9)
    rs = [float(r) for r, _, _ in chans]
    assert abs(rs[0] - 10 / 3) < 1e-12 and 5.0 not in rs and 7.0 not in rs


@pytest.mark.slow
def test_legendre_constant_beta100():
    mp.dps = 200
    C = ln(mpf(15) / 7) ** 2 / (4 * ln(2) ** 2)
    b = mpf(100)

    def tri():
        a = 0
        while True:
            n = 4 ** a * 7
            k = n
            while True:
                yield k
                k += 8 * 4 ** a
                if k > 10 ** 9:
                    break
            a += 1

    def forb():
        n = 0
        while True:
            n += 1
            m, v = n, 0
            while m % 2 == 0:
                m //= 2
                v += 1
            if v % 2 == 0 and m % 8 == 7:
                yield n

    K = S.K_marked_set(b, forb(), 7, rel_target=mpf(10) ** -45)
    assert abs((mpf(1) / 4 - K) / (mpf(14) / 15) ** b - C) < mpf(10) ** -16
