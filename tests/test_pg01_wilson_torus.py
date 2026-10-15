"""PG-6 and CC-66 (2026-10-14): the whole Wilson-line torus — four critical points, two saddles, the exact zero, the metric curvatures.

Two routes (E2). Route 1 is the package's Newton census (`mtft.matter.vacuum_atlas.wilson_line_torus`, on PA-01's box lattice). Route 2 is
written here and uses neither the package's lattice nor its tau: a lattice sum over the disc |w| <= R, which the pairing (m, n) -> (-m - n, n)
preserves, so the real-structure zeros are exact at every cutoff; a census without Newton's method (the gradient's winding count on a
1024 x 1024 grid, a Lipschitz bound that excludes every cell away from the half-periods, and a contraction bound that leaves one zero at
each); and the curvatures in the kinetic metric from the disc's own Hessian."""
import math

import numpy as np
import pytest

from mtft.matter import ledger as L
from mtft.matter import puzzle as Z
from mtft.matter import vacuum_atlas as VA

TAU = complex(0.5, 1.0232745926964612)            # 143a1 (route 2 does not use the package's lattice or tau)
EPS = ((0, 1), (1, 0), (1, 1))
HALF_PERIODS = ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5))


def _disc(eps, R=40.0, box=64):
    m, n = np.meshgrid(np.arange(-box, box + 1), np.arange(-box, box + 1), indexing="ij")
    m, n = m.ravel().astype(float), n.ravel().astype(float)
    r = np.abs(m + n * TAU) / math.sqrt(TAU.imag)
    keep = (r > 0) & (r <= R) & (((eps[0] * m + eps[1] * n) % 2) == 1)
    assert np.abs(m[keep]).max() < box and np.abs(n[keep]).max() < box     # the disc lies inside the box
    return m[keep], n[keep], r[keep] ** -8


def _V(eps, th):
    m, n, w8 = _disc(eps)
    return float(-2 * np.sum(w8 * np.cos(2 * math.pi * (th[0] * m + th[1] * n))))


def _fields(eps, N):
    """The gradient of the disc sum at the nodes ((j + 1/2)/N, (k + 1/2)/N), by FFT, so that every half-period sits at a cell centre, never on
    a node or an edge; and the Lipschitz constants of the gradient (L2) and of the Hessian (L3) in theta."""
    m, n, w8 = _disc(eps)
    c, tp = -2 * w8, 2 * math.pi
    shift = np.exp(1j * math.pi * (m + n) / N)

    def at_nodes(weights):
        C = np.zeros((N, N), dtype=complex)
        np.add.at(C, (m.astype(int) % N, n.astype(int) % N), weights * shift)
        return np.real(np.fft.ifft2(C) * N * N)

    L2 = float(np.sum(np.abs(c) * tp ** 2 * (m * m + n * n)))           # |grad V(x) - grad V(y)| <= L2 |x - y|
    L3 = float(np.sum(np.abs(c) * tp ** 3 * (m * m + n * n) ** 1.5))    # ||H(x) - H(y)|| <= L3 |x - y|
    return at_nodes(c * 1j * tp * m), at_nodes(c * 1j * tp * n), L2, L3, (m, n, c)


def _hess(mnc, t):
    m, n, c = mnc
    k = -(2 * math.pi) ** 2 * c * np.cos(2 * math.pi * (t[0] * m + t[1] * n))
    return np.array([[np.sum(k * m * m), np.sum(k * m * n)], [np.sum(k * m * n), np.sum(k * n * n)]])


@pytest.fixture(scope="module")
def census():
    return {eps: VA.wilson_line_torus(eps) for eps in EPS}


def test_four_critical_points_two_saddles(census):
    for eps, w in census.items():
        assert w["counts"] == {"sink": 1, "saddle": 2, "source": 1} and w["morse_sum"] == 0 and w["all_at_half_periods"]
        types = {r["theta"]: r["type"] for r in w["critical_points"]}
        assert types[(0.0, 0.0)] == "sink" and types[(eps[0] / 2, eps[1] / 2)] == "source"
        assert all("saddle (1 in, 1 out)" in r["flow_type"] for r in w["critical_points"] if r["type"] == "saddle")


def test_values_by_two_routes(census):
    for eps, w in census.items():
        for r in w["critical_points"]:
            assert r["V"] == pytest.approx(_V(eps, r["theta"]), abs=1e-9)
    real = {r["theta"]: r["V"] for r in census[(0, 1)]["critical_points"]}
    assert real[(0.0, 0.0)] == pytest.approx(-3.178063, abs=1e-6) and real[(0.0, 0.5)] == pytest.approx(3.178063, abs=1e-6)
    assert abs(_V((0, 1), (0.5, 0.0))) < 1e-13 and abs(_V((0, 1), (0.5, 0.5))) < 1e-13          # exact on the disc: the pairing
    assert abs(real[(0.5, 0.0)]) < 1e-9 and abs(real[(0.5, 0.5)]) < 1e-9                          # the box route, to its truncation
    tau_loop = {r["theta"]: r["V"] for r in census[(1, 0)]["critical_points"]}
    assert tau_loop[(0.0, 0.5)] == pytest.approx(-2.813625, abs=1e-6) and tau_loop[(0.5, 0.5)] == pytest.approx(2.813625, abs=1e-6)
    assert abs(_V((1, 0), (0, 0)) - (_V((1, 0), (0, 0.5)) + _V((0, 1), (0, 0)))) < 1e-12           # V_(1,0)(0) = V_(1,0)(0,1/2) + V_(0,1)(0)
    assert sorted(r["V"] for r in census[(1, 0)]["critical_points"]) == pytest.approx(sorted(r["V"] for r in census[(1, 1)]["critical_points"]), abs=1e-9)


def test_pairing_is_exact():
    p = VA.saddle_pairing()
    assert p["re_tau_is_half"] and p["norm_preserved"] and p["parity_of_m_flipped"] and p["windings_checked"] > 1000


def test_curvatures_belong_to_the_metric(census):
    G = VA.kinetic_metric()
    x, y = TAU.real, TAU.imag
    assert np.allclose(G, [[abs(TAU) ** 2 / y, -x / y], [-x / y, 1 / y]], atol=1e-12) and np.linalg.det(G) == pytest.approx(1.0, abs=1e-12)
    w = census[(0, 1)]
    assert w["jet_curvatures_theta"] == pytest.approx([28.893, 166.935], abs=2e-3)                # PG-01's 28.9 and 166.9
    assert w["jet_curvatures_metric"] == pytest.approx([37.394, 128.985], abs=2e-3)
    assert np.prod(w["jet_curvatures_theta"]) == pytest.approx(np.prod(w["jet_curvatures_metric"]), rel=1e-9)   # det G = 1
    assert w["mass_ratio_metric"] == pytest.approx(1.8572, abs=1e-4) and w["mass_ratio_if_theta_coordinates"] == pytest.approx(2.4037, abs=1e-4)
    a, b = census[(1, 0)], census[(1, 1)]                                                          # mirror images under the real structure
    assert a["jet_curvatures_metric"] == pytest.approx(b["jet_curvatures_metric"], rel=1e-6)                 # equal; the box truncation is not mirror-symmetric
    assert a["jet_curvatures_metric"] == pytest.approx([60.834, 194.868], abs=2e-3)
    assert abs(a["jet_curvatures_theta"][0] - b["jet_curvatures_theta"][0]) > 10                    # the theta chart tells them apart


@pytest.mark.parametrize("eps", EPS)
def test_census_by_index_and_exclusion(eps, N=1024):
    """Route 2 for the census, without Newton's method. The winding count of grad V puts index +1 at the sink and the source and -1 at the
    other two half-periods, and nowhere else. Then a proof for the disc sum that these four are the only critical points: the Lipschitz bound
    excludes every cell but four small clusters, each within R < 0.02 of a half-period (critical because V(-theta) = V(theta)), and there
    ||H^-1|| L3 R < 1 makes theta -> theta - H^-1 grad V a contraction, so the cluster holds no second zero."""
    g1, g2, L2, L3, mnc = _fields(eps, N)
    a = np.arctan2(g2, g1)

    def wrap(d):
        return (d + np.pi) % (2 * np.pi) - np.pi

    a10, a01 = np.roll(a, -1, 0), np.roll(a, -1, 1)
    a11 = np.roll(a10, -1, 1)
    idx = np.rint((wrap(a10 - a) + wrap(a11 - a10) + wrap(a01 - a11) + wrap(a - a01)) / (2 * np.pi)).astype(int)
    assert idx.sum() == 0                                                                           # Poincare-Hopf, chi(T^2) = 0
    cells = {((j + 1) % N / N, (k + 1) % N / N): int(idx[j, k]) for j, k in zip(*np.nonzero(idx))}  # keyed by cell centre
    ends = {(0.0, 0.0), (eps[0] / 2, eps[1] / 2)}                                                   # sink and source
    assert cells == {hp: (1 if hp in ends else -1) for hp in HALF_PERIODS}
    gn = np.hypot(g1, g2)
    low = np.minimum.reduce([gn, np.roll(gn, -1, 0), np.roll(gn, -1, 1), np.roll(np.roll(gn, -1, 0), -1, 1)])
    j, k = np.nonzero(low <= L2 * math.sqrt(2) / N)                       # the cells the Lipschitz bound cannot exclude
    t = np.stack([(j + 1.0) / N, (k + 1.0) / N], 1) % 1                     # their centres
    h = np.round(2 * t) / 2 % 1                                             # and the nearest half-period
    r = np.hypot(*(((t - h + 0.5) % 1) - 0.5).T) + math.sqrt(2) / N          # bounds the distance from it to every point of the cell
    clusters = {tuple(x) for x in h}
    assert clusters == set(HALF_PERIODS)
    for hp in clusters:
        R = r[(h == hp).all(1)].max()
        assert R < 0.02 and np.linalg.norm(np.linalg.inv(_hess(mnc, hp)), 2) * L3 * R < 1


def test_metric_curvatures_route2(census):
    x, y = TAU.real, TAU.imag
    G = np.array([[abs(TAU) ** 2 / y, -x / y], [-x / y, 1 / y]])
    for eps in EPS:
        m, n, w8 = _disc(eps)
        ev = np.sort(np.linalg.eigvals(np.linalg.solve(G, _hess((m, n, -2 * w8), (0.0, 0.0)))).real)
        assert ev == pytest.approx(census[eps]["jet_curvatures_metric"], rel=1e-6)


def test_ledger_atlas_and_map():
    e = {x["id"]: x for x in L.E}
    assert e["PG-6"]["status"] == "active" and e["PG-6"]["module"] == "mtft.matter.vacuum_atlas" and e["PG-6"]["test"] == "test_pg01_wilson_torus.py"
    assert e["CC-66"]["status"] == "active" and "PG-4" in e["CC-66"]["statement"] and "PG-6" in e["CC-66"]["statement"]
    assert e["PG-4"]["status"] == "caveat" and "CC-66" in e["PG-4"]["note"] and "No saddle appears anywhere" in e["PG-4"]["statement"]   # append-only
    a = Z.assignment(); assert a["PG-6"] == "IV" and a["CC-66"] == "IV"
    assert {"PG-6", "CC-66"} <= set(Z.subs()["IV.7"][1].claims) and all(Z.check().values())
    A = VA.atlas()
    assert A["saddles"] == 2 and [r[3] for r in A["rows"] if "saddle at" in r[0]] == ["PG-6", "PG-6"]
