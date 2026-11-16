r"""mtft.gravity.fields — the labels, the knotted light fields and the squared field, in space (MODEL; v1.1.0, 2026-11-16).

Conventions (the lab's MP08; α and β are the Bateman variables of Kedia, Bialynicki-Birula, Peralta-Salas and
Irvine, Phys. Rev. Lett. 111, 150404 (2013), who built these fields)::

    a = x − iy,   b = t − i − z,   d = r² − (t − i)²,
    α = 1 − 2ib/d,   β = 2a/d,          |α|² + |β|² = 1   (the label sphere),
    n = (4(b² − a²), −4i(a² + b²), 8ab),          ∇α × ∇β = n/d³,
    F_pq = E + iB = ∇α^p × ∇β^q = pq α^{p−1} β^{q−1} n/d³     (a null solution of Maxwell's equations),
    W    = n ⊗ n / d⁵                                          (the squared hopfion: a linearised vacuum field of type N),
    G    = √d · F₁₁                                            (its lines of force are the lines of Re G and Im G).

The labels are constant along the light rays, whose velocity is ``flow``.  Every field here is a weak wave on
flat space; nothing in this module is a solution of the full Einstein equations.

Lines are integrated by ``step``/``integrate``, an extrapolation of the modified midpoint rule (Gragg, Bulirsch and
Stoer) written here so that the package needs no optional dependency.  It is for smooth, non-stiff fields.
"""
from __future__ import annotations

import math

import numpy as np

SEQUENCE = (2, 4, 6, 8, 10, 12)


# ----------------------------------------------------------------------------------------------------------------- the fields
def ab(t, X):
    """a = x − iy and b = t − i − z at the points X[..., 3]."""
    X = np.asarray(X, float)
    return X[..., 0] - 1j * X[..., 1], t - 1j - X[..., 2]


def dd(t, X):
    """d = r² − (t − i)²."""
    X = np.asarray(X, float)
    return np.sum(X * X, -1) - (t - 1j) ** 2


def labels(t, X):
    """The labels (α, β) of the points X at the time t."""
    a, b = ab(t, X)
    d = dd(t, X)
    return 1 - 2j * b / d, 2 * a / d


def grad_labels(t, X):
    """(∇α, ∇β) at the points X, exactly: arrays [..., 3]."""
    X = np.asarray(X, float)
    a, b = ab(t, X)
    d = dd(t, X)[..., None]
    gd = 2 * X
    ga = np.array([1, -1j, 0]) + 0 * gd
    gb = np.array([0, 0, -1]) + 0j * gd
    return -2j * gb / d + 2j * b[..., None] * gd / d ** 2, 2 * ga / d - 2 * a[..., None] * gd / d ** 2


def nvec(t, X):
    """n = (4(b² − a²), −4i(a² + b²), 8ab), so that ∇α × ∇β = n/d³."""
    a, b = ab(t, X)
    return np.stack([4 * (b * b - a * a), -4j * (a * a + b * b), 8 * a * b], -1)


def F(t, X, p: int = 1, q: int = 1):
    """F_pq = E + iB = pq α^{p−1} β^{q−1} n/d³ (Kedia et al. 2013); (1, 1) is the hopfion."""
    al, be = labels(t, X)
    return (p * q * al ** (p - 1) * be ** (q - 1) / dd(t, X) ** 3)[..., None] * nvec(t, X)


def G(t, X):
    """G = √d · F₁₁, principal root (smooth for t ≥ 0).  The tendex lines of the squared hopfion are the lines of Re G and
    Im G: the light's own lines turned by the local duality angle arg(d)/2 (Swearngin et al., arXiv:1302.1431, eq. 29)."""
    d = dd(t, X)
    return (np.sqrt(d + 0j) / d ** 3)[..., None] * nvec(t, X)


def squared(t, X):
    """W = n ⊗ n / d⁵, the squared hopfion, as a complex 3 × 3 array (E = Re W, B = −Im W in ``mtft.gravity.tetrad``)."""
    n = nvec(t, X)
    return n[..., :, None] * n[..., None, :] / dd(t, X)[..., None, None] ** 5


def flow(t, X):
    """The velocity of the light rays, Re(i n × n̄)/(n · n̄): a unit vector, the common null direction of F_pq and of W."""
    n = nvec(t, X)
    return np.real(1j * np.cross(n, np.conj(n))) / np.real(np.sum(n * np.conj(n), -1))[..., None]


def chart(t, alpha, beta):
    """The point of space that carries the labels (α, β) at the time t (α ≠ 1, which is the label of infinity)."""
    alpha, beta = np.asarray(alpha, complex), np.asarray(beta, complex)
    w = 1 / (1 - alpha)
    tau = (t - w.imag) / w.real
    xy = beta * (1 + 1j * tau) * w
    return np.stack([xy.real, -xy.imag, t - tau], -1)


def unit(V):
    """V scaled to unit length along its last axis."""
    V = np.asarray(V)
    return V / np.linalg.norm(V, axis=-1)[..., None]


def line_field(t, kind: str = "B", p: int = 1, q: int = 1, squared_field: bool = False, label_speed: bool = True):
    """The right-hand side dX/ds for a line of force at the time t.

    ``kind`` is "B" (the imaginary part) or "E" (the real part); the field is F_pq, or G when ``squared_field``.
    With ``label_speed`` the parameter s is length on the label sphere, which keeps a line that runs far out in
    space at a bounded speed; otherwise s is length in space."""
    def f(X):
        V = G(t, X) if squared_field else F(t, X, p, q)
        u = unit(V.imag if kind == "B" else V.real)
        if not label_speed:
            return u
        ga, gb = grad_labels(t, X)
        return u / math.sqrt(abs(ga @ u) ** 2 + abs(gb @ u) ** 2)
    return f


# ----------------------------------------------------------------------------------------------------------------- lines
def _midpoint(f, y, H, n):
    h = H / n
    z0, z1 = y, y + h * f(y)
    for _ in range(n - 1):
        z0, z1 = z1, z0 + 2 * h * f(z1)
    return 0.5 * (z0 + z1 + h * f(z1))


def step(f, y, H, sequence=SEQUENCE):
    """One step of length H of y' = f(y) by extrapolating the modified midpoint rule; returns (y(H), error estimate)."""
    y = np.asarray(y, float)
    table = []
    for k, n in enumerate(sequence):
        row = [_midpoint(f, y, H, n)]
        for j in range(1, k + 1):
            row.append(row[j - 1] + (row[j - 1] - table[k - 1][j - 1]) / ((n / sequence[k - j]) ** 2 - 1))
        table.append(row)
    return table[-1][-1], float(np.max(np.abs(table[-1][-1] - table[-1][-2])))


def integrate(f, y0, length, h, event=None, max_events=None, min_s: float = 1e-9, sequence=SEQUENCE):
    """Integrate y' = f(y) from y0 over ``length`` in steps of h.

    Returns a dict: ``s`` and ``y`` (the points at the step ends), ``error`` (the largest one-step estimate) and, when
    ``event`` is a function g(y), ``events``: the (s, y) at which g crosses zero going upward, each located by solving
    for the length of the step that lands on it, so to the accuracy of the integrator itself.  A crossing before
    ``min_s`` is not counted: a line that starts on the surface has not crossed it."""
    y = np.asarray(y0, float)
    n = max(1, round(length / h))
    h = length / n
    S, Y, worst, events = [0.0], [y], 0.0, []
    g0 = event(y) if event else None
    for i in range(n):
        y1, err = step(f, y, h, sequence)
        worst = max(worst, err)
        if event:
            g1 = event(y1)
            if g0 < 0 <= g1:
                lo, hi, glo, ghi, side, last = 0.0, h, g0, g1, 0, None
                for _ in range(60):                                      # the Illinois rule: a bracketed secant
                    mid = (lo * ghi - hi * glo) / (ghi - glo)
                    ym, _ = step(f, y, mid, sequence)
                    gm = event(ym)
                    if gm == 0 or (last is not None and abs(mid - last) < 1e-14 * h):
                        break
                    last = mid
                    if gm < 0:
                        lo, glo = mid, gm
                        ghi = ghi / 2 if side < 0 else ghi
                        side = -1
                    else:
                        hi, ghi = mid, gm
                        glo = glo / 2 if side > 0 else glo
                        side = 1
                if i * h + mid > min_s:                                  # a start that lies on the surface is not a crossing
                    events.append((i * h + mid, ym))
                    if max_events and len(events) >= max_events:
                        S.append(i * h + mid)
                        Y.append(ym)
                        break
            g0 = g1
        y = y1
        S.append((i + 1) * h)
        Y.append(y)
    out = {"s": np.array(S), "y": np.array(Y), "error": worst}
    if event:
        out["events"] = events
    return out


def trace(t, X0, length, kind: str = "B", p: int = 1, q: int = 1, squared_field: bool = False, h: float = 0.05, **kw):
    """A line of force through X0 at the time t, over ``length`` on the label sphere (see ``line_field``, ``integrate``)."""
    return integrate(line_field(t, kind, p, q, squared_field), X0, length, h, **kw)


def gauss_linking(C1, C2) -> float:
    """Gauss's linking integral of two closed polygons [n, 3], by the midpoint rule."""
    C1, C2 = np.asarray(C1, float), np.asarray(C2, float)
    d1, d2 = np.roll(C1, -1, 0) - C1, np.roll(C2, -1, 0) - C2
    R = (C1 + d1 / 2)[:, None, :] - (C2 + d2 / 2)[None, :, :]
    tri = np.einsum("ijk,ijk->ij", R, np.cross(d1[:, None, :], d2[None, :, :]))
    return float(np.sum(tri / np.linalg.norm(R, axis=2) ** 3) / (4 * math.pi))
