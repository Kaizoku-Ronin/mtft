r"""mtft.plane.kernel — one real 2 × 2 matrix read eight ways (EXACT structure; v1.1.0, 2026-11-16).

The Modular Plane lab reads every matrix it meets on one strip, eight ways.  This module is that
strip as functions.  Write

    M = αI + s₁S₁ + s₂S₂ + wJ,     S₁ = diag(1, −1),  S₂ = [[0, 1], [1, 0]],  J = [[0, −1], [1, 0]],

so α = tr M / 2 is the dilation, s = √(s₁² + s₂²) the strain, w the rotation, and
β = s² − w² = Δ/4 with Δ = tr² − 4 det.  The sign of Δ is the *kind*: elliptic (Δ < 0), parabolic
(Δ = 0), hyperbolic (Δ > 0).  The eight readings (``READINGS``, in the lab's order):

  motion   τ ↦ (aτ + b)/(cτ + d): a rotation about a point, a translation to a cusp, or a boost
           along a geodesic, by 2η with η = arccosh(|tr| / 2√det)
  form     its axis as a binary quadratic form (c, d − a, −b) of discriminant Δ, and that form
           as a 2 + 1 momentum: timelike, lightlike or spacelike
  number   the number a + εb it is in ℝ[M] = ℝ[ε]/(ε² = κ): complex, dual or split-complex, with
           norm det M and an angle that winds, is a slope, or has a horizon
  flow     as a deformation rate: which of rotation, shear and strain wins, and its place on the
           (α, β) chart (centre, focus, saddle, node)
  field    the mode equation ẍ = κx it generates: wave, ultralocal or Laplace; as a Bogoliubov
           map a phase, free evolution or a squeeze
  vacuum   for det M = 1: the distance d it moves the vacuum i, and the quanta N = sinh²(d/2)
  two      two levels: the gap 2s of its symmetric part and the rotation w
  range    its numerical range W(M) = {⟨x, Mx⟩ : |x| = 1}: an ellipse with the eigenvalues as foci
           (the elliptical range theorem); numerical radius, spectral norm, their ratio ≤ 2 (Crouzeix)

Two identities tie the last three readings together (EXACT; gated on random matrices against a
route with no formula).  For every real M the numerical range is the ellipse centred at α with
semi-axis s along the real axis and |w| along the imaginary one: the seventh and eighth readings
are the same two numbers.  The spectral norm is |α + iw| + s.  And for det M = 1, s² = sinh²(d/2):
the strain squared is the number of quanta the sixth reading counts.

``read(M)`` reads a matrix as a motion; ``read(M, generator=True)`` reads it as a generator L,
whose motion is exp(tL) (``expm``).  Entries may be integers, Fractions or floats; with exact
entries the kind and every rational quantity are exact, and square roots are floats.  For an
integer matrix of determinant 1, ``read`` adds Rademacher's Ψ and, given a level, the cycle type
on P¹(ℤ/N), both through ``mtft.lines``.

Provenance: the Modular Plane lab's core (version 8.3), the strip of Room 0; the numerical range
was pre-registered in the record MP06 (8 November 2026).  The lab's page and glossary are teaching
material and are not part of this package.
"""
from __future__ import annotations

import math
from fractions import Fraction as Fr

READINGS = ("motion", "form", "number", "flow", "field", "vacuum", "two", "range")
LABELS = {"motion": "a motion of τ", "form": "a form, a momentum", "number": "a number a + εb", "flow": "a deformation, a flow",
          "field": "a field equation", "vacuum": "one mode's vacuum", "two": "two levels", "range": "its numerical range"}
KIND_NAME = {"ell": "elliptic", "par": "parabolic", "hyp": "hyperbolic"}
WORDS = {
    "ell": {"number": "complex, ε² = −1", "mode": "rotation", "equation": "the wave equation", "vacuum": "the vacuum stands still",
            "flow": "rotation wins (β < 0)", "time": "real time, a clock", "bogoliubov": "a phase on a mode", "momentum": "timelike"},
    "par": {"number": "dual, ε² = 0", "mode": "shear", "equation": "the ultralocal equation", "vacuum": "a horocycle",
            "flow": "a shear (β = 0)", "time": "the degenerate middle", "bogoliubov": "free evolution of a mode", "momentum": "lightlike"},
    "hyp": {"number": "split-complex, ε² = +1", "mode": "boost", "equation": "Laplace's equation", "vacuum": "a geodesic",
            "flow": "strain wins (β > 0)", "time": "modular time, a thermometer", "bogoliubov": "a squeeze of a mode", "momentum": "spacelike"},
}
"""The lab's words for the three kinds, reading by reading."""

_ORDER = {1: 6, 0: 4, -1: 3}          # the order in SL(2, ℤ) of an elliptic integer matrix, by its trace


def _entries(M):
    (a, b), (c, d) = M
    exact = all(isinstance(v, (int, Fr)) and not isinstance(v, bool) for v in (a, b, c, d))
    if exact:
        return Fr(a), Fr(b), Fr(c), Fr(d), True
    return float(a), float(b), float(c), float(d), False


def _tol(a, b, c, d, exact):
    return 0 if exact else 1e-12 * max(1.0, a * a + b * b + c * c + d * d)


def _num(v):
    """An exact value as an int where it is one, else as it is."""
    if isinstance(v, Fr) and v.denominator == 1:
        return int(v)
    return v


def decompose(M) -> dict:
    """α, s₁, s₂, w, the strain s and β = s₁² + s₂² − w² of M = αI + s₁S₁ + s₂S₂ + wJ."""
    a, b, c, d, exact = _entries(M)
    alpha, s1, s2, w = (a + d) / 2, (a - d) / 2, (b + c) / 2, (c - b) / 2
    strain_sq = s1 * s1 + s2 * s2
    return {"alpha": _num(alpha), "s1": _num(s1), "s2": _num(s2), "w": _num(w), "strain_sq": _num(strain_sq),
            "strain": math.sqrt(strain_sq), "beta": _num(strain_sq - w * w), "exact": exact}


def invariants(M) -> dict:
    """Trace, determinant, discriminant Δ = tr² − 4 det, the kind and the eigenvalues."""
    a, b, c, d, exact = _entries(M)
    t, dt = a + d, a * d - b * c
    D = t * t - 4 * dt
    tol = _tol(a, b, c, d, exact)
    if abs(D) <= tol:
        D = D * 0
    kind = "ell" if D < 0 else ("hyp" if D > 0 else "par")
    scalar = abs(b) <= tol and abs(c) <= tol and abs(a - d) <= tol
    root = math.sqrt(abs(D))
    eig = ((float(t) - root) / 2, (float(t) + root) / 2) if D >= 0 else (complex(float(t) / 2, -root / 2), complex(float(t) / 2, root / 2))
    integer = exact and all(v.denominator == 1 for v in (a, b, c, d))
    return {"tr": _num(t), "det": _num(dt), "disc": _num(D), "kind": kind, "scalar": bool(scalar), "eigenvalues": eig,
            "exact": exact, "integer": integer, "unimodular": bool(integer and dt == 1)}


def kind(M) -> str:
    """'ell', 'par' or 'hyp' by the sign of Δ."""
    return invariants(M)["kind"]


def region(alpha, beta, tol: float = 1e-6) -> str:
    """The region of the (α, β) chart for the flow ẋ = Lx with eigenvalues α ± √β."""
    z = abs(alpha) < tol
    if beta < -tol:
        return "centre" if z else ("attracting focus" if alpha < 0 else "repelling focus")
    if beta > tol:
        if abs(beta - alpha * alpha) < tol:
            return "det L = 0: one direction frozen"
        return "saddle" if beta > alpha * alpha else ("attracting node" if alpha < 0 else "repelling node")
    return "pure shear" if z else ("attracting star or improper node" if alpha < 0 else "repelling star or improper node")


def expm(L, t: float):
    """exp(tL) in closed form: e^{αt}(c·I + s·X) with X = L − αI, (c, s) = (cosh, sinh/√β), (cos, sin/√−β) or (1, t)."""
    dec = decompose(L)
    alpha, beta = float(dec["alpha"]), float(dec["beta"])
    X = ((float(dec["s1"]), float(dec["s2"]) - float(dec["w"])), (float(dec["s2"]) + float(dec["w"]), -float(dec["s1"])))
    a, b, c, d, exact = _entries(L)
    tol = _tol(a, b, c, d, exact) if not exact else 0.0
    e = math.exp(alpha * t)
    if beta < -tol:
        om = math.sqrt(-beta)
        co, si = math.cos(om * t), math.sin(om * t) / om
    elif beta > tol:
        om = math.sqrt(beta)
        co, si = math.cosh(om * t), math.sinh(om * t) / om
    else:
        co, si = 1.0, t
    return ((e * (co + si * X[0][0]), e * si * X[0][1]), (e * si * X[1][0], e * (co + si * X[1][1])))


def fixed_points(M) -> dict:
    """The fixed points of τ ↦ (aτ + b)/(cτ + d) on the closed upper half-plane."""
    a, b, c, d, exact = _entries(M)
    inv = invariants(M)
    D, tol = inv["disc"], _tol(a, b, c, d, exact)
    if inv["scalar"]:
        return {"kind": "par", "every_point": True}
    if abs(c) <= tol:
        if inv["kind"] == "par":
            return {"kind": "par", "cusp": math.inf}
        return {"kind": "hyp", "ends": (float(b / (d - a)), math.inf)}
    root = math.sqrt(abs(D))
    if inv["kind"] == "ell":
        return {"kind": "ell", "tau": (float((a - d) / (2 * c)), root / (2 * abs(float(c))))}
    if inv["kind"] == "hyp":
        e1, e2 = (float(a - d) - root) / (2 * float(c)), (float(a - d) + root) / (2 * float(c))
        return {"kind": "hyp", "ends": (min(e1, e2), max(e1, e2))}
    return {"kind": "par", "cusp": _num((a - d) / (2 * c)) if exact else float((a - d) / (2 * c))}


def numerical_radius(alpha: float, p: float, q: float) -> float:
    """The farthest point from 0 of the ellipse centred at α with semi-axis p along the reals and q along the imaginaries."""
    alpha, p, q = abs(float(alpha)), abs(float(p)), abs(float(q))
    if q * q - p * p <= alpha * p:
        return alpha + p
    return q * math.sqrt(1.0 + alpha * alpha / (q * q - p * p))


def numerical_range(M) -> dict:
    """W(M) for a real 2 × 2 matrix: the ellipse centred at α with semi-axes s (real) and |w| (imaginary)."""
    a, b, c, d, exact = _entries(M)
    dec, inv = decompose(M), invariants(M)
    tol = _tol(a, b, c, d, exact)
    D = inv["disc"]
    p_sq, q_sq = dec["strain_sq"], dec["w"] * dec["w"]
    if abs(p_sq) <= tol:
        p_sq = p_sq * 0
    if abs(q_sq) <= tol:
        q_sq = q_sq * 0
    p, q = math.sqrt(p_sq), math.sqrt(q_sq)
    real = D >= 0
    minor_sq = q_sq if real else p_sq                       # the squared semi-axis across the line of the foci
    focal = math.sqrt(abs(D)) / 2
    alpha = float(dec["alpha"])
    foci = ((alpha - focal, 0.0), (alpha + focal, 0.0)) if real else ((alpha, -focal), (alpha, focal))
    shape = ("point" if D == 0 else "segment") if minor_sq == 0 else ("disk" if D == 0 else "ellipse")
    radius = numerical_radius(alpha, p, q)
    norm = math.hypot(alpha, float(dec["w"])) + p           # the spectral norm: |α + iw| + s
    return {"shape": shape, "centre": dec["alpha"], "foci": foci, "focal": focal, "real_foci": bool(real),
            "semi_real": p, "semi_imag": q, "semi_real_sq": _num(p_sq), "semi_imag_sq": _num(q_sq),
            "semi_major": max(p, q), "semi_minor": min(p, q), "minor_axis": 2 * min(p, q),
            "radius": radius, "norm": norm, "ratio": (norm / radius) if radius > 1e-300 else 1.0,
            "departure": _num(4 * minor_sq), "normal": minor_sq == 0}


def range_boundary(M, n: int = 72) -> list:
    """n + 1 points (x, y) round the boundary of W(M), for drawing."""
    R = numerical_range(M)
    c = float(R["centre"])
    return [(c + R["semi_real"] * math.cos(2 * math.pi * i / n), R["semi_imag"] * math.sin(2 * math.pi * i / n)) for i in range(n + 1)]


def jacobi_monodromy(M) -> dict:
    """Round the closed geodesic of a hyperbolic M of determinant 1: the Jacobi equation ξ″ = ξ has monodromy M², of trace 2 cosh L."""
    a, b, c, d, _ = _entries(M)
    inv = invariants(M)
    if inv["kind"] != "hyp" or abs(float(inv["det"]) - 1.0) > 1e-9:
        return {"kind": inv["kind"], "defined": False}
    lam = (abs(float(inv["tr"])) + math.sqrt(float(inv["disc"]))) / 2
    L = 2 * math.log(lam)
    M2 = ((a * a + b * c, a * b + b * d), (c * a + d * c, c * b + d * d))
    return {"kind": "hyp", "defined": True, "length": L, "exp_length": lam * lam, "two_cosh_length": 2 * math.cosh(L),
            "trace_M2": _num(M2[0][0] + M2[1][1]), "M2": tuple(tuple(_num(v) for v in row) for row in M2)}


def _vacuum(a, b, c, d, dt, exact):
    if (dt != 1) if exact else (abs(dt - 1.0) > 1e-9):
        return {"label": LABELS["vacuum"], "symplectic": False}
    nf = a * a + b * b + c * c + d * d
    quanta = (nf - 2) / 4
    den = c * c + d * d
    return {"label": LABELS["vacuum"], "symplectic": True, "image_of_i": (float((a * c + b * d) / den), float(1 / den)),
            "distance": math.acosh(max(1.0, float(nf) / 2)), "quanta": _num(quanta)}


def _number(inv, dec):
    D, k = inv["disc"], inv["kind"]
    a = dec["alpha"]
    b = math.sqrt(abs(D)) / 2
    out = {"label": LABELS["number"], "norm": inv["det"], "a": a, "b": b, "kappa": {"ell": -1, "par": 0, "hyp": 1}[k]}
    if inv["scalar"]:
        out.update(system="real", words="a real number: M is a scalar", angle=0.0)
        return out
    out.update(system={"ell": "complex", "par": "dual", "hyp": "split-complex"}[k], words=WORDS[k]["number"])
    af, dt = float(a), float(inv["det"])
    if k == "ell":
        out.update(angle=math.atan2(b, af), angle_is="an angle that winds: every + 2π is another log")
    elif k == "par":
        out.update(b=1.0, angle=(1.0 / af if af != 0 else math.inf), angle_is="a slope")
    elif dt > 0:
        out.update(angle=math.acosh(max(1.0, abs(af) / math.sqrt(dt))), angle_is="a rapidity, with a horizon at |b| = |a|")
    else:
        out.update(angle=None, angle_is="beyond the horizon: no log (the norm is not positive)")
    return out


def read(M, generator: bool = False, level: int | None = None) -> dict:
    """The eight readings of a real 2 × 2 matrix, as a motion of the plane or (``generator=True``) as a generator.

    Returns {'matrix', 'invariants', 'kind', 'kind_name', 'as', 'readings'}; ``readings`` has the
    eight keys of ``READINGS``.  ``level`` adds the cycle type of an integer matrix of determinant 1
    on P¹(ℤ/level) to the motion reading.
    """
    a, b, c, d, exact = _entries(M)
    inv, dec = invariants(M), decompose(M)
    k = inv["kind"]
    W = WORDS[k]
    dt, D = inv["det"], inv["disc"]
    R = numerical_range(M)
    two = {"label": LABELS["two"], "gap": 2 * dec["strain"], "gap_sq": _num(4 * dec["strain_sq"]), "rotation": dec["w"],
           "beta": dec["beta"], "levels": (-dec["strain"], dec["strain"])}
    rng = dict(R, label=LABELS["range"])
    number = _number(inv, dec)
    if generator:
        alpha, beta = float(dec["alpha"]), float(dec["beta"])
        growth = "" if abs(alpha) < 1e-12 else (", damped" if alpha < 0 else ", growing")
        readings = {
            "motion": {"label": LABELS["motion"], "mode": W["mode"], "words": f"exp(tL): a {W['mode']}{growth}",
                       "rate": math.sqrt(abs(beta)), "returns": k == "ell"},
            "form": {"label": LABELS["form"], "preserves": "p² − κx² with κ = β" + ("" if abs(alpha) < 1e-12 else " only at α = 0"),
                     "kappa": dec["beta"], "type": {"ell": "definite", "par": "degenerate", "hyp": "indefinite"}[k], "momentum": W["momentum"]},
            "number": dict(number, words=W["number"]),
            "flow": {"label": LABELS["flow"], "region": region(alpha, beta), "wins": W["flow"], "alpha": dec["alpha"], "beta": dec["beta"]},
            "field": {"label": LABELS["field"], "equation": W["equation"], "time": W["time"], "kappa": dec["beta"]},
            "vacuum": {"label": LABELS["vacuum"], "path": W["vacuum"], "on_the_plane": abs(alpha) < 1e-12},
            "two": two, "range": rng}
        return {"matrix": tuple(tuple(_num(v) for v in row) for row in ((a, b), (c, d))), "invariants": inv, "kind": k,
                "kind_name": KIND_NAME[k], "as": "generator", "decomposition": dec, "readings": readings}
    fix = fixed_points(M)
    motion = {"label": LABELS["motion"], "fixed": fix, "acts_on_the_plane": bool(dt > 0)}
    if inv["scalar"]:
        motion["words"] = "every τ fixed: a scalar matrix"
    elif k == "ell":
        motion["words"] = "a rotation about a point"
        motion["angle"] = 2 * math.acos(max(-1.0, min(1.0, float(inv["tr"]) / (2 * math.sqrt(float(dt))))))
    elif k == "par":
        motion["words"] = "a translation to a cusp"
    else:
        motion["words"] = "a boost along a geodesic"
        if dt > 0:
            eta = math.acosh(max(1.0, abs(float(inv["tr"])) / (2 * math.sqrt(float(dt)))))
            motion.update(rapidity=eta, length=2 * eta, multiplier=(abs(float(inv["tr"])) + math.sqrt(float(D))) / 2)
    if inv["unimodular"]:
        Mi = ((int(a), int(b)), (int(c), int(d)))
        if k == "ell":
            motion["order"] = _ORDER[int(inv["tr"])]
        elif k == "par":
            motion["order"] = 1 if Mi == ((1, 0), (0, 1)) else (2 if Mi == ((-1, 0), (0, -1)) else math.inf)
        else:
            from mtft.lines import lorenz as _lz
            from mtft.lines import rademacher as _rm
            motion["psi"] = _rm.Psi(Mi)
            motion["lorenz_word"] = _lz.lorenz_word(Mi)
        if level is not None:
            from mtft.lines import darts as _dt
            motion["level"] = int(level)
            motion["cycle_type"] = _dt.darts(int(level)).cycle_type(Mi)
    elif level is not None:
        raise ValueError("a cycle type needs an integer matrix of determinant 1")
    A, B, C = c, d - a, -b                                  # the axis: A τ² + B τ + C = 0
    form = {"label": LABELS["form"], "form": (_num(A), _num(B), _num(C)), "disc": D,
            "type": {"ell": "definite", "par": "degenerate", "hyp": "indefinite"}[k],
            "momentum": {"E": _num((A + C) / 2), "px": _num((A - C) / 2), "py": _num(B / 2), "kind": W["momentum"]}}
    if inv["integer"]:
        g = math.gcd(math.gcd(abs(int(A)), abs(int(B))), abs(int(C)))
        form["content"] = g
        if g:
            form["primitive"] = (int(A) // g, int(B) // g, int(C) // g)
            form["primitive_disc"] = int(D) // (g * g)
    readings = {
        "motion": motion, "form": form, "number": number,
        "flow": {"label": LABELS["flow"], "wins": W["flow"], "generated_by": {"ell": "a rotation (centre)", "par": "a shear", "hyp": "a strain (saddle)"}[k],
                 "alpha": dec["alpha"], "strain": dec["strain"], "rotation": dec["w"], "beta": dec["beta"]},
        "field": {"label": LABELS["field"], "bogoliubov": W["bogoliubov"], "equation": W["equation"], "time": W["time"]},
        "vacuum": _vacuum(a, b, c, d, dt, exact), "two": two, "range": rng}
    return {"matrix": tuple(tuple(_num(v) for v in row) for row in ((a, b), (c, d))), "invariants": inv, "kind": k,
            "kind_name": KIND_NAME[k], "as": "motion", "decomposition": dec, "readings": readings}


def lab_numerical_range(M) -> dict:
    """The lab core's own computation of W(M) (version 7), kept as the second route: foci from Δ, the minor axis from
    |b − c| or √((b − c)² + Δ), the numerical radius by a golden-section search of the support function, the norm from
    the Frobenius norm and the determinant."""
    (a, b), (c, d) = ((float(v) for v in row) for row in M)
    t, det, nf2 = a + d, a * d - b * c, a * a + b * b + c * c + d * d
    tol = 1e-12 * max(1.0, nf2)
    D = t * t - 4 * det
    if abs(D) < tol:
        D = 0.0
    cx, real, focal = t / 2, D >= 0, math.sqrt(abs(D)) / 2
    minor2 = (b - c) * (b - c) if real else (b - c) * (b - c) + D
    if minor2 < tol:
        minor2 = 0.0
    semi_minor = math.sqrt(minor2) / 2
    semi_major = math.sqrt(focal * focal + semi_minor * semi_minor)
    angle = 0.0 if real else math.pi / 2

    def h(phi):
        return cx * math.cos(phi) + math.sqrt(semi_major ** 2 * math.cos(phi - angle) ** 2 + semi_minor ** 2 * math.sin(phi - angle) ** 2)

    best, bphi = 0.0, 0.0
    for i in range(180):
        phi = 2 * math.pi * i / 180
        if h(phi) > best:
            best, bphi = h(phi), phi
    lo, hi, gr = bphi - 2 * math.pi / 180, bphi + 2 * math.pi / 180, 0.618033988749895
    for _ in range(60):
        m1, m2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
        if h(m1) > h(m2):
            hi = m2
        else:
            lo = m1
    w = h((lo + hi) / 2)
    norm2 = math.sqrt((nf2 + math.sqrt(max(0.0, nf2 * nf2 - 4 * det * det))) / 2)
    shape = ("point" if D == 0 else "segment") if minor2 == 0 else ("disk" if D == 0 else "ellipse")
    return {"shape": shape, "centre": cx, "focal": focal, "semi_major": semi_major, "semi_minor": semi_minor, "real_foci": real,
            "radius": w, "norm": norm2, "ratio": norm2 / w if w > 1e-300 else 1.0}
