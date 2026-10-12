"""TR-01 (2026-10-10): time reversal in the programme — the two mirror classes of the modular tiling, the period lattices of the thirteen
star-quotient curves, the t -> -t parity of their L-functions, and the pre-registered strong-CP test for piece I.4.

Time reversal acts on a gauge coupling tau = theta/(2 pi) + 4 pi i/g^2 as tau -> -conj(tau).  Its fixed lines are theta = 0 (Re tau = 0) and
theta = pi (Re tau = 1/2).  Under the duality group SL(2, Z) each line sweeps out one class of mirrors of the modular (2, 3, infinity) tiling,
and the two classes never meet: an integer matrix of determinant -1 and trace 0 acts on the upper half-plane by tau -> (a conj(tau) + b)/(c conj(tau) + d)
as the reflection in the geodesic c(x^2 + y^2) - 2 a x - b = 0, whose endpoints (a -+ 1)/c are cusps (x = -b/2a and infinity when c = 0).  The
binary form (c, -2a, -b) has discriminant 4 and content gcd(c, 2a, b) = 2 or 1; equivalently the two endpoint cusps p/q, r/s have |ps - qr| = 1
(Farey neighbours) or 2.  Content 2 is the class of theta = 0 (the Farey edges), content 1 the class of theta = pi (the Farey medians); these are
the two conjugacy classes of reflections in PGL(2, Z) — in the (2, 3, infinity) triangle group only the label 3 is odd, joining the reflections
in |tau| = 1 and Re tau = 1/2, while Re tau = 0 meets them at the labels 2 and infinity — so no duality carries theta = 0 to theta = pi.  The oloid move q -> -q/(1 + q) and the triangle mirror n -> -1 - n (the Holes and
Charges note, 2026-10-10) are both of the theta = pi class; so is the period lattice of 143a1.

A point of the upper half-plane is time-reversal symmetric when a mirror passes through it.  For a CM point this is exact: tau = (-B + sqrt D)/(2A)
lies on a mirror iff its reduced form (A, B, C) has B = 0 (class theta = 0, a rectangular lattice), |B| = A or A = C (class theta = pi, a rhombic
lattice); these are Gauss's ambiguous classes, the 2-torsion of the class group, the CM points with real j.  This is the level-1 picture: a
fixed point of an Atkin–Lehner involution on X0(143) lies over such a class, two points of X0(143) over each (W_13's four fixed points lie over
the two classes of discriminant -52, one of each mirror class), and the level-143 real structure is a finer question, open.

An elliptic curve over Q has a real period lattice: rectangular with two real circles when its discriminant is positive, rhombic with one real
circle (Re tau = 1/2 exactly) when it is negative.  Its completed L-function obeys Lambda(s) = epsilon Lambda(2 - s) with real coefficients, hence
Lambda(1 + it) = epsilon conj(Lambda(1 + it)) exactly: the root number is the parity of Lambda under t -> -t, and epsilon = -1 makes Lambda purely
imaginary and odd, forcing the zero at t = 0 (SQ-2's odd rank, restated).  PARI's values are consistent with this but assume the functional
equation, so they are not an independent check of it.

The thirteen star quotients X0(N)* of SQ-01 are identified as curves: in each isogeny class exactly one curve has modular degree 4 = |<W_p, W_q>|,
and it is the quotient.  A degree-4 Atkin–Lehner-invariant parametrisation of a curve E' in the class pulls the Neron differential back to
m times 2 pi i f(tau) dtau with m a nonzero integer (m = lambda c, lambda an integer from the Neron models and c the Manin constant, an integer by
Edixhoven and 1 here by Cremona), and its degree is m^2 times E''s modular degree as PARI defines it (degree over the Manin constant squared).
So 4 = m^2 d with d in {1, 4}; d = 1 does not occur in these classes (every curve has degree 2, 4 or 8), so d = 4 and m = 1.  Nine quotients have
one real circle (theta = pi) and four have two (82, 141, 142, 145).

Nothing here derives the level or a physical constant; the pre-registered test (STRONG_CP_RULE, sealed in mtft.matter.prereg as PR-S2; PR-S1,
its first form, rejected theta = pi points outright and was superseded the same day after independent review, CC-63) waits for a candidate.  Every count is asserted in tests/test_tr01_time_reversal.py; the curve table is reproduced from PARI's database when cypari2 is present."""
from fractions import Fraction
from functools import lru_cache
from math import gcd, sqrt

# ------------------------------------------------------------------------------------------------------------------------------------ reflections
THETA_0, THETA_PI = "theta_0", "theta_pi"


def _check(m):
    if len(m) != 2 or any(len(r) != 2 for r in m) or not all(isinstance(x, int) for r in m for x in r):
        raise ValueError("a reflection is a 2 x 2 integer matrix")
    (a, b), (c, d) = m
    if a * d - b * c != -1 or a + d != 0:
        raise ValueError(f"{m} is not a reflection of the upper half-plane (determinant -1 and trace 0 needed)")
    return a, b, c, d


def is_reflection(m):
    try:
        _check(m); return True
    except ValueError:
        return False


def reflect(m, tau):
    """The action tau -> (a conj(tau) + b)/(c conj(tau) + d) of an integer matrix of determinant -1 on the upper half-plane."""
    (a, b), (c, d) = m
    if a * d - b * c != -1:
        raise ValueError("determinant -1 needed")
    z = complex(tau).conjugate()
    return (a * z + b) / (c * z + d)


def fixed_form(m):
    """Binary form (c, -2a, -b) whose zeros are the endpoints of the mirror c(x^2 + y^2) - 2 a x - b = 0; its discriminant is 4."""
    a, b, c, d = _check(m)
    return (c, -2 * a, -b)


def content(form):
    g = 0
    for x in form:
        g = gcd(g, abs(x))
    return g


def mirror_endpoints(m):
    """The mirror's endpoints (a - 1)/c and (a + 1)/c as cusps: Fractions, or None for infinity."""
    a, b, c, d = _check(m)
    if c == 0:
        return (Fraction(-b, 2 * a), None)                   # the vertical line x = -b/(2a)
    return (Fraction(a - 1, c), Fraction(a + 1, c))


def farey_det(m):
    """|p1 q2 - p2 q1| for the endpoint cusps p1/q1, p2/q2 (infinity = 1/0): 1 for Farey neighbours, 2 for the medians."""
    e1, e2 = mirror_endpoints(m)
    p1, q1 = (e1.numerator, e1.denominator) if e1 is not None else (1, 0)
    p2, q2 = (e2.numerator, e2.denominator) if e2 is not None else (1, 0)
    return abs(p1 * q2 - p2 * q1)


def mirror_class(m, route="content"):
    """theta_0 (content 2, Farey edge: the duality images of theta = 0) or theta_pi (content 1, Farey median: the images of theta = pi)."""
    if route == "content":
        return THETA_0 if content(fixed_form(m)) == 2 else THETA_PI
    if route == "farey":
        return THETA_0 if farey_det(m) == 1 else THETA_PI
    raise ValueError("route is 'content' or 'farey'")


def conjugate(m, g):
    """g m g^-1 for g in GL(2, Z)."""
    (p, q), (r, s) = g
    det = p * s - q * r
    if det not in (1, -1):
        raise ValueError("g must be in GL(2, Z)")
    ginv = ((s * det, -q * det), (-r * det, p * det))
    def mul(x, y):
        return ((x[0][0] * y[0][0] + x[0][1] * y[1][0], x[0][0] * y[0][1] + x[0][1] * y[1][1]),
                (x[1][0] * y[0][0] + x[1][1] * y[1][0], x[1][0] * y[0][1] + x[1][1] * y[1][1]))
    return mul(mul(g, m), ginv)


def reflections(bound):
    """Every reflection with entries of absolute value at most bound."""
    out = []
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            for c in range(-bound, bound + 1):
                if a * a + b * c == 1:
                    out.append(((a, b), (c, -a)))
    return out


STANDARD = {
    "Re tau = 0 (theta = 0)": ((1, 0), (0, -1)),
    "Re tau = 1/2 (theta = pi)": ((-1, 1), (0, 1)),
    "|tau| = 1": ((0, 1), (1, 0)),
    "Re tau = -1/2 (the triangle mirror n -> -1 - n)": ((-1, -1), (0, 1)),
    "|tau + 1| = 1 (the oloid move q -> -q/(1 + q))": ((-1, 0), (1, 1)),
}


def physics_dictionary():
    return {"coupling": "tau = theta/(2 pi) + 4 pi i/g^2", "time_reversal": "tau -> -conj(tau)",
            "fixed_lines": {"theta = 0": "Re tau = 0", "theta = pi": "Re tau = 1/2"},
            "classes": {THETA_0: "SL(2,Z) images of theta = 0: the Farey edges (form content 2)",
                        THETA_PI: "SL(2,Z) images of theta = pi: the Farey medians (form content 1)"},
            "no_duality_between_them": "the two conjugacy classes of reflections in PGL(2, Z)"}


# ------------------------------------------------------------------------------------------------------------------------------------ CM points
def _form_check(form):
    A, B, C = form
    D = B * B - 4 * A * C
    if A <= 0 or D >= 0:
        raise ValueError("a positive definite form (A > 0, B^2 - 4AC < 0) is needed")
    return A, B, C, D


def cm_point(form):
    A, B, C, D = _form_check(form)
    return complex(-B / (2 * A), sqrt(-D) / (2 * A))


def reduce_form(form):
    """Gauss reduction of a positive definite form: the reduced form (|B| <= A <= C, B >= 0 if |B| = A or A = C) and g in SL(2, Z) with tau_red = g tau."""
    A, B, C, D = _form_check(form)
    g = ((1, 0), (0, 1))
    def act(h):
        nonlocal g
        g = ((h[0][0] * g[0][0] + h[0][1] * g[1][0], h[0][0] * g[0][1] + h[0][1] * g[1][1]),
             (h[1][0] * g[0][0] + h[1][1] * g[1][0], h[1][0] * g[0][1] + h[1][1] * g[1][1]))
    for _ in range(10_000):
        if not (-A < B <= A):                                   # substitute x -> x + k y: B -> B + 2Ak, and the point moves to tau - k
            k = -((B + A - 1) // (2 * A))
            B, C = B + 2 * A * k, A * k * k + B * k + C
            act(((1, -k), (0, 1)))
            continue
        if A > C:                                               # invert: tau -> -1/tau swaps A and C and negates B
            A, B, C = C, -B, A
            act(((0, -1), (1, 0)))
            continue
        if A == C and B < 0:
            A, B, C = C, -B, A
            act(((0, -1), (1, 0)))
            continue
        return (A, B, C), g
    raise RuntimeError("reduction did not terminate")


def cm_mirror_class(form):
    """Which mirror class passes through the CM point: 'theta_0', 'theta_pi', 'both' (only tau = i) or None (no mirror: j is not real)."""
    (A, B, C), _ = reduce_form(form)
    on0, onpi = B == 0, (abs(B) == A or A == C)
    if on0 and onpi:
        return "both"
    return THETA_0 if on0 else THETA_PI if onpi else None


def reduced_forms(D):
    """The primitive reduced forms of discriminant D < 0 (one per class)."""
    if D >= 0 or D % 4 not in (0, 1):
        raise ValueError("D must be a negative discriminant")
    out = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (c == a and b < 0) or gcd(gcd(a, abs(b)), c) != 1:
                continue
            out.append((a, b, c))
        a += 1
    return out


def ambiguous_classes(D):
    """The reduced forms on a mirror (the 2-torsion classes, real j), each with its class."""
    return [(f, cm_mirror_class(f)) for f in reduced_forms(D) if cm_mirror_class(f)]


def reflections_fixing(form, bound=None):
    """Route 2 for cm_mirror_class: the reflections (a, b; c, -a) with c C + a B - b A = 0 and a^2 + b c = 1, by bounded search."""
    A, B, C, D = _form_check(form)
    bound = bound or 4 * (A + abs(B) + C)
    out = []
    for a in range(-bound, bound + 1):
        for c in range(-bound, bound + 1):
            if (c * C + a * B) % A:
                continue
            b = (c * C + a * B) // A
            if a * a + b * c == 1:
                out.append(((a, b), (c, -a)))
    return out


def fricke_fixed_point_classes():
    """The W_143 fixed points lie over the classes of discriminant -143 and -572 and W_13's four over those of -52 (SQ-01): which classes are real on X(1), and in which mirror class."""
    return {D: {"classes": len(reduced_forms(D)), "real": ambiguous_classes(D)} for D in (-143, -572, -52, -44, -11)}


# ------------------------------------------------------------------------------------------------------------------------------------ curves
def discriminant(a):
    a1, a2, a3, a4, a6 = a
    b2, b4, b6 = a1 * a1 + 4 * a2, 2 * a4 + a1 * a3, a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6


def real_circles(a):
    """E(R) is a circle (discriminant < 0) or two (discriminant > 0)."""
    return 2 if discriminant(a) > 0 else 1


def lattice_class(a):
    """Rectangular lattice (two real circles, the theta = 0 mirror Re tau = 0) or rhombic (one circle, Re tau = 1/2: the theta = pi mirror)."""
    return THETA_0 if discriminant(a) > 0 else THETA_PI


def l_parity(root_number):
    """Lambda(1 + it) = epsilon conj(Lambda(1 + it)): even and real for epsilon = +1, odd and purely imaginary for epsilon = -1."""
    if root_number not in (1, -1):
        raise ValueError("root number +-1")
    return "even (real on the critical line)" if root_number == 1 else "odd (purely imaginary on the critical line; a zero at t = 0 is forced)"


# PARI/GP 2.17.2 with elldata (scripts/time_reversal_curves.py): the (+,+) isogeny class at each torus level of SQ-01 and the curve of modular degree 4 in it.
STAR_QUOTIENT_CURVES = {
    57: {"curve": "57a1", "a": (0, -1, 1, -2, 2), "class_moddeg": {"57a1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": -2},
    58: {"curve": "58a1", "a": (1, -1, 0, -1, 1), "class_moddeg": {"58a1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": None},
    65: {"curve": "65a2", "a": (1, 0, 0, 4, 1), "class_moddeg": {"65a1": 2, "65a2": 4}, "torsion": 2, "root_number": -1, "rank": 1, "a2": -1},
    77: {"curve": "77a1", "a": (0, 0, 1, 2, 0), "class_moddeg": {"77a1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": 0},
    82: {"curve": "82a1", "a": (1, 0, 1, -2, 0), "class_moddeg": {"82a1": 4, "82a2": 8}, "torsion": 2, "root_number": -1, "rank": 1, "a2": None},
    91: {"curve": "91a1", "a": (0, 0, 1, 1, 0), "class_moddeg": {"91a1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": -2},
    118: {"curve": "118a1", "a": (1, 1, 0, 1, 1), "class_moddeg": {"118a1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": None},
    123: {"curve": "123b1", "a": (0, -1, 1, 1, -1), "class_moddeg": {"123b1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": 0},
    141: {"curve": "141d1", "a": (0, -1, 1, -1, 0), "class_moddeg": {"141d1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": 0},
    142: {"curve": "142b1", "a": (1, 1, 0, -1, -1), "class_moddeg": {"142b1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": None},
    143: {"curve": "143a1", "a": (0, -1, 1, -1, -2), "class_moddeg": {"143a1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": 0},
    145: {"curve": "145a1", "a": (1, -1, 1, -3, 2), "class_moddeg": {"145a1": 4, "145a2": 8}, "torsion": 2, "root_number": -1, "rank": 1, "a2": -1},
    155: {"curve": "155c1", "a": (0, -1, 1, -1, 1), "class_moddeg": {"155c1": 4}, "torsion": 1, "root_number": -1, "rank": 1, "a2": 0},
}
COMPARISON_CURVES = {"11a1": {"a": (0, -1, 1, -10, -20), "root_number": 1, "rank": 0}, "37a1": {"a": (0, 0, 1, -1, 0), "root_number": -1, "rank": 1}}


def curve_table():
    rows = {}
    for N, c in STAR_QUOTIENT_CURVES.items():
        rows[N] = dict(c, discriminant=discriminant(c["a"]), real_circles=real_circles(c["a"]), lattice=lattice_class(c["a"]),
                       parity=l_parity(c["root_number"]), supersingular_at_2=(c["a2"] == 0))
    return rows


def circles_census():
    t = curve_table()
    return {"one_circle": [N for N, r in t.items() if r["real_circles"] == 1], "two_circles": [N for N, r in t.items() if r["real_circles"] == 2],
            "all_odd": all(r["root_number"] == -1 for r in t.values()),
            "moddeg_4_unique": all(list(r["class_moddeg"].values()).count(4) == 1 for r in t.values())}


# ------------------------------------------------------------------------------------------------------------------------------------ the fold
def fold_preimages(value, r=4.0):
    """The two preimages of a value under F(n) = 2 r T(n) = r n (n + 1): mirror images under n -> -1 - n, so the iteration cannot be run backwards."""
    disc = 1 + 4 * value / r
    if disc < 0:
        return ()
    s = sqrt(disc)
    return ((-1 - s) / 2, (-1 + s) / 2)


# ------------------------------------------------------------------------------------------------------------------------------------ strong CP (I.4)
STRONG_CP_RULE = (
    "Pre-registered test for I.4 (sealed as PR-S2 on 2026-10-10, before any candidate exists; PR-S1, sealed earlier the same day, rejected theta = pi "
    "points outright and was superseded after independent review). The effective angle is theta-bar = theta + arg det(M_u M_d). A lift that fixes "
    "the strong coupling at an arithmetic point tau_s (a CM point given by its form) passes without an axion if and only if tau_s lies on a mirror of "
    "the modular tiling and the mirror's class matches the quark-mass phase: class theta = 0 with arg det = 0 mod 2 pi, or class theta = pi with "
    "arg det = pi mod 2 pi (a point on both classes, such as tau = i, passes with either phase). Any other combination gives theta-bar = pi, a "
    "time-reversal-symmetric vacuum that is not ours (QCD at theta-bar = pi is theta-bar = 0 with a negative light-quark mass, excluded by meson "
    "spectroscopy and the lattice, not by the neutron EDM, which vanishes there): FAIL. A tau_s on no mirror has no antiunitary symmetry in any "
    "duality frame: FAIL unless an axion relaxes theta-bar. A cusp is the g -> 0 limit, not a coupling. A numerical candidate is classified to the "
    "stated tolerance and cannot pass exactly. The verdict is read off; nothing is tuned after the candidate is seen.")
PI = 3.141592653589793


def strong_cp_test(candidate, det_phase=None, tolerance=None):
    """candidate: ('cm', (A, B, C)) | ('numeric', complex) | ('cusp', anything).  det_phase: arg det(M_u M_d) in radians, or None if unknown."""
    kind, value = candidate
    if kind == "cusp":
        return {"mirror_class": "both", "exact": True, "det_phase_ok": None, "theta_bar": None, "verdict": "not a candidate: a cusp is the g -> 0 limit, not a coupling"}
    if kind == "cm":
        cls, exact = cm_mirror_class(value), True
    elif kind == "numeric":
        if tolerance is None:
            raise ValueError("a numerical candidate needs a tolerance")
        z = complex(value)
        for _ in range(10_000):
            z = complex(z.real - round(z.real), z.imag)
            if abs(z) < 1 - 1e-15:
                z = -1 / z
            else:
                break
        on0 = abs(z.real) < tolerance
        onpi = abs(abs(z.real) - 0.5) < tolerance or abs(abs(z) - 1) < tolerance
        cls, exact = ("both" if on0 and onpi else THETA_0 if on0 else THETA_PI if onpi else None), False
    else:
        raise ValueError("candidate kind is 'cm', 'numeric' or 'cusp'")
    angles = {THETA_0: (0.0,), THETA_PI: (PI,), "both": (0.0, PI), None: ()}[cls]
    if cls is None:
        return {"mirror_class": None, "exact": exact, "det_phase_ok": None, "theta_bar": None, "verdict": "FAIL: off every mirror (needs an axion)"}
    if det_phase is None:
        need = " or ".join("0" if a == 0 else "pi" for a in angles)
        return {"mirror_class": cls, "exact": exact, "det_phase_ok": None, "theta_bar": None, "verdict": f"pending: class {cls} needs arg det(M_u M_d) = {need}"}
    tol = tolerance or 1e-12
    wrap = lambda x: (x + PI) % (2 * PI) - PI
    bars = [wrap(a + det_phase) for a in angles]
    ok = any(abs(b) < tol for b in bars)
    theta_bar = min(bars, key=abs)
    if ok:
        verdict = "PASS" if exact else "pass to tolerance only"
    else:
        verdict = "FAIL: theta-bar = pi" if abs(abs(theta_bar) - PI) < tol else f"FAIL: theta-bar = {theta_bar:.4g}"
    return {"mirror_class": cls, "exact": exact, "det_phase_ok": ok, "theta_bar": theta_bar, "verdict": verdict}


# ------------------------------------------------------------------------------------------------------------------------------------ report
@lru_cache(maxsize=None)
def _sweep(bound):
    refl = reflections(bound)
    return {"reflections": len(refl), "routes_agree": all(mirror_class(m) == mirror_class(m, "farey") for m in refl),
            "theta_0": sum(mirror_class(m) == THETA_0 for m in refl), "theta_pi": sum(mirror_class(m) == THETA_PI for m in refl)}


def report(bound=12):
    return {"dictionary": physics_dictionary(), "standard": {k: mirror_class(m) for k, m in STANDARD.items()}, "sweep": _sweep(bound),
            "fricke_points": fricke_fixed_point_classes(), "curves": curve_table(), "census": circles_census(), "rule": STRONG_CP_RULE}


def text():
    r = report()
    L = ["TR-01 — time reversal: tau -> -conj(tau) on tau = theta/(2 pi) + 4 pi i/g^2; theta = 0 and theta = pi are the two mirror classes of the tiling", ""]
    for k, c in r["standard"].items():
        L.append(f"  {k:52s} {c}")
    s = r["sweep"]
    L.append(f"  all {s['reflections']} matrices with entries <= 12 ({s['reflections'] // 2} reflections, M and -M acting alike): content and Farey routes agree ({s['routes_agree']}); theta_0 {s['theta_0']}, theta_pi {s['theta_pi']}")
    L.append("")
    L.append("CM points on a mirror (ambiguous classes, real j):")
    for D, v in r["fricke_points"].items():
        L.append(f"  D = {D}: {v['classes']} classes, real: " + ", ".join(f"{f} {c}" for f, c in v["real"]))
    L.append("")
    L.append("The thirteen star quotients as curves (modular degree 4 in the class):")
    for N, c in r["curves"].items():
        L.append(f"  {N:4d} {c['curve']:6s} a = {c['a']}  disc {c['discriminant']:7d}  circles {c['real_circles']}  {c['lattice']:8s}  eps {c['root_number']:+d}  rank {c['rank']}")
    cz = r["census"]
    L.append(f"  one circle: {cz['one_circle']}; two: {cz['two_circles']}; all root numbers -1 (Lambda odd under t -> -t): {cz['all_odd']}")
    L.append("")
    L.append(r["rule"])
    return "\n".join(L)
