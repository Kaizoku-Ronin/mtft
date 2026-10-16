"""PG-01 (2026-10-12): the vacuum as a dynamical system — the critical-point atlas of the potentials the programme has, the shape modulus's cosh well,
what Hubble friction does to each critical point, a labelled toy run, and the gate a big-bang-to-matter simulation would have to pass.

The eight types of three-dimensional critical point are the sign and reality patterns of a linearisation's eigenvalues: three real eigenvalues give the
source (all out), the sink (all in) and the two saddles (one in two out; two in one out); a complex pair with a real eigenvalue gives the four spirals
(attracting or repelling in the plane, in or out along the axis).  A potential's critical point is hyperbolic when no Hessian eigenvalue vanishes.

What the atlas finds (EXACT unless marked):
  1. the shape modulus r = A_X/A_E sits in a cosh well: with phi = log r, E(r)/E* = cosh(phi - phi*), floor at the self-dual point r* = sqrt(A/B)
     (PV-01's S1/S2: A = 19, B = 11), and the well's curvature at the floor is exactly 1 in units of the floor energy — a sink (the reading of phi as the
     canonical field is the natural log metric of a ratio modulus, not derived here);
  2. the volume has no finite critical point: with both coefficients positive, both terms of the Einstein-frame potential fall with Vol (PV-01), so
     the gradient flow runs to decompactification — the only critical point is at infinity, a sink there;
  3. the condensate's D-term potential depends on u = v^2 - vbar^2 only, so its Hessian in the field variables (v, vbar) at the minimum (v*, 0) has a
     zero eigenvalue along vbar (the flat family u = const): a non-hyperbolic line of critical points; a mass term m^2 (v^2 + vbar^2) with
     0 < m^2 < xi moves the minimum to v^2 = (8/13)(xi - m^2) and lifts that eigenvalue to exactly 4 m^2 — a sink (CO-1 made quantitative);
  4. the relative-U(1) Wilson line on the torus of flat connections has a two-dimensional sink at the jet vacuum theta = 0 and a two-dimensional source
     at the untwisted point theta_eta (the full Hessians, extending PA1-1's statement along the straight path), for all three even spin structures —
     EXACT in PA-01's product approximation to the potential, as PA1-1 is; on the whole torus two nondegenerate saddles at the other two
     half-periods complete the picture, as a Morse function on the torus requires (#min - #saddle + #max = chi(T^2) = 0), and in the real spin
     structure both sit at V = 0 exactly (PG-6, 2026-10-14; CC-66 caveats the first filing's "no saddle", which BL-3 had already contradicted);
     curvatures, and so the Wilson-line masses, are those of G^-1 H in the kinetic metric G of 143a1 (`wilson_line_torus`, `kinetic_metric`,
     `saddle_pairing`);
  5. the bosonic adjoints have positive one-loop masses (M05-1: sinks) and the one-loop Casimir energy has no minimum in A_X (MS-6: no critical point).
  Two saddles appear, both on the Wilson-line torus (PG-6; the first filing said none, CC-66), and no spiral: spirals need time, and with Hubble
  friction 3H a sink of mass m is a spiral iff m > 3H/2 (otherwise a node).  The type of every sink in a cosmological run is therefore set by
  m/H, which the programme has not fixed (STRUCTURAL).
  6. a toy run (DIAGNOSTIC, normalisation not derived): phi'' + 3H phi' + sinh(phi - phi*) = 0 with constant H and in a radiation era H = 1/(2t) shows
     the spiral capture of the shape at the self-dual point (18 floor crossings in a run of length 60/m at H = 0.2 m; the count grows with the run);
     the energy identity dE/dt = -3H phi'^2 <= 0 for E = phi'^2/2 + cosh(phi - phi*) - 1 is checked on every run (`energy_argument`), and capture
     follows from it when H is bounded, of bounded variation and not integrable — constant H and the radiation and matter eras, the cases run here
     — but not for every H > 0: an integrable H leaves energy behind and a fast-growing H can freeze the field off the floor (CC-67, 2026-10-16;
     `capture_conditions`); the "captured" flag is only a smoke test of the envelope; the volume rolls away alongside, like (2/3) log t at late
     times with constant H.
  7. the simulation gate: what a big-bang-to-matter run of this vacuum needs — a declared parent (I.1), a vacuum energy from the parent action with a
     volume stationary point (I.6; GS-1's fixed Kaehler moduli), the couplings as functions of the moduli (GC-01), a cosmological background (the
     package's Friedmann module is pheno-tier), and then standard physics at confinement and nucleosynthesis, to be used, not re-derived.
Tests: tests/test_pg01_energy.py, tests/test_pg01_wilson_torus.py, tests/test_pg01_corrections.py.  Report: docs/SM/PG01_ENERGY_TWO_CLOCKS.md
(section 5: the 2026-10-14 addendum; section 6: the 2026-10-16 addendum, CC-67 and CC-68)."""
from __future__ import annotations

import math

import numpy as np
import sympy as sp

# ------------------------------------------------------------------------------------------------------------------------------ the eight types
TYPES = {1: "source (all out)", 2: "sink (all in)", 3: "saddle, 1 in 2 out", 4: "saddle, 2 in 1 out",
         5: "spiral attract, in", 6: "spiral repel, out", 7: "spiral repel, in", 8: "spiral attract, out"}


def classify(eigs, tol=1e-12):
    """Type of a hyperbolic critical point of a three-dimensional flow from its eigenvalues (real or complex); lower dimensions pad with 'absent'.
    Returns (type number or None, name)."""
    eigs = [complex(e) for e in eigs]
    if any(abs(e.real) <= tol for e in eigs):
        return None, "non-hyperbolic (a zero real part)"
    cplx = [e for e in eigs if abs(e.imag) > tol]; real = [e.real for e in eigs if abs(e.imag) <= tol]
    if not cplx:
        n_in = sum(1 for x in real if x < 0); n = len(real)
        if n_in == n: return 2, TYPES[2]
        if n_in == 0: return 1, TYPES[1]
        if n == 3: return (4, TYPES[4]) if n_in == 2 else (3, TYPES[3])
        return None, f"saddle ({n_in} in, {n - n_in} out) in {n} dimensions"
    sigma = cplx[0].real; axis = real[0] if real else None
    if axis is None:
        return (5, "spiral attract (plane only)") if sigma < 0 else (6, "spiral repel (plane only)")
    if sigma < 0: return (5, TYPES[5]) if axis < 0 else (8, TYPES[8])
    return (7, TYPES[7]) if axis < 0 else (6, TYPES[6])


def gradient_flow_eigs(hessian_eigs):
    """A gradient flow x' = -grad V has eigenvalues minus the Hessian's: a minimum (all positive) is a sink."""
    return [-float(h) for h in hessian_eigs]


# ------------------------------------------------------------------------------------------------------------------------------ 1. the shape well
def shape_well(model="S2"):
    """E(r)/E* = cosh(log r - log r*) with r* = sqrt(A/B): the identity, the floor and the curvature (EXACT)."""
    from mtft.research.vacuum_energy import flux_energy
    fe = flux_energy(model=model); A, B = int(fe["A"]), int(fe["B"])
    phi, phis = sp.symbols("phi phi_star", real=True)
    E = A * sp.exp(-phi) + B * sp.exp(phi); Estar = 2 * sp.sqrt(A * B); r_star = sp.sqrt(sp.Rational(A, B))
    ident = sp.simplify(E / Estar - sp.cosh(phi - sp.log(r_star)))
    curvature_at_floor = sp.simplify(sp.diff(E / Estar, phi, 2).subs(phi, sp.log(r_star)))
    hess_r = sp.simplify(sp.diff(A / sp.Symbol("r", positive=True) + B * sp.Symbol("r", positive=True), sp.Symbol("r", positive=True), 2).subs(sp.Symbol("r", positive=True), r_star))
    return {"model": model, "A": A, "B": B, "r_star": r_star, "E_star_over_4pi2": Estar, "cosh_identity": ident == 0, "curvature_at_floor_phi": curvature_at_floor,
            "hessian_r": hess_r, "type": classify(gradient_flow_eigs([float(hess_r)])), "field": "phi = log r (the log metric of a ratio modulus; not derived)"}


# ------------------------------------------------------------------------------------------------------------------------------ 2. the volume
def volume_runaway():
    """The Einstein-frame potential V = c_f E(r) Vol^-2 + c_c r^-1/2 Vol^-3/2 (PV-01) with both coefficients positive: dV/dVol < 0 everywhere, no finite
    critical point; along the valley r_valley(Vol) the potential still falls, so the (phi, log Vol) gradient flow has its only critical point at
    infinity (EXACT, given the signs)."""
    from mtft.research.vacuum_energy import einstein_frame_scaling
    es = einstein_frame_scaling()
    vol, r, cf, cc, A, B = sp.symbols("Vol r c_f c_c A B", positive=True)
    V = cf * (A / r + B * r) * vol ** -2 + cc * r ** sp.Rational(-1, 2) * vol ** sp.Rational(-3, 2)
    dV = sp.diff(V, vol)
    negative = sp.simplify(dV * vol ** 3) == sp.simplify(-2 * cf * (A / r + B * r) - sp.Rational(3, 2) * cc * r ** sp.Rational(-1, 2) * vol ** sp.Rational(1, 2))
    return {"V": str(V), "dV_dVol_negative": negative and es["volume_stationary_point"] is False, "finite_critical_point": False,
            "type": (None, "no finite critical point: runaway to Vol -> infinity (a sink at infinity)"), "pv01": es["statement"]}


# ------------------------------------------------------------------------------------------------------------------------------ 3. the condensate
def condensate_dterm(xi=1, m2=None):
    """V = (1/2)[(u - xi)^2 + (5/8) u^2], u = v^2 - vbar^2 (CO-01), in the field variables (v, vbar): at the minimum (v*, 0), v*^2 = 8 xi/13, the Hessian
    has a zero eigenvalue along vbar (the flat family u = const) — non-hyperbolic; with m2 (v^2 + vbar^2) and 0 < m2 < xi the minimum moves to
    v^2 = (8/13)(xi - m2) and the vbar eigenvalue becomes exactly 4 m2 > 0: a sink (EXACT).  For m2 >= xi the lifted minimum is at the origin."""
    v, vb = sp.symbols("v vbar", real=True); xi = sp.nsimplify(xi)
    u = v ** 2 - vb ** 2; V = sp.Rational(1, 2) * ((u - xi) ** 2 + sp.Rational(5, 8) * u ** 2)
    vstar = sp.sqrt(sp.Rational(8, 13) * xi)
    H = sp.hessian(V, (v, vb)).subs({v: vstar, vb: 0}); eig0 = sorted(float(e) for e in np.linalg.eigvalsh(np.array(H.evalf(), dtype=float)))
    out = {"v_star_squared_over_xi": sp.Rational(8, 13), "gradient_zero": [sp.simplify(sp.diff(V, x).subs({v: vstar, vb: 0})) for x in (v, vb)] == [0, 0],
           "hessian_eigenvalues": eig0, "degenerate": any(abs(e) < 1e-12 for e in eig0), "type_unlifted": classify(gradient_flow_eigs(eig0))}
    if m2 is not None:
        m2 = sp.nsimplify(m2); Vl = V + m2 * (v ** 2 + vb ** 2)
        # the lifted minimum: v^2 = (8/13)(xi - m2) exactly (dV/dv = 2v[(13/8) v^2 - xi + m2] = 0)
        vl2 = sp.Rational(8, 13) * (xi - m2); vl = sp.sqrt(vl2)
        Hl = sp.hessian(Vl, (v, vb)).subs({v: vl, vb: 0}); eigl = sorted(float(e) for e in np.linalg.eigvalsh(np.array(Hl.evalf(), dtype=float)))
        out.update({"m2": float(m2), "v_lifted": float(vl), "hessian_lifted": eigl, "vbar_eigenvalue_is_4m2": abs(eigl[0] - 4 * float(m2)) < 1e-9,
                    "type_lifted": classify(gradient_flow_eigs(eigl))})
    return out


# ------------------------------------------------------------------------------------------------------------------------------ 4. the Wilson line
def wilson_line_hessians(eps=(0, 1), N=60):
    """V(theta) = -2 sum_{odd w} |w|^-8 cos(2 pi theta . (m, n)) on the torus of the relative U(1) (PA-01's product approximation): the 2 x 2 Hessians
    at theta = 0 (jet vacuum) and theta_eta = eps/2 (untwisted point), their eigenvalues and types — EXACT in PA-01's product approximation to the
    potential (the same standing as PA1-1), not in the full Casimir sum."""
    from mtft.matter.path_a import _lattice
    L = _lattice(eps, N)
    def hess(theta):
        H = np.zeros((2, 2))
        for m, n, r2, odd in L:
            if not odd: continue
            c = 2 * r2 ** -4 * (2 * math.pi) ** 2 * math.cos(2 * math.pi * (theta[0] * m + theta[1] * n))
            H += c * np.array([[m * m, m * n], [m * n, n * n]])
        return H
    H0 = hess((0.0, 0.0)); He = hess((eps[0] / 2, eps[1] / 2))
    e0 = sorted(np.linalg.eigvalsh(H0)); ee = sorted(np.linalg.eigvalsh(He))
    return {"eps": eps, "hessian_jet": H0.tolist(), "eigs_jet": [float(x) for x in e0], "type_jet": classify(gradient_flow_eigs(e0)),
            "hessian_untwisted": He.tolist(), "eigs_untwisted": [float(x) for x in ee], "type_untwisted": classify(gradient_flow_eigs(ee)),
            "minus_each_other": bool(np.allclose(H0, -He))}


# ------------------------------------------------------------------------------------------------------------------------------ 4b. the whole torus (CC-66)
def kinetic_metric():
    """The kinetic metric of the Wilson-line moduli in the holonomy coordinates theta = (theta_1, theta_2) around the 1-loop and the tau-loop of 143a1:
    its Hodge metric G = [[|tau|^2/y, -x/y], [-x/y, 1/y]] with tau = x + i y (det G = 1; the overall normalisation of the kinetic term is not fixed here).
    The gradient flow of V is -G^-1 grad V, and small oscillations about a minimum have frequencies^2 equal to the eigenvalues of G^-1 H."""
    from mtft.research.theta_torus import tau_143a1
    tau = complex(tau_143a1()); x, y = tau.real, tau.imag
    return np.array([[abs(tau) ** 2 / y, -x / y], [-x / y, 1 / y]])


def _odd_windings(eps, N):
    from mtft.matter.path_a import _lattice
    a = np.array([(m, n, r2) for m, n, r2, odd in _lattice(eps, N) if odd], dtype=float)
    return a[:, 0], a[:, 1], a[:, 2] ** -4                                            # m, n, |w|^-8


def wilson_line_torus(eps=(0, 1), N=40, grid=12):
    """Every critical point of V(theta) = -2 sum_{odd w} |w|^-8 cos(2 pi theta . w) on the whole torus R^2/Z^2, not only along PA-01's straight path
    (PG-6, CC-66): Newton's method from a grid x grid set of starts; each point typed by the eigenvalues of G^-1 H in the kinetic metric
    (`kinetic_metric`); the Morse count #min - #saddle + #max, which must be chi(T^2) = 0 for a Morse function; the curvatures at the jet vacuum in
    theta-coordinates and in the metric (the masses follow the latter; the product agrees because det G = 1).  In PA-01's product approximation it
    is EXACT that the four half-periods are critical points (V(-theta) = V(theta); V(theta_eta + d) = -V(d) maps the sink to the source) and
    NUMERICAL that there are no others (this census, and in tests/test_pg01_wilson_torus.py a winding-number count with a Lipschitz exclusion
    and a contraction bound); all four are nondegenerate, so the torus forces exactly two saddles.  The values are lattice sums over the box
    |m|, |n| <= N."""
    m, n, w8 = _odd_windings(eps, N)
    Gi = np.linalg.inv(kinetic_metric())
    ph = lambda t: 2 * math.pi * (t[0] * m + t[1] * n)
    V = lambda t: float(-2 * np.sum(w8 * np.cos(ph(t))))
    def grad(t):
        s = np.sin(ph(t)) * w8
        return np.array([4 * math.pi * np.sum(s * m), 4 * math.pi * np.sum(s * n)])
    def hess(t):
        c = np.cos(ph(t)) * w8 * 8 * math.pi ** 2
        return np.array([[np.sum(c * m * m), np.sum(c * m * n)], [np.sum(c * m * n), np.sum(c * n * n)]])
    found = []
    for i in range(grid):
        for j in range(grid):
            t = np.array([(i + 0.5) / grid, (j + 0.5) / grid])
            for _ in range(80):
                g = grad(t)
                if np.linalg.norm(g) < 1e-11:
                    break
                try:
                    step = np.linalg.solve(hess(t), g)
                except np.linalg.LinAlgError:
                    break
                if np.linalg.norm(step) > 0.05:
                    step *= 0.05 / np.linalg.norm(step)
                t = t - step
            if np.linalg.norm(grad(t)) < 1e-9:
                t = np.round(t % 1.0, 9) % 1.0
                if not any(np.allclose(((t - f + 0.5) % 1.0) - 0.5, 0, atol=1e-6) for f in found):
                    found.append(t)
    rows = []
    for t in sorted(found, key=lambda u: (round(float(u[0]), 6), round(float(u[1]), 6))):
        ev = sorted(float(x) for x in np.linalg.eigvals(Gi @ hess(t)).real)
        kind = "sink" if ev[0] > 0 else "source" if ev[1] < 0 else "saddle"
        rows.append({"theta": (round(float(t[0]), 6), round(float(t[1]), 6)), "V": V(t), "type": kind,
                     "flow_type": classify(gradient_flow_eigs(ev))[1], "metric_curvatures": ev})
    kinds = [r["type"] for r in rows]
    H0 = hess((0.0, 0.0))
    e_theta = sorted(float(x) for x in np.linalg.eigvalsh(H0)); e_metric = sorted(float(x) for x in np.linalg.eigvals(Gi @ H0).real)
    return {"eps": tuple(eps), "critical_points": rows, "counts": {k: kinds.count(k) for k in ("sink", "saddle", "source")},
            "morse_sum": kinds.count("sink") - kinds.count("saddle") + kinds.count("source"),
            "all_at_half_periods": all(abs(2 * x - round(2 * x)) < 1e-6 for r in rows for x in r["theta"]),
            "jet_curvatures_theta": e_theta, "jet_curvatures_metric": e_metric,
            "mass_ratio_metric": math.sqrt(e_metric[1] / e_metric[0]), "mass_ratio_if_theta_coordinates": math.sqrt(e_theta[1] / e_theta[0])}


def saddle_pairing(N=40):
    """The argument for the exact zeros at the two saddles of the real spin structure eps = (0, 1) (odd n): (m, n) -> (-m - n, n) keeps
    |m + n tau| because Re tau = 1/2, and for odd n it flips the parity of m, so at theta = (1/2, 0) and (1/2, 1/2) the odd windings cancel in
    pairs.  The two facts are algebraic identities once Re tau = 1/2, so this is a regression guard on tau (it fails only if Re tau moves off 1/2);
    the zeros themselves are evaluated in `wilson_line_torus` and by the disc sum in tests/test_pg01_wilson_torus.py."""
    eps = (0, 1)
    from mtft.research.theta_torus import tau_143a1
    tau = complex(tau_143a1()); s = math.sqrt(tau.imag)
    norm = parity = True; pairs = 0
    for mm in range(-N, N + 1):
        for nn in range(-N, N + 1):
            if (eps[0] * mm + eps[1] * nn) % 2 != 1 or abs(mm + nn) > N:
                continue
            pairs += 1
            norm &= abs(abs(mm + nn * tau) - abs(-mm - nn + nn * tau)) / s < 1e-9
            parity &= (-mm - nn - mm) % 2 == 1
    return {"re_tau_is_half": abs(tau.real - 0.5) < 1e-12, "windings_checked": pairs, "norm_preserved": bool(norm),
            "parity_of_m_flipped": bool(parity)}


# ------------------------------------------------------------------------------------------------------------------------------ 5. the rest, cited
def cited_directions():
    from mtft.matter import one_loop_adjoints as OL, moduli_casimir as MC
    K = OL.mass_coefficients()
    return {"adjoints": {"claim": "M05-1", "positive_mass_coefficients": {k: v > 0 for k, v in K.items()}, "type": (2, TYPES[2])},
            "casimir_A_X": {"claim": "MS-6", "minimum_possible": MC.minimum_possible(-24), "type": (None, "no critical point in A_X at one loop (N10)")}}


# ------------------------------------------------------------------------------------------------------------------------------ the atlas
def atlas():
    sw = shape_well(); vr = volume_runaway(); cd = condensate_dterm(m2=sp.Rational(1, 20)); wl = {e: wilson_line_hessians(e) for e in ((0, 1), (1, 0), (1, 1))}
    wt = wilson_line_torus((0, 1))
    cited = cited_directions()
    rows = [("shape r = A_X/A_E (flux landscape, PV-01)", "cosh well, floor at r* = sqrt(19/11)", sw["type"][1], "PG-1"),
            ("volume Vol (Einstein frame, PV-01)", "no finite critical point", vr["type"][1], "PG-4"),
            ("condensate (v^2, vbar^2) (CO-01)", "degenerate line, lifted by m^2 > 0", f"{cd['type_unlifted'][1]} -> {cd['type_lifted'][1]}", "PG-4"),
            ("Wilson line theta on the torus (PA-01), jet vacuum", "2D minimum", wl[(0, 1)]["type_jet"][1], "PG-4"),
            ("Wilson line theta, untwisted point theta_eta", "2D maximum", wl[(0, 1)]["type_untwisted"][1], "PG-4"),
            *[(f"Wilson line theta, saddle at {r['theta']}", "2D saddle (V = 0 exactly, real spin structure)", r["flow_type"], "PG-6")
              for r in wt["critical_points"] if r["type"] == "saddle"],
            ("27 bosonic adjoints per generator (M05-1)", "positive one-loop masses", cited["adjoints"]["type"][1], "PG-4"),
            ("A_X at one loop (MS-01)", "no Casimir minimum", cited["casimir_A_X"]["type"][1], "PG-4")]
    return {"rows": rows, "shape": sw, "volume": vr, "condensate": cd, "wilson": wl, "wilson_torus": wt, "cited": cited,
            "saddles": wt["counts"]["saddle"], "spirals_without_friction": 0,
            "reading": "sinks, one source, two saddles (the Wilson-line torus), one degenerate line, two runaways; spirals need time"}


# ------------------------------------------------------------------------------------------------------------------------------ friction and the toy
def friction_type(m, H):
    """phi'' + 3H phi' + m^2 phi = 0: eigenvalues (-3H +- sqrt(9H^2 - 4m^2))/2; spiral iff m > 3H/2 (EXACT)."""
    disc = 9 * H * H - 4 * m * m
    if disc < 0:
        eigs = [complex(-1.5 * H, math.sqrt(-disc) / 2), complex(-1.5 * H, -math.sqrt(-disc) / 2)]; kind = "spiral attract (underdamped)"
    elif disc == 0:
        eigs = [-1.5 * H, -1.5 * H]; kind = "critically damped node"
    else:
        eigs = [(-3 * H + math.sqrt(disc)) / 2, (-3 * H - math.sqrt(disc)) / 2]; kind = "node (overdamped sink)"
    return {"m": m, "H": H, "threshold_H": 2 * m / 3, "eigenvalues": eigs, "kind": kind, "type": classify(eigs)}


def toy_run(H=0.2, mode="constant", phi0=1.5, v0=0.0, t_end=60.0, n=6000, model="S2"):
    """DIAGNOSTIC: the shape modulus rolling in its cosh well with Hubble friction, phi'' + 3H(t) phi' + sinh(phi - phi*) = 0, in units where the
    well's curvature at the floor is 1 (time in units of 1/m); mode 'constant' or 'radiation' (H = 1/(2t), started at t0 = 1/(2H)).  Returns the
    trajectory, the crossings of the floor — a count that grows with the run length (18 in a run of length 60 at H = 0.2, 15 in one of length 50) —
    and the final distance.  Capture: E = phi'^2/2 + cosh(phi - phi*) - 1 has dE/dt = -3H phi'^2 <= 0 (`energy_argument`), and with H bounded,
    of bounded variation and not integrable — both modes here — E -> 0 (CC-67; `capture_conditions` has the theorem and the cases where it fails);
    the "captured" flag is only a smoke test that the envelope has shrunk within the run.  The normalisation of H against the modulus mass is not
    derived (the point of the gate)."""
    from scipy.integrate import solve_ivp
    sw = shape_well(model); phis = float(sp.log(sw["r_star"]))
    t0 = 0.0 if mode == "constant" else 1 / (2 * H)
    Hf = (lambda t: H) if mode == "constant" else (lambda t: 1 / (2 * t))
    def rhs(t, y):
        return [y[1], -3 * Hf(t) * y[1] - math.sinh(y[0] - phis)]
    ts = np.linspace(t0, t0 + t_end, n)
    sol = solve_ivp(rhs, (t0, t0 + t_end), [phis + phi0, v0], t_eval=ts, rtol=1e-9, atol=1e-12)
    x = sol.y[0] - phis; crossings = int(np.sum(np.sign(x[1:]) != np.sign(x[:-1])))
    k = max(1, n // 5); amp_first = float(np.max(np.abs(x[:k]))); amp_last = float(np.max(np.abs(x[-k:])))
    run = {"mode": mode, "H": H, "phi_star": phis, "t": sol.t, "phi": sol.y[0], "phidot": sol.y[1], "floor_crossings": crossings,
           "run_length": t_end, "final_distance": float(abs(x[-1])), "amplitude_first_fifth": amp_first, "amplitude_last_fifth": amp_last,
           "captured": amp_last < 0.25 * amp_first and amp_last < 0.2,                      # smoke test: the envelope shrinks within the run
           "friction": friction_type(1.0, H if mode == "constant" else 1 / (2 * (t0 + t_end))), "class": "DIAGNOSTIC"}
    run["energy"] = energy_argument(run)
    return run


def energy_argument(run):
    """The energy identity, checked on a run: E = phi'^2/2 + cosh(phi - phi*) - 1 satisfies dE/dt = -3H phi'^2 <= 0 (the friction only removes energy).
    Capture at the floor follows when H is bounded, of bounded variation and not integrable, as in both modes of `toy_run` — not for every
    H(t) > 0 (CC-67 corrected that; see `capture_conditions`).  Returns E along the run, whether it is non-increasing, and the fraction of the
    initial energy that remains."""
    E = 0.5 * np.asarray(run["phidot"]) ** 2 + np.cosh(np.asarray(run["phi"]) - run["phi_star"]) - 1
    return {"E_initial": float(E[0]), "E_final": float(E[-1]), "non_increasing": bool(np.all(np.diff(E) <= 1e-9 * max(1.0, float(E[0])))),
            "fraction_remaining": float(E[-1] / E[0]) if E[0] > 0 else 0.0,
            "statement": "dE/dt = -3H phi'^2 <= 0; capture follows for H bounded, of bounded variation and not integrable (CC-67)"}


def capture_conditions(phi0=1.5, t_late=(400.0, 1500.0)):
    """When Hubble friction captures the shape, and when it does not (CC-67, 2026-10-16).  In the well, x = phi - phi*, E = x'^2/2 + cosh x - 1 has
    dE/dt = -3H x'^2 <= 0 for any H >= 0, but capture needs more.
    EXACT (proof in the PG-01 report, section 6): in a convex coercive single well, with H >= 0 of bounded variation (so bounded) and
    int H dt = infinity, every trajectory is captured, E -> 0 — with W = x x', the convexity inequality x V'(x) >= V(x) and two
    Riemann–Stieltjes integrations by parts give int H E dt < infinity (the nonincreasing case is Cabot–Engler–Gadat, Trans. AMS 361 (2009)
    5983).  Neither hypothesis can simply be dropped, though neither is necessary (H = t also captures): an integrable H leaves
    E(t) >= E(0) exp(-6 int_0^t H) (Gronwall, since x'^2 <= 2E), and a fast-growing H overdamps the field to rest off the floor:
    x = a + b e^-t with a, b > 0 solves x'' + 3H x' + sinh x = 0 exactly with H = 1/3 + e^t sinh(x)/(3b) > 0, whose integral diverges (Codex's
    example).  NUMERICAL: the integrable case H = 0.1 e^-t, the overdamped example, and captured cases — constant H, radiation H = 1/(2t) and
    matter H = 2/(3t) from t0 > 0 (the toy runs use the first two) — whose energies decay like t^-1.5 and t^-2 at late times (E ~ t^(-3p) for
    H = p/t, from the Bessel asymptotics of x'' + (3p/t) x' + x = 0)."""
    from scipy.integrate import solve_ivp

    def run(Hf, t0, t1, x0, v0=0.0, method="DOP853"):
        sol = solve_ivp(lambda t, y: [y[1], -3.0 * Hf(t) * y[1] - math.sinh(y[0])], (t0, t1), [x0, v0], method=method, rtol=1e-10, atol=1e-12,
                        dense_output=True)
        def E(t):
            x, v = sol.sol(t)
            return 0.5 * v * v + math.cosh(x) - 1.0
        return sol, E

    E0 = math.cosh(phi0) - 1.0
    sol, E = run(lambda t: 0.1 * math.exp(-t), 0.0, 200.0, phi0)
    ts = np.linspace(0.0, 200.0, 2001); Es = np.array([E(t) for t in ts])
    integrable = {"H": "0.1 exp(-t)", "fraction_left": float(Es[-1] / E0), "gronwall_floor": math.exp(-0.6),
                  "bound_holds": bool(np.all(Es >= E0 * np.exp(-0.6 * (1 - np.exp(-ts))) * (1 - 1e-9)))}
    a, b = 1.0, 0.5
    t_, a_, b_ = sp.symbols("t a b", positive=True)
    xs = a_ + b_ * sp.exp(-t_); Hs = sp.Rational(1, 3) + sp.exp(t_) * sp.sinh(xs) / (3 * b_)
    exact = sp.simplify(sp.diff(xs, t_, 2) + 3 * Hs * sp.diff(xs, t_) + sp.sinh(xs)) == 0
    sol, E = run(lambda t: 1 / 3 + math.exp(t) * math.sinh(a + b * math.exp(-t)) / (3 * b), 0.0, 12.0, a + b, -b, method="Radau")
    overdamped = {"exact_solution": bool(exact), "x_at_12": float(sol.y[0][-1]), "x_limit": a, "H_integral_diverges": True,
                  "energy_left": float(E(12.0) / (0.5 * b * b + math.cosh(a + b) - 1.0))}
    captured = {}
    for name, Hf, t0, t1, p in (("constant H = 0.2", lambda t: 0.2, 0.0, 200.0, None),
                                ("radiation H = 1/(2t)", lambda t: 1 / (2 * t), 1.25, t_late[1], 0.5),
                                ("matter H = 2/(3t)", lambda t: 2 / (3 * t), 1.0, t_late[1], 2 / 3)):
        sol, E = run(Hf, t0, t1, phi0)
        row = {"E_end_over_E0": float(E(t1) / E0)}
        if p is not None:
            row.update({"decay_exponent": math.log(E(t_late[1]) / E(t_late[0])) / math.log(t_late[1] / t_late[0]), "predicted": -3 * p})
        captured[name] = row
    return {"integrable": integrable, "overdamped": overdamped, "captured": captured,
            "theorem": "convex coercive single well, H >= 0 bounded, of bounded variation, int H = infinity => E -> 0",
            "class": "EXACT (theorem and identities) / NUMERICAL (runs)"}


def volume_companion(t_end=60.0, n=600, c_f=1.0, c_c=1.0, H=0.2):
    """DIAGNOSTIC: log Vol rolling in V = c_f Vol^-2 + c_c Vol^-3/2 (both coefficients positive) with the same constant friction: it never stops (the
    runaway), and at late times, once the slower term dominates and the motion is overdamped, log Vol ~ (2/3) log t — without bound, but only
    logarithmically in this toy."""
    from scipy.integrate import solve_ivp
    def rhs(t, y):
        L = y[0]; dV = -2 * c_f * math.exp(-2 * L) - 1.5 * c_c * math.exp(-1.5 * L)     # dV/dlogVol
        return [y[1], -3 * H * y[1] - dV]
    sol = solve_ivp(rhs, (0, t_end), [0.0, 0.0], t_eval=np.linspace(0, t_end, n), rtol=1e-9, atol=1e-12)
    late = sol.t > t_end / 2
    slope = float(np.polyfit(np.log(sol.t[late]), sol.y[0][late], 1)[0]) if late.sum() > 2 else float("nan")
    return {"t": sol.t, "log_vol": sol.y[0], "monotone_increasing": bool(np.all(np.diff(sol.y[0]) > -1e-12)), "final_log_vol": float(sol.y[0][-1]),
            "late_time_slope_dlogVol_dlogt": slope, "asymptotic_slope": 2 / 3, "late_time_law": "log Vol ~ (2/3) log t (overdamped, Vol^-3/2 term dominant)"}


# ------------------------------------------------------------------------------------------------------------------------------ the gate
GATE = (("I.1", "a declared parent: supersymmetric U(8) (particles on the difference lines, GF-6) or PV-01's non-supersymmetric U(8)"),
        ("I.6", "a vacuum energy computed from the parent action with a stationary point in the volume — today the classical terms run away (PG-2) and the one-loop Casimir has no minimum (MS-6)"),
        ("GS-1", "fixed Kaehler moduli, which the FI-driven condensates need"),
        ("GC-01", "the couplings as functions of the moduli at the stabilised point (volume and kappa), with the matching scale (VI.1)"),
        ("cosmology", "a cosmological background to supply H(t): the package's Friedmann module (mtft.cosmology) is pheno-tier and would have to be re-derived or replaced"),
        ("standard physics", "after the moduli settle and the couplings are fixed: confinement, hadronisation and nucleosynthesis are standard and should be used, not re-derived"))


def simulation_gate():
    return {"needs": GATE, "status": "open: none of the first three is settled", "what_runs_today": "the shape modulus (toy_run) and the volume's runaway (volume_companion)"}


def report():
    A = atlas()
    runs = {"constant_H_0.2": toy_run(0.2), "constant_H_1.0": toy_run(1.0), "radiation": toy_run(0.4, mode="radiation", t_end=80.0)}
    A["wilson"] = {str(k): v for k, v in A["wilson"].items()}
    return {"atlas": A, "friction": {H: friction_type(1.0, H) for H in (0.2, 2 / 3, 1.0)},
            "runs": {k: {kk: v for kk, v in r.items() if kk not in ("t", "phi", "phidot")} for k, r in runs.items()},
            "volume": {k: v for k, v in volume_companion().items() if k not in ("t", "log_vol")}, "gate": simulation_gate()}


def text():
    r = report()
    L = ["PG-01 — the vacuum as a dynamical system (2026-10-12)", "", "critical-point atlas:"]
    for d, what, typ, claim in r["atlas"]["rows"]:
        L.append(f"  {d:56s} {what:38s} -> {typ}  [{claim}]")
    sw = r["atlas"]["shape"]
    L.append(f"shape well: E/E* = cosh(phi - phi*) {sw['cosh_identity']}, floor r* = {sw['r_star']}, curvature at the floor {sw['curvature_at_floor_phi']} (units of E*)")
    wl = r["atlas"]["wilson"]["(0, 1)"]
    L.append(f"Wilson line (eps = (0,1)): Hessian eigenvalues jet {[round(x, 3) for x in wl['eigs_jet']]}, untwisted {[round(x, 3) for x in wl['eigs_untwisted']]}, minus each other {wl['minus_each_other']}")
    wt = r["atlas"]["wilson_torus"]
    residue = ", ".join(f"{c['V']:.1e}" for c in wt["critical_points"] if c["type"] == "saddle")
    L.append(f"Wilson-line torus (eps = (0,1)): {wt['counts']['sink']} sink, {wt['counts']['saddle']} saddles, {wt['counts']['source']} source at the half-periods, "
             f"Morse sum {wt['morse_sum']}; saddle values exactly 0 by the pairing (box-truncation residue {residue}); "
             f"jet curvatures in the metric {[round(x, 3) for x in wt['jet_curvatures_metric']]} (theta-coordinates {[round(x, 3) for x in wt['jet_curvatures_theta']]}) [PG-6, CC-66]")
    for H, f in r["friction"].items():
        L.append(f"friction: m = 1, H = {H:.3f}: {f['kind']} (threshold H = 2m/3 = {f['threshold_H']:.3f})")
    for k, run in r["runs"].items():
        L.append(f"toy run {k} (length {run['run_length']:.0f}/m): floor crossings {run['floor_crossings']}, envelope shrunk {run['captured']}, "
                 f"energy non-increasing {run['energy']['non_increasing']} ({run['friction']['kind']}) [DIAGNOSTIC]")
    L.append(f"volume companion: monotone increasing {r['volume']['monotone_increasing']}, final log Vol {r['volume']['final_log_vol']:.2f}, "
             f"late slope d log Vol/d log t = {r['volume']['late_time_slope_dlogVol_dlogt']:.2f} -> 2/3 [DIAGNOSTIC]")
    L.append("simulation gate: " + "; ".join(f"{k}: {v[:60]}…" for k, v in r["gate"]["needs"][:3]))
    return "\n".join(L)
