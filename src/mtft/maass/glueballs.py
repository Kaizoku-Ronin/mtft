r"""mtft.maass.glueballs — H24 against the lattice: the dictionary, every registered forecast and its outcome (PHENO; GB-01, 2026-11-16).

Hypothesis H24 (R&D notes, 16 October 2026): if MTFT's configuration space is X₀(143) with the
hyperbolic Laplacian as Hamiltonian, the glueball mass ratios of Yang–Mills theory are ratios of
the Maass spectrum of Γ₀(143).  The lab tested it on 12 to 14 November 2026 (records MP15, MP18,
MP19, MP20).  This module holds that test: the lattice numbers, the dictionary, what was predicted
before each run, and what happened.  **H24 as written fails**, in both of its readings, and the
status of everything here is PHENO: a comparison under a stated dictionary, not a derivation.

The lattice (``LATTICE``).  Athenodorou and Teper, JHEP 11 (2020) 172, arXiv:2007.06422: SU(3),
continuum limit, masses in units of √σ.  Six values, as the lab transcribed them on 12 November
(``TRANSCRIPTION`` says what that means).

The dictionary (``DICTIONARY``).  D1: a state's energy is its eigenvalue λ = ¼ + R², and only
ratios are tested.  D2: weight 0 means J = 0 only.  D3: parity is the mirror.  D4: C is the Fricke
sign ε₁₁ε₁₃; or D4′: C = −1 for the η-twisted forms.  D5: reading a takes all of L² with the 0⁺⁺ at
the edge ¼ of the continuum; reading b takes the cusp forms only.

What happened (``REGISTER``; the numbers are recomputed here from the frozen spectrum):

  reading b      all five predictions miss (``reading_b``): the lowest cusp form is odd, and the
                 three ratios are 0.870, 2.236, 1.247 against windows near 1.55, 1.72, 2.14
  reading a      one prediction passed, the 0⁻⁺ at 1.559 against 1.549(16); the reading fails on
                 structure, since its 0⁺⁺ is the edge of a continuum and no isolated state
  the η-twist    a registered miss: C = −1 states at 1.86 and 1.92 times the edge against the
                 lattice's 2.87(14) and 2.74(8) (``comparison``)
  an interaction ``self_interaction_fit``: χ² = 10.2 on 2; ``one_coupling_fit``: 170 on 4
  a twist energy ``twist_energy``: μ = 0.218(18) fits two numbers with one, χ² = 1.57 on 1

Two forecasts are open, and sealed (``SEALS``, ``twist_test``).  With μ fixed, a second odd
C = −1 state at 2.87(7) and a second even one at 3.24(7), in units of the 0⁺⁺.  The lab registered
them on 14 November before any lookup; the criterion was fixed and hash-sealed on 16 November,
still before any lookup.  ``evaluate_open`` runs the sealed rule and refuses if the seal no longer
verifies.

A new dictionary is tested forecast first: ``register_forecast`` seals the predicted ratios, and
``compare_forecast`` takes only such a record and answers in the form of ``mtft.falsify``.
"""
from __future__ import annotations

import hashlib
import json
import math

from . import spectrum as _sp

SOURCE = ("A. Athenodorou and M. Teper, JHEP 11 (2020) 172, arXiv:2007.06422: SU(3), continuum limit of M/√σ; "
          "Table 17 (states with J = 0 identified) and Table 14 (the lightest states of A1+- and A1--)")
CITATION = "Athenodorou and Teper, JHEP 11 (2020) 172"
LATTICE = {"0++": (3.405, 0.021, 17), "0++*": (5.855, 0.041, 17), "0-+": (5.276, 0.045, 17), "0-+*": (7.29, 0.13, 17),
           "A1+-": (9.32, 0.28, 14), "A1--": (9.78, 0.46, 14)}
"""state → (M/√σ, error, the source's table)."""
TRANSCRIPTION = ("as transcribed in the lab record MP15 on 12 November 2026; not re-read against the paper when this module was "
                 "written, because the same tables would show the excited C = −1 states that the sealed forecasts are about; "
                 "to be checked when those forecasts are evaluated")

DICTIONARY = {
    "D1": "energy: the eigenvalue λ = 1/4 + R²; the unit is not fixed, so only ratios are tested",
    "D2": "spin: Maass forms of weight 0 are scalars, so only J = 0 is tested; the 2++ has no partner",
    "D3": "parity: the mirror τ ↦ −τ̄; even forms P = +, odd forms P = −",
    "D4": "charge conjugation (MP15): the Fricke sign, C = ε₁₁ε₁₃",
    "D4'": "charge conjugation (MP18): C = +1 for functions on X₀(143), C = −1 for η-twisted ones",
    "D5a": "states, reading a: all of L², the continuum included; the 0++ is the edge 1/4",
    "D5b": "states, reading b: the cusp forms only",
}
EDGE = 0.25
READING_A = {"0-+": "u.o.mm.1", "0++*": "u.e.mm.1", "0-+*": "u.o.pp.1", "A1+-": "t.e.pp.1", "A1--": "t.o.mm.1"}
"""Reading a with D4′: lattice state → form of the frozen table (``mtft.maass.spectrum``)."""
OTHER_TWISTED = {"second odd": "t.o.pm.1", "second even": "t.e.pm.1"}
PARTICIPATION = {"u.o.mm.1": 1.8905032283795062, "u.e.mm.1": 2.3287582173942587, "u.o.pp.1": 2.2821076294022764,
                 "t.e.pp.1": 2.191984862548154, "t.o.mm.1": 1.9034732534165089, "t.e.pm.1": 1.8883844544772608,
                 "t.o.pm.1": 2.4722362026688436}
"""The participation ratios I = area·∫u⁴/(∫u²)² of the seven forms, as the lab record MP20 found them."""


# ------------------------------------------------------------------ small tools
def chi2_sf(x: float, dof: int) -> float:
    """The probability that a χ² variable with ``dof`` degrees of freedom exceeds x."""
    import mpmath as mp
    return float(mp.gammainc(mp.mpf(dof) / 2, mp.mpf(x) / 2, mp.inf, regularized=True))


def _solve(A, b):
    n = len(A)
    M = [list(map(float, A[i])) + [float(b[i])] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(M[i][c]))
        M[c], M[p] = M[p], M[c]
        for i in range(n):
            if i != c:
                f = M[i][c] / M[c][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def _quad(u, C, v) -> float:
    """u · C⁻¹ · v."""
    return sum(a * b for a, b in zip(u, _solve(C, v)))


def lattice_ratio(state: str) -> tuple:
    """(M/M(0++), error) for a lattice state, the two errors combined in quadrature."""
    m, s, _ = LATTICE[state]
    m0, s0, _ = LATTICE["0++"]
    r = m / m0
    return r, r * math.hypot(s / m, s0 / m0)


def curve_ratio(form_id: str) -> float:
    """λ/(1/4) = 1 + 4R² for a form of the frozen table: its energy in units of the edge."""
    return _sp.form(form_id)["lambda"] / EDGE


def covariance(states, correlated: bool = True) -> list:
    """The covariance of the needed shifts of λ for some lattice states (ratios over 4).

    The ratios share the error of the 0++; ``correlated=False`` keeps the diagonal only."""
    m0, s0, _ = LATTICE["0++"]
    r = [LATTICE[k][0] / m0 for k in states]
    rel = [LATTICE[k][1] / LATTICE[k][0] for k in states]
    n = len(states)
    return [[(r[i] * r[j] * ((s0 / m0) ** 2 + (rel[i] ** 2 if i == j else 0.0)) / 16) if (correlated or i == j) else 0.0
             for j in range(n)] for i in range(n)]


def needed_shifts() -> dict:
    """For each state of reading a: the lattice ratio, the curve's, and the shift of λ that would reconcile them."""
    out = {}
    for state, fid in READING_A.items():
        lr, le = lattice_ratio(state)
        cr = curve_ratio(fid)
        out[state] = {"form": fid, "lattice_ratio": lr, "error": le, "curve_ratio": cr, "shift": (lr - cr) / 4,
                      "shift_error": le / 4, "sigma": (cr - lr) / le, "percent": 100 * (lr - cr) / cr}
    return out


def comparison() -> list:
    """Reading a with D4′, state by state, in the form of ``mtft.falsify`` (PASS within 2σ, TENSION to 3σ, FAIL beyond)."""
    from mtft.falsify import Prediction
    notes = {"0-+": "the one registered pass (S2); levels 133 and 145 have no form there: a coincidence without a mechanism",
             "0++*": "first looked at after MP15's run: not a registered test", "0-+*": "found after looking at several ratios",
             "A1+-": "MP18's P1: a registered miss", "A1--": "MP18's P1: a registered miss"}
    rows = []
    for n, (state, d) in enumerate(needed_shifts().items(), start=1):
        sig = abs(d["sigma"])
        rows.append(Prediction(number=n, relation=f"lambda({d['form']})/edge vs M({state})/M(0++)", predicted=d["curve_ratio"],
                               observed=d["lattice_ratio"], obs_error=d["error"],
                               error_percent=abs(d["curve_ratio"] - d["lattice_ratio"]) / d["lattice_ratio"] * 100, sigma=sig,
                               sector="glueball, C = " + ("-1" if state.startswith("A1") else "+1"), source="MP15, MP18; " + CITATION,
                               status="PASS" if sig <= 2 else ("TENSION" if sig <= 3 else "FAIL"), rating="",
                               deviation_ppm=abs(d["curve_ratio"] - d["lattice_ratio"]) / d["lattice_ratio"] * 1e6,
                               multiplicity=notes[state], group="H24 reading a"))
    return rows


# ------------------------------------------------------------------ reading b (MP15)
def channel(parity: int, c: int) -> list:
    """The untwisted cusp forms of a channel under D3 and D4, by increasing R: parity 0 or 1, C = ε₁₁ε₁₃ = ±1."""
    return [f for f in _sp.forms(twist=False, parity=parity) if f["eps11"] * f["eps13"] == c]


def reading_b() -> dict:
    """MP15's five predictions for reading b, recomputed: each with its value, its window and whether it holds."""
    even_p, odd_p = channel(0, 1), channel(1, 1)
    lowest = min(_sp.forms(twist=False), key=lambda f: f["R"])
    l0 = even_p[0]["lambda"]
    minus = min(channel(0, -1) + channel(1, -1), key=lambda f: f["R"])
    ratios = {"N1": (odd_p[0]["lambda"] / l0, (1.472, 1.627)), "N2": (even_p[1]["lambda"] / l0, (1.634, 1.806)),
              "N3": (odd_p[1]["lambda"] / l0, (2.020, 2.262))}
    out = {"N0": {"statement": "the lowest cusp form is even with Fricke sign +", "value": (lowest["id"], lowest["parity"]),
                  "holds": lowest["parity"] == 0 and lowest["eps11"] * lowest["eps13"] == 1}}
    for k, (v, w) in ratios.items():
        out[k] = {"value": v, "window": w, "holds": w[0] <= v <= w[1]}
    out["N4"] = {"statement": "no cusp form with Fricke sign − below 2.60 λ(0++)", "value": minus["lambda"] / l0, "form": minus["id"],
                 "bound": 2.60, "holds": minus["lambda"] / l0 >= 2.60}
    out["forms"] = {"0++": even_p[0]["id"], "0++*": even_p[1]["id"], "0-+": odd_p[0]["id"], "0-+*": odd_p[1]["id"]}
    out["all_miss"] = not any(out[k]["holds"] for k in ("N0", "N1", "N2", "N3", "N4"))
    return out


def reading_a_checks() -> dict:
    """What can be recomputed of reading a and the η-twist: S2's window, MP18's bound, MP19's empty window."""
    r0 = _sp.form(READING_A["0-+"])["R"]
    lr, le = lattice_ratio("0-+")
    twisted_low = [f["id"] for f in _sp.forms(twist=True) if f["R"] < 0.6325]
    window = [f["id"] for f in channel(0, 1) if 0.410 <= f["R"] <= 0.438]
    return {"S2": {"R": r0, "window": (0.344, 0.396), "holds": 0.344 <= r0 <= 0.396, "sigma": (curve_ratio(READING_A["0-+"]) - lr) / le},
            "MP18-P1": {"bound": 0.6325, "twisted_forms_below": twisted_low, "holds": not twisted_low},
            "MP19-P4": {"window": (0.410, 0.438), "even_forms_with_fricke_plus": window, "holds": not window}}


# ------------------------------------------------------------------ MP20: could it be a hidden self-interaction?
def participation(recompute: bool = False) -> dict:
    """The seven participation ratios: the lab's recorded values, or (``recompute=True``, about 12 s) the package's own
    quadrature over the Ford domain (``mtft.maass.ford``)."""
    if not recompute:
        return dict(PARTICIPATION)
    from . import ford as _fd
    ids = list(PARTICIPATION)
    vals = _fd.participation_ratios([_sp.form(i) for i in ids])
    return {i: float(v) for i, v in zip(ids, vals)}


def self_interaction_fit(I: dict | None = None) -> dict:
    """MP20's P1: one coupling κ₊ for the three C = +1 states, shifts κ₊·I, with the ratios' covariance."""
    I = I or PARTICIPATION
    states = ["0-+", "0++*", "0-+*"]
    need = needed_shifts()
    d = [need[k]["shift"] for k in states]
    iv = [I[READING_A[k]] for k in states]
    C = covariance(states)
    kap = _quad(iv, C, d) / _quad(iv, C, iv)
    res = [x - kap * y for x, y in zip(d, iv)]
    chi = _quad(res, C, res)
    return {"kappa": kap, "chi2": chi, "dof": 2, "p": chi2_sf(chi, 2), "passes": chi2_sf(chi, 2) > 0.05, "needed": d,
            "fitted": [kap * y for y in iv], "chi2_without": _quad(d, C, d), "needed_ratio": d[0] / d[1], "actual_ratio": iv[0] / iv[1]}


def one_coupling_fit(I: dict | None = None, correlated: bool = False) -> dict:
    """One coupling for all five states.  The lab's recorded 170 on 4 used independent errors (``correlated=False``);
    with the shared error of the 0++ carried through all five it is 178 on 4."""
    I = I or PARTICIPATION
    states = list(READING_A)
    need = needed_shifts()
    d = [need[k]["shift"] for k in states]
    iv = [I[READING_A[k]] for k in states]
    C = covariance(states, correlated=correlated)
    kap = _quad(iv, C, d) / _quad(iv, C, iv)
    res = [x - kap * y for x, y in zip(d, iv)]
    chi = _quad(res, C, res)
    return {"kappa": kap, "chi2": chi, "dof": 4, "p": chi2_sf(chi, 4), "chi2_without": _quad(d, C, d), "dof_without": 5,
            "signs_needed": [(-1 if x < 0 else 1) for x in d]}


def minus_sector_fit(I: dict | None = None) -> dict:
    """MP20's P2: a coupling κ₋ for the C = −1 states, fixed by the lowest even twisted form and A1+-; where the others land."""
    I = I or PARTICIPATION
    need = needed_shifts()
    ie, io = I[READING_A["A1+-"]], I[READING_A["A1--"]]
    km, skm = need["A1+-"]["shift"] / ie, need["A1+-"]["shift_error"] / ie
    lam = {k: _sp.form(k)["lambda"] for k in ("t.o.mm.1", "t.o.pm.1", "t.e.pm.1")}
    pred = {k: 4 * (lam[k] + km * I[k]) for k in lam}
    low = min(("t.o.mm.1", "t.o.pm.1"), key=lambda k: pred[k])
    a_odd, a_even = lattice_ratio("A1--"), lattice_ratio("A1+-")
    sig = (pred[low] - a_odd[0]) / math.hypot(a_odd[1], 4 * skm * I[low])
    sig_bare = (pred[low] - a_odd[0]) / a_odd[1]
    n_even2 = (pred["t.e.pm.1"] - a_even[0]) / math.hypot(a_even[1], 4 * skm * I["t.e.pm.1"])
    w = [1 / need[k]["shift_error"] ** 2 for k in ("A1+-", "A1--")]
    d, iv = [need[k]["shift"] for k in ("A1+-", "A1--")], [ie, io]
    kfit = sum(wi * i * x for wi, i, x in zip(w, iv, d)) / sum(wi * i * i for wi, i in zip(w, iv))
    chi = sum(wi * (x - kfit * i) ** 2 for wi, i, x in zip(w, iv, d))
    return {"kappa": km, "kappa_error": skm, "landed": pred, "lowest_odd": low, "sigma": sig, "sigma_lattice_error_alone": sig_bare,
            "second_even_sigma": n_even2, "passes_by_the_letter": abs(sig) <= 2 and n_even2 >= -2,
            "as_a_fit": {"kappa": kfit, "chi2": chi, "dof": 1, "p": chi2_sf(chi, 1)}}


def twist_energy() -> dict:
    """A constant energy μ on every η-twisted form, fitted to the two C = −1 ground states with independent errors, and
    where the four lowest twisted states then sit, in units of the 0++ (each ± 4 times the error of μ)."""
    need = needed_shifts()
    d = [need[k]["shift"] for k in ("A1+-", "A1--")]
    w = [1 / need[k]["shift_error"] ** 2 for k in ("A1+-", "A1--")]
    mu = sum(wi * x for wi, x in zip(w, d)) / sum(w)
    err = 1 / math.sqrt(sum(w))
    chi = sum(wi * (x - mu) ** 2 for wi, x in zip(w, d))
    ids = ["t.o.mm.1", "t.e.pp.1", "t.o.pm.1", "t.e.pm.1"]
    states = {i: 4 * (_sp.form(i)["lambda"] + mu) for i in ids}
    return {"mu": mu, "error": err, "mu_over_edge": mu / EDGE, "chi2": chi, "dof": 1, "p": chi2_sf(chi, 1), "states": states,
            "state_error": 4 * err, "second_odd": states["t.o.pm.1"], "second_even": states["t.e.pm.1"]}


# ------------------------------------------------------------------ the two open forecasts, sealed
TWIST_RULE = (
    "TWIST ENERGY, OPEN FORECASTS. Registered in the lab record MP20 at 08:21 New York time on 14 November 2026, before any lookup: "
    "with a constant energy mu = 0.2176(179) in lambda on every eta-twisted form, fixed by the lattice's two C = -1 ground states "
    "A1+- 9.32(28) and A1-- 9.78(46) over 0++ 3.405(21), a second odd twisted state at 2.87(7) and a second even one at 3.24(7), in "
    "units of the 0++ mass. Criterion, fixed on 16 November 2026 at 06:54:50 New York time (the build record of mtft 1.1.0) and "
    "hash-sealed the same morning, before any lookup: under D4' "
    "(odd forms P = -, even forms P = +, twisted forms C = -1) the second odd twisted state is the first excited state of the "
    "lattice's A1-- channel and the second even one the first excited state of A1+-, as listed in the source's continuum table for "
    "SU(3). For each channel sigma = |forecast - lattice ratio| / sqrt(forecast error^2 + e^2), the lattice ratio being the state's "
    "M/sqrt(sigma) over the 0++ value and e its error with the 0++ error in quadrature. sigma <= 2 is PASS, 2 < sigma <= 3 TENSION, "
    "sigma > 3 FAIL. A channel with no excited state listed has no verdict. A FAIL in either channel closes the constant twist "
    "energy as a description; a PASS is weak support (one fitted number; the spins of the A1 states are not identified).")
TWIST_PARAMS = {"second_odd": 2.8658, "second_even": 3.2399, "forecast_error": 0.0715, "unit": 3.405, "unit_error": 0.021,
                "mu": 0.2176, "mu_error": 0.0179}


def twist_test(excited_odd=None, excited_even=None) -> dict:
    """The sealed rule.  Each argument is (M/√σ, error) of the first excited state of A1-- (odd) or A1+- (even), or None."""
    p = TWIST_PARAMS
    out = {}
    for name, cand, fc in (("A1-- first excited (second odd)", excited_odd, p["second_odd"]),
                           ("A1+- first excited (second even)", excited_even, p["second_even"])):
        if cand is None:
            out[name] = {"forecast": fc, "verdict": "NO DATA"}
            continue
        m, s = float(cand[0]), float(cand[1])
        r = m / p["unit"]
        e = r * math.hypot(s / m, p["unit_error"] / p["unit"])
        sig = abs(fc - r) / math.hypot(p["forecast_error"], e)
        out[name] = {"forecast": fc, "lattice_ratio": r, "lattice_error": e, "sigma": sig,
                     "verdict": "PASS" if sig <= 2 else ("TENSION" if sig <= 3 else "FAIL")}
    verdicts = [v["verdict"] for v in out.values()]
    out["overall"] = ("FAIL" if "FAIL" in verdicts else "NO DATA" if set(verdicts) == {"NO DATA"} else
                      "TENSION" if "TENSION" in verdicts else "PASS" if "NO DATA" not in verdicts else "PASS (one channel)")
    return out


SEALS = (
    {"id": "GB-S1", "sealed": "2026-11-16", "registered": "2026-11-14 08:21 America/New_York (lab record MP20)", "status": "sealed",
     "statement": "The constant twist energy's two open forecasts: a second odd C = -1 state at 2.87(7) and a second even one at "
                  "3.24(7), in units of the 0++ mass, against the first excited states of the lattice's A1-- and A1+- channels.",
     "rule": "mtft.maass.glueballs.twist_test", "constants": ("mtft.maass.glueballs.TWIST_RULE", "mtft.maass.glueballs.TWIST_PARAMS"),
     "params": {"source": "arXiv:2007.06422, SU(3), continuum", "channels": ["A1--", "A1+-"], "thresholds": [2, 3]},
     "window": "sigma <= 2 PASS; 2 < sigma <= 3 TENSION; sigma > 3 FAIL; no excited state listed: NO DATA",
     "sha256": "e75106afbfc3e3bd92c29ebcfeea81db48a0c67ab715c50fdc29dad1735500e7",
     "outcome": "pending", "evaluated": None, "result": "not looked up: the plan's experiment E5"},
)


def seal_hash(entry: dict) -> str:
    """SHA-256 of a seal's id, statement, parameters, window, its rule's code as it stands, and its sealed constants
    (``mtft.matter.prereg``'s recipe)."""
    from mtft.matter.prereg import seal_hash as _h
    return _h(entry)


def verify_seals() -> dict:
    """Recompute each seal from the code as it stands: True where it still matches."""
    return {e["id"]: seal_hash(e) == e["sha256"] for e in SEALS}


def evaluate_open(excited_odd=None, excited_even=None) -> dict:
    """Run the sealed rule on lattice values for the first excited A1-- and A1+- states; refuses a seal that no longer verifies."""
    e = SEALS[0]
    if seal_hash(e) != e["sha256"]:
        raise RuntimeError("GB-S1 no longer verifies: the rule or its text changed after sealing; re-seal under a new id")
    return {"seal": e["id"], "sealed": e["sealed"], "registered": e["registered"], "result": twist_test(excited_odd, excited_even)}


# ------------------------------------------------------------------ a new dictionary: forecast first
def _forecast_hash(rec: dict) -> str:
    body = {k: rec[k] for k in ("name", "statement", "predictions", "registered")}
    return hashlib.sha256(json.dumps(body, sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest()


def register_forecast(name: str, statement: str, predictions: dict, registered: str) -> dict:
    """Seal the predicted ratios M/M(0++) of a new dictionary before they are compared.

    ``predictions`` maps lattice states of ``LATTICE`` to a ratio or to (ratio, error); ``registered`` is the time the
    forecast was made, read from a clock.  Keep the returned record (commit it) and pass it to ``compare_forecast``."""
    if not registered or not str(registered).strip():
        raise ValueError("a forecast needs the time it was registered")
    preds = {}
    for state, v in predictions.items():
        if state not in LATTICE or state == "0++":
            raise ValueError(f"{state!r} is not a lattice state a ratio can be tested against")
        pair = (float(v[0]), float(v[1])) if isinstance(v, (tuple, list)) else (float(v), 0.0)
        preds[state] = [pair[0], pair[1]]
    if not preds:
        raise ValueError("a forecast needs at least one prediction")
    rec = {"name": str(name), "statement": str(statement), "predictions": preds, "registered": str(registered)}
    rec["sha256"] = _forecast_hash(rec)
    return rec


def compare_forecast(record: dict) -> list:
    """Compare a sealed forecast with the lattice, in the form of ``mtft.falsify``.  Refuses a record that was not made by
    ``register_forecast`` or was changed since."""
    from mtft.falsify import Prediction
    if not isinstance(record, dict) or "sha256" not in record or _forecast_hash(record) != record["sha256"]:
        raise ValueError("not a sealed forecast: register it first, with register_forecast, and do not change it")
    rows = []
    for n, (state, (pred, perr)) in enumerate(record["predictions"].items(), start=1):
        lr, le = lattice_ratio(state)
        sig = abs(pred - lr) / math.hypot(le, perr)
        rows.append(Prediction(number=n, relation=f"{record['name']}: M({state})/M(0++)", predicted=pred, observed=lr,
                               obs_error=math.hypot(le, perr), error_percent=abs(pred - lr) / lr * 100, sigma=sig, sector="glueball",
                               source=f"registered {record['registered']}", status="PASS" if sig <= 2 else ("TENSION" if sig <= 3 else "FAIL"),
                               rating="", deviation_ppm=abs(pred - lr) / lr * 1e6, multiplicity="", group=record["name"]))
    return rows


# ------------------------------------------------------------------ the register of forecasts and outcomes
REGISTER = (
    {"id": "MP15-N0", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "H24, reading b", "probability": None,
     "statement": "the lowest cusp form of Gamma_0(143) is even with Fricke sign +: the 0++ is the mass gap", "outcome": "missed",
     "result": "the lowest cusp form is odd"},
    {"id": "MP15-N1", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "H24, reading b", "probability": None,
     "statement": "lambda(0-+)/lambda(0++) in [1.472, 1.627]", "outcome": "missed", "result": "0.870"},
    {"id": "MP15-N2", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "H24, reading b", "probability": None,
     "statement": "lambda(0++*)/lambda(0++) in [1.634, 1.806]", "outcome": "missed", "result": "2.236"},
    {"id": "MP15-N3", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "H24, reading b", "probability": None,
     "statement": "lambda(0-+*)/lambda(0++) in [2.020, 2.262]", "outcome": "missed", "result": "1.247"},
    {"id": "MP15-N4", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "H24, reading b", "probability": None,
     "statement": "no cusp form with Fricke sign - below 2.60 lambda(0++)", "outcome": "missed", "result": "one at 2.152"},
    {"id": "MP15-S1", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "derived", "probability": None,
     "statement": "under reading a the 0++ is the edge of a continuum: no isolated state, no plateau", "outcome": "stands",
     "result": "reading a fails on structure; tested against the paper's description of its plateaux, not its correlators"},
    {"id": "MP15-S2", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "H24, reading a", "probability": 0.1,
     "statement": "the lowest odd cusp form with Fricke sign + lies in R in [0.344, 0.396]", "outcome": "passed",
     "result": "R = 0.37385: 1.559 against the lattice's 1.549(16), 0.6 sigma; the forecast that it would miss (0.9) was wrong"},
    {"id": "MP15-F1", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "forecast", "probability": 0.7,
     "statement": "in every one of the eight sectors the lowest cusp form has R < 1.5", "outcome": "held", "result": "0.37 to 1.41"},
    {"id": "MP15-F2", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "forecast", "probability": 0.9,
     "statement": "reading b misses at least one of N0 to N4", "outcome": "held", "result": "all five"},
    {"id": "MP15-F3", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "forecast", "probability": 0.85,
     "statement": "N4 fails", "outcome": "held", "result": ""},
    {"id": "MP15-F4", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "forecast", "probability": 0.9,
     "statement": "S2 fails: the window is narrow", "outcome": "missed", "result": "S2 passed"},
    {"id": "MP15-F5", "record": "MP15", "registered": "2026-11-12 13:31", "whose": "forecast", "probability": 0.55,
     "statement": "N0 fails: the lowest cusp form overall is odd", "outcome": "held", "result": ""},
    {"id": "MP15-control", "record": "MP15", "registered": "2026-11-12 13:58", "whose": "control", "probability": None,
     "statement": "at levels 133 and 145 no odd cusp form with Fricke sign + lies below R = 0.8", "outcome": "held",
     "result": "the S2 match is specific to 143: a coincidence without a mechanism; a recorded result of the lab's scripts, not recomputed here"},
    {"id": "MP18-P1", "record": "MP18", "registered": "2026-11-13 08:52", "whose": "D4', the test", "probability": 0.1,
     "statement": "no eta-twisted form has R < 0.6325", "outcome": "missed",
     "result": "four do; the lightest C = -1 states at 1.857 and 1.922 times the edge, 7.5 and 9.7 sigma below the lattice; the forecast of a miss (0.9) held"},
    {"id": "MP18-P2", "record": "MP18", "registered": "2026-11-13 08:52", "whose": "two routes", "probability": 0.7,
     "statement": "the twisted trace formula closes on the first full list", "outcome": "missed, then held",
     "result": "three forms were missing from the first list of 39; with all 42 it closes"},
    {"id": "MP18-F1", "record": "MP18", "registered": "2026-11-13 08:52", "whose": "forecast", "probability": 0.75,
     "statement": "the lowest twisted form has R < 0.5", "outcome": "held", "result": "0.4627"},
    {"id": "MP19-P4", "record": "MP19", "registered": "2026-11-13 13:39", "whose": "forecast", "probability": 0.97,
     "statement": "no even form with Fricke sign + lies in R in [0.410, 0.438], where the lattice's 1.720(16) would need one",
     "outcome": "held", "result": "reading a's 0++* stays at 1.792, 4.5 sigma from the lattice"},
    {"id": "MP20-V13", "record": "MP20", "registered": "2026-11-14 07:59", "whose": "check", "probability": None,
     "statement": "the quadrature gives the area to 1e-4 and orthogonality to 1e-5", "outcome": "missed, then held",
     "result": "the first rule gave orthogonality 1.4e-4; on panels between the domain's breakpoints 2e-10"},
    {"id": "MP20-SIU", "record": "MP20", "registered": "2026-11-14 07:59", "whose": "decided by sign", "probability": None,
     "statement": "one self-coupling for every state", "outcome": "excluded",
     "result": "the needed shifts are negative for C = +1 and positive for C = -1; chi2 = 170 on 4"},
    {"id": "MP20-P1", "record": "MP20", "registered": "2026-11-14 07:59", "whose": "test", "probability": 0.25,
     "statement": "one coupling for the three C = +1 states passes at p > 0.05", "outcome": "failed",
     "result": "chi2 = 10.2 on 2, p = 0.006; the forecast of a fail (0.75) held"},
    {"id": "MP20-P2", "record": "MP20", "registered": "2026-11-14 07:59", "whose": "test", "probability": 0.65,
     "statement": "a coupling for C = -1 fixed by A1+- puts the lowest odd twisted state within 2 sigma of A1--",
     "outcome": "passed by the letter", "result": "-1.99 sigma with the coupling's error, -2.26 with the lattice's alone; no support"},
    {"id": "MP20-F1", "record": "MP20", "registered": "2026-11-14 07:59", "whose": "forecast", "probability": 0.8,
     "statement": "all seven participation ratios lie between 1.5 and 4.5", "outcome": "held", "result": "1.89 to 2.47"},
    {"id": "GB-S1", "record": "MP20", "registered": "2026-11-14 08:21", "whose": "the twist energy", "probability": None,
     "statement": "a second odd C = -1 state at 2.87(7) and a second even one at 3.24(7), in units of the 0++", "outcome": "open",
     "result": "sealed on 2026-11-16 before any lookup"},
)


def register() -> list:
    """Every registered forecast of the lattice comparison with its outcome (copies)."""
    return [dict(r) for r in REGISTER]


def summary() -> dict:
    """The verdicts in one place, recomputed."""
    rb, ra, tw = reading_b(), reading_a_checks(), twist_energy()
    comp = {p.relation.split("vs ")[1]: (round(p.sigma, 2), p.status) for p in comparison()}
    return {"reading_b_all_miss": rb["all_miss"], "S2_holds": ra["S2"]["holds"], "eta_twist_bound_holds": ra["MP18-P1"]["holds"],
            "comparison": comp, "self_interaction_chi2": self_interaction_fit()["chi2"], "one_coupling_chi2": one_coupling_fit()["chi2"],
            "twist_energy": (tw["mu"], tw["error"], tw["chi2"]), "open": (tw["second_odd"], tw["second_even"], tw["state_error"]),
            "seals": verify_seals(), "status": "H24 as written fails in both readings; no claim about Yang–Mills follows"}
