"""mtft.maass.spectrum — the frozen Maass spectrum of Γ₀(143) (CERTIFIED; data with SHA-1 provenance; MA-01, 2026-11-15).

The table holds the 41 Maass cusp forms of Γ₀(143) with R < 2.0 (λ = ¼ + R²) and the 42 forms of
the η-twist with R < 1.8, each with its parity, its Atkin–Lehner signs (ε₁₁, ε₁₃) and the first
400 coefficients of its expansion at the cusp ∞ (``mtft.maass.hejhal``).  Completeness: by the
package's own trace formulas (``mtft.maass.trace``, lengths to 17) no form is missing below
R = 1.95 in any untwisted sector, nor below R = 1.60 in the twisted lists together; the lab's
records MP17 to MP19, with lengths to 23.5, give 1.98 and 1.8, which stay recorded results.

Every eigenvalue in the table was found again by the package's own solver when the table was
built (``studies/maass_2026nov/build_frozen_spectrum.py``): the lab's value is kept beside it as
``R_lab``.  ``rederive`` repeats that for one form at call time.

The first form.  The lowest eigenvalue is λ₁ = ¼ + 0.3738518157² = 0.3897652, an odd form with
signs (−, −).  Nothing lies in (0, ¼): Selberg's conjecture is a theorem at this level.

Form ids read ``u.o.mm.1``: untwisted or twisted (u, t); even or odd (e, o); the signs of ε₁₁ and
ε₁₃ (p, m); the place in its sector by increasing R.
"""
from __future__ import annotations

import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

from .domain import reduce_point
from .kbessel import k_bessel

_DATA = Path(__file__).resolve().parent / "_data"
N = 143


@lru_cache(maxsize=1)
def _load():
    meta = json.loads((_DATA / "maass_143.json").read_text(encoding="utf-8"))
    npz = np.load(_DATA / "maass_143_coefficients.npz")
    ids = [str(s) for s in npz["ids"]]
    if ids != [f["id"] for f in meta["forms"]]:
        raise RuntimeError("the coefficient file and the table list different forms")
    return meta, {i: npz["c"][k] for k, i in enumerate(ids)}


def metadata() -> dict:
    """Everything in the table's header: the group, the expansion, the solver's setup, the controls."""
    meta = _load()[0]
    return {k: v for k, v in meta.items() if k != "forms"}


def verify_provenance() -> dict:
    """SHA-1 of the two data files against PROVENANCE.txt."""
    want = {}
    for line in (_DATA / "PROVENANCE.txt").read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 2 and len(parts[0]) == 40:
            want[parts[1]] = parts[0]
    return {name: hashlib.sha1((_DATA / name).read_bytes()).hexdigest() == want.get(name)
            for name in ("maass_143.json", "maass_143_coefficients.npz")}


def _attach(f: dict) -> dict:
    out = dict(f)
    out["eps"] = (f["eps11"], f["eps13"])
    out["c"] = _load()[1][f["id"]]
    return out


def forms(twist=None, parity=None, eps=None, r_max: float | None = None, newform=None) -> list:
    """The forms of the table, by increasing R, filtered.

    ``twist``: False for the untwisted forms, True for the η-twisted ones, None for both.
    ``parity``: 0 (even) or 1 (odd).  ``eps``: (ε₁₁, ε₁₃).  Each form is a dictionary carrying its
    coefficient vector under ``"c"``, ready for ``mtft.maass.ford``.
    """
    out = []
    for f in _load()[0]["forms"]:
        if twist is not None and f["twist"] != twist:
            continue
        if parity is not None and f["parity"] != parity:
            continue
        if eps is not None and (f["eps11"], f["eps13"]) != tuple(eps):
            continue
        if r_max is not None and f["R"] >= r_max:
            continue
        if newform is not None and f["newform"] != newform:
            continue
        out.append(_attach(f))
    return sorted(out, key=lambda f: f["R"])


def form(fid: str) -> dict:
    """One form by its id."""
    for f in _load()[0]["forms"]:
        if f["id"] == fid:
            return _attach(f)
    raise KeyError(f"no form {fid!r} in the table")


def first(twist: bool = False) -> dict:
    """The lowest form: at twist = False the bottom of the discrete spectrum of X₀(143)."""
    return forms(twist=twist)[0]


def gap() -> dict:
    """The spectral gap of X₀(143): the continuum's edge ¼ and the first cusp form above it."""
    f = first(False)
    return {"edge": 0.25, "lambda_1": f["lambda"], "R_1": f["R"], "id": f["id"], "parity": "odd" if f["parity"] else "even",
            "eps": f["eps"], "class": "CERTIFIED: two horocycles, Hecke relations; no form missing below R = 1.95 by mtft.maass.trace"}


def sectors(twist: bool = False) -> dict:
    """The number of forms in each sector (parity, ε₁₁, ε₁₃)."""
    out = {}
    for f in forms(twist=twist):
        key = ("odd" if f["parity"] else "even", f["eps11"], f["eps13"])
        out[key] = out.get(key, 0) + 1
    return out


def _trig(parity):
    return np.cos if parity == 0 else np.sin


def series(f: dict, x, y) -> np.ndarray:
    """The form's Fourier series at the points (x, y), with no reduction.

    Use it inside the Ford domain, or a little below it: with 400 coefficients the series is
    good to about 1e-9 down to y = 0.0088.  For any point of the plane use ``evaluate``.
    """
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    out = np.empty(x.shape)
    c, R, trig = f["c"], f["R"], _trig(f["parity"])
    shift = 0.5 if f["twist"] else 0.0
    for i in np.ndindex(x.shape):
        M = min(len(c), int(36.0 / (2 * math.pi * y[i])) + 2)
        fr = np.arange(1, M + 1) - shift
        out[i] = math.sqrt(y[i]) * float((c[:M] * k_bessel(R, 2 * np.pi * fr * y[i])) @ trig(2 * np.pi * fr * x[i]))
    return out


def evaluate(f: dict, x, y) -> np.ndarray:
    """The form's value at any points of the upper half-plane.

    Each point is carried into the Ford domain of Γ₀*(143); the value there is multiplied by the
    character of that motion (the Atkin–Lehner signs, and the η-twist's sign for a twisted form).
    """
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    out = np.empty(x.shape)
    for i in np.ndindex(x.shape):
        xs, ys, par, sign = reduce_point(float(x[i]), float(y[i]), N, True)
        chi = (f["eps11"] if par[0] else 1) * (f["eps13"] if par[1] else 1) * (sign if f["twist"] else 1)
        out[i] = chi * series(f, xs, ys)
    return out


def rederive(f, half_width: float = 1e-6) -> dict:
    """Find the form's eigenvalue again with the solver, on the records' two horocycles.

    Returns the solver's result with ``difference`` = R found − R frozen.  About three seconds.
    """
    from . import hejhal
    if isinstance(f, str):
        f = form(f)
    r = hejhal.refine(f["R"], f["parity"], f["eps"], f["twist"], half_width=half_width)
    r["difference"] = r["R"] - f["R"]
    r["id"] = f["id"]
    return r
