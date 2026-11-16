r"""mtft.lines.table — the lab's frozen record of the closed lines of Γ₀(143), with SHA-1 provenance (LN-01, 2026-11-16).

``_data/lines_143.json`` holds what the Modular Plane lab found on 12 November 2026 (records MP12,
MP13 and MP14): the 29 generators of Γ₀(143) from the graph on the darts, the thirty integers of
each of the 28 closed lines, the crossing matrices of the 24 geodesics they trace and of 14 more.
The package computes all of it again with its own code (``mtft.lines.closed``, ``periods``,
``crossings``); this record is what the gates compare with, value for value.

One part is used and not only compared: ``coord``, the integer coordinates of the 29 generators on
the four Hecke blocks in the lab's bases (``mtft.lines.periods.coordinates``).  A basis of a
lattice is a choice; freezing the lab's keeps its table of periods reproducible number for number,
and ``mtft.lines.periods.compare_bases`` checks the choice against the package's own Manin symbols.
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path

_DATA = Path(__file__).resolve().parent / "_data"
FILES = ("lines_143.json",)
BLOCKS = ("ell", "old", "q4", "q6")
BLOCK_DIMS = {"ell": 2, "old": 4, "q4": 8, "q6": 12}
BLOCK_NAMES = {"ell": "143a1", "old": "the level-11 oldspace (the ghost of X₀(11))", "q4": "f2 (quartic)", "q6": "f3 (sextic)"}


@lru_cache(maxsize=1)
def _load() -> dict:
    return json.loads((_DATA / "lines_143.json").read_text(encoding="utf-8"))


def verify_provenance() -> dict:
    """SHA-1 of the data file against PROVENANCE.txt."""
    want = {}
    for line in (_DATA / "PROVENANCE.txt").read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 2 and len(parts[0]) == 40:
            want[parts[1]] = parts[0]
    return {name: hashlib.sha1((_DATA / name).read_bytes()).hexdigest() == want.get(name) for name in FILES}


def record(key: str | None = None):
    """The frozen record, or one of its parts (a fresh copy: change it freely)."""
    rec = _load()
    return json.loads(json.dumps(rec if key is None else rec[key]))


def metadata() -> dict:
    return record("meta")


def lines() -> list:
    """The lab's 28 closed lines: motion, ℓ, dart, Ψ, windings, class, periods by block, orbit lengths modulo 11 and 13."""
    return record("lines")


def generators() -> list:
    """The lab's 29 generators of Γ₀(143), as tuples of rows."""
    return [tuple(tuple(int(v) for v in row) for row in g) for g in _load()["generators"]]


@lru_cache(maxsize=1)
def coord() -> dict:
    """Block name → 29 rows of integer coordinates, one row per generator."""
    c = _load()["coord"]
    out = {name: tuple(tuple(int(v) for v in row) for row in c[name]) for name in BLOCKS}
    for name in BLOCKS:
        if len(out[name]) != 29 or any(len(row) != BLOCK_DIMS[name] for row in out[name]):
            raise RuntimeError("the frozen block coordinates have the wrong shape")
    return out


def motions() -> dict:
    """The six motions of the records, by name."""
    return {k: tuple(tuple(int(v) for v in row) for row in m) for k, m in _load()["motions"].items()}
