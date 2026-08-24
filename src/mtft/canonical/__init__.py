"""mtft.canonical — the canonical ideal of X0(143) and the Atkin-Lehner descent.

Frozen certified data from the 2026-08-16 canonical-ideal arc, with the
pre-registered predictions re-derived from that data at call time rather
than stored as assertions.  See `mtft.canonical.gates`.

The curve.  X0(143) has genus 13 and is non-hyperelliptic, so the canonical
map embeds it in P^12 as a curve of degree 2g-2 = 24.  By Petri the ideal
is cut by (g-2)(g-3)/2 = 55 quadrics, and those quadrics generate it.

The descent.  W_11 and W_13 are automorphisms, so they act on the ideal
(Hecke operators do not).  In the adapted basis both are diagonal, the 55
quadrics split 26/5/4/20 across the four Atkin-Lehner sectors, and
pi_* K_X descends to four eigen-line-bundles on the genus-1 quotient
E = X0(143)* = 143a1 of degrees (0, 6, 5, 1).

SECTOR ORDERING — PINNED.  Two orderings are in circulation across the
corpus.  This module uses `SECTOR_ORDER` throughout:

    ("(+,+)", "(+,-)", "(-,+)", "(-,-)")   ->  S_2 dims (1, 6, 5, 1)

The v0.16.0 CHANGELOG quotes the same data in the order
((+,+), (-,+), (+,-), (-,-)) -> (1, 5, 6, 1).  Both are correct; they are
not independent results.  Anything comparing the two must transpose the
middle pair.  Use `reorder_sectors()` rather than reindexing by hand.

Provenance: `_data/PROVENANCE.txt` carries SHA-1s matching the byte-
identical wave artifacts in `studies/ci_2026aug/`.
"""

from __future__ import annotations

from itertools import combinations_with_replacement
from pathlib import Path

__all__ = [
    "SECTOR_ORDER",
    "SECTOR_ORDER_CHANGELOG",
    "GENUS",
    "LEVEL",
    "COORDINATE_LABELS",
    "MONOMIALS",
    "DESCENT",
    "PREDICTIONS",
    "reorder_sectors",
    "monomial_sector",
    "s2_qexpansions",
    "adapted_basis",
    "adapted_qexpansions",
    "ideal_basis",
    "ideal_by_sector",
    "ci_a_quadric",
    "data_path",
    "IDEAL_BASIS_FRAME",
    "ideal_basis_adapted",
]

LEVEL = 143
GENUS = 13

#: FRAME NOTE (v0.19.0).  `ideal_basis()` loads X0_143_I2_quadric_basis.txt,
#: whose coordinates are the *s2* basis (as that file's header states);
#: COORDINATE_LABELS describes the *adapted* basis.  Mixing them was the
#: v0.18.0 pitfall documented in certificate v6.  For adapted-frame
#: quadrics use `ideal_basis_adapted()`; for mod-p geometry use
#: `mtft.canonical.integral`.
IDEAL_BASIS_FRAME = "s2"

#: Pinned sector ordering for this module.  S_2 dims (1, 6, 5, 1).
SECTOR_ORDER = ("(+,+)", "(+,-)", "(-,+)", "(-,-)")

#: The ordering used in the v0.16.0 CHANGELOG.  S_2 dims (1, 5, 6, 1).
SECTOR_ORDER_CHANGELOG = ("(+,+)", "(-,+)", "(+,-)", "(-,-)")

#: (Galois block, sector) for the adapted coordinates y_1 .. y_13.
COORDINATE_LABELS = (
    ("f1", "(+,+)"),
    *(("f3", "(+,-)"),) * 6,
    ("old+", "(-,+)"),
    *(("f2", "(-,+)"),) * 4,
    ("old-", "(-,-)"),
)

#: The 91 degree-2 monomials as 0-based coordinate pairs (i, j), i <= j.
MONOMIALS = tuple(combinations_with_replacement(range(GENUS), 2))

_SECTOR_INDEX = {s: i for i, s in enumerate(SECTOR_ORDER)}
# character multiplication on (Z/2)^2, indexed by SECTOR_ORDER
_MULT = ((0, 1, 2, 3), (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0))

#: Everything the four sessions established about the descent, in one table.
DESCENT = {
    "quotient": "E = X0(143)* = 143a1",
    "quotient_genus": 1,
    "quotient_rank": 1,
    "quotient_torsion": 1,
    "weierstrass": {"a1": 0, "a2": -1, "a3": 1, "a4": -1, "a6": -2},
    "modular_degree": 4,
    "origin": "common image of the four cusps (a single free AL orbit)",
    "sector_dims_S2": {"(+,+)": 1, "(+,-)": 6, "(-,+)": 5, "(-,-)": 1},
    "bundle_degrees": {"(+,+)": 0, "(+,-)": 6, "(-,+)": 5, "(-,-)": 1},
    "fixed_points": {"W11": 0, "W13": 4, "W143": 20},
    "branch_degrees": {"W11": 0, "W13": 2, "W143": 10},
    "quotient_genera": {"W11": 7, "W13": 6, "W143": 2},
    "W11_acts_freely": True,
    "cusp_action": {"W11": "(1 11)(13 143)", "W13": "(1 13)(11 143)"},
    "branch_points": {"Q1": "(2i, -3+2i)", "Q2": "(-2i, -3-2i)", "disc": -52},
    "L_minusminus": {
        "point": (4, -7),
        "note": "a generator of E(Q) = Z; unique rational square root of Q1+Q2",
        "section": "g(q) - 13*g(q^13), g = 11a1 = eta(t)^2 eta(11t)^2",
    },
    "L_minusplus_partner_section": "g(q) + 13*g(q^13)",
    "C2": {
        "curve": "X0(143)/W143, genus 2",
        "model": "w^2 = (y - x + 3)/(y + 3x - 5) over 143a1",
        "jacobian": "~ 143a1 x 11a1",
    },
    "ideal_grading": {"(+,+)": 26, "(+,-)": 5, "(-,+)": 4, "(-,-)": 20},
    "h0_2K_grading": {"(+,+)": 12, "(+,-)": 6, "(-,+)": 7, "(-,-)": 11},
    "monomials_per_sector": {"(+,+)": 38, "(+,-)": 11, "(-,+)": 11, "(-,-)": 31},
    "corrections": ("CC-08", "CC-09"),
}

#: Pre-registered predictions P1-P9, re-derived by `gates`, not asserted.
PREDICTIONS = {
    "P1_dim_S2": 13,
    "P2_dim_Sym2": 91,
    "P3_h0_2K": 36,
    "P4_dim_I2": 55,
    "P5_dim_Sym3": 455,
    "P6_h0_3K": 60,
    "P7_dim_I3": 395,
    "P8_rank_V_I2": 395,
    "P9_sturm_weight4": 56,
}

_DATA = Path(__file__).resolve().parent / "_data"


def data_path(name: str) -> Path:
    """Absolute path to a frozen data file."""
    p = _DATA / name
    if not p.exists():
        raise FileNotFoundError(
            f"{name} missing from mtft/canonical/_data. If this is an sdist, "
            "check that MANIFEST.in ships package data (see the v0.16.1 "
            "hotfix: studies/tests globs previously dropped .txt and .json)."
        )
    return p


def _rows(name: str):
    """Yield (label, [int, ...]) for a shipped labelled-CSV data file."""
    header_seen = False
    for line in data_path(name).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("##"):
            continue
        parts = line.split(",")
        if not header_seen:
            header_seen = True
            continue
        yield parts[0], [int(v) for v in parts[1:]]


def _matrix(name: str, nrows: int, ncols: int):
    out = [r for _, r in _rows(name)]
    if len(out) != nrows or any(len(r) != ncols for r in out):
        raise ValueError(f"{name}: expected {nrows}x{ncols}, got "
                         f"{len(out)}x{len(out[0]) if out else 0}")
    return out


def s2_qexpansions():
    """PARI mfbasis of S_2(Gamma_0(143)), q^0..q^140.  141 x 13 integers."""
    return _matrix("X0_143_S2_qexpansions.txt", 141, GENUS)


def adapted_basis():
    """B, 13 x 13 integers.  e_k = sum_i B[i][k] f_i.  det B = -1078272."""
    return _matrix("X0_143_AL_adapted_basis.txt", GENUS, GENUS)


def adapted_qexpansions():
    """q-expansions of e_1..e_13, q^0..q^140.  141 x 13 integers."""
    return _matrix("X0_143_AL_adapted_qexpansions.txt", 141, GENUS)


def ideal_basis():
    """I_2 in the original x-coordinates.  91 monomials x 55 quadrics."""
    return _matrix("X0_143_I2_quadric_basis.txt", 91, 55)


def ideal_by_sector():
    """I_2 split by Atkin-Lehner class, in adapted y-coordinates.

    Returns {sector: (monomial_indices, columns)} where `columns` is a list
    of dicts mapping a monomial index (into MONOMIALS) to its coefficient.
    Keys follow SECTOR_ORDER.
    """
    text = data_path("X0_143_I2_by_AL_sector.txt").read_text().splitlines()
    out, sector, idx, cols = {}, None, None, None
    for line in text:
        line = line.strip()
        if line.startswith("## class"):
            if sector is not None:
                out[sector] = (idx, cols)
            sector = line.split()[2]
            idx, cols = [], None
            continue
        if not line or line.startswith("#") or sector is None:
            continue
        parts = line.split(",")
        if parts[0] == "monomial":
            cols = [dict() for _ in parts[1:]]
            continue
        m = _parse_monomial(parts[0])
        idx.append(m)
        for c, v in enumerate(parts[1:]):
            v = int(v)
            if v:
                cols[c][m] = v
    if sector is not None:
        out[sector] = (idx, cols)
    return {s: out[s] for s in SECTOR_ORDER}


def ci_a_quadric():
    """Q*, the unique quadric on the f1+f2 projection.

    Returns {monomial_index: coefficient}.  Not f2-specific: see
    studies/ci_2026aug/X0_143_CI_AB_REPORT.md, the decoy result.
    """
    out = {}
    for label, vals in _rows("X0_143_CI_A_quadric.txt"):
        out[_parse_monomial(label)] = vals[0]
    return out


def _parse_monomial(label: str) -> int:
    """'y1^2' or 'x1x2' -> index into MONOMIALS."""
    s = label.replace("^2", "")
    lead = s[0]
    parts = [int(p) for p in s.split(lead) if p]
    if len(parts) == 1:
        parts = parts * 2
    i, j = sorted(p - 1 for p in parts)
    return MONOMIALS.index((i, j))


def monomial_sector(m: int) -> str:
    """Atkin-Lehner sector of monomial index `m` in adapted coordinates."""
    i, j = MONOMIALS[m]
    a = _SECTOR_INDEX[COORDINATE_LABELS[i][1]]
    b = _SECTOR_INDEX[COORDINATE_LABELS[j][1]]
    return SECTOR_ORDER[_MULT[a][b]]


def reorder_sectors(values, frm=SECTOR_ORDER, to=SECTOR_ORDER_CHANGELOG):
    """Move a 4-tuple between the two sector orderings in the corpus.

    >>> reorder_sectors((1, 6, 5, 1))
    (1, 5, 6, 1)
    """
    lookup = dict(zip(frm, values))
    return tuple(lookup[s] for s in to)


def ideal_basis_adapted():
    """The 55 quadrics in the ADAPTED frame, integrally saturated.

    91 x 55 rows-by-columns like `ideal_basis()`, but expressed in the
    coordinates that COORDINATE_LABELS describes, and saturated as an
    integer lattice.  See the data-file header for the mixed-model
    warning; use `mtft.canonical.integral.count_points_modp` for counts.
    """
    out = []
    for _label, vals in _rows("X0_143_I2_adapted_saturated.txt"):
        out.append(vals)
    return out
