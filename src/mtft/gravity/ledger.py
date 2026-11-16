r"""mtft.gravity.ledger — what the gravity work in this package claims, what the gravity map records against itself, and every miss (v1.1.0, 2026-11-16).

The matter programme has a ledger (``mtft.matter.ledger``) checked against its puzzle map.  The gravity map was
written by hand and had none; this is it, in the same format: an entry is (id, wave, class, status, statement,
module, test, note), with the matter ledger's classes and statuses.

  * ``E``          the entries for ``mtft.gravity`` (waves GR-01 the tube, GR-02 the lines of force, GR-03 the tetrad),
                   three of them retracted statements that stay exposed
  * ``PIECES``     the map's eleven pieces with their status
  * ``TESTS``      its pre-registered tests T-G1 to T-G10, in the map's own words
  * ``NEGATIVES``  its sixteen honest negatives N-G1 to N-G16, in the map's own words
  * ``MISSES``     the recorded misses that are not among the negatives, with the ones this release's build found
  * ``OPEN``       what is open

The negatives and the tests are records with their own ids, not rows of ``E``: a row needs a class, and the map gives
them none; to class them here would be this module's reading of the map, not the map.

The map's text is carried, not retyped: ``_data/gravity_map.json`` is parsed from the lab's map of 14 November 2026
(``studies/gravity_2026nov``), whose SHA-1 is checked by ``check``.  Statements about X₀(143) in the map (T-G4, T-G9,
T-G10, N-G13, N-G14) are tested in ``mtft.lines`` and entered in the matter ledger as LN-1 to LN-8; they are listed
here because the map lists them.

Nothing here is a claim that gravity is derived.  N-G1 is the first negative: Newton's constant is not.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from functools import lru_cache
from pathlib import Path

CLASSES = ("EXACT", "NUMERICAL", "STRUCTURAL", "ESTIMATE", "MODEL")
STATUSES = ("active", "caveat", "superseded", "retracted")
WAVES = {"GR-01": ("2026-11-16", "the tube round a core line (lab records MP09, MP21)"),
         "GR-02": ("2026-11-16", "the squared hopfion's lines of force after t = 0 (lab record MP22)"),
         "GR-03": ("2026-11-16", "a tetrad frozen into the null flow, and its helicity (lab record MP23)")}
REPORT = "docs/gravity/KNOTTED_FIELDS.md"
MAP_SHA1 = "8730ed451d2fc3f6c85a5c281fa95b2c53ccab11"            # of the lab's map as parsed (14 November 2026)
DATA_SHA1 = "1b61fb857070f007106e680b1cb104e707f212bf"           # of _data/gravity_map.json

E: list = []


def _e(i, wave, cls, status, text, module, test, note=""):
    E.append({"id": i, "wave": wave, "class": cls, "status": status, "statement": text, "module": module, "test": test, "note": note})


_K, _H, _T = "mtft.gravity.knotted", "mtft.gravity.hopfion", "mtft.gravity.tetrad"
_e("GR-1", "GR-01", "EXACT", "active", "Next to a core line of the (p, q) field the lines turn by 2 pi nu per circuit in the torus frame, nu = sqrt(2pq(p + q)): the cross-section is a Hamiltonian system for 2 ln H and its linear part at the core line has determinant 2pq/(p + q). The circuit matrix is the identity exactly when 2pq(p + q) is a square: among coprime p <= q <= 60, (1, 1), (1, 8), (1, 49), (32, 49).", _K, "test_gravity_knotted.py", "the fields, their invariant Re(alpha^p beta^q) and its tori are Kedia, Bialynicki-Birula, Peralta-Salas and Irvine's (2013); nu was not found in the sources read on 2026-11-16, which is not the same as new; derived in MP09, unreviewed")
_e("GR-2", "GR-01", "EXACT", "active", "The shear at the core line, with u = 1 - h/h_max: d iota/du = nu [3pq(p + q) - 2(p^2 + pq + q^2)]/(12pq(p + q)), zero for the hopfion only; 35/36 for (1, 8).", _K, "test_gravity_knotted.py", "three routes: the closed form (exact for symbolic p, q), the profile integral at 30 digits (4.9e-10 of nu), and the section integrated as a differential equation; derived in MP21, unreviewed")
_e("GR-3", "GR-01", "NUMERICAL", "active", "On the torus H = h the winding iota(h) rises from nu at the core line to 2pq at the edge.", _K, "test_gravity_knotted.py", "MP09's integral, valid on every torus")
_e("GR-4", "GR-01", "EXACT", "active", "The helicity of the (p, q) field on the label sphere is 2 pi^2 pq p! q!/(p + q)! and its energy is (p + q) times that.", _K, "test_gravity_knotted.py", "by hand and by quadrature to 1e-14; the ratio 1/(p + q) is Kedia et al.'s published value, the one outside anchor of the module's normalisation; not in the lab's records")
_e("GR-5", "GR-02", "EXACT", "active", "At a time t > 0 the lines of force of the squared hopfion are the orbits of H = |P|^2 + |s^2 - sigma^2|^2/t^2 on H = 1, with the second invariant I_2 = L^2 + 4m p_x^2 + (4m/t^2) x^2 (r^2 - 2m): {H, I_2} = 0 for symbolic t. They are never chaotic.", _H, "test_gravity_hopfion.py", "that the lines of force are the light's lines turned by arg(d)/2 is Swearngin et al.'s (2013); the potential reached is Garnier's integrable quartic; the reduction itself was not found in the sources read; derived in MP22, unreviewed")
_e("GR-6", "GR-02", "NUMERICAL", "active", "In space, on the field itself: the reduction holds pointwise to 1e-14 in angle; I_2 is kept to 1e-13 of its range along lines of force over 100 circuits; the axial great circle and the cross line are closed lines of force, linked -1 (magnetic type) and +1 (electric).", _H, "test_gravity_hopfion.py", "integrated by the package's own routine; the lab's runs with scipy kept I_2 to 2.7e-10")
_e("GR-7", "GR-02", "NUMERICAL", "active", "The winding on a torus is W = T_nu/T_mu, two periods of Y^2 = X(X - 1)Pi(X) depending on t only through eps = t/(2 sqrt(1 + t^2)); to leading order W - 1 = (5/12) eps^2 on every torus (exact on the axial line to that order). For the line through (0.6, -0.3, 0.5) at t = 0.4, W - 1 = 0.015660.", _H, "test_gravity_hopfion.py", "MP22's table to 1e-16 of the lab's file; measured along an orbit to 1e-15 of the integrals")
_e("GR-8", "GR-03", "NUMERICAL", "active", "The tetrad e = dx + B, B_mu = Re[(i a^2 (d_mu b) d gamma - b^2 (d_mu a) d beta)/d^2], is frozen into the null flow of the squared hopfion; the curvature of h = B + B^T is n (x) n/d^5 and its Ricci tensor is zero.", _T, "test_gravity_tetrad.py", "exactly zero differences over the Gaussian rationals at the events tried; a randomized certificate (Schwartz-Zippel, degree at most 22) bounds a false pass by 3e-21; not a symbolic proof; the degree count is unreviewed. Frozen-in gravitational fields for a timelike flow are Asenjo, Winkler and Comisso's (2026); this null-flow frame was not found in the sources read; derived in MP23, unreviewed")
_e("GR-9", "GR-03", "NUMERICAL", "caveat", "In that frame the helicity is H_g = (2/3) pi^2, beside H_em = pi^2 for the light field, on the label sphere and on slices of space, constant in time.", _T, "test_gravity_tetrad.py", "the number belongs to the frame (N-G5): in transverse-traceless gauge Huang, Long, Zhao and Zhou (arXiv:2608.30154) find a helicity of the same kind that changes in time, and the two normalisations have not been compared; the frame is unique only within the lab's holomorphic class up to degree 4; 2/3 has no mechanism (AG-D5)")
_e("GR-10", "GR-02", "NUMERICAL", "retracted", "At t = 0.4 the line through (0.6, -0.3, 0.5) makes one turn in about 29 circuits (MP10).", _H, "test_gravity_hopfion.py", "a miscount of returns; one turn in 63.9 circuits (MP22; N-G15)")
_e("GR-11", "GR-02", "ESTIMATE", "retracted", "At t = 0.8 the line of force through (0.6, -0.3, 0.5) is lost (MP11's reading); it reaches more than 10 from the origin (MP22's forecast, given 0.6).", _H, "test_gravity_hopfion.py", "it was not lost: the search plane missed it for an arc of 69.8; it keeps I_2 on a torus at 0.43 of the range, returns, and reaches 4.7 (MP22; N-G15)")
_e("GR-12", "GR-03", "ESTIMATE", "retracted", "H_g/H_em = 3/4 (a hand value, registered in MP23 before the integral was computed).", _T, "test_gravity_tetrad.py", "a slip of algebra; the ratio is 2/3 (N-G16)")

MISSES = (
    {"id": "M-G1", "record": "MP09", "kind": "hand value", "what": "the plan printed the trace 0.9868 for (2, 5); the formula gives 0.9871991, as the record's own result has it"},
    {"id": "M-G2", "record": "MP21", "kind": "hand value", "what": "three hand values of d iota/dh in the plan were wrong (-15.2404 for -15.0920, -4.4244 for -4.6719, -6.018516 for -6.018519); the closed form was right"},
    {"id": "M-G3", "record": "MP22", "kind": "fault of a test", "what": "a count of crossings included the starting point, so a check printed MISS on a line that closes; corrected in the record"},
    {"id": "M-G4", "record": "MP09", "kind": "found in the 1.1.0 build", "what": "the line field on the label sphere is said to hold 'up to a positive factor'; with B = Im F_pq the factor is negative; no result depends on the sense"},
    {"id": "M-G5", "record": "MP23", "kind": "found in the 1.1.0 build", "what": "the number file holds the two helicities on the label sphere with a minus sign and the test took absolute values: the order (eta, xi_1, xi_2) is the orientation opposite to space's; the package uses (xi_1, eta, xi_2)"},
    {"id": "M-G6", "record": "the release plan of 1.1.0", "kind": "found in the 1.1.0 build", "what": "the plan's gate for the whole-number pairs left out 'coprime'; without it there are 22 pairs for p <= q <= 60, not four"},
)
MISSES_AMONG_NEGATIVES = ("N-G7", "N-G10", "N-G13", "N-G14", "N-G15", "N-G16")
OPEN = {"T-G3": "a three-dimensional example with resistivity: not attempted",
        "H-G7": "the core line of (1, 8) answers a perturbation by its cube root: a hypothesis of MP21, not tested",
        "the Kerr-Schild comparison": "whether the frozen tetrad's metric is the Kerr-Schild potential of the same field in another gauge: not checked",
        "the helicity comparison": "the frozen frame's (2/3) pi^2 against Huang et al.'s transverse-traceless value in one normalisation: not done",
        "the other gravitational (p, q) fields": "no invariant is known for their later lines of force",
        "the frame's uniqueness": "a proof that the frozen frame is unique among all frozen frames: not attempted"}


@lru_cache(maxsize=1)
def _map() -> dict:
    return json.loads((Path(__file__).resolve().parent / "_data" / "gravity_map.json").read_text(encoding="utf-8"))


def pieces() -> list:
    """The map's eleven pieces with their status (copies)."""
    return [dict(p) for p in _map()["pieces"]]


def tests() -> list:
    """The map's pre-registered tests T-G1 to T-G10, in its own words (copies)."""
    return [dict(t) for t in _map()["tests"]]


def negatives() -> list:
    """The map's sixteen honest negatives, in its own words (copies)."""
    return [dict(n) for n in _map()["negatives"]]


def entries(status=None, wave=None) -> list:
    return [dict(e) for e in E if (status is None or e["status"] == status) and (wave is None or e["wave"] == wave)]


def _module_exists(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):                    # a parent package that does not exist raises
        return False


def check(tests_dir=None) -> dict:
    """The ledger's own consistency: every value should be True (``tests_exist`` is None where no tests folder is at hand)."""
    data_file = Path(__file__).resolve().parent / "_data" / "gravity_map.json"
    m = _map()
    ids = [e["id"] for e in E]
    root = Path(__file__).resolve().parents[3]
    tdir = Path(tests_dir) if tests_dir else root / "tests"
    return {
        "ids_unique": len(set(ids)) == len(ids) and ids == [f"GR-{k}" for k in range(1, len(ids) + 1)],
        "classes_and_statuses": all(e["class"] in CLASSES and e["status"] in STATUSES and e["wave"] in WAVES for e in E),
        "modules_exist": all(_module_exists(e["module"]) for e in E),
        "tests_exist": all((tdir / e["test"]).exists() for e in E) if tdir.is_dir() else None,
        "retracted_have_pointers": all(e["note"] for e in E if e["status"] in ("retracted", "superseded", "caveat")),
        "sixteen_negatives": [n["id"] for n in m["negatives"]] == [f"N-G{k}" for k in range(1, 17)],
        "ten_tests_one_pending": [t["id"] for t in m["tests"]] == [f"T-G{k}" for k in range(1, 11)] and [t["id"] for t in m["tests"] if t["status"] != "done"] == ["T-G3"],
        "eleven_pieces": len(m["pieces"]) == 11,
        "misses_point_to_negatives": set(MISSES_AMONG_NEGATIVES) <= {n["id"] for n in m["negatives"]},
        "map_as_parsed": m["source_sha1"] == MAP_SHA1 and hashlib.sha1(data_file.read_bytes()).hexdigest() == DATA_SHA1,
        "open_is_listed": "T-G3" in OPEN and "H-G7" in OPEN,
    }


def to_markdown() -> str:
    """The ledger as a document (``docs/gravity/GRAVITY_LEDGER.md``)."""
    m = _map()

    def cell(s):
        return str(s).replace("|", "\\|").replace("\n", " ")
    out = ["# The gravity ledger", "",
           "Generated from `mtft.gravity.ledger` (`to_markdown()`); do not edit by hand. Entries are in the matter ledger's format. The pieces, the",
           f"pre-registered tests and the negatives are the lab's gravity map's own words (last changed {m['last_changed']}; SHA-1 `{m['source_sha1']}`).",
           "Nothing here is a claim that gravity is derived: N-G1 is the first negative.", "",
           (f"{len(E)} entries ({sum(e['status'] == 'active' for e in E)} active, {sum(e['status'] == 'caveat' for e in E)} with a caveat, "
            f"{sum(e['status'] == 'retracted' for e in E)} retracted and kept); {len(m['tests'])} pre-registered tests "
            f"({sum(t['status'] == 'done' for t in m['tests'])} done, {sum(t['status'] != 'done' for t in m['tests'])} pending); {len(m['negatives'])} negatives; "
            f"{len(MISSES)} further recorded misses."), "",
           "## Entries", "", "| id | wave | class | status | statement | module | test | note |", "|---|---|---|---|---|---|---|---|"]
    out += [f"| {e['id']} | {e['wave']} | {e['class']} | {e['status']} | {cell(e['statement'])} | `{e['module']}` | `{e['test']}` | {cell(e['note'])} |" for e in E]
    out += ["", "## Waves", ""] + [f"- **{w}** ({d}): {t}. Report: `{REPORT}`." for w, (d, t) in WAVES.items()]
    out += ["", "## The map's pieces", ""] + [f"- **{p['piece']}** {p['title']}: {p['status']}" for p in m["pieces"]]
    out += ["", "## Pre-registered tests (the map's words)", "", "| test | registered | criterion | outcome |", "|---|---|---|---|"]
    out += [f"| {t['id']} {cell(t['title'])} | {cell(t['registered'])} | {cell(t['criterion'])} | {cell(t['outcome'])} |" for t in m["tests"]]
    out += ["", "## Honest negatives (the map's words)", "", "| id | pieces | what was closed | record |", "|---|---|---|---|"]
    out += [f"| {n['id']} | {n['pieces']} | {cell(n['statement'])} | {cell(n['record'])} |" for n in m["negatives"]]
    out += ["", "## Further recorded misses", "", f"Among the negatives, {', '.join(MISSES_AMONG_NEGATIVES)} are misses of a forecast or corrections of an earlier statement. Beside them:", "",
            "| id | record | kind | what |", "|---|---|---|---|"]
    out += [f"| {x['id']} | {x['record']} | {x['kind']} | {cell(x['what'])} |" for x in MISSES]
    out += ["", "## Open", ""] + [f"- **{k}**: {v}" for k, v in OPEN.items()]
    return "\n".join(out) + "\n"
