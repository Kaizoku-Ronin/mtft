"""PR-01 (2026-10-10): sealed pre-registrations and the trials counter.

Pre-registration in the programme has been a sentence in a report (T0–T5 in mtft.matter.puzzle.PREREG).  This module makes it mechanical.
A rule is sealed by the SHA-256 of its id, its statement, the source code of the callable that evaluates it, its parameters and its pass
window, dated before any candidate is evaluated; verify() recomputes the hash from the code as it stands, so a rule edited after sealing fails
the test suite until it is re-sealed under a new id with the old entry kept as superseded (append-only, like the CC corrections).  The legacy
tests T0–T5 were stated in prose before evaluation but not hash-sealed; they are listed here as such and not re-sealed after the fact.

The trials counter is the other half of SQ-8 and AG-D5 (proximity is not evidence).  Every property tried as a reason for the level 143 is
logged against the arena of IV.6, the 13 torus levels of SQ-01, pre-registered or not, with the subset it selects.  The report states the honest
denominator: n properties tried, k that single out 143, and the chance baseline — a property selecting m of 13 exchangeable levels includes 143
with probability m/13, so sum(m_i)/13 trials are expected to include 143 and n_1/13 of the single-level ones to be it.  A self-referential
property (one defined from 143's own value, such as 'has the orbit pattern (1, 4, 6)') is scored by how many of the 13 levels pass their own
version.  The list is append-only; the test suite refuses a distinction in star_quotients.distinctions() that is not logged here."""
import hashlib
import inspect
import json
from importlib import import_module

TORI = (57, 58, 65, 77, 82, 91, 118, 123, 141, 142, 143, 145, 155)
LEVEL = 143


# ------------------------------------------------------------------------------------------------------------------------------------ seals
def _callable(path):
    mod, name = path.rsplit(".", 1)
    return getattr(import_module(mod), name)


def seal_hash(entry):
    """SHA-256 of (id, statement, parameters, window, the evaluating code as it stands now, any sealed text constant)."""
    parts = [entry["id"], entry["statement"], json.dumps(entry["params"], sort_keys=True), entry["window"], inspect.getsource(_callable(entry["rule"]))]
    for const in entry.get("constants", ()):
        parts.append(str(_callable(const)))
    return hashlib.sha256("\n".join(parts).encode("utf-8")).hexdigest()


SEALS = (
    {"id": "PR-S1", "sealed": "2026-10-10", "piece": "I", "sub": "I.4", "status": "superseded by PR-S2 (CC-63)",
     "statement": "Strong CP after the lift: the lifted strong coupling tau_s passes without an axion iff it lies on a theta = 0 mirror of the "
                  "modular tiling and arg det(M_u M_d) = 0 mod 2 pi; a theta = pi mirror or no mirror fails.",
     "rule": "mtft.matter.time_reversal.strong_cp_test", "constants": ("mtft.matter.time_reversal.STRONG_CP_RULE",),
     "params": {"tolerance": None, "candidates": "arithmetic (cusp or CM form); numerical candidates classify to tolerance only"},
     "window": "mirror_class in {theta_0, both} and det_phase 0 mod 2 pi -> PASS; theta_pi -> FAIL (theta-bar = pi); None -> FAIL (axion needed)",
     "sha256": "fb959ccad7fea62c602ffb74951e8e2dffb161408bf578f56894a89e797d8444",
     "outcome": "superseded", "evaluated": None,
     "result": "never run; the rule rejected theta = pi points outright, but theta-bar = theta + arg det(M_u M_d), so a theta = pi point with "
               "arg det = pi passes — found by independent review on the day of sealing (CC-63); the code it hashed is kept in its history"},
    {"id": "PR-S2", "sealed": "2026-10-10", "piece": "I", "sub": "I.4", "status": "sealed",
     "statement": "Strong CP after the lift: with theta-bar = theta + arg det(M_u M_d), the lifted strong coupling tau_s passes without an axion "
                  "iff it lies on a mirror of the modular tiling whose class matches the quark-mass phase (theta = 0 with arg det = 0, or "
                  "theta = pi with arg det = pi, mod 2 pi); any other combination is theta-bar = pi and fails; no mirror fails unless an axion exists.",
     "rule": "mtft.matter.time_reversal.strong_cp_test", "constants": ("mtft.matter.time_reversal.STRONG_CP_RULE",),
     "params": {"tolerance": None, "candidates": "arithmetic (CM form); numerical candidates classify to tolerance only; a cusp is not a candidate"},
     "window": "theta-bar = 0 mod 2 pi for some mirror through tau_s -> PASS; theta-bar = pi -> FAIL; no mirror -> FAIL (axion needed)",
     "sha256": "c95f415adff1f6bac3cec9a5aaccb5687e945b27ac35797524f3f40ad22b9b7d",
     "outcome": "pending", "evaluated": None, "result": "no candidate tau_s exists yet (PV-01's lift is open)"},
    {"id": "PR-S3", "sealed": "2026-10-11", "piece": "IX", "sub": "IX.10", "status": "sealed",
     "statement": "The touch rule (DS-01, 2026-10-07): a condensate on a block with charge vector q_s reaches a field with charge vector q_f at "
                  "first order iff q_s . q_f != 0.  A parent's block menu passes iff some Standard-Model-singlet block of the menu (colour and weak "
                  "singlet, Y = 0, q != 0) has q . q_f != 0 for q_f the d^c line (B <- c) or the H_d line (B -> L): the down sector is reachable.",
     "rule": "mtft.matter.forces.touch_test", "constants": ("mtft.matter.forces.TOUCH_RULE",),
     "params": {"menus": ["U8", "SO16", "E8"], "hypercharge": "the stack assignment (1/6, -1/2, 0, 1/2, -1/2)", "orientation": "PS-01 (Hom(j -> i) = e_i - e_j)"},
     "window": "a singlet block with nonzero dot product against d^c or H_d -> PASS (reachable); none -> FAIL (unreachable at first order)",
     "sha256": "b14d0c86d7035e3dd9cb6bd6d24366f85c3fbfd97bea9ad65ca04c9223b34c7a",
     "outcome": "mixed", "evaluated": "2026-10-11",
     "result": "U8 FAIL (the only singlet line is nu^c on A1-A2, DS-1 reproduced); SO16 PASS (B+A2, A1+B, L+L); E8 PASS (the sum lines and six spinor "
               "sheets) — but neither the SO(16) nor the E8 column gives the exact Standard Model (GF-11); PV-01's non-supersymmetric U(8) carries the "
               "sum lines as fermions only and was not scanned"},
    {"id": "PR-S4", "sealed": "2026-10-11", "piece": "V", "sub": "V.7", "status": "superseded by PR-S5 (CC-64)",
     "statement": "Charge quantisation: a candidate (a block menu with stack hypercharges, Y(Q) = 1/6) passes iff the generator (omega, -1, "
                  "e^{2 pi i/6}) of the centre Z3 x Z2 x U(1)_Y acts trivially on every block that carries net chirality: t/3 + d/2 + Y is an integer "
                  "for each (t triality, d weak duality).  Then the gauge group is [SU(3) x SU(2) x U(1)]/Z6 or a further quotient and every free "
                  "particle carries a whole electric charge, Q = (2/3) t mod 1.",
     "rule": "mtft.matter.forces.z6_test", "constants": ("mtft.matter.forces.Z6_RULE",),
     "params": {"blocks": "the Standard-Model lines of the SM grouping, or any list of charge vectors", "hypercharge": "stack assignment with Y(Q) = 1/6"},
     "window": "every block Z6-trivial -> PASS; any block with a nonzero phase -> FAIL",
     "sha256": "83448da1318801953ce2bd8b54116dec76f579224f9bb1dbb7965513c588d201",
     "outcome": "pass", "evaluated": "2026-10-11",
     "result": "run once (PASS on the Standard-Model lines) and superseded the same day: the statement's reading of the global form was backwards — "
               "Z6 trivial on the matter allows G/Z6 or any of its covers, not a further quotient (independent review, CC-64); the criterion itself "
               "is unchanged and re-sealed as PR-S5; the code it hashed is kept in its history"},
    {"id": "PR-S5", "sealed": "2026-10-11", "piece": "V", "sub": "V.7", "status": "sealed",
     "statement": "Charge quantisation: a candidate (a block menu with stack hypercharges, Y(Q) = 1/6) passes iff the generator (omega, -1, "
                  "e^{2 pi i/6}) of the centre Z3 x Z2 x U(1)_Y, U(1)_Y acting by e^{2 pi i Y}, acts trivially on every block that carries net "
                  "chirality: t/3 + d/2 + Y is an integer for each (t triality, d weak duality).  Then the matter allows the gauge group "
                  "[SU(3) x SU(2) x U(1)]/Z6 or any of its covers, and every free particle carries a whole electric charge, Q = (2/3) t mod 1.",
     "rule": "mtft.matter.forces.z6_test", "constants": ("mtft.matter.forces.Z6_RULE",),
     "params": {"blocks": "the Standard-Model lines of the SM grouping, or any list of charge vectors", "hypercharge": "stack assignment with Y(Q) = 1/6"},
     "window": "every block Z6-trivial -> PASS; any block with a nonzero phase -> FAIL",
     "sha256": "54da8fcd80627492d592fbf8606c98a8bb56b62ae4626c842793374a571c184d",
     "outcome": "pass", "evaluated": "2026-10-11",
     "result": "PASS on the Standard-Model lines; passes for every odd T_3R coefficient a of the anomaly-free family, so it is necessary and not what "
               "fixes a = +-1 (GF-4)"},
)

LEGACY = ("T0", "T1", "T2", "T3", "T4", "T5")          # stated before evaluation in reports, not hash-sealed: see puzzle.PREREG


def verify(entry_id=None):
    """Recompute each live seal from the code as it stands; a changed rule shows as False.  A superseded seal is reported as 'superseded' (its code is history)."""
    out = {}
    for e in SEALS:
        if entry_id is None or e["id"] == entry_id:
            out[e["id"]] = "superseded" if e["status"].startswith("superseded") else seal_hash(e) == e["sha256"]
    if entry_id is not None and entry_id not in out:
        raise ValueError(f"no seal {entry_id!r}")
    return out


def evaluate(entry_id, *args, **kwargs):
    """Run a sealed rule on a candidate; refuses if the seal no longer verifies.  Outcomes are recorded in the ledger, never here."""
    e = next((e for e in SEALS if e["id"] == entry_id), None)
    if e is None:
        raise ValueError(f"no seal {entry_id!r}")
    if e["status"].startswith("superseded"):
        raise RuntimeError(f"{entry_id} is {e['status']}; run the seal that supersedes it")
    if verify(entry_id)[entry_id] is not True:
        raise RuntimeError(f"{entry_id} no longer verifies: the rule changed after sealing; re-seal under a new id")
    return {"seal": entry_id, "sealed": e["sealed"], "result": _callable(e["rule"])(*args, **kwargs)}


# ------------------------------------------------------------------------------------------------------------------------------------ trials
def _rows():
    from . import star_quotients as S
    return {r["N"]: r for r in S.census()}


def _curves():
    from . import time_reversal as T
    return T.curve_table()


def _orders():
    from . import star_quotients as S
    return set(S.MONSTER_ELEMENT_ORDERS)


def _unique_among_105(key):
    """The torus levels whose own value of key is unique among the 105 two-prime levels (the null for a self-referential property)."""
    rows = _rows()
    return [N for N in TORI if sum(1 for r in rows.values() if r[key] == rows[N][key]) == 1]


# id, date, source, description, kind ('property' | 'self'), pre-registered seal id or None, predicate on (N, census row, curve row)
TRIALS = (
    ("X01", "2026-10-08", "SQ-3", "twin primes", "property", None, lambda N, r, c: r["twin"]),
    ("X02", "2026-10-08", "SQ-3", "consecutive supersingular primes", "property", None, lambda N, r, c: r["consecutive"]),
    ("X03", "2026-10-08", "SQ-3", "new dimension equals the smaller prime", "property", None, lambda N, r, c: r["new_dim"] == r["p"]),
    ("X04", "2026-10-08", "SQ-3", "new dimension equals the smaller prime and the primes are consecutive", "property", None,
     lambda N, r, c: r["new_dim"] == r["p"] and r["consecutive"]),
    ("X05", "2026-10-08", "SQ-3", "new dimension equals the larger prime (the book's reading)", "property", None, lambda N, r, c: r["new_dim"] == r["q"]),
    ("X06", "2026-10-08", "SQ-3", "Galois-orbit pattern (1, 4, 6)", "self:orbits", None, lambda N, r, c: r["orbits"] == (1, 4, 6)),
    ("X07", "2026-10-08", "SQ-4", "index 168", "self:index", None, lambda N, r, c: r["index"] == 168),
    ("X08", "2026-10-08", "SQ-4", "genus 13", "self:genus", None, lambda N, r, c: r["genus"] == 13),
    ("X09", "2026-10-08", "SQ-4", "no elliptic points (torsion-free)", "property", None, lambda N, r, c: r["torsion_free"]),
    ("X10", "2026-10-08", "SQ-5", "an element order of the Monster", "property", None, lambda N, r, c: N in _orders()),
    ("X11", "2026-10-10", "Holes and Charges note", "N = 7 mod 8: Legendre's exception, h(-4N) = h(-N)", "property", None, lambda N, r, c: N % 8 == 7),
    ("X12", "2026-10-10", "Holes and Charges note", "the star quotient is supersingular at 2 (a_2 even; the y^2 + y = cubic shape)", "property", None,
     lambda N, r, c: c["a2"] is not None and c["a2"] % 2 == 0),
    ("X13", "2026-10-10", "TR-3", "the star quotient has one real circle (rhombic lattice, theta = pi)", "property", None, lambda N, r, c: c["real_circles"] == 1),
    ("X14", "2026-10-10", "TR-3", "the star quotient has no rational torsion", "property", None, lambda N, r, c: c["torsion"] == 1),
    ("X15", "2026-10-10", "TR-3", "the optimal curve already has modular degree 4", "property", None, lambda N, r, c: c["curve"].endswith("1")),
    ("X16", "2026-10-10", "TR-3", "odd level (good reduction of the quotient at 2)", "property", None, lambda N, r, c: N % 2 == 1),
    ("X17", "2026-10-11", "GF-01", "the shortest closed loop of X0(N) is the twelve loop (trace 4, discriminant 12)", "property", None,
     lambda N, r, c: _systole_trace(N) == 4),
)


def _systole_trace(N):
    from .forces import systole_trace
    return systole_trace(N)
DISTINCTION_TRIALS = {"only_twin_torus": "X01", "only_new_dim_smaller_and_consecutive": "X04", "only_orbit_pattern": "X06",
                      "only_consecutive_torus": "X02", "only_index": "X07", "only_torus": None}       # 'only_torus' is the arena itself


def trials():
    rows, curves = _rows(), _curves()
    out = []
    for tid, date, source, desc, kind, seal, pred in TRIALS:
        selected = [N for N in TORI if pred(N, rows[N], curves[N])]
        entry = {"id": tid, "date": date, "source": source, "property": desc, "kind": kind, "sealed": seal, "selects": selected, "m": len(selected),
                 "includes_143": LEVEL in selected, "singles_out_143": selected == [LEVEL]}
        if kind.startswith("self:"):
            entry["levels_unique_in_their_own_version"] = _unique_among_105(kind.split(":", 1)[1])
        out.append(entry)
    return out


def counter():
    t = trials()
    n = len(t)
    n1 = sum(1 for e in t if e["m"] == 1)
    exp_include = sum(e["m"] for e in t) / 13
    return {"arena": list(TORI), "trials": n, "pre_registered": sum(1 for e in t if e["sealed"]),
            "include_143": sum(e["includes_143"] for e in t), "expected_include_143_by_chance": round(exp_include, 2),
            "single_level": n1, "single_out_143": sum(e["singles_out_143"] for e in t), "expected_single_out_143_by_chance": round(n1 / 13, 2),
            "self_referential": [(e["id"], len(e["levels_unique_in_their_own_version"])) for e in t if e["kind"].startswith("self:")],
            "reading": "every property was chosen by looking at 143, so 'includes 143' is near-certain by construction and the chance columns apply "
                       "only to rules sealed blind to the level — none has been run; the record is a look-elsewhere account, not evidence (AG-D5, SQ-8)"}


def report():
    v = verify()
    return {"seals": [dict(e, verifies=v[e["id"]]) for e in SEALS], "legacy": LEGACY, "trials": trials(), "counter": counter()}


def text():
    r = report()
    L = ["PR-01 — sealed pre-registrations", ""]
    for e in r["seals"]:
        L.append(f"  {e['id']} ({e['sealed']}, {e['sub']}): {e['statement']}")
        L.append(f"      rule {e['rule']}  sha256 {e['sha256'][:16]}…  verifies {e['verifies']}  outcome {e['outcome']}: {e['result']}")
    L.append(f"  legacy (stated before evaluation, not hash-sealed): {', '.join(r['legacy'])} — see the puzzle map")
    L.append("")
    L.append("The trials counter — properties tried as a reason for 143, on the 13 torus levels:")
    for e in r["trials"]:
        mark = "143 alone" if e["singles_out_143"] else ("includes 143" if e["includes_143"] else "not 143")
        extra = f"; levels unique in their own version: {len(e['levels_unique_in_their_own_version'])}/13" if e["kind"].startswith("self:") else ""
        L.append(f"  {e['id']} {e['date']} {e['source']:22s} {e['property']:78s} -> {e['m']:2d} level(s): {mark}{extra}")
    c = r["counter"]
    L.append(f"  {c['trials']} tried ({c['pre_registered']} pre-registered): {c['include_143']} include 143 (chance {c['expected_include_143_by_chance']}); "
             f"{c['single_out_143']} of {c['single_level']} single-level properties are 143 (chance {c['expected_single_out_143_by_chance']})")
    L.append("  " + c["reading"])
    L.append("  (X12 was logged with a_2 = 0 at first; the independent review corrected it to a_2 even, the actual supersingular condition, before release;")
    L.append("   X17 was computed in the GF-01 reply before it was logged: tried after looking, like the rest)")
    return "\n".join(L)
