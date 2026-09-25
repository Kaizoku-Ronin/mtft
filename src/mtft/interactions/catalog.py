"""Reference catalog queries and evidence-ledger validation.

Validation checks bookkeeping and recorded numerical comparisons. It does not
certify a physical derivation, inspect evidence documents, or compute amplitudes.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter
from importlib.resources import files
from pathlib import Path

SECTORS = {
    "qcd": "Strong interactions", "ew_fermion": "Electroweak fermions",
    "gauge_self": "Electroweak gauge self-interactions", "higgs_yukawa": "Higgs Yukawas",
    "higgs_gauge": "Higgs–gauge interactions", "higgs_self": "Higgs self-interactions",
    "goldstone": "Goldstone interactions", "ghost": "Ghost interactions",
}
STATUSES = ("unmapped", "symmetry-compatible", "vertex-derived", "amplitude-tested")


def load_catalog(path=None):
    """Load a fresh bundled snapshot or a supplied normalized JSON catalog."""
    source = Path(path) if path else files(__package__).joinpath("_data/sm.json")
    result = json.loads(source.read_text(encoding="utf-8"))
    validate_catalog(result)
    return result


def catalog_digest(catalog):
    """Content identity binds a ledger to definitions, not merely vertex numbers."""
    encoded = json.dumps(catalog, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def validate_catalog(catalog):
    """Reject unresolved references, inconsistent tensor indices and duplicate IDs."""
    if not isinstance(catalog, dict) or catalog.get("schema_version") != 1:
        raise ValueError("unsupported catalog schema_version")
    provenance = catalog.get("provenance", {})
    for key in ("source_url", "revision", "model_path", "gauge", "restrictions"):
        if not isinstance(provenance.get(key), str) or not provenance[key].strip():
            raise ValueError(f"missing provenance.{key}")
    if not re.fullmatch(r"[0-9a-f]{40}", provenance["revision"]):
        raise ValueError("unversioned source revision")
    if provenance.get("perturbative_order") != "tree-level":
        raise ValueError("only tree-level catalogs are supported")
    if not provenance.get("source_files") or not all(
        re.fullmatch(r"[0-9a-f]{64}", h) for h in provenance["source_files"].values()
    ):
        raise ValueError("missing or invalid source file hashes")
    for key in ("particles", "couplings", "lorentz", "parameters", "coupling_orders", "functions"):
        if not isinstance(catalog.get(key), dict) or not catalog[key]:
            raise ValueError(f"missing catalog table {key}")
    for symbol, particle in catalog["particles"].items():
        for key in ("pdg_code", "spin", "color"):
            if type(particle.get(key)) is not int:
                raise ValueError(f"{symbol}: {key} must be an integer")
        if not _finite(particle.get("charge")) or not _text(particle.get("name")):
            raise ValueError(f"{symbol}: invalid charge or particle name")
        for key in ("mass", "width"):
            if particle.get(key) not in catalog["parameters"]:
                raise ValueError(f"{symbol}: unresolved {key} parameter")
    for symbol, coupling in catalog["couplings"].items():
        if not _text(coupling.get("value")) or not isinstance(coupling.get("order"), dict):
            raise ValueError(f"{symbol}: invalid symbolic coupling or order")
        for name, power in coupling["order"].items():
            if name not in catalog["coupling_orders"] or type(power) is not int or power < 0:
                raise ValueError(f"{symbol}: invalid coupling order")
    vertices = catalog.get("vertices")
    if not isinstance(vertices, list) or not vertices:
        raise ValueError("empty vertex catalog")
    ids = set()
    for v in vertices:
        vid = v["id"]
        if vid in ids or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", vid):
            raise ValueError(f"duplicate or invalid vertex ID: {vid}")
        ids.add(vid)
        if v["sector"] not in SECTORS or not 3 <= len(v["particles"]) <= 4:
            raise ValueError(f"{vid}: unsupported sector or valence")
        for p in v["particles"]:
            if p not in catalog["particles"]:
                raise ValueError(f"{vid}: unknown particle {p}")
        spins = [catalog["particles"][p]["spin"] for p in v["particles"]]
        for lorentz in v["lorentz"]:
            if lorentz not in catalog["lorentz"] or catalog["lorentz"][lorentz]["spins"] != spins:
                raise ValueError(f"{vid}: unresolved or inconsistent Lorentz structure {lorentz}")
        if not v["color"] or not v["couplings"]:
            raise ValueError(f"{vid}: missing color or coupling tensors")
        slots = set()
        for c in v["couplings"]:
            slot = (c["color_index"], c["lorentz_index"])
            if (slot in slots or any(type(i) is not int for i in slot)
                    or not 0 <= slot[0] < len(v["color"])
                    or not 0 <= slot[1] < len(v["lorentz"])):
                raise ValueError(f"{vid}: invalid coupling tensor slot {slot}")
            slots.add(slot)
            if c["coupling"] not in catalog["couplings"]:
                raise ValueError(f"{vid}: unknown coupling {c['coupling']}")
    # JSON encoding also rejects non-finite numerical data.
    catalog_digest(catalog)


def select_vertices(catalog, *, sector=None, particle=None, physical_only=False):
    """Filter by sector or an exact UFO symbol/name/PDG code (signed, per particle).

    ``physical_only`` hides ghost/Goldstone vertices. It is a display filter,
    not a gauge transformation or sufficient input for a loop calculation.
    """
    if sector is not None and sector not in SECTORS:
        raise ValueError(f"unknown sector {sector!r}; choose from {', '.join(SECTORS)}")
    particle_ids = set()
    if particle is not None:
        token = str(particle)
        particle_ids = {k for k, p in catalog["particles"].items()
                        if token in (k, p["name"], str(p["pdg_code"]))}
        if not particle_ids:
            raise ValueError(f"unknown particle {token!r}")
    return [v for v in catalog["vertices"]
            if (sector is None or sector == v["sector"])
            and (not physical_only or v["sector"] not in ("ghost", "goldstone"))
            and (particle is None or particle_ids.intersection(v["particles"]))]


def summary(catalog=None):
    catalog = load_catalog() if catalog is None else catalog
    return {
        "model_id": catalog["model_id"], "catalog_sha256": catalog_digest(catalog),
        "vertices": len(catalog["vertices"]), "particles": len(catalog["particles"]),
        "physical_vertices": len(select_vertices(catalog, physical_only=True)),
        "sectors": dict(sorted(Counter(v["sector"] for v in catalog["vertices"]).items())),
        "provenance": catalog["provenance"], "scope_notes": catalog.get("scope_notes", []),
    }


def new_ledger(catalog=None):
    """Make an explicit unmapped row for every reference vertex, without promotions."""
    catalog = load_catalog() if catalog is None else catalog
    return {
        "schema_version": 1, "catalog_sha256": catalog_digest(catalog),
        "model_id": catalog["model_id"],
        "entries": [{"vertex_id": v["id"], "status": "unmapped", "evidence": [],
                     "notes": "", "benchmarks": []} for v in catalog["vertices"]],
    }


def _text(value):
    return isinstance(value, str) and bool(value.strip())


def _finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def validate_ledger(ledger, catalog=None):
    """Return status counts; raise on stale catalog, missing evidence or failed gates.

    Evidence strings are citations/artifact paths, not proof validation. A tested
    amplitude needs a documented process, energy, order, observable, units, input
    conventions and two result artifacts, plus a passing numerical comparison.
    """
    catalog = load_catalog() if catalog is None else catalog
    if (not isinstance(ledger, dict) or ledger.get("schema_version") != 1
            or ledger.get("catalog_sha256") != catalog_digest(catalog)):
        raise ValueError("ledger schema/catalog mismatch; regenerate or explicitly migrate it")
    if ledger.get("model_id") != catalog["model_id"]:
        raise ValueError("ledger model_id mismatch")
    expected, seen = {v["id"] for v in catalog["vertices"]}, set()
    counts = Counter({s: 0 for s in STATUSES})
    for row in ledger.get("entries", []):
        vid, status = row.get("vertex_id"), row.get("status")
        if vid not in expected or vid in seen:
            raise ValueError(f"duplicate or unknown ledger vertex {vid}")
        seen.add(vid)
        if status not in STATUSES:
            raise ValueError(f"{vid}: unknown mapping status {status!r}")
        evidence = row.get("evidence")
        if not isinstance(evidence, list) or any(not _text(e) for e in evidence):
            raise ValueError(f"{vid}: evidence must be a list of nonempty references")
        if status != "unmapped" and not evidence:
            raise ValueError(f"{vid}: promoted status requires evidence references")
        if not isinstance(row.get("notes"), str) or not isinstance(row.get("benchmarks"), list):
            raise TypeError(f"{vid}: notes must be text and benchmarks must be a list")
        if status == "amplitude-tested":
            if not row["benchmarks"]:
                raise ValueError(f"{vid}: amplitude-tested requires at least one benchmark")
            for b in row["benchmarks"]:
                fields = ("process", "perturbative_order", "observable", "unit", "conventions",
                          "reference_artifact", "mtft_artifact", "mtft_revision")
                if any(not _text(b.get(k)) for k in fields):
                    raise ValueError(f"{vid}: incomplete benchmark context/artifacts")
                fields = ("energy_gev", "reference_value", "mtft_value",
                          "absolute_tolerance", "relative_tolerance")
                if any(not _finite(b.get(k)) for k in fields):
                    raise ValueError(f"{vid}: benchmark values must be finite numbers")
                if b["energy_gev"] <= 0 or min(b["absolute_tolerance"], b["relative_tolerance"]) < 0:
                    raise ValueError(f"{vid}: invalid benchmark energy/tolerance")
                error = abs(b["mtft_value"] - b["reference_value"])
                tolerance = b["absolute_tolerance"] + b["relative_tolerance"] * abs(b["reference_value"])
                if error > tolerance:
                    raise ValueError(f"{vid}: recorded benchmark fails its tolerance")
        counts[status] += 1
    if seen != expected:
        raise ValueError(f"ledger is missing {len(expected - seen)} vertices")
    return dict(counts)
