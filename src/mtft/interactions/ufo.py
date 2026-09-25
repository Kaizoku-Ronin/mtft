"""Read a conservative, declarative subset of tree-level UFO without executing it.

Only literal constructor assignments, particle ``anti()`` declarations, metadata,
and imports are accepted. Imports are inspected, never performed. Expressions in
couplings, parameters and Lorentz structures remain strings. Dynamic Python and
counterterm models fail explicitly rather than yielding a partial catalog.
"""
from __future__ import annotations

import ast
import hashlib
import operator
import re
from copy import deepcopy
from pathlib import Path

FILES = {
    "particles.py": ("Particle", "particles"),
    "vertices.py": ("Vertex", "vertices"),
    "couplings.py": ("Coupling", "couplings"),
    "lorentz.py": ("Lorentz", "lorentz"),
    "parameters.py": ("Parameter", "parameters"),
    "coupling_orders.py": ("CouplingOrder", "coupling_orders"),
    "function_library.py": ("Function", "functions"),
}


class UFOError(ValueError):
    """Input requires semantics outside the supported static UFO subset."""


def _literal(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, (str, int, float, bool)):
        return node.value
    if isinstance(node, (ast.List, ast.Tuple)):
        items = [_literal(n) for n in node.elts]
        return tuple(items) if isinstance(node, ast.Tuple) else items
    if isinstance(node, ast.Dict):
        return {_literal(k): _literal(v) for k, v in zip(node.keys, node.values)}
    if (isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)
            and node.value.id in ("P", "C", "L", "Param")):
        return node.attr
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = _literal(node.operand)
        if type(value) in (int, float):
            return -value if isinstance(node.op, ast.USub) else value
    if isinstance(node, ast.BinOp):
        operations = {ast.Add: operator.add, ast.Sub: operator.sub,
                      ast.Mult: operator.mul, ast.Div: operator.truediv}
        left, right = _literal(node.left), _literal(node.right)
        if type(node.op) in operations and type(left) in (int, float) and type(right) in (int, float):
            return operations[type(node.op)](left, right)
    raise UFOError(f"unsupported expression at line {getattr(node, 'lineno', '?')}")


def _imports_only(node):
    if isinstance(node, (ast.Import, ast.ImportFrom, ast.Pass)):
        return True
    if isinstance(node, ast.Try):
        bodies = node.body + node.orelse + node.finalbody
        bodies += [n for h in node.handlers for n in h.body]
        return all(_imports_only(n) for n in bodies)
    return False


def _anti(particle):
    result = deepcopy(particle)
    if particle["name"] == particle["antiname"]:
        raise UFOError("anti() on a self-conjugate particle")
    for a, b in (("name", "antiname"), ("texname", "antitexname")):
        result[a], result[b] = particle[b], particle[a]
    for key in ("pdg_code", "charge", "GhostNumber", "LeptonNumber", "Y"):
        if key in result:
            result[key] = -result[key]
    if result["color"] not in (1, 8):
        result["color"] *= -1
    return result


def _read(path, constructor):
    records = {}
    for node in ast.parse(path.read_text(encoding="utf-8"), filename=str(path)).body:
        if _imports_only(node):
            continue
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
            continue  # module docstrings
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            raise UFOError(f"{path.name}:{node.lineno}: expected a declarative assignment")
        if not isinstance(node.targets[0], ast.Name):
            raise UFOError(f"{path.name}:{node.lineno}: unsupported assignment target")
        symbol = node.targets[0].id
        call = node.value
        if symbol.startswith("__"):
            _literal(call)
            continue
        if symbol in records:
            raise UFOError(f"duplicate symbol {symbol}")
        if isinstance(call, ast.Call) and isinstance(call.func, ast.Name):
            if call.func.id != constructor or call.args or any(k.arg is None for k in call.keywords):
                raise UFOError(f"{path.name}:{node.lineno}: unsupported constructor")
            record = {k.arg: _literal(k.value) for k in call.keywords}
            if constructor == "Particle":
                supported = {"pdg_code", "name", "antiname", "spin", "color", "mass", "width",
                             "texname", "antitexname", "charge", "GhostNumber", "LeptonNumber",
                             "Y", "goldstoneboson"}
                if set(record) - supported:
                    raise UFOError("unsupported particle fields; their anti() semantics need review")
        elif (constructor == "Particle" and isinstance(call, ast.Call)
              and isinstance(call.func, ast.Attribute) and call.func.attr == "anti"
              and isinstance(call.func.value, ast.Name) and not call.args and not call.keywords):
            parent = call.func.value.id
            if parent not in records:
                raise UFOError(f"unknown anti() source {parent}")
            record = _anti(records[parent])
        else:
            raise UFOError(f"{path.name}:{node.lineno}: dynamic UFO declarations are not supported")
        record["source_line"] = node.lineno
        records[symbol] = record
    return records


def _sector(particles):
    if any(p.get("GhostNumber", 0) or p["spin"] == -1 for p in particles):
        return "ghost"
    if any(p.get("goldstoneboson", False) for p in particles):
        return "goldstone"
    pdgs = {abs(p["pdg_code"]) for p in particles}
    if 21 in pdgs:
        return "qcd"
    if 25 in pdgs:
        if pdgs == {25}:
            return "higgs_self"
        return "higgs_yukawa" if any(p["spin"] == 2 for p in particles) else "higgs_gauge"
    return "ew_fermion" if any(p["spin"] == 2 for p in particles) else "gauge_self"


def import_ufo(directory, *, source_url, revision, model_path, gauge, restrictions):
    """Normalize a local, tree-level SM UFO directory; no imported code runs.

    ``revision`` must be a full Git commit SHA. Gauge and restrictions are supplied
    explicitly: they cannot be inferred reliably from the folder's name. This
    parser does not apply parameter cards or claim support for arbitrary UFOs.
    """
    from .catalog import validate_catalog

    root = Path(directory)
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise UFOError("revision must be a full lowercase 40-character Git commit SHA")
    if not source_url.startswith("https://") or not all((model_path, gauge, restrictions)):
        raise UFOError("provide an HTTPS source URL, model path, gauge, and restrictions")
    if list(root.glob("CT_*.py")):
        raise UFOError("counterterm UFOs are not supported; use an explicit tree-level model")
    result = {
        "schema_version": 1,
        "model_id": "sm-" + revision[:12],
        "provenance": {
            "source_url": source_url, "revision": revision, "model_path": model_path,
            "gauge": gauge, "restrictions": restrictions, "perturbative_order": "tree-level",
            "source_files": {},
        },
        "scope_notes": [],
    }
    for filename, (constructor, key) in FILES.items():
        path = root / filename
        if not path.is_file():
            raise UFOError(f"missing required model file: {filename}")
        result["provenance"]["source_files"][filename] = hashlib.sha256(path.read_bytes()).hexdigest()
        result[key] = _read(path, constructor)
    vertices = []
    for symbol, record in result["vertices"].items():
        if symbol != record["name"]:
            raise UFOError(f"vertex symbol/name mismatch: {symbol}")
        record["id"] = symbol
        record["couplings"] = [
            {"color_index": i, "lorentz_index": j, "coupling": c}
            for (i, j), c in record["couplings"].items()
        ]
        try:
            record["sector"] = _sector([result["particles"][p] for p in record["particles"]])
        except KeyError as error:
            raise UFOError(f"{symbol}: unknown particle {error}") from error
        vertices.append(record)
    result["vertices"] = vertices
    validate_catalog(result)
    return result
