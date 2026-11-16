"""Every library the package imports is declared (v1.1.0, 2026-11-16).

Up to v1.0.0 scipy was imported at module level by four modules and declared only in the extra ``full``, so a plain
``pip install mtft`` gave a package in which those modules could not be imported.  These tests read the declared
dependencies from ``pyproject.toml`` and every import in ``src/mtft`` with ``ast``, so they do not depend on what
happens to be installed: an import must be the standard library, the package itself, a declared dependency, or an
optional library that is named in an extra and imported only inside a function or under a guard.
"""
import ast
import re
import sys
import sysconfig
from importlib.util import find_spec
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "mtft"


def _toml():
    try:
        import tomllib
    except ModuleNotFoundError:                     # Python < 3.11: pytest itself depends on tomli there
        import tomli as tomllib
    return tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))


def _names(requirements) -> set:
    return {re.split(r"[<>=!~\[;\s]", r.strip(), maxsplit=1)[0].lower().replace("-", "_") for r in requirements}


def _is_stdlib(name: str) -> bool:
    if hasattr(sys, "stdlib_module_names"):         # Python >= 3.10
        return name in sys.stdlib_module_names
    spec = find_spec(name)                          # Python 3.9: by where the module lives
    if spec is None:
        return False
    if spec.origin in ("built-in", "frozen"):
        return True
    paths = sysconfig.get_paths()
    origin = str(spec.origin or (spec.submodule_search_locations or [""])[0])
    return "site-packages" not in origin and any(origin.startswith(paths[k]) for k in ("stdlib", "platstdlib"))


_GUARDS = ("ImportError", "ModuleNotFoundError", "Exception")
_SCOPES = (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)


def _guarded(node: ast.Try) -> bool:
    """A try statement whose handlers catch a failed import."""
    for h in node.handlers:
        if h.type is None:
            return True
        types = h.type.elts if isinstance(h.type, ast.Tuple) else [h.type]
        if any(isinstance(t, ast.Name) and t.id in _GUARDS for t in types):
            return True
    return False


def _walk(node, lazy: bool, out: list) -> None:
    """Collect (line, root name, lazy) below node; lazy = inside a function or the body of an import guard."""
    for child in ast.iter_child_nodes(node):
        inner = lazy or isinstance(child, _SCOPES)
        if isinstance(child, ast.Import):
            out.extend((child.lineno, a.name.split(".")[0], lazy) for a in child.names)
        elif isinstance(child, ast.ImportFrom) and child.level == 0 and child.module:
            out.append((child.lineno, child.module.split(".")[0], lazy))
        if isinstance(child, ast.Try) and _guarded(child):
            for stmt in child.body:
                _walk(ast.Module(body=[stmt], type_ignores=[]), True, out)
            for stmt in child.handlers + child.orelse + child.finalbody:
                _walk(ast.Module(body=[stmt], type_ignores=[]), inner, out)
            continue
        _walk(child, inner, out)


def _imports():
    """(file, line, root name, lazy) for every absolute import in src/mtft."""
    found = []
    for path in sorted(SRC.rglob("*.py")):
        out = []
        _walk(ast.parse(path.read_text(encoding="utf-8"), filename=str(path)), False, out)
        found.extend((path, line, name, lazy) for line, name, lazy in out)
    return found


def test_every_import_is_declared():
    project = _toml()["project"]
    core = _names(project["dependencies"])
    extras = set().union(*(_names(v) for v in project.get("optional-dependencies", {}).values())) - core
    undeclared, eager_optional = [], []
    for path, line, name, lazy in _imports():
        where = f"{path.relative_to(ROOT)}:{line} {name}"
        if name == "mtft" or name in core or _is_stdlib(name):
            continue
        if name in extras:
            if not lazy:
                eager_optional.append(where)            # an optional library imported when the module loads
            continue
        undeclared.append(where)
    assert not undeclared, f"imported but declared nowhere: {undeclared}"
    assert not eager_optional, f"optional libraries imported at module level: {eager_optional}"


def test_the_declared_set():
    project = _toml()["project"]
    core = _names(project["dependencies"])
    assert core == {"numpy", "scipy", "sympy", "mpmath"}
    used = {name for _, _, name, _ in _imports()}
    assert core <= used                                 # nothing declared that nothing imports
    floors = {re.split(r"[<>=!~\s]", r, maxsplit=1)[0]: r for r in project["dependencies"]}
    # numpy.trapezoid needs numpy 2.0; scipy 1.13 is the first scipy built for numpy 2; matplotlib 3.9 likewise
    assert floors["numpy"] == "numpy>=2.0" and floors["scipy"] == "scipy>=1.13"
    assert all("matplotlib>=3.9" in v for k, v in project["optional-dependencies"].items() if any("matplotlib" in r for r in v))
    calls = [p for p in SRC.rglob("*.py") if "np.trapezoid(" in p.read_text(encoding="utf-8")]
    assert calls                                        # the reason for the numpy floor is still there


def test_the_guard_sees_what_it_guards(tmp_path, monkeypatch):
    """The import walk finds module-level, in-function and guarded imports, and an undeclared import fails the check."""
    (tmp_path / "pyproject.toml").write_text((ROOT / "pyproject.toml").read_text(encoding="utf-8"), encoding="utf-8")
    pkg = tmp_path / "src" / "mtft"
    pkg.mkdir(parents=True)
    (pkg / "a.py").write_text("import numpy\nimport no_such_lib_xyz\n\ndef f():\n    import matplotlib.pyplot\n"
                              "try:\n    import plotly\nexcept ImportError:\n    plotly = None\n", encoding="utf-8")
    monkeypatch.setattr(sys.modules[__name__], "SRC", pkg)
    monkeypatch.setattr(sys.modules[__name__], "ROOT", tmp_path)
    got = {(name, lazy) for _, _, name, lazy in _imports()}
    assert got == {("numpy", False), ("no_such_lib_xyz", False), ("matplotlib", True), ("plotly", True)}
    try:
        test_every_import_is_declared()
    except AssertionError as e:
        assert "no_such_lib_xyz" in str(e)
    else:
        raise AssertionError("an undeclared import passed the check")
