"""Generate a source-linked module index without importing research engines."""
import argparse
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def render():
    lines = [
        "# Capability index", "",
        "This index covers every public Python module in `src/mtft`, including the",
        "early arithmetic/phenomenology tools and current geometry/research engines.",
        "Descriptions come from module docstrings; they describe software scope and",
        "do not independently certify mathematical or physical claims. Check the",
        "[SM correction registers](SM/) and each engine's assumptions before use.", "",
        "Names in the last column are a few defined entry points, not an exhaustive",
        "API listing. Follow the source link for signatures, exports and limitations.",
        "Use `help(module)` or `python -m mtft.legend search TERM` for more detail.", "",
        "Regenerate with `python scripts/build_capability_index.py`; use `--check` in",
        "reviews to detect modules missing from this index.", "",
    ]
    groups = {}
    for path in sorted((ROOT / "src/mtft").rglob("*.py")):
        relative = path.relative_to(ROOT / "src/mtft")
        if any(part.startswith("_") for part in relative.parts[:-1]):
            continue
        if path.name.startswith("_") and path.name not in ("__init__.py", "__main__.py"):
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        doc = ast.get_docstring(tree) or "See module source for its public interface."
        description = doc.splitlines()[0].strip().replace("|", "\\|")
        names = [n.name for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))
                 and not n.name.startswith("_")]
        module = "mtft." + str(relative.with_suffix("")).replace("/", ".")
        module = module.removesuffix(".__init__")
        if relative.name == "__init__.py" and len(relative.parts) == 1:
            module = "mtft"
        link = "../" + str(path.relative_to(ROOT))
        group = relative.parts[0] if len(relative.parts) > 1 else "Core modules"
        points = ", ".join(f"`{n}`" for n in names[:4]) or "Package exports / data"
        if len(names) > 4:
            points += f" (+{len(names)-4} more in source)"
        groups.setdefault(group, []).append(f"| [`{module}`]({link}) | {description} | {points} |")
    for group, rows in groups.items():
        lines.extend([f"## {group}", "", "| Module | Scope | Entry points |",
                      "|---|---|---|", *rows, ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    path = ROOT / "docs/CAPABILITIES.md"
    result = render()
    if args.check:
        if not path.exists() or path.read_text(encoding="utf-8") != result:
            raise SystemExit("Capability index is stale; run python scripts/build_capability_index.py")
        print("Capability index matches the source tree.")
    else:
        path.write_text(result, encoding="utf-8")
        print(path)


if __name__ == "__main__":
    main()
