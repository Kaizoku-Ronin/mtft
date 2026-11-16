"""Build src/mtft/gravity/_data/gravity_map.json from the lab's gravity puzzle map (lab/MTFT_Gravity_Puzzle_Map.md).

The map was written by hand; its pieces, its pre-registered tests and its honest negatives are parsed from its headings and
its two tables and carried in its own words.  Nothing is retyped.

    python studies/gravity_2026nov/build_gravity_map.py
"""
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "lab" / "MTFT_Gravity_Puzzle_Map.md"
OUT = HERE.parents[1] / "src" / "mtft" / "gravity" / "_data"


def cells(line: str) -> list:
    """The cells of a table row.  A cell separator is ' | '; a bare '|' inside a cell (|W|, |P|²) is part of the cell."""
    return [c.strip() for c in line.strip().strip("|").split(" | ")]


def build() -> dict:
    text = SOURCE.read_text(encoding="utf-8")
    lines = text.splitlines()
    pieces = [{"piece": m.group(1), "title": m.group(2), "status": m.group(3)}
              for m in (re.match(r"^## ([IVX]+) — (.+) \((open|partial|answered|closed)\)$", ln) for ln in lines) if m]
    tests, negatives = [], []
    for ln in lines:
        if re.match(r"^\| T-G\d+ ", ln):
            c = cells(ln)
            if len(c) != 4:
                raise ValueError(f"a test row with {len(c)} cells: {ln[:60]}")
            ident, title = c[0].split(" ", 1)
            done = re.search(r"\*\*done (\d{4}-\d{2}-\d{2})\*\*", c[3])
            tests.append({"id": ident, "title": title, "registered": c[1], "criterion": c[2], "outcome": c[3],
                          "status": "done" if done else ("pending" if c[3].startswith("pending") else "other"), "done": done.group(1) if done else None})
        elif re.match(r"^\| N-G\d+ \|", ln):
            c = cells(ln)
            if len(c) != 4:
                raise ValueError(f"a negative row with {len(c)} cells: {ln[:60]}")
            negatives.append({"id": c[0], "pieces": c[1], "statement": c[2], "record": c[3]})
    head = next(ln for ln in lines if "pre-registered tests" in ln and "recorded negatives" in ln)
    return {"source": SOURCE.name, "source_sha1": hashlib.sha1(SOURCE.read_bytes()).hexdigest(), "last_changed": re.search(r"Last changed (\d+ \w+ \d{4})", head).group(1),
            "summary_line": head.split(". Last changed")[0] + ".", "pieces": pieces, "tests": tests, "negatives": negatives}


if __name__ == "__main__":
    data = build()
    OUT.mkdir(parents=True, exist_ok=True)
    body = json.dumps(data, ensure_ascii=False, indent=1) + "\n"
    (OUT / "gravity_map.json").write_text(body, encoding="utf-8")
    (OUT / "PROVENANCE.txt").write_text(
        "gravity_map.json: the pieces, pre-registered tests and honest negatives of the lab's gravity puzzle map, in the map's own words.\n"
        f"Source: studies/gravity_2026nov/lab/{SOURCE.name} (last changed {data['last_changed']}), SHA-1 {data['source_sha1']}.\n"
        "Built by studies/gravity_2026nov/build_gravity_map.py on 2026-11-16.\n"
        f"gravity_map.json SHA-1 {hashlib.sha1(body.encode('utf-8')).hexdigest()}\n", encoding="utf-8")
    print(len(data["pieces"]), "pieces;", len(data["tests"]), "tests;", len(data["negatives"]), "negatives;", data["summary_line"][:90])
    print([t["id"] for t in data["tests"] if t["status"] != "done"], sorted({t["status"] for t in data["tests"]}))
