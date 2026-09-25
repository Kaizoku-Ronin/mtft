"""Command-line interface to the SM reference catalog and mapping ledger."""
import argparse
import json
import sys
from pathlib import Path

from .catalog import (
    SECTORS,
    load_catalog,
    new_ledger,
    select_vertices,
    summary,
    validate_ledger,
)
from .render import render_html, vertex_svg
from .ufo import import_ufo


def _json(value):
    return json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n"


def _write(path, content, force=False):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Fail rather than destroy a user's evidence ledger or previous export.
    with destination.open("w" if force else "x", encoding="utf-8") as stream:
        stream.write(content)
    print(destination)


def main(argv=None):
    parser = argparse.ArgumentParser(description="SM interaction atlas: reference vertices and MTFT evidence")
    parser.add_argument("--catalog", help="normalized catalog JSON (default: bundled pinned SM)")
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("summary", help="print counts, provenance and scope as JSON")
    listing = sub.add_parser("list", help="list/filter vertices")
    listing.add_argument("--sector", choices=SECTORS)
    listing.add_argument("--particle", help="exact UFO symbol, name, or signed PDG code")
    listing.add_argument("--physical-only", action="store_true")
    listing.add_argument("--json", action="store_true")
    show = sub.add_parser("show", help="show one vertex and its tensor definitions")
    show.add_argument("vertex")
    for command, description in (
        ("export", "export the full normalized reference JSON"),
        ("render", "write the searchable offline HTML atlas"),
        ("init-ledger", "write an all-unmapped evidence ledger"),
        ("svg", "write an individual vertex diagram"),
        ("import-ufo", "statically read a declarative tree-level SM UFO directory"),
    ):
        cmd = sub.add_parser(command, help=description)
        cmd.add_argument("-o", "--output", required=True)
        cmd.add_argument("--force", action="store_true", help="replace an existing output")
        if command == "render":
            cmd.add_argument("--ledger", help="validated evidence ledger to display")
        if command == "svg":
            cmd.add_argument("vertex")
        if command == "import-ufo":
            cmd.add_argument("directory")
            for option in ("source-url", "revision", "model-path", "gauge", "restrictions"):
                cmd.add_argument("--" + option, required=True)
    validate = sub.add_parser("validate-ledger", help="check coverage, evidence and recorded comparisons")
    validate.add_argument("ledger")
    args = parser.parse_args(argv)
    try:
        if args.command == "import-ufo":
            catalog = import_ufo(args.directory, source_url=args.source_url, revision=args.revision,
                                 model_path=args.model_path, gauge=args.gauge, restrictions=args.restrictions)
            _write(args.output, _json(catalog), args.force)
            return 0
        catalog = load_catalog(args.catalog)
        if args.command in (None, "summary"):
            print(_json(summary(catalog)), end="")
        elif args.command == "list":
            vertices = select_vertices(catalog, sector=args.sector, particle=args.particle,
                                       physical_only=args.physical_only)
            if args.json:
                print(_json(vertices), end="")
            else:
                for v in vertices:
                    names = ", ".join(catalog["particles"][p]["name"] for p in v["particles"])
                    print(f"{v['id']:7s} {v['sector']:14s} {names}")
        elif args.command in ("show", "svg"):
            matches = [v for v in catalog["vertices"] if v["id"] == args.vertex]
            if not matches:
                raise ValueError(f"unknown vertex {args.vertex}")
            v = matches[0]
            if args.command == "svg":
                _write(args.output, vertex_svg(v, catalog), args.force)
            else:
                print(_json({"vertex": v,
                             "particles": {p: catalog["particles"][p] for p in v["particles"]},
                             "lorentz": {p: catalog["lorentz"][p] for p in v["lorentz"]},
                             "couplings": {c["coupling"]: catalog["couplings"][c["coupling"]]
                                           for c in v["couplings"]}}), end="")
        elif args.command == "export":
            _write(args.output, _json(catalog), args.force)
        elif args.command == "init-ledger":
            _write(args.output, _json(new_ledger(catalog)), args.force)
        elif args.command == "validate-ledger":
            counts = validate_ledger(json.loads(Path(args.ledger).read_text(encoding="utf-8")), catalog)
            print(_json({"valid_bookkeeping": True, "statuses": counts,
                         "note": "Evidence references are not a certified physical derivation."}), end="")
        elif args.command == "render":
            ledger = json.loads(Path(args.ledger).read_text(encoding="utf-8")) if args.ledger else None
            _write(args.output, render_html(catalog, ledger), args.force)
    except (ValueError, OSError, KeyError, TypeError, SyntaxError) as error:
        print(f"atlas: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
