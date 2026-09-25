#!/usr/bin/env python3
"""Rebuild the bundled SM snapshot from a pinned MadGraph commit.

Run from a checkout with mtft installed. Downloaded Python is parsed as data,
never imported. No restrictions card is applied. Review the diff and scope
notes before accepting an updated source revision.
"""
import argparse
import json
import tempfile
import urllib.request
from pathlib import Path

from mtft.interactions.ufo import FILES, import_ufo

PIN = "abfd2c92873b3cfe4baa65c661b580e5cc07d45f"
SOURCE = "https://github.com/mg5amcnlo/mg5amcnlo"
ROOT = Path(__file__).resolve().parents[2]


def build(directory, revision=PIN):
    catalog = import_ufo(
        directory, source_url=SOURCE, revision=revision, model_path="models/sm",
        gauge="Feynman-gauge vertices, including ghosts and Goldstone bosons",
        restrictions="No restrict_*.dat applied; intrinsic stock-model assumptions remain",
    )
    catalog["scope_notes"] = [
        "Complete vertices.py inventory of this pinned stock model, not all possible SM diagrams.",
        "The source fixes u, d, s and neutrino masses to zero. It omits the u, d and s Higgs Yukawa vertices.",
        "CKM entries use the source's Wolfenstein parametrization; parameter defaults are inherited reference inputs, not new measurements or MTFT predictions.",
        "Tree-level only: no counterterms, loop amplitudes, or effective hgg/hgamma-gamma vertices.",
        "Hiding ghosts and Goldstones is a display filter, not a change of gauge.",
        "The coupling powers labelled QED also include weak and Higgs interactions.",
        "No source expressions are numerically evaluated and no MTFT vertex is claimed derived.",
    ]
    catalog["coverage_gaps"] = [
        {"target": "h ubar u", "reason": "u mass and Yukawa absent from stock source"},
        {"target": "h dbar d", "reason": "d mass and Yukawa absent from stock source"},
        {"target": "h sbar s", "reason": "s mass and Yukawa absent from stock source"},
    ]
    catalog["provenance"]["license_file"] = "MADGRAPH_LICENSE.txt"
    catalog["provenance"]["license_text"] = (
        ROOT / "src/mtft/interactions/_data/MADGRAPH_LICENSE.txt"
    ).read_text(encoding="utf-8")
    return catalog


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revision", default=PIN, help="full upstream Git SHA")
    parser.add_argument("--local", type=Path, help="use already downloaded models/sm files")
    parser.add_argument("--check", action="store_true", help="compare without writing")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        directory = args.local or Path(temp)
        if not args.local:
            for filename in FILES:
                url = f"https://raw.githubusercontent.com/mg5amcnlo/mg5amcnlo/{args.revision}/models/sm/{filename}"
                with urllib.request.urlopen(url, timeout=30) as response:
                    (directory / filename).write_bytes(response.read())
        catalog = build(directory, args.revision)
    payload = json.dumps(catalog, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    destination = ROOT / "src/mtft/interactions/_data/sm.json"
    if args.check:
        if destination.read_text(encoding="utf-8") != payload:
            raise SystemExit("Reference differs. Review the source, model scope, license and diff.")
        print("Pinned reference reproduced exactly.")
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(payload, encoding="utf-8")
        print(f"Wrote {len(catalog['vertices'])} vertices to {destination}")


if __name__ == "__main__":
    main()
