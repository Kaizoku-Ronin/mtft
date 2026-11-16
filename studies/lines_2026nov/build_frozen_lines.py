"""Build the frozen record of the closed lines of Gamma_0(143) for mtft.lines (2026-11-16).

Input, as the Modular Plane lab left them on 12 November 2026:
  lab/mp12/mp12_numbers.json                      (MP12, test T-G4: the 28 closed lines and their thirty integers)
  lab/mp13/mp13_numbers.json                      (MP13, test T-G9: every crossing of the 24 geodesics)
  ../maass_2026nov/lab/mp14/mp14_numbers.json     (MP14, test T-G10: 14 more geodesics over three new motions)

Output (in src/mtft/lines/_data): lines_143.json and PROVENANCE.txt.

What the file is for.  The package computes the graph, the classes, the windings, the closed lines
and the crossings by its own code; this file holds the lab's values for the gates to compare with.
One part of it is used and not only compared: ``coord``, the integer coordinates of the 29
generators on the four Hecke blocks in the lab's bases, which make the lab's table of periods
reproducible number for number.  The large matrices of the lines (up to 199 digits) are not
copied; the package builds them and the gate compares their classes.

Run from the repository root:  python studies/lines_2026nov/build_frozen_lines.py
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "src" / "mtft" / "lines" / "_data"
SRC = {"mp12": HERE / "lab" / "mp12" / "mp12_numbers.json",
       "mp13": HERE / "lab" / "mp13" / "mp13_numbers.json",
       "mp14": HERE.parent / "maass_2026nov" / "lab" / "mp14" / "mp14_numbers.json"}


def sha1(path):
    return hashlib.sha1(Path(path).read_bytes()).hexdigest()


def main():
    n12 = json.loads(SRC["mp12"].read_text(encoding="utf-8"))["numbers"]
    n13 = json.loads(SRC["mp13"].read_text(encoding="utf-8"))["numbers"]
    n14 = json.loads(SRC["mp14"].read_text(encoding="utf-8"))["numbers"]["A"]
    lines = []
    for l, w in zip(n12["lines"], n12["P13"]["lines"]):
        assert (l["motion"], l["ell"]) == (w[0], w[1])
        lines.append({"motion": l["motion"], "ell": l["ell"], "dart": l["dart"], "psi": l["psi"], "phi": l["phi"], "cls": l["cls"],
                      "blocks": l["blocks"], "m11": l["m11"], "m13": l["m13"], "mirror": l["mirror"], "al": w[2:]})
    out = {
        "meta": {"level": 143, "built": "2026-11-16",
                 "source": "Modular Plane lab, records MP12 (T-G4), MP13 (T-G9) and MP14 (T-G10), 12 November 2026",
                 "note": "the lab's values, for the gates of mtft.lines to compare with; 'coord' is also used, as the bases of the periods"},
        "motions": {"cat": [[2, 1], [1, 1]], "s10": [[84, -47], [143, -80]], "pell": [[12, 1], [143, 12]],
                    "uv3": [[4, 1], [3, 1]], "u2v2": [[5, 2], [2, 1]], "uv5": [[6, 1], [5, 1]]},
        "graph": n12["P1"],
        "generators": n12["P2"]["generators"],
        "face_cycles": n12["P2"]["face_cycles"],
        "phi_on_generators": n12["P3"]["phi_on_generators"],
        "cusp_loops": n12["P3"]["cusp_loops"],
        "chi_generators": n12["P4"]["chi_generators"],
        "psi_motions": n12["P4"]["psi_motions"],
        "letters": n12["P4"]["letters"],
        "coord": n12["coord"],
        "lines": lines,
        "block_signs": n12["P13"]["block_signs"],
        "al_images": n12["post"]["al_images"],
        "found": {"rank28": n12["found"]["rank28"], "rank28_with_psi": n12["found"]["rank28_with_psi"],
                  "rank_by_motion": n12["found"]["rank_by_motion"], "lattice": n12["P10"],
                  "index_cuspidal": n12["P8"]["index_cuspidal"], "index_product_lattice": n12["P8"]["index_product_lattice"],
                  "relations_H1": n12["post"]["relations_H1"], "kept": n12["post"]["kept"]},
        "periods": {"Omega1": n12["P9"]["Omega1"], "Omega2": n12["P9"]["Omega2"], "j_curve": n12["P9"]["j_curve"],
                    "small": n12["small"], "ell_small": n12["P9"]["ell_small"]},
        "crossings_24": {"walks": n13["P1"]["walks"], "geodesics": n13["geodesics"], "I": n13["I"], "iota": n13["iota"],
                         "self": n13["self"], "totals": n13["totals"], "mirror": n13["mirror"]},
        "crossings_38": {"geodesics": n14["geodesics"], "I": n14["I"], "iota": n14["iota"], "pell": n14["pell"], "new": n14["new"]},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "lines_143.json"
    path.write_text(json.dumps(out, separators=(",", ":"), sort_keys=True) + "\n", encoding="utf-8")
    prov = ["mtft.lines frozen data: provenance (SHA-1)", "", f"{sha1(path)}  lines_143.json", "",
            "built by studies/lines_2026nov/build_frozen_lines.py from",
            f"{sha1(SRC['mp12'])}  studies/lines_2026nov/lab/mp12/mp12_numbers.json",
            f"{sha1(SRC['mp13'])}  studies/lines_2026nov/lab/mp13/mp13_numbers.json",
            f"{sha1(SRC['mp14'])}  studies/maass_2026nov/lab/mp14/mp14_numbers.json", ""]
    (OUT / "PROVENANCE.txt").write_text("\n".join(prov), encoding="utf-8")
    print(f"wrote {path} ({path.stat().st_size} bytes): {len(lines)} lines, {len(out['generators'])} generators, "
          f"{len(out['crossings_24']['I'])} and {len(out['crossings_38']['I'])} geodesics")


if __name__ == "__main__":
    main()
