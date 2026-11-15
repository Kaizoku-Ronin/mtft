"""Build the frozen spectrum of Gamma_0(143) for mtft.maass from the lab's final lists (2026-11-15).

Input: lab/mp19/mp19_numbers.json (final_list_2p0.json: the 41 untwisted forms below R = 2.0) and
lab/mp18/mp18_numbers.json (master_v3.json: the 42 eta-twisted forms below R = 1.8), as the Modular
Plane lab left them on 13 November 2026.

For every form the package's own solver (mtft.maass.hejhal.refine) looks again for the eigenvalue
within 1e-5 of the lab's value, on the records' two horocycles, and the coefficient vector of the
first horocycle is kept (the first NCOEF coefficients).  Nothing is copied from the lab but the
sector and the starting value of R.

Output (in src/mtft/maass/_data): maass_143.json, maass_143_coefficients.npz, PROVENANCE.txt.
Run from the repository root:  python studies/maass_2026nov/build_frozen_spectrum.py   (about five minutes)
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from mtft.maass import hejhal as H

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "src" / "mtft" / "maass" / "_data"
NCOEF = 400
SIGN = {1: "p", -1: "m"}


def sha1(path):
    return hashlib.sha1(Path(path).read_bytes()).hexdigest()


def main():
    t0 = time.time()
    lab19 = json.loads((HERE / "lab" / "mp19" / "mp19_numbers.json").read_text(encoding="utf-8"))["final_list_2p0.json"]["found"]
    lab18 = json.loads((HERE / "lab" / "mp18" / "mp18_numbers.json").read_text(encoding="utf-8"))["master_v3.json"]["found"]
    st = H.setup()
    forms, vectors, count = [], [], {}
    for twist, lab in ((False, lab19), (True, lab18)):
        for e in sorted(lab, key=lambda e: e["R"]):
            par = 0 if e["parity"] == "even" else 1
            eps = (int(e["eps11"]), int(e["eps13"]))
            key = (twist, par, eps)
            count[key] = count.get(key, 0) + 1
            fid = f"{'t' if twist else 'u'}.{'e' if par == 0 else 'o'}.{SIGN[eps[0]]}{SIGN[eps[1]]}.{count[key]}"
            note = ""
            try:
                r = H.refine(e["R"], par, eps, twist, st)
            except ValueError as err:                       # no sign change of c_2: try the next coefficients
                r, note = None, str(err)[:60]
                for index in (2, 3, 4):
                    try:
                        r = H.refine(e["R"], par, eps, twist, st, index=index)
                        note = f"refined on coefficient index {index}"
                        break
                    except ValueError:
                        continue
                if r is None:
                    raise RuntimeError(f"{fid}: no refinement near R = {e['R']}")
            c = r["coefficients"]
            p11 = float(c[10 if not twist else 5] * np.sqrt(11))
            p13 = float(c[12 if not twist else 6] * np.sqrt(13))
            forms.append({"id": fid, "twist": twist, "parity": par, "eps11": eps[0], "eps13": eps[1], "R": r["R"],
                          "lambda": 0.25 + r["R"] ** 2, "R_lab": e["R"], "newform": bool(e["new"]), "agreement": r["agreement"],
                          "hecke": r["hecke"], "c11_sqrt11": p11, "c13_sqrt13": p13, "note": note})
            vectors.append(c[:NCOEF])
            print(f"{fid:10s} R = {r['R']:.13f}  lab {e['R']:.12f}  diff {r['R'] - e['R']:+.1e}  agree {r['agreement']:.1e}  hecke {r['hecke']:.1e}  "
                  f"{note}  [{time.time() - t0:.0f} s]", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / "maass_143_coefficients.npz", ids=np.array([f["id"] for f in forms]), c=np.array(vectors))
    meta = {"level": 143,
            "group": "Gamma_0*(143) = Gamma_0(143) with its Atkin-Lehner matrices; characters (eps11, eps13); eta-twist chi = nu_eta^12",
            "expansion": "f = sum_n c_n sqrt(y) K_iR(2 pi nu_n y) trig(2 pi nu_n x), c_1 = 1, nu_n = n (untwisted) or n - 1/2 (twisted), trig = cos (parity 0) or sin (parity 1)",
            "setup": {"Y0": H.Y0_143, "heights": list(st.heights), "M": st.M, "samples": st.horocycles[0].Q},
            "coefficients_kept": NCOEF,
            "complete_below_R": {"untwisted": 1.98, "twisted": 1.8},
            "searched_below_R": {"untwisted": 2.0, "twisted": 1.8},
            "controls": {"level 1": {"odd": 9.5336952613535575, "even": 13.779751351890738},
                         "level 11": {"even (+)": 2.48359105955, "odd (+)": 2.96820576382, "odd (-)": 2.03309099386},
                         "level 13": {"LMFDB 13.1, odd, Fricke -1": 0.97081541}},
            "records": "Modular Plane lab, MP15 to MP19 (12 and 13 November 2026); this file built 2026-11-15 by studies/maass_2026nov/build_frozen_spectrum.py",
            "forms": forms}
    (OUT / "maass_143.json").write_text(json.dumps(meta, indent=1) + "\n", encoding="utf-8")
    prov = ["mtft.maass frozen data: provenance (SHA-1)", "",
            f"{sha1(OUT / 'maass_143.json')}  maass_143.json",
            f"{sha1(OUT / 'maass_143_coefficients.npz')}  maass_143_coefficients.npz", "",
            "built by studies/maass_2026nov/build_frozen_spectrum.py from",
            f"{sha1(HERE / 'lab' / 'mp19_numbers.json')}  studies/maass_2026nov/lab/mp19/mp19_numbers.json",
            f"{sha1(HERE / 'lab' / 'mp18_numbers.json')}  studies/maass_2026nov/lab/mp18/mp18_numbers.json", ""]
    (OUT / "PROVENANCE.txt").write_text("\n".join(prov), encoding="utf-8")
    worst = max(abs(f["R"] - f["R_lab"]) for f in forms)
    print(f"{len(forms)} forms; largest |R - R_lab| = {worst:.2e}; largest disagreement {max(f['agreement'] for f in forms):.2e}; "
          f"largest Hecke residual {max(f['hecke'] for f in forms):.2e}; {time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
