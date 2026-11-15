"""MP18: merge every accepted form (find2, find3 on both horocycle pairs, and direct probes) into one list, deduplicated by
subsector and R (1e-7); each entry keeps the pairs and searches that found it.  Usage: python3 -I mp18_merge.py out.json"""
import glob, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
forms = {}
def add(e, how):
    key = (e["subsector"], round(e["R"], 7))
    if key not in forms:
        forms[key] = dict(e, found_by=[])
    if how not in forms[key]["found_by"]:
        forms[key]["found_by"].append(how)
for f in sorted(glob.glob(os.path.join(HERE, "find[23][AB]*.json"))):
    tag = os.path.basename(f)[:-5]
    for e in json.load(open(f))["found"]:
        add(e, tag)
for f in sorted(glob.glob(os.path.join(HERE, "probe_*.json"))):
    for row in json.load(open(f)):
        sub = row["subsector"]
        for pair in ("A", "B", "C"):
            v = row.get(pair)
            if v and "R" in v and v["agree"] <= 1e-6 and v["hecke"] <= 1e-4:
                par = sub.split()[0]
                e11 = int(sub.split("e11=")[1][:2]); e13 = int(sub.split("e13=")[1][:2])
                new = abs(abs(v["b11_sqrt11"]) - 1) < 1e-4 and abs(abs(v["b13_sqrt13"]) - 1) < 1e-4
                add({"subsector": sub, "parity": par, "eps11": e11, "eps13": e13, "R": v["R"], "lambda": 0.25 + v["R"] ** 2, "new": bool(new),
                     "agree": v["agree"], "hecke": v["hecke"], "b11_sqrt11": v["b11_sqrt11"], "b13_sqrt13": v["b13_sqrt13"]}, "probe " + pair)
out = sorted(forms.values(), key=lambda e: e["R"])
json.dump({"found": out}, open(sys.argv[1], "w"), indent=1)
print(len(out), "forms;", sum(1 for e in out if not e["new"]), "old;", "found by one search only:",
      [(e["subsector"], round(e["R"], 6)) for e in out if len(e["found_by"]) == 1])
