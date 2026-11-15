#!/usr/bin/env python3
"""MP16: gather the results of every run into mp16_numbers.json (the file that is stored; the inputs below are not).
Usage: python3 -I mp16_collect.py"""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    p = os.path.join(HERE, name)
    return json.load(open(p)) if os.path.exists(p) else None


out = {"plan": "MP16_H24_Mechanism_Plan.md", "vchecks": load("mp16_vchecks.json"), "arith": load("mp16_arith.json"),
       "level143 identities against MP15's list": load("mp16_level143.json"),
       "residuals against MP15's list": load("mp16_residual_mp15.json"),
       "prime levels 11 and 13": load("mp16_prime.json"), "identities with the four new forms": load("mp16_close.json"),
       "cross-level": load("mp16_cross.json"), "searches at 143": {}, "prime-level searches": {}, "C4": {}}
for f in sorted(glob.glob(os.path.join(HERE, "find_*.json"))):
    out["searches at 143"][os.path.basename(f)] = load(os.path.basename(f))
for f in ("old_11.json", "old_13.json"):
    out["prime-level searches"][f] = load(f)
for f in sorted(glob.glob(os.path.join(HERE, "c4_*.json"))):
    out["C4"][os.path.basename(f)] = load(os.path.basename(f))
json.dump(out, open(os.path.join(HERE, "mp16_numbers.json"), "w"), indent=1)
print(f"mp16_numbers.json: {os.path.getsize(os.path.join(HERE, 'mp16_numbers.json'))} bytes; C4 files {sorted(out['C4'])}")
