#!/usr/bin/env python3
"""Dependency-free exact certificate replay; exits nonzero on any mismatch."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from dyadic_bounds import certify

ROOT=Path(__file__).resolve().parent


def main():
    stored=json.loads((ROOT/'BOUND_CERTIFICATE.json').read_text())
    models=json.loads((ROOT/'FROZEN_MODELS.json').read_text())
    for name, expected in stored['provenance'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise RuntimeError('Provenance mismatch: '+name)
    actual=certify(models['H'],models['Hr'],8,[F(1,8),F(1,4),F(1,2),F(3,4),F(1)],
                   [4,6,8,10,12],F(1,100))
    expected={k:v for k,v in stored.items() if k!='provenance'}
    if actual!=expected:
        raise RuntimeError('Certificate differs on exact replay')
    print(json.dumps({'exact_certificate_replay':'PASS','rows_verified':len(actual['rows']),
                      'largest_certified_listed_horizon':actual['largest_certified_listed_horizon']}))


if __name__=='__main__':
    main()
