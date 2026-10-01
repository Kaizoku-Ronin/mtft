"""Bind Matter-01 to the SM atlas without promoting MTFT vertex claims."""
from __future__ import annotations

import json
from pathlib import Path
from mtft.interactions.catalog import load_catalog, new_ledger, validate_ledger, catalog_digest

ROOT = Path(__file__).resolve().parents[1]


def main():
    catalog = load_catalog()
    ledger = new_ledger(catalog)
    electron_photon = next(v for v in catalog['vertices']
                           if v['particles']==['e__plus__','e__minus__','a'])
    row = next(r for r in ledger['entries'] if r['vertex_id']==electron_photon['id'])
    row['notes'] = ('Matter-01 validates an ordinary supplied-parameter QED reference amplitude. '
                    'No selected-vacuum MTFT electron-photon vertex or coefficient has been derived; '
                    'status remains unmapped. See matter01_evidence_ledger.json.')
    counts = validate_ledger(ledger,catalog)
    assert counts['unmapped']==len(catalog['vertices'])
    report = {
        'milestone':'Matter-01', 'catalog_sha256':catalog_digest(catalog),
        'reference_vertex':electron_photon,
        'reference_kernel':{'status':'DIAGNOSTIC: ordinary tree-level QED validated on declared grids',
                            'evidence':['qed_reference_report.json','qed_oracle_checks.json'],
                            'mass':'INPUT','alpha':'INPUT','photon_massless':'DEFINED reference assumption'},
        'conditional_reduction':{'status':'EXACT conditional identities; sampled numerical evaluation',
                                 'evidence':['u1_reduction_checks.json'],
                                 'assumptions':['constant unbroken photon zero mode','positive kinetic Gram',
                                                'uniform charge within each mode block','supplied parent coupling']},
        'selected_MTFT_vacuum':{'status':'OPEN', 'evidence':['vacuum_audit.json'],
                               'missing':['full cohomology and field dictionary','normalized electron Yukawa',
                                          'photon axion/Stueckelberg kernel','physical scales and gauge matching',
                                          'unwanted-mode masses']},
        'atlas_MTFT_mapping_counts':counts,
        'atlas_promotions':0,
    }
    output=ROOT/'results'; output.mkdir(exist_ok=True)
    for name,value in [('sm_atlas_mapping_ledger.json',ledger),('matter01_evidence_ledger.json',report)]:
        (output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'reference_vertex':electron_photon['id'],'mapping_counts':counts,'promotions':0}))


if __name__=='__main__':main()
