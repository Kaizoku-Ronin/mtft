"""Export exact integer coefficients for the standalone JavaScript engine."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
data=json.loads((ROOT/'boundary_coefficients.json').read_text())
keys=['partition','contacts_in_boundary_order','strict_interior','local_counts',
      'strict_interior_up_counts','local_minimum_cut','exterior_counts']
out={k:data[k] for k in keys}
out['models']={k:{'defect_mask':v['defect_mask']} for k,v in data['models'].items()}
out['positions']=np.round(np.load(ROOT/'graph143_drawing_positions.npy'),7).tolist()
(ROOT/'browser_coefficients.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
print('Exported browser_coefficients.json')
