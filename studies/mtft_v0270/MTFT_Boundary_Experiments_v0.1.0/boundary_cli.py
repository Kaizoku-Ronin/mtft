"""Query the boundary model and averaged sensor information.

python boundary_cli.py --beta .64 --observe 6=+ 18=- 24=+
python boundary_cli.py --beta .64 --model one_reversed_bond --observe all=+
python boundary_cli.py --beta .64 --best-sensors 2
"""
import argparse
import json
from boundary_model import BoundaryModel


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--beta',type=float,default=.64)
    p.add_argument('--model',choices=['ferro','one_reversed_bond'],default='ferro')
    p.add_argument('--observe',nargs='*',default=[],help='Boundary observations such as 6=+ 18=-, or all=+; omitted vertices are hidden.')
    p.add_argument('--best-sensors',type=int,choices=range(9))
    args=p.parse_args();m=BoundaryModel();mask=values=0
    for item in args.observe:
        try:vertex,value=item.split('=');assert value in ('+','-')
        except (ValueError,AssertionError):p.error('observations use vertex=+ or vertex=-')
        if vertex=='all':mask=255;values=255 if value=='+' else 0
        else:
            try:bit=m.B.index(int(vertex))
            except ValueError:p.error('boundary vertices are '+str(m.B))
            mask|=1<<bit
            if value=='+':values|=1<<bit
            else:values&=255^(1<<bit)
    try:s=m.evaluate(args.beta,args.model)
    except ValueError as exc:p.error(str(exc))
    rows=m.all_sensor_masks(s)
    if args.best_sensors is not None:
        candidates=[r for r in rows if r['observed_count']==args.best_sensors]
        best=max(r['mutual_information'] for r in candidates)
        out={'beta':args.beta,'model':args.model,'sensor_count':args.best_sensors,'maximum_MI_bits':best,
             'masks_within_1e_minus9_bits':[{'mask':r['mask'],'vertices':[v for j,v in enumerate(m.B) if r['mask']&(1<<j)]} for r in candidates if abs(r['mutual_information']-best)<1e-9]}
    else:
        out={'observation':m.observation(args.beta,mask,values,args.model),'averaged_over_all_outcomes_for_this_sensor_set':rows[mask],
             'unobserved_interior_entropy':s['H_A'],'interior_entropy_given_contact_spins':s['H_A_given_contacts']}
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=='__main__':main()
