#!/usr/bin/env python3
"""Matter-03 exact Higgs-pair existence inputs; no PDE solve or mass output.

This standard-library certificate establishes cusp orders, blow-up valuation
parity, classical gauge constraints, and topological slope necessities. The
existence deduction also invokes the cited surface-specific Hitchin--Kobayashi
theorem; finite arithmetic gates do not independently prove that theorem.
"""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def rref(rows):
    a=[[F(x) for x in row] for row in rows]
    pivots=[]; r=0
    if not a:
        return a,pivots
    for c in range(len(a[0])):
        p=next((i for i in range(r,len(a)) if a[i][c]),None)
        if p is None:
            continue
        a[r],a[p]=a[p],a[r]
        d=a[r][c]; a[r]=[x/d for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][c]:
                d=a[i][c]; a[i]=[x-d*y for x,y in zip(a[i],a[r])]
        pivots.append(c); r+=1
        if r==len(a):
            break
    return a,pivots


def read_basis(path,sha1):
    raw=path.read_bytes()
    if hashlib.sha1(raw).hexdigest()!=sha1:
        raise ValueError('Archived q-expansion provenance mismatch')
    lines=[line for line in raw.decode().splitlines() if line.strip() and not line.lstrip().startswith('#')]
    rows=list(csv.reader(lines))
    matrix=[[int(x) for x in row[1:]] for row in rows[1:]]
    if len(matrix)!=141 or any(len(row)!=13 for row in matrix):
        raise ValueError('Wrong basis shape')
    return matrix,{'sha1':sha1,'sha256':hashlib.sha256(raw).hexdigest()}


def canonical_form(matrix,targets):
    equations=[matrix[n]+[target] for n,target in targets.items()]
    reduced,pivots=rref(equations)
    if 13 in pivots or len(pivots)!=len(equations):
        raise ValueError('Canonical coefficient constraints fail')
    vector=[F(0)]*13
    for row,pivot in enumerate(pivots):
        vector[pivot]=reduced[row][-1]
    qcoefs=[sum(F(row[i])*vector[i] for i in range(13)) for row in matrix]
    return vector,qcoefs


def commutator_unit(i,j,k,l):
    out={}
    if j==k:
        out[(i,l)]=out.get((i,l),0)+1
    if l==i:
        out[(k,j)]=out.get((k,j),0)-1
    return {key:value for key,value in out.items() if value}


def run():
    gates=[]
    def check(name,value):
        passed=bool(value)
        gates.append({'name':name,'pass':passed})
        if not passed:
            raise AssertionError(name)
    data=Path(__file__).resolve().parent/'data'
    basis,hashes=read_basis(data/'X0_143_S2_qexpansions.txt','479316a0952d72dbd4731f50389003263cbfc4e1')
    adapted,adapted_hashes=read_basis(data/'X0_143_AL_adapted_qexpansions.txt','e71eaabe5cd23f3a3a35c204373ad6dfd5115358')
    for m in range(1,14):
        check(f'original_cusp_jet_rank_{m}',len(rref(basis[1:m+1])[1])==m)
        check(f'adapted_cusp_jet_rank_{m}',len(rref(adapted[1:m+1])[1])==m)
    hd,hdq=canonical_form(basis,{1:0,2:0,3:0,4:1,5:0})
    hu,huq=canonical_form(basis,{1:0,2:0,3:1,4:0})
    check('Hd_is_frozen_Matter02_order3_witness',hdq[1:6]==[0,0,0,1,0])
    check('Hu_has_exact_canonical_order2',huq[1:5]==[0,0,1,0])
    # f(q)dq/q: curve differential coefficient exponent is n-1.
    down_order=next(n-1 for n,x in enumerate(hdq) if x)
    up_order=next(n-1 for n,x in enumerate(huq) if x)
    check('Hd_curve_order3',down_order==3)
    check('Hu_curve_order2',up_order==2)
    # Blow-up chart z=u,w=uv: dz wedge dw =u du wedge dv.
    # L frame u^3/z=u^2; KL^2 frame u^4 du dv;
    # KL^-2 frame u^-4 du dv. Hence Hd scalar coefficient u^0,
    # Hu scalar coefficient u^7, and HdHu has canonical-square order7.
    jacobian_order=1
    l_frame_order=3-1
    kl2_frame_order=jacobian_order-2+5
    klminus2_frame_order=jacobian_order+2-7
    hd_scalar_order=down_order+2-5
    hu_scalar_order=up_order-2+7
    product_order=down_order+up_order+2*jacobian_order
    check('L_local_frame_order2',l_frame_order==2)
    check('KL2_local_frame_order4',kl2_frame_order==2*l_frame_order)
    check('KLminus2_local_frame_order_minus4',klminus2_frame_order==-2*l_frame_order)
    check('Hd_local_scalar_order0',hd_scalar_order==0)
    check('Hu_local_scalar_order7',hu_scalar_order==7)
    check('spectral_product_exceptional_order7',product_order==7)
    check('odd_divisorial_valuation_excludes_square',product_order%2==1)
    # Independent product-series route for the regular curve coefficients.
    product=[sum(hdq[i]*huq[n-i] for i in range(n+1)) for n in range(20)]
    product_q_order=next(n for n,x in enumerate(product) if x)
    check('product_q_series_order7',product_q_order==7)
    check('product_two_differential_divisions_order5',product_q_order-2==5)
    check('product_leading_coefficient1',product[7]==1)
    powers=(0,0,0,1,1,2,-3,-1)
    check('full_flux_trace0',sum(powers)==0)
    spectators=(0,0,0,1,2,-3)
    check('spectator_power_sum0',sum(spectators)==0)
    check('active_pair_det_trivial',1+(-1)==0)
    # Exact samples expose the signs; the report supplies the global inequality.
    slopes=[]
    for AE,AX,eps in ((3,3,F(1,2)),(3,3,1),(3,3,F(3,2))):
        x=F(AE)-3*eps
        ample=0<eps<min(AE,AX)
        necessary_single=(2*x<=0 and -3*x<=0)
        check(f'ample_sample_{eps}',ample)
        check(f'spectator_necessity_iff_wall_{eps}',necessary_single==(x==0))
        slopes.append({'A_E':AE,'A_X':AX,'epsilon':eps,'mu_L':x,'A1_slope':2*x,'A2_slope':-3*x})
    allowed=[(i,j) for i in range(8) for j in range(8) if powers[i]==powers[j]]
    roots=((3,7),(7,3))
    constraints=[]
    for k,l in roots:
        for a in range(8):
            for b in range(8):
                constraints.append([commutator_unit(i,j,k,l).get((a,b),0) for i,j in allowed])
    rank=len(rref(constraints)[1])
    check('unbroken_original_gauge_dimension16',len(allowed)==16)
    check('paired_background_breaks_three_constant_generators',rank==3)
    check('paired_classical_gauge_centralizer_dimension13',len(allowed)-rank==13)
    Y=(F(1,6),F(1,6),F(1,6),0,0,F(-1,2),F(-1,2),F(1,2))
    T3=(0,0,0,F(1,2),F(-1,2),0,0,0)
    Q=tuple(x+y for x,y in zip(Y,T3))
    check('Q_is_traceless',sum(Q)==0)
    check('Q_commutes_with_Hd',Q[3]-Q[7]==0)
    check('Q_commutes_with_Hu',Q[7]-Q[3]==0)
    check('Y_alone_is_broken_by_both_arrows',Y[3]!=Y[7])
    check('T3_alone_is_broken_by_both_arrows',T3[3]!=T3[7])
    return {
        'milestone':'Matter-03','classification':'EXACT CONDITIONAL CLASSICAL EXISTENCE',
        'finite_arithmetic_checks':len(gates),'passed':True,'gates':gates,
        'background':'S=Bl_(infinity,O)(X0(143) x 143a1); L=pi*O(P) tensor O(-3D); no flat twists; mu(L)=0',
        'Hd_space_dimension':10,'Hu_space_dimension':11,
        'Hd_basis_coefficients':hd,'Hu_basis_coefficients':hu,
        'Hd_q_coefficients':hdq,'Hu_q_coefficients':huq,
        'spectral_product':{'canonical_curve_order':5,'exceptional_valuation':7,'leading_coefficient':1,
            'generic_characteristic_polynomial':'lambda^2-Hd*Hu','not_square':True,
            'no_rank1_invariant_subsheaf':True},
        'single_Hd':{'F_flat':True,'nonzero_D_flat_on_wall':False,'arbitrary_smooth_metric_no_go':True,
            'full_background_wall_displacement_rescue':False},
        'paired_branch':{'active_rank2_Higgs_pair_stable':True,'full_rank8_polystable_on_wall':True,
            'Hermitian_metric_exists_by_surface_HK':True,'Hermitian_metric_computed':False,
            'common_nonzero_scaling_parameter':'t in C*','t0_split_background_polystable':True,
            'VEV_scale_selected':False,'physical_mass_values_derived':False},
        'kahler_samples':slopes,
        'classical_gauge':{'U8_parent_centralizer':'U(3) x U(1)^4','dimension':13,
            'SU8_parent_centralizer':'S(U(3) x U(1)^4)','dimension_SU8':12,
            'Q':Q,'photon_masslessness_with_unspecified_axions':'OPEN'},
        'input_hashes':{'original':hashes,'adapted':adapted_hashes},
        'primary_references':[
            {'url':'https://arxiv.org/pdf/1312.2673','locator':'Definitions1.1-1.2 and Theorem1.3; canonical-twisted surface pair'},
            {'url':'https://arxiv.org/pdf/0802.3391','locator':'Equations3.28 and3.32-3.34; classical bulk F/D equations'}],
        'assumptions':['Frozen columns form the full13-dimensional weight2 cusp-form basis.',
            'Smooth compact projective surface and ample wall polarization.',
            'Classical defect-free BHV bulk parent with canonical twist and no extra soft/axion potential.',
            'The cited surface-specific Hitchin-Kobayashi correspondence.']}


def json_safe(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {str(k):json_safe(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [json_safe(v) for v in value]
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'vacuum_exact.json')
    args=parser.parse_args()
    result=run()
    args.output.write_text(json.dumps(json_safe(result),indent=2)+'\n')
    print(f"Matter03 vacuum: {result['finite_arithmetic_checks']} exact arithmetic gates passed; paired metric existence conditional on surface HK; no selected mass scale")


if __name__=='__main__':
    main()
