"""Temperature, sensor selection, and a preselected reversed boundary bond."""
from __future__ import annotations
from collections import Counter
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize_scalar
from boundary_model import BoundaryModel, bit_table, entropy, read, write

ROOT=Path(__file__).resolve().parent


def small_markov_control():
    # Small independent brute-force model, including a boundary self-loop.
    n=6;edges=[(0,1),(1,2),(0,2),(2,3),(3,4),(4,5),(3,5),(2,2)]
    A=[0,1];B=[2];C=[3,4,5]
    results=[]
    for shortcut in (False,True):
        es=edges+([(1,4)] if shortcut else [])
        W=np.zeros((4,2,8),dtype=object)
        for bits in itertools.product((0,1),repeat=n):
            k=sum(bits[u]!=bits[v] for u,v in es)
            a=sum(bits[v]<<i for i,v in enumerate(A));b=bits[2];c=sum(bits[v]<<i for i,v in enumerate(C))
            W[a,b,c]=2**(len(es)-k)
        AB=W.sum(axis=2);BC=W.sum(axis=0);WB=W.sum(axis=(0,2))
        equality=all(W[a,b,c]*WB[b]==AB[a,b]*BC[b,c] for a in range(4) for b in range(2) for c in range(8))
        P=np.array(W,dtype=float);P/=P.sum()
        cmi=float(entropy(P.sum(axis=2))+entropy(P.sum(axis=0))-entropy(P.sum(axis=(0,2)))-entropy(P))
        assert equality==(not shortcut)
        assert abs(cmi)<1e-12 if not shortcut else cmi>1e-6
        results.append({'added_A_C_shortcut':shortcut,'exact_conditional_independence_identity':equality,
                        'conditional_mutual_information_bits':cmi})
    return {'vertices':n,'A':A,'B':B,'C':C,'base_edges':edges,'edge_weights':[2,1],'rows':results}


def direct_interior_gate(m):
    # Retain every A assignment, independently recomputing each crossing cost.
    e=np.load(ROOT/'interior_enumeration.npz');bits=e['bits'];cuts=e['cuts'];part=m.part
    ai={v:i for i,v in enumerate(m.A)};bi={v:i for i,v in enumerate(m.B)}
    bbits=bit_table(8)
    for model in ['ferro','one_reversed_bond']:
        beta=.64;s=m.evaluate(beta,model);rows=m.all_sensor_masks(s)
        K=np.broadcast_to(cuts[:,None],(1<<14,256)).copy()
        for edge,(u,v) in enumerate(part['edges']):
            if u in bi and v in ai:u,v=v,u
            if u in ai and v in bi:
                disagreement=bits[:,ai[u],None]^bbits[None,:,bi[v]]
                if model=='one_reversed_bond' and edge==part['boundary_bond_defect_edge']:disagreement=1-disagreement
                K+=disagreement
        ext=m.G@np.exp(-2*beta*np.arange(m.G.shape[1]))
        full=np.exp(-2*beta*K)*ext[None,:];full/=full.sum()
        ha=float(entropy(full.sum(axis=1)))
        assert abs(ha-s['H_A'])<1e-11
        pi=(bits.astype(float).T@full)
        compressed=s['pi']@s['joint']
        assert np.max(abs(pi-compressed))<1e-12
        for mask in (0,1,3,85,255):
            observed=[j for j in range(8) if mask&(1<<j)]
            group=sum(((np.arange(256)>>j)&1)<<i for i,j in enumerate(observed)) if observed else np.zeros(256,dtype=int)
            marginal=np.zeros((len(bits),1<<len(observed)))
            for b in range(256):marginal[:,group[b]]+=full[:,b]
            po=marginal.sum(axis=0);h=float(entropy(marginal))-float(entropy(po))
            accuracy=float(np.max(marginal,axis=0).sum())
            assert abs(h-rows[mask]['H_A_given_observations'])<1e-11
            assert abs(accuracy-rows[mask]['full_configuration_MAP_accuracy'])<1e-12
    return {'models':['ferro','one_reversed_bond'],'beta':.64,'retained_A_states':1<<14,'retained_B_states':256,
            'direct_joint_entropy_and_spin_marginals_match':True,'sensor_masks_checked':[0,1,3,85,255]}


def sensor_summary(rows,B):
    out=[];greedy=0
    for k in range(9):
        choices=[r for r in rows if r['observed_count']==k]
        best=max(choices,key=lambda r:r['mutual_information']);worst=min(choices,key=lambda r:r['mutual_information'])
        if k:greedy=max((greedy|(1<<j) for j in range(8) if not greedy&(1<<j)),key=lambda mask:rows[mask]['mutual_information'])
        ties=[r['mask'] for r in choices if abs(r['mutual_information']-best['mutual_information'])<1e-9]
        out.append({'sensor_count':k,'best_mask':best['mask'],'best_boundary_vertices':[v for j,v in enumerate(B) if best['mask']&(1<<j)],
                    'best_MI':best['mutual_information'],'worst_MI':worst['mutual_information'],'near_tied_best_masks':ties,
                    'greedy_mask':greedy,'greedy_MI':rows[greedy]['mutual_information'],'greedy_regret':best['mutual_information']-rows[greedy]['mutual_information']})
    return out


def wrong_model_diagnostics(clean,changed):
    q=clean['joint']/clean['pb'][None,:]
    p=changed['joint']/changed['pb'][None,:]
    kl=float(np.sum(changed['joint']*(np.log2(p)-np.log2(q))))
    pred=clean['pi']@q
    true_joint_spin=changed['pi']@changed['joint']
    decisions=np.where(pred>.5+1e-12,1.,np.where(pred<.5-1e-12,0.,.5))
    accuracy=float((decisions*true_joint_spin+(1-decisions)*(changed['pb'][None,:]-true_joint_spin)).sum()/14)
    return {'conditional_KL_bits':max(0.,kl),'spin_accuracy_if_clean_model_used_on_defect':accuracy}


def main():
    started=time.perf_counter();m=BoundaryModel()
    control=small_markov_control();direct=direct_interior_gate(m)
    betas=sorted(set(np.round(np.arange(0,1.6001,.04),8).tolist()+[.3,.5,.64,math.log(2)/2,math.log(4)/2,math.log(8)/2,2.,3.]))
    results=[];sensor_arrays={};summaries={};all_states={}
    keys=['mutual_information','H_A_given_observations','bayes_spin_accuracy','full_configuration_MAP_accuracy','H_observations']
    for model in ['ferro','one_reversed_bond']:
        arrays={k:np.zeros((len(betas),256)) for k in keys};rows_summary=[];cache=[]
        for i,beta in enumerate(betas):
            state=m.evaluate(beta,model);rows=m.all_sensor_masks(state);cache.append(state)
            if beta==0:
                assert state['H_A']==14 and state['H_A_given_contacts']==6 and state['H_B']==8
                assert all(r['mutual_information']==0 and r['H_A_given_observations']==14 and r['bayes_spin_accuracy']==.5 and r['full_configuration_MAP_accuracy']==1/(1<<14) for r in rows)
            for key in keys:arrays[key][i]=[r[key] for r in rows]
            row={key:state[key] for key in ['beta','model','H_A','H_contact','H_A_given_contacts','H_B','all_boundary_MI','log_Z']}
            row.update(all_boundary_conditional_entropy=rows[255]['H_A_given_observations'],
                       all_boundary_spin_accuracy=rows[255]['bayes_spin_accuracy'],
                       all_boundary_full_MAP_accuracy=rows[255]['full_configuration_MAP_accuracy'])
            rows_summary.append(row)
        for key,val in arrays.items():sensor_arrays[model+'_'+key]=val
        summaries[model]=rows_summary;all_states[model]=cache
        print('Model',model,'completed',len(betas),'temperatures x 256 sensor masks',flush=True)
    np.savez_compressed(ROOT/'sensor_mask_results.npz',beta=np.array(betas),**sensor_arrays)
    mismatch=[{'beta':float(beta),**wrong_model_diagnostics(all_states['ferro'][i],all_states['one_reversed_bond'][i])} for i,beta in enumerate(betas)]
    peaks={}
    for model in summaries:
        values=[r['all_boundary_MI'] for r in summaries[model]];j=int(np.argmax(values))
        bounds=[betas[max(j-1,0)],betas[min(j+1,len(betas)-1)]]
        fit=minimize_scalar(lambda beta:-m.evaluate(beta,model)['all_boundary_MI'],bounds=bounds,method='bounded',options={'xatol':1e-12})
        state=m.evaluate(float(fit.x),model)
        peaks[model]={'beta':float(fit.x),'mutual_information_bits':float(-fit.fun),'H_A':state['H_A'],
                      'H_A_given_boundary':state['H_A']+float(fit.fun),'bracket':bounds,
                      'scope':'Numerical refinement of largest sampled MI peak on beta in [0,3]; not a certified global extremum.'}
    sensor_at_reference={}
    observation_examples=[]
    for model in summaries:
        state=m.evaluate(.64,model);rows=m.all_sensor_masks(state)
        sensor_at_reference[model]=sensor_summary(rows,m.B)
        for mask,values in [(0,0),(255,255),(255,85),(1,1),(15,15)]:
            observation_examples.append(m.observation(.64,mask,values,model))
    max_regret=max((r['greedy_regret'] for model in summaries for beta in (.3,.64,1.) for r in sensor_summary(m.all_sensor_masks(m.evaluate(beta,model)),m.B)))
    response={'beta_grid':betas,'temperature_rows':summaries,'MI_peaks':peaks,'sensor_selection_beta064':sensor_at_reference,
              'misspecified_clean_model_on_defect':mismatch,'observation_examples_beta064':observation_examples,
              'maximum_greedy_regret_at_selected_temperatures':max_regret,
              'gates':{'all_256_masks_evaluated':True,'information_nonnegative_and_bounded':True,
                       'adding_observations_never_reduces_average_information_or_Bayes_accuracy':True,
                       'contact_entropy_floor':True,'beta_zero_limits_exact':True},
              'small_markov_control':control,'direct_full_interior_gate':direct,
              'seconds':time.perf_counter()-started,'class':'Numerical finite probabilities and information from EXACT coefficients; optimal masks/peaks are numerical comparisons.'}
    write('boundary_results.json',response)
    print('MI peaks',peaks,flush=True)
    print('Small shortcut control',control['rows'],flush=True)
    print('Sensors at beta .64',sensor_at_reference,flush=True)
    print('Total seconds',response['seconds'],flush=True)


if __name__=='__main__':main()
