"""MTFT boundary experiment: exact coefficients and numerical information.

A: 14 interior spins, S: 8 contact spins within A, B: 8 boundary spins,
C: 34 exterior spins. The eight A-B edges form a matching S<->B.
There are no A-C or B-B edges in the declared partition.

D_s[k] counts interior states with contact pattern s and internal cut k.
G_b[k] counts exterior states with boundary b and exterior/BC cut k.
The full Gibbs model reduces to P(s,b) proportional to
  D_s(exp(-2 beta)) G_b(exp(-2 beta)) exp(-2 beta * Hamming(s,b xor defect)).
All coefficients are exact Python integers; information uses numerical logs.
"""
from __future__ import annotations
from collections import Counter
import json
import math
from pathlib import Path
import time
import numpy as np
from ising_exact import min_fill_order
from ising_extensions import eliminate_sum, signed_counts, IntegerSampler

ROOT=Path(__file__).resolve().parent


def read(name):return json.loads((ROOT/name).read_text())


def write(name,obj):
    (ROOT/name).write_text(json.dumps(obj,indent=2)+'\n')


def bit_table(n):
    return ((np.arange(1<<n,dtype=np.uint32)[:,None]>>np.arange(n,dtype=np.uint32))&1).astype(np.uint8)


def entropy(p,axis=None):
    p=np.asarray(p,dtype=float)
    # Zero probabilities contribute zero, with no invented small-probability floor.
    logp=np.zeros_like(p);np.log2(p,out=logp,where=p>0)
    return -np.sum(p*logp,axis=axis)


def classify(part):
    A,B,C=map(set,(part['A'],part['B'],part['C']))
    assert not(A&B or A&C or B&C) and A|B|C==set(range(56))
    groups={key:[] for key in ('AA','AB','AC','BB','BC','CC')}
    label={v:k for k,S in [('A',A),('B',B),('C',C)] for v in S}
    for e,(u,v) in enumerate(part['edges']):
        key=''.join(sorted((label[u],label[v])));groups[key].append((e,u,v))
    assert not groups['AC'] and not groups['BB']
    contacts=[]
    for b in part['B']:
        nbr=[v if u==b else u for _,u,v in groups['AB'] if u==b or v==b]
        assert len(nbr)==1;contacts.append(nbr[0])
    assert len(set(contacts))==8 and len(groups['AB'])==8
    return groups,contacts


def conditioned_exterior_counts(part,groups,bits,reverse=False):
    C=part['C'];ci={v:i for i,v in enumerate(C)};B=part['B'];bi={v:i for i,v in enumerate(B)}
    n=len(C);degree=len(groups['CC'])+len(groups['BC']);digit_bits=n+1;base=1<<digit_bits
    c_edges=[(ci[u],ci[v]) for _,u,v in groups['CC']]
    costs=np.zeros((n,2),dtype=int)
    for _,u,v in groups['BC']:
        if u in bi:u,v=v,u
        b=int(bits[bi[v]]);costs[ci[u]]+=[b,1-b]
    factors=[((i,),np.array([base**int(costs[i,0]),base**int(costs[i,1])],dtype=object)) for i in range(n)]
    for u,v in c_edges:
        if u!=v:factors.append((tuple(sorted((u,v))),np.array([[1,base],[base,1]],dtype=object)))
    packed,cert,_=eliminate_sum(n,c_edges,factors,reverse_ties=reverse)
    counts=[(packed>>(digit_bits*k))&(base-1) for k in range(degree+1)]
    assert packed>>(digit_bits*(degree+1))==0 and sum(counts)==1<<n
    return counts,cert


def construct():
    part=read('partition.json');groups,contacts=classify(part)
    A,B,C=part['A'],part['B'],part['C'];ai={v:i for i,v in enumerate(A)}
    abits=bit_table(len(A));bbits=bit_table(len(B))
    local_cut=np.zeros(len(abits),dtype=np.uint8)
    for _,u,v in groups['AA']:local_cut+=abits[:,ai[u]]^abits[:,ai[v]]
    contact_pattern=np.sum(abits[:,[ai[v] for v in contacts]].astype(np.uint32)*(1<<np.arange(8)),axis=1).astype(int)
    local_degree=len(groups['AA']);inner=sorted(set(A)-set(contacts))
    D=np.bincount(contact_pattern*(local_degree+1)+local_cut,minlength=256*(local_degree+1)).reshape(256,local_degree+1)
    up={}
    for v in inner:
        chosen=abits[:,ai[v]]==1
        up[str(v)]=np.bincount((contact_pattern*(local_degree+1)+local_cut)[chosen],minlength=D.size).reshape(D.shape).tolist()
    assert all(int(row.sum())==1<<len(inner) for row in D)
    assert np.array_equal(D,D[::-1])
    minimum=[]
    for s in range(256):minimum.append(int(np.flatnonzero(D[s])[0]))
    G=[];certs=[];started=time.perf_counter()
    for b in range(256):
        counts,cert=conditioned_exterior_counts(part,groups,bbits[b]);G.append(counts);certs.append(cert)
    exterior_seconds=time.perf_counter()-started
    assert G==G[::-1]
    for b in (0,13,55,127,128,255):
        other,_=conditioned_exterior_counts(part,groups,bbits[b],reverse=True);assert other==G[b]
    # Independently enumerate all interior configurations for each fixed boundary.
    Acounts=[]
    for b in range(256):
        crossing=np.sum(abits[:,[ai[v] for v in contacts]]^bbits[b],axis=1)
        Acounts.append(np.bincount(local_cut+crossing,minlength=local_degree+9).tolist())
    e_def=part['boundary_bond_defect_edge'];u,v=part['edges'][e_def]
    boundary_vertex=u if u in B else v;defect_mask=1<<B.index(boundary_vertex)
    gates={'partition_separates_A_C':True,'matching_AB_interface':True,'all_interior_fibers_have_64_states':True,
           'all_exterior_fibers_have_2pow34_states':all(sum(row)==1<<34 for row in G),
           'exterior_spin_reversal':G==G[::-1],'six_alternate_exterior_orders':True}
    baselines={}
    for model,dm in [('ferro',0),('one_reversed_bond',defect_mask)]:
        total=[0]*85
        for b in range(256):
            for k,x in enumerate(Acounts[b^dm]):
                if x:
                    for j,y in enumerate(G[b]):
                        if y:total[k+j]+=x*y
        preferred=[0]*84
        if dm:preferred[e_def]=1
        independent,_=signed_counts(56,part['edges'],preferred)
        assert total==independent
        gates[model+'_full_DOS_matches_independent_elimination']=True
        assert sum(total)==1<<56
        for a in (1,2,4,8):
            Z=sum(c*a**(84-k) for k,c in enumerate(total))
            other=IntegerSampler(56,part['edges'],a,1,preferred=preferred).Z_integer
            assert Z==other
            gates[f'{model}_integer_partition_a{a}']=True
        baselines[model]={'defect_mask':dm,'density_of_states':total,
                          'minimum_frustration':next(k for k,c in enumerate(total) if c),
                          'ground_degeneracy':next(c for c in total if c)}
    result={'partition':part,'contacts_in_boundary_order':contacts,'strict_interior':inner,
            'edge_counts':{key:len(value) for key,value in groups.items()},
            'local_degree':local_degree,'local_counts':D.tolist(),'strict_interior_up_counts':up,
            'local_minimum_cut':minimum,'exterior_counts':G,
            'exterior_degree':len(G[0])-1,'models':baselines,'gates':gates,
            'exterior_construction':{'seconds':exterior_seconds,'largest_width_upper_bound':max(c['width_upper_bound'] for c in certs),
                                     'largest_spin_table':max(c['peak_spin_table_entries'] for c in certs)},
            'class':'EXACT finite integer coefficients. Entropies and probabilities are evaluated numerically.'}
    write('boundary_coefficients.json',result)
    np.savez_compressed(ROOT/'interior_enumeration.npz',bits=abits,cuts=local_cut,contact_patterns=contact_pattern)
    print('Constructed exact model',result['edge_counts'],result['exterior_construction'],gates,flush=True)
    return result


class BoundaryModel:
    def __init__(self,data=None):
        self.data=read('boundary_coefficients.json') if data is None else data
        self.part=self.data['partition'];self.A=self.part['A'];self.B=self.part['B']
        self.contact=self.data['contacts_in_boundary_order'];self.inner=self.data['strict_interior']
        self.D=np.array(self.data['local_counts'],dtype=float)
        self.G=np.array(self.data['exterior_counts'],dtype=float)
        self.up={int(v):np.array(a,dtype=float) for v,a in self.data['strict_interior_up_counts'].items()}
        self.bits=bit_table(8);self.popcount=np.array([s.bit_count() for s in range(256)])
        self.xor=np.bitwise_xor(np.arange(256)[:,None],np.arange(256)[None,:])

    def evaluate(self,beta,model='ferro'):
        if not math.isfinite(beta) or not 0<=beta<=3:raise ValueError('beta must be finite and in [0,3]')
        dm=self.data['models'][model]['defect_mask']
        power=np.exp(-2*beta*np.arange(self.D.shape[1]))
        local=self.D@power
        mean_k=(self.D@(power*np.arange(len(power))))/local
        h_local=(np.log(local)+2*beta*mean_k)/math.log(2)
        exterior=self.G@np.exp(-2*beta*np.arange(self.G.shape[1]))
        kernel=np.exp(-2*beta*self.popcount[self.xor^dm])
        joint=local[:,None]*exterior[None,:]*kernel
        Z_reduced=float(joint.sum());joint/=Z_reduced
        ps=joint.sum(axis=1);pb=joint.sum(axis=0)
        h_contact=float(entropy(ps));floor=float(ps@h_local);h_A=h_contact+floor
        pi=np.empty((14,256))
        for i,v in enumerate(self.A):
            if v in self.contact:pi[i]=self.bits[:,self.contact.index(v)]
            else:pi[i]=(self.up[v]@power)/local
        best_local=np.exp(-2*beta*np.array(self.data['local_minimum_cut']))/local
        return {'beta':float(beta),'model':model,'joint':joint,'ps':ps,'pb':pb,'pi':pi,'h_local':h_local,
                'best_local':best_local,'H_A':h_A,'H_contact':h_contact,'H_A_given_contacts':floor,
                'H_B':float(entropy(pb)),'log_Z':math.log(Z_reduced)+84*beta,
                'all_boundary_MI':h_contact+float(entropy(pb))-float(entropy(joint))}

    def observation(self,beta,mask=255,values=255,model='ferro'):
        if not 0<=mask<256 or not 0<=values<256:raise ValueError('8-bit masks required')
        state=self.evaluate(beta,model)
        selected=(np.arange(256)&mask)==(values&mask)
        weight=state['joint'][:,selected].sum(axis=1);probability=float(weight.sum());posterior=weight/probability
        h_contacts=float(entropy(posterior));h_hidden=float(posterior@state['h_local'])
        spin_p=state['pi']@posterior
        return {'beta':float(beta),'model':model,'mask':mask,'values':values&mask,'event_probability':probability,
                'posterior_spin_plus':spin_p.tolist(),'H_A_given_this_observation':h_contacts+h_hidden,
                'contact_uncertainty':h_contacts,'strict_interior_uncertainty_given_contacts':h_hidden,
                'bayes_spin_accuracy_given_this_observation':float(np.maximum(spin_p,1-spin_p).mean()),
                'best_full_configuration_probability':float(np.max(posterior*state['best_local']))}

    def all_sensor_masks(self,state):
        rows=[None]*256
        nbits=8
        initial=state['joint'].reshape((256,)+(2,)*nbits)
        initial_spin=(state['pi']@state['joint']).reshape((14,)+(2,)*nbits)
        def recurse(J,spin,bits_remaining,mask,last_removed):
            matrix=J.reshape(256,-1);p_o=matrix.sum(axis=0)
            h_o=float(entropy(p_o));h_so=float(entropy(matrix))
            mi=state['H_contact']+h_o-h_so
            if abs(mi)<1e-13:mi=0.
            spin1=spin.reshape(14,-1)
            accuracy=float(np.maximum(spin1,p_o[None,:]-spin1).sum()/14)
            map_accuracy=float(np.max(matrix*state['best_local'][:,None],axis=0).sum())
            rows[mask]={'mask':mask,'observed_count':mask.bit_count(),'mutual_information':mi,
                        'H_A_given_observations':state['H_A']-mi,'H_observations':h_o,
                        'bayes_spin_accuracy':accuracy,'full_configuration_MAP_accuracy':map_accuracy}
            for bit in range(last_removed+1,8):
                if mask&(1<<bit):
                    axis=1+bits_remaining.index(bit)
                    recurse(J.sum(axis=axis),spin.sum(axis=axis),[b for b in bits_remaining if b!=bit],mask^(1<<bit),bit)
        # C-order reshaping puts highest pattern bit on the first boundary axis.
        recurse(initial,initial_spin,list(range(7,-1,-1)),255,-1)
        assert all(r is not None for r in rows)
        for mask,row in enumerate(rows):
            assert -1e-10<=row['mutual_information']<=min(row['H_observations'],state['H_A'])+1e-10
            assert row['H_A_given_observations']>=state['H_A_given_contacts']-1e-10
            for bit in range(8):
                if not mask&(1<<bit):
                    assert rows[mask|(1<<bit)]['mutual_information']>=row['mutual_information']-1e-10
                    assert rows[mask|(1<<bit)]['bayes_spin_accuracy']>=row['bayes_spin_accuracy']-1e-10
                    assert rows[mask|(1<<bit)]['full_configuration_MAP_accuracy']>=row['full_configuration_MAP_accuracy']-1e-10
        return rows


if __name__=='__main__':construct()
