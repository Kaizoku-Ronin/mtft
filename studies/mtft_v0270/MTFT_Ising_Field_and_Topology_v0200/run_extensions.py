"""Reproduce the field, topology, defect, bilayer, and sampler experiments."""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np

from ising_extensions import (graph, signed_counts, brute_joint, IntegerSampler,
                              signed_minimum, bilayer_counts, joint_observables)
from ising_exact import thermodynamics

ROOT = Path(__file__).resolve().parent


def write(name, data):
    (ROOT/name).write_text(json.dumps(data, indent=2)+'\n')


def run_joint():
    baseline={r['N']:r for r in json.loads((ROOT/'baseline_results.json').read_text())['levels']}
    for N in (6,11,15,35,143):
        cx,n,es=graph(N)
        counts,cert=signed_counts(n,es,joint=True)
        gates={'zero_field_marginal':list(map(sum,counts))==baseline[N]['density_of_states'],
               'binomial_magnetization_marginal':[sum(row[m] for row in counts) for m in range(n+1)]==[math.comb(n,m) for m in range(n+1)],
               'spin_flip_symmetry':all(row==row[::-1] for row in counts),
               'cubic_cut_parity':all(not c or (k-m)%2==0 for k,row in enumerate(counts) for m,c in enumerate(row))}
        if n<=16:
            gates['full_brute_joint']=counts==brute_joint(n,es)
        else:
            alt,_=signed_counts(n,es,joint=True,reverse_ties=True)
            gates['alternate_order_identical']=counts==alt
            for m in range(4):
                direct=[0]*(len(es)+1)
                for up in itertools.combinations(range(n),m):
                    up=set(up);k=sum((u in up)!=(v in up) for u,v in es);direct[k]+=1
                gates[f'low_m_brute_{m}']=direct==[row[m] for row in counts]
        for a,d,fw in ((1,1,(1,1)),(2,1,(1,2)),(4,1,(3,2))):
            sampler=IntegerSampler(n,es,a,d,fw)
            z=sum(c*a**(len(es)-k)*d**k*fw[0]**(n-m)*fw[1]**m for k,row in enumerate(counts) for m,c in enumerate(row))
            gates[f'integer_weight_Z_{a}_{d}_{fw}']=z==sampler.Z_integer
        assert all(gates.values()),gates
        write(f'joint_{N}.json',{'N':N,'spins':n,'edges':es,'E':len(es),'genus':cx.inv.genus,'counts':counts,'certificate':cert,'gates':gates})
        print('joint',N,cert['seconds'],gates,flush=True)


def signed_even_transform(counts, n, E):
    out = [0]*(E+1)
    for k,dk in enumerate(counts):
        if not dk: continue
        prev,cur=0,1
        for j in range(E+1):
            out[j]+=dk*cur
            if j<E:
                num=(E-2*k)*cur-(E-j+1)*prev
                assert num%(j+1)==0
                prev,cur=cur,num//(j+1)
    assert all(a%(1<<n)==0 for a in out)
    return [a//(1<<n) for a in out]


def run_twists():
    from mtft.surface.cycles import tree_cotree
    from mtft.surface.ising import gf2_rref
    results=[]
    for N in (6,11,35,143):
        cx,n,es=graph(N); E=len(es)
        cb=tree_cotree(cx)
        P=(cb.basis_matrix.T%2).astype(np.uint8)
        b1=(cx.boundary_1%2).astype(np.uint8)
        b2=(cx.boundary_2%2).astype(np.uint8)
        rank=lambda a:len(gf2_rref(a)[1])
        gates={'flat_cocycles':not np.any(b1@P.T%2),
               'correct_quotient_rank':rank(np.vstack((b2.T,P)))-rank(b2.T)==2*cx.inv.genus,
               'basis_size':len(P)==2*cx.inv.genus}
        baseline,_=signed_counts(n,es)
        rng=np.random.default_rng(N)
        rows=[]
        for i,eta in enumerate(P):
            counts,cert=signed_counts(n,es,eta)
            spin_change=rng.integers(0,2,n,dtype=np.uint8)
            gauged=eta^(b2@spin_change%2)
            other,_=signed_counts(n,es,gauged)
            assert counts==other
            minimum=next(k for k,d in enumerate(counts) if d)
            assert minimum>0
            temps=[]
            for beta in (.2,.4,.64,1.,1.4):
                z=thermodynamics(counts,n,E,beta)['log_Z']
                z0=thermodynamics(baseline,n,E,beta)['log_Z']
                temps.append({'beta':beta,'dimensionless_twist_cost':z0-z,
                              'free_energy_cost_J1':(z0-z)/beta})
            rows.append({'basis_index':i,'negative_edge_indices':np.flatnonzero(eta).tolist(),
                         'counts':counts,'min_frustrated':minimum,'ground_degeneracy':counts[minimum],
                         'gauge_spin_bits':spin_change.tolist(),'gauge_counts_identical':counts==other,
                         'thermal_costs':temps})
        gates['all_gauge_controls']=all(row['gauge_counts_identical'] for row in rows)
        result={'N':N,'genus':cx.inv.genus,'spins':n,'E':E,'basis_representatives':P.tolist(),
                'basis_warning':'Lexicographic tree/cotree basis; individual basis statistics are not intrinsic.',
                'gates':gates,'baseline_counts':baseline,'basis_twists':rows}
        if len(P)>=2:
            signed=[]
            for mask in range(4):
                eta=((mask&1)*P[0])^(((mask>>1)&1)*P[1])
                counts,_=signed_counts(n,es,eta)
                signed.append(signed_even_transform(counts,n,E))
            classes=[]
            for h in range(4):
                nums=[sum((-1)**((h&mask).bit_count())*signed[mask][j] for mask in range(4)) for j in range(E+1)]
                assert all(v%4==0 and v>=0 for v in nums)
                classes.append([v//4 for v in nums])
            assert all(sum(row)==1<<(E-n-1) for row in classes)
            result['two_cycle_projection']={'basis_indices':[0,1],'class_even_counts':classes,
                 'class_bits':'low bit pairs with basis 0; high bit with basis 1; only two homology coordinates retained'}
            if N==143:
                input_text=f'{n} {E}\n'+'\n'.join(f'{u} {v}' for u,v in es)+'\n'
                input_text+='\n'.join(' '.join(str(int(x)) for x in row) for row in P[:2])+'\n'
                (ROOT/'twist_cycle_input143.txt').write_text(input_text)
        assert all(gates.values())
        results.append(result)
        write('surface_twists.json',{'class':'EXACT integer census; numerical free-energy evaluations','levels':results})
        print('twists',N,'basis',len(P),'supports',Counter(len(x['negative_edge_indices']) for x in rows),
              'minimums',Counter(x['min_frustrated'] for x in rows),flush=True)


def run_defects():
    cx,n,es=graph(143);E=len(es)
    preferred=[1]*E
    base,cert=signed_counts(n,es,preferred)
    ground=signed_minimum(n,es,preferred)
    k0=ground['min_frustrated'];ref=ground['spin_bits']
    assert k0==9 and base[k0]==2
    rows=[]
    for e,(u,v) in enumerate(es):
        b=preferred.copy();b[e]=0
        counts,_=signed_counts(n,es,b)
        mn=next(k for k,d in enumerate(counts) if d)
        witness=signed_minimum(n,es,b,reference=ref)
        assert witness['min_frustrated']==mn
        rows.append({'edge':e,'endpoints':[u,v],'self_loop':u==v,
                     'was_frustrated':ref[u]==ref[v],'min_frustrated':mn,
                     'ground_energy':2*mn-E,'ground_degeneracy':counts[mn],
                     'minimum_rearrangement':witness['minimum_hamming_distance'],
                     'witness_bits':witness['spin_bits'],'counts':counts})
    out={'N':143,'spins':n,'E':E,'baseline_ground':ground,'baseline_degeneracy':base[k0],
         'intervention':'Reverse one AF bond to ferro; keep all magnitudes at one.',
         'distance_definition':'Minimum Hamming distance to a fixed baseline coloring, over all new ground states; global spin-flip symmetry is included.',
         'rows':rows,'gates':{'all_witnesses_match_counts':True,'census_complete':len(rows)==E}}
    # Independent counting certificate for the largest forced rearrangement:
    # gauge the reference to all-zero spins, so m becomes Hamming distance.
    worst=max(rows,key=lambda row:row['minimum_rearrangement'])
    b=preferred.copy();b[worst['edge']]=0
    gauged=[bit^ref[u]^ref[v] for (u,v),bit in zip(es,b)]
    distance_counts,distance_cert=signed_counts(n,es,gauged,joint=True)
    ground_distances=distance_counts[worst['min_frustrated']]
    assert list(map(sum,distance_counts))==worst['counts']
    assert next(m for m,c in enumerate(ground_distances) if c)==worst['minimum_rearrangement']
    assert sum(ground_distances)==worst['ground_degeneracy']
    out['largest_rearrangement_counting_certificate']={
        'edge':worst['edge'],'min_frustrated':worst['min_frustrated'],
        'ground_state_counts_by_hamming_distance':ground_distances,
        'gauge_preferred_bits':gauged,'certificate':distance_cert}
    out['gates']['largest_rearrangement_independently_counted']=True
    write('defect_census.json',out)
    print('defects','minimums',Counter(r['min_frustrated'] for r in rows),'degeneracies',Counter(r['ground_degeneracy'] for r in rows),
          'rearrangements',Counter(r['minimum_rearrangement'] for r in rows),flush=True)


def run_bilayer():
    rows=[]
    for N in (6,11,15,35):
        cx,n,es=graph(N);E=len(es)
        counts,cert=bilayer_counts(n,es)
        single,_=signed_counts(n,es)
        convolution=[sum(single[i]*single[k-i] for i in range(len(single)) if 0<=k-i<len(single)) for k in range(2*E+1)]
        gates={'zero_coupling_factorization':list(map(sum,counts))==convolution,
               'overlap_binomial_marginal':[sum(row[j] for row in counts) for j in range(n+1)]==[(1<<n)*math.comb(n,j) for j in range(n+1)],
               'opposite_coupling_symmetry':all(row==row[::-1] for row in counts),
               'locked_layer_limit':all(row[0]==(single[k//2] if k%2==0 else 0) for k,row in enumerate(counts))}
        if n<=8:
            brute=np.zeros((2*E+1,n+1),dtype=np.int64)
            for cfg in range(1<<(2*n)):
                s=[(cfg>>i)&1 for i in range(n)]
                t=[(cfg>>(n+i))&1 for i in range(n)]
                k=sum((s[u]^s[v])+(t[u]^t[v]) for u,v in es)
                j=sum(a^b for a,b in zip(s,t));brute[k,j]+=1
            gates['full_brute_bilayer']=counts==brute.tolist()
        assert all(gates.values())
        rows.append({'N':N,'spins_per_layer':n,'E_per_layer':E,'counts':counts,'certificate':cert,'gates':gates})
        write('bilayer.json',{'interaction':'p proportional to exp(beta*S_intra + kappa*Q), Q=sum_i sigma_i*tau_i. Two fields on one graph; kappa is an explicit new coupling.','levels':rows})
        print('bilayer',N,cert['seconds'],gates,flush=True)


def run_sampling():
    # Complete exact probability check on the smallest graph, including a field.
    _,small_n,small_es=graph(6)
    small=IntegerSampler(small_n,small_es,3,2,(2,3))
    probs=[]
    for bits in itertools.product((0,1),repeat=small_n):
        p=small.probability(bits)
        assert p==Fraction(small.weight(bits),small.Z_integer)
        probs.append(p)
    assert sum(probs)==1
    joint=json.loads((ROOT/'joint_143.json').read_text())
    n,es,E=joint['spins'],joint['edges'],joint['E']
    metadata=[];archive={}
    for q in (1,2,4,8):
        sam=IntegerSampler(n,es,q,1)
        samples=sam.sample(20000,seed=14300+q)
        cuts=np.zeros(len(samples),dtype=int)
        for u,v in es:cuts+=samples[:,u]^samples[:,v]
        M=2*samples.astype(int).sum(axis=1)-n
        obs=joint_observables(joint['counts'],n,E,sam.beta)
        exactS=-n*obs['interaction_energy_per_spin']
        exactM=n*obs['magnetization_per_spin']
        g=np.array(obs['fisher_metric'])
        observedS=float(np.mean(E-2*cuts));observedM=float(M.mean())
        seS=math.sqrt(g[0,0]/len(samples));seM=math.sqrt(g[1,1]/len(samples))
        zS=abs(observedS-exactS)/seS;zM=abs(observedM-exactM)/seM
        assert max(zS,zM)<6
        archive[f'q{q}_spin_bits']=samples[:256]
        archive[f'q{q}_cut_histogram']=np.bincount(cuts,minlength=E+1)
        archive[f'q{q}_magnetization_histogram']=np.bincount((M+n)//2,minlength=n+1)
        metadata.append({'a':q,'d':1,'beta':sam.beta,'sample_count':len(samples),'seed':14300+q,
                         'partition_integer':str(sam.Z_integer),'S_mean':observedS,'S_exact_expectation':exactS,
                         'M_mean':observedM,'M_exact_expectation':exactM,'S_standardized_error':zS,'M_standardized_error':zM})
    np.savez_compressed(ROOT/'equilibrium_samples.npz',**archive)
    write('sampling.json',{'method':'Exact integer conditional sampling; each configuration is an independent equilibrium draw. Pseudorandom integer generator, fixed recorded seeds; no physical time interpretation.',
         'small_graph_all_configuration_probabilities_exact':True,'rows':metadata})
    print('sampling',[(r['a'],r['S_standardized_error'],r['M_standardized_error']) for r in metadata],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stages',nargs='+',choices=['joint','twists','defects','bilayer','sampling'],default=['joint','twists','defects','bilayer','sampling'])
    args=p.parse_args()
    for name in args.stages:
        globals()['run_'+name]()
