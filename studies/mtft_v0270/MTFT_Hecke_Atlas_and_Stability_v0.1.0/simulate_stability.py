"""Fixed, dimensionless conservative perturbations in the Hodge metric."""
import numpy as np
from scipy.linalg import expm
from common import ROOT, BLOCKS, DIMS, CLASSES, dump, numeric_stage

EPSILONS = np.array([0.,.05,.1,.2,.4])
TIMES = np.linspace(0.,80.,161)
SEEDS = np.arange(2026090600,2026090616)

def sym(a):
    return (a+a.T)/2

def comm(a,b):
    return a@b-b@a

def perturbations(stage,seed):
    J,Ps,W11,W13 = (stage[k] for k in ['J','P','W11','W13'])
    x = sym(np.random.default_rng(int(seed)).normal(size=(26,26)))
    x = sym((x-J@x@J)/2)
    offsets = sum(np.trace(p@x)/d*p for p,d in zip(Ps,DIMS))
    within = sum(p@x@p for p in Ps)
    both = W11@W13
    al = (x+W11@x@W11+W13@x@W13+both@x@both)/4
    balanced = Ps[1]@al@Ps[2]+Ps[2]@al@Ps[1]
    return np.array([sym(v)/np.linalg.norm(sym(v),'fro') for v in [offsets,within,al,balanced,x]])

def population_matrix(u, slices):
    result = np.empty((len(u),4,4))
    for i,a in enumerate(slices):
        for j,b in enumerate(slices):
            result[:,j,i] = np.sum(np.abs(u[:,b,a])**2,axis=(1,2))/DIMS[i]
    return result

def main():
    st = numeric_stage()
    J,Ps = st['J'],st['P']
    bases=[]
    for P in Ps:
        w,v=np.linalg.eigh(sym(P)); bases.append(v[:,w>.5])
    B=np.hstack(bases)
    assert np.linalg.norm(B.T@B-np.eye(26))<1e-10
    cuts=np.r_[0,np.cumsum(DIMS)]; slices=[slice(cuts[i],cuts[i+1]) for i in range(4)]
    T=sym(st['T2']); tnorm=np.linalg.norm(T,2)
    A0=np.eye(26)+T/(4*tnorm)
    base_freq=np.linalg.eigvalsh(A0)[::2]
    transfers=np.empty((len(SEEDS),len(CLASSES),len(EPSILONS),len(TIMES),4,4))
    Vs=np.empty((len(SEEDS),len(CLASSES),26,26))
    stats=np.empty((len(SEEDS),len(CLASSES),7))
    frequencies=np.empty((len(SEEDS),len(CLASSES),len(EPSILONS),13))
    residuals={'J_squared':float(np.linalg.norm(J@J+np.eye(26))),
               'J_skew':float(np.linalg.norm(J+J.T)),
               'projector_symmetry':float(max(np.linalg.norm(p-p.T) for p in Ps)),
               'projector_completeness':float(np.linalg.norm(Ps.sum(axis=0)-np.eye(26))),
               'perturbation_J_commutator':0., 'generator_skew_projection':0.,
               'generator_symplectic':0., 'population_column_sum':0.,
               'orthogonal_flow':0., 'flow_imaginary_part':0.,
               'expm_crosscheck':0., 'paired_frequency_split':0.,
               'frequency_shift_bound_excess':0., 'balanced_diagonal_block_norm':0.}
    min_eigenvalue=float('inf')
    for si,seed in enumerate(SEEDS):
        variations=perturbations(st,seed); Vs[si]=variations
        for ci,V in enumerate(variations):
            off=V-sum(p@V@p for p in Ps)
            stats[si,ci]=[np.linalg.norm(off),np.linalg.norm(comm(V,T)),np.linalg.norm(comm(V,st['W11'])),
                          np.linalg.norm(comm(V,st['W13'])),np.linalg.norm(V,2),np.linalg.norm(V,'fro'),np.trace(V)]
            residuals['perturbation_J_commutator']=max(residuals['perturbation_J_commutator'],float(np.linalg.norm(comm(V,J))))
            if CLASSES[ci]=='balanced_old_q4':
                residuals['balanced_diagonal_block_norm']=max(residuals['balanced_diagonal_block_norm'],float(max(np.linalg.norm(p@V@p) for p in Ps)))
            for ei,epsilon in enumerate(EPSILONS):
                A=sym(A0+epsilon*V)
                values=np.linalg.eigvalsh(A); frequencies[si,ci,ei]=values[::2]
                min_eigenvalue=min(min_eigenvalue,float(values[0]))
                residuals['paired_frequency_split']=max(residuals['paired_frequency_split'],float(np.max(np.abs(values[::2]-values[1::2]))))
                excess=np.max(np.abs(values[::2]-base_freq))-epsilon*np.linalg.norm(V,2)
                residuals['frequency_shift_bound_excess']=max(residuals['frequency_shift_bound_excess'],float(excess))
                F=J@A; K=(F-F.T)/2
                residuals['generator_skew_projection']=max(residuals['generator_skew_projection'],float(np.linalg.norm(F-K)))
                residuals['generator_symplectic']=max(residuals['generator_symplectic'],float(np.linalg.norm(F.T@st['E']+st['E']@F)))
                eig,vec=np.linalg.eigh(1j*K)
                q=B.T@vec; phase=np.exp(-1j*TIMES[:,None]*eig[None,:])
                u=np.einsum('ak,tk,bk->tab',q,phase,q.conj(),optimize=True)
                transfers[si,ci,ei]=population_matrix(u,slices)
                residuals['flow_imaginary_part']=max(residuals['flow_imaginary_part'],float(np.max(np.abs(u.imag))))
                ortho=np.einsum('tki,tkj->tij',u.conj(),u,optimize=True)-np.eye(26)
                residuals['orthogonal_flow']=max(residuals['orthogonal_flow'],float(np.max(np.linalg.norm(ortho,axis=(1,2)))))
                residuals['population_column_sum']=max(residuals['population_column_sum'],float(np.max(np.abs(transfers[si,ci,ei].sum(axis=1)-1))))
                if si==0 and ei==len(EPSILONS)-1:
                    for ti in [17,len(TIMES)-1]:
                        other=B.T@expm(TIMES[ti]*F)@B
                        residuals['expm_crosscheck']=max(residuals['expm_crosscheck'],float(np.linalg.norm(u[ti]-other)))
        if (si+1)%4==0:
            print(f'Completed {si+1}/{len(SEEDS)} paired seeds',flush=True)
    assert min_eigenvalue>.35-1e-12
    assert all(v<1e-9 for v in residuals.values())
    assert np.max(np.abs(transfers[:,:2,:,:,:,:]-np.eye(4)))<1e-9
    # Under both AL symmetries, only the old–quartic off-diagonal entries survive.
    forbidden=np.ones((4,4),bool); np.fill_diagonal(forbidden,False); forbidden[1,2]=forbidden[2,1]=False
    forbidden_max=float(transfers[:,2:4,:,:,:, :][...,forbidden].max())
    assert forbidden_max<1e-9
    balanced_off_error=float(np.max(np.abs(stats[:,3,0]-1)))
    assert balanced_off_error<1e-9
    diagonal=np.diagonal(transfers,axis1=-2,axis2=-1)
    leakage=1-diagonal
    aggregate=np.sum(leakage*(DIMS/26),axis=-1)
    summaries=[]
    for ci,key in enumerate(CLASSES):
        peak=leakage[:,ci,-1].max(axis=1)
        aggregate_peak=aggregate[:,ci,-1].max(axis=1)
        summaries.append({'class':key,'off_block_norm_mean':float(stats[:,ci,0].mean()),
                          'commutator_T2_mean':float(stats[:,ci,1].mean()),
                          'commutator_W11_mean':float(stats[:,ci,2].mean()),
                          'commutator_W13_mean':float(stats[:,ci,3].mean()),
                          'mean_peak_leakage_by_source':peak.mean(axis=0).tolist(),
                          'min_peak_leakage_by_source':peak.min(axis=0).tolist(),
                          'max_peak_leakage_by_source':peak.max(axis=0).tolist(),
                          'mean_peak_aggregate_leakage':float(aggregate_peak.mean()),
                          'time_and_seed_mean_transfer':transfers[:,ci,-1].mean(axis=(0,1)).tolist(),
                          'mean_peak_aggregate_leakage_by_epsilon':aggregate[:,ci].max(axis=2).mean(axis=0).tolist(),
                          'largest_frequency_shift':float(np.max(np.abs(frequencies[:,ci,-1]-base_freq)))})
    np.savez_compressed(ROOT/'stability_arrays.npz',transfer=transfers,perturbations=Vs,
                        perturbation_stats=stats,frequencies=frequencies,baseline_frequencies=base_freq,
                        epsilons=EPSILONS,times=TIMES,seeds=SEEDS,block_dimensions=DIMS)
    dump('stability_summary.json',{'block_order':BLOCKS,'classes':CLASSES,'seeds':SEEDS.tolist(),
         'epsilon_values':EPSILONS.tolist(),'time_range':[0,80],'time_samples':len(TIMES),
         'flows':int(len(SEEDS)*len(CLASSES)*len(EPSILONS)),'baseline_min_eigenvalue':float(np.linalg.eigvalsh(A0)[0]),
         'minimum_perturbed_eigenvalue':min_eigenvalue,'T2_operator_norm':float(tnorm),
         'normalization':'Frobenius norm(V) = 1 in Hodge-orthonormal coordinates',
         'frequency_units':'dimensionless, chosen energy model',
         'peak_definition':'For each seed, maximum over the fixed 161-point time grid; then arithmetic mean over seeds. Not the continuous-time supremum.',
         'forbidden_AL_transfer_max':forbidden_max,'balanced_offblock_fraction_error':balanced_off_error,
         'numerical_tolerance':1e-9,'residuals':residuals,'class_results':summaries})
    print('Saved 400 flows; largest numerical residual:',max(residuals.values()),flush=True)
    print('Mean peak aggregate leakage:',[(c['class'],c['mean_peak_aggregate_leakage']) for c in summaries],flush=True)

if __name__=='__main__':
    main()
