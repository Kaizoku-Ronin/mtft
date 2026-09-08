"""Check saved results, original-coordinate dynamics, and inspector behavior."""
import hashlib
import json
import subprocess
import sys
from fractions import Fraction
import numpy as np
import sympy as sp
from scipy.linalg import expm
from common import ROOT, BLOCKS, DIMS, frozen, read, dump, numeric_stage

def main():
    manifest=read('inputs/manifest.json'); raw=frozen()
    assert hashlib.sha256((ROOT/'inputs/x0143_certified.json').read_bytes()).hexdigest()==manifest['frozen_json_sha256']
    assert hashlib.sha256((ROOT/'EXPERIMENT_PLAN.md').read_bytes()).hexdigest()==manifest['plan_sha256']
    assert all(manifest['frozen_package_gates'].values())
    rows=read('prime_atlas.json'); types=read('splitting_types.json'); fields=read('field_certificate.json')
    sector=read('sector_certificate.json'); summary=read('stability_summary.json')
    assert all(sector['checks'].values())
    assert sum(Fraction(c['density']) for c in types['joint'])==1
    assert sum(c['count'] for c in types['joint'])==1248
    assert len(types['joint'])==55
    for row in rows:
        p,n=row['p'],row['square_n']
        assert n*n<p<(n+1)*(n+1)
        assert row['wieferich_base2']==(p!=2 and pow(2,p-1,p*p)==1)
        assert (row['irregularity_index']==-1)==(p==2)
    assert [r['p'] for r in rows if not r['joint_unramified']]==[5,7,19,103]
    assert all(r['joint_unramified'] for r in rows if r['p'] in [11,13])
    st=numeric_stage(); G=st['G']; R,Ri=st['R'],st['Ri']
    # Reconstruct J from the period coordinates, rather than its serialized J field.
    Q=raw['period_Q']; standard=np.block([[np.zeros((13,13)),-np.eye(13)],[np.eye(13),np.zeros((13,13))]])
    JfromQ=np.linalg.solve(Q,standard@Q)
    if np.linalg.eigvalsh((raw['intersection_cycles']@JfromQ+(raw['intersection_cycles']@JfromQ).T)/2)[0]<0:
        JfromQ=-JfromQ
    jerror=float(np.linalg.norm(JfromQ-raw['J_true']))
    assert jerror<1e-9
    with np.load(ROOT/'stability_arrays.npz') as data:
        transfer=data['transfer']; Vs=data['perturbations']
        assert transfer.shape==(16,5,5,161,4,4)
        assert np.min(transfer)>=0 and np.max(np.abs(transfer.sum(axis=-2)-1))<1e-9
        assert np.max(np.abs(np.linalg.norm(Vs,axis=(-2,-1))-1))<1e-12
        # Compare the saved spectral solution with expm in the original integral cycle frame.
        # The Hodge Hilbert–Schmidt norm there is Tr(G^-1 M^T G M), not Tr(M^T M).
        Porig=[np.array(sp.Matrix(sector['projectors'][k]),float) for k in BLOCKS]
        A0=np.eye(26)+(st['T2']+st['T2'].T)/2/(4*summary['T2_operator_norm'])
        frame_error=0.; metric_flow_error=0.; symplectic_flow_error=0.
        E=raw['intersection_cycles'].astype(float)
        for ci in [2,3,4]:
            Fcycle=Ri@st['J']@(A0+.4*Vs[0,ci])@R
            U=expm(80*Fcycle)
            metric_flow_error=max(metric_flow_error,float(np.linalg.norm(U.T@G@U-G)/np.linalg.norm(G)))
            symplectic_flow_error=max(symplectic_flow_error,float(np.linalg.norm(U.T@E@U-E)/np.linalg.norm(E)))
            for i in range(4):
                for j in range(4):
                    M=Porig[j]@U@Porig[i]
                    weight=float(np.trace(np.linalg.solve(G,M.T@G@M))/DIMS[i])
                    frame_error=max(frame_error,abs(weight-transfer[0,ci,-1,-1,j,i]))
        assert max(frame_error,metric_flow_error,symplectic_flow_error)<1e-9
        # Half the oldspace has character (-,-), so it cannot reach the quartic block
        # under AL-preserving flows. This gives a sharp algebraic 1/2 upper bound.
        old_q4=transfer[:,2:4,:,:,2,1]
        assert old_q4.max()<=.5+1e-9
        old_q4_reciprocity=float(np.max(np.abs(DIMS[1]*transfer[:,2:4,:,:,2,1]-DIMS[2]*transfer[:,2:4,:,:,1,2])))
        assert old_q4_reciprocity<1e-9
        for ci,c in enumerate(summary['class_results']):
            leakage=1-np.diagonal(transfer[:,ci,-1],axis1=-2,axis2=-1)
            assert np.allclose(leakage.max(axis=1).mean(axis=0),c['mean_peak_leakage_by_source'],atol=1e-13,rtol=0)
    cases=[(['--prime','1093'],0),(['--prime','3511'],0),(['--prime','103'],0),
           (['--interval','9'],0),(['--sectors'],0),(['--prime','100'],2),(['--interval','101'],2)]
    for args,code in cases:
        proc=subprocess.run([sys.executable,str(ROOT/'inspect_hecke.py'),*args],capture_output=True,text=True)
        assert proc.returncode==code,(args,proc.stderr)
        if code==0:
            json.loads(proc.stdout)
    dump('validation_results.json',{
        'status':'passed','exact_field_checks':fields['checks'],
        'exact_sector_checks':sector['checks'],'prime_label_and_interval_checks':len(rows),
        'joint_density_sum':'1','all_55_categories_included':True,
        'frozen_input_and_plan_hashes_match':True,'frozen_package_gates_passed':len(manifest['frozen_package_gates']),
        'period_reconstruction_J_residual':jerror,'original_cycle_frame_transfer_error':frame_error,
        'original_cycle_frame_Hodge_conservation_error':metric_flow_error,
        'original_cycle_frame_symplectic_conservation_error':symplectic_flow_error,
        'old_quartic_reciprocity_error':old_q4_reciprocity,
        'AL_oldspace_leakage_half_bound':True,'maximum_simulation_residual':max(summary['residuals'].values()),
        'numerical_tolerance':1e-9,'inspector_cases':len(cases),
        'scope':'Exact rational/finite-field certificates plus double-precision numerical checks. Irregularity labels reuse the preceding independently verified Bernoulli census. No Lean proof or physical calibration.'})
    print(json.dumps(read('validation_results.json'),indent=2))

if __name__=='__main__':
    main()
