/* Numerical browser engine for the exact MTFT boundary coefficients.
   Natural inverse temperature beta; all information quantities in bits.
   Contact-pattern bit j and boundary bit j are connected by one cut edge. */
function createBoundaryEngine(data) {
  const pop = Array.from({length:256},(_,n)=>{let k=0;for(;n;n&=n-1)k++;return k;});
  const contactIndex = new Map(data.contacts_in_boundary_order.map((v,i)=>[v,i]));
  const sum = a => a.reduce((s,v)=>s+v,0);
  const H = a => {let h=0;for(const p of a)if(p>0)h-=p*Math.log2(p);return h;};
  const poly = (a,p) => {let s=0;for(let i=0;i<a.length;i++)s+=a[i]*p[i];return s;};
  function evaluate(beta,model='ferro') {
    if(!Number.isFinite(beta)||beta<0||beta>3)throw Error('beta outside [0,3]');
    const dm=data.models[model].defect_mask;
    const power=Array.from({length:60},(_,k)=>Math.exp(-2*beta*k));
    const local=new Float64Array(256), ext=new Float64Array(256), hlocal=new Float64Array(256), best=new Float64Array(256);
    for(let s=0;s<256;s++) {
      local[s]=poly(data.local_counts[s],power);ext[s]=poly(data.exterior_counts[s],power);
      let ek=0;for(let k=0;k<data.local_counts[s].length;k++)ek+=k*data.local_counts[s][k]*power[k];
      hlocal[s]=(Math.log(local[s])+2*beta*ek/local[s])/Math.LN2;
      best[s]=power[data.local_minimum_cut[s]]/local[s];
    }
    const joint=new Float64Array(65536), ps=new Float64Array(256), pb=new Float64Array(256);let Z=0;
    for(let s=0;s<256;s++)for(let b=0;b<256;b++){const w=local[s]*ext[b]*power[pop[s^b^dm]];joint[s*256+b]=w;Z+=w;}
    for(let s=0;s<256;s++)for(let b=0;b<256;b++){const p=joint[s*256+b]/Z;joint[s*256+b]=p;ps[s]+=p;pb[b]+=p;}
    const pi=data.partition.A.map(v=>{
      const out=new Float64Array(256);const j=contactIndex.get(v);
      for(let s=0;s<256;s++)out[s]=j!==undefined?(s>>j)&1:poly(data.strict_interior_up_counts[String(v)][s],power)/local[s];
      return out;
    });
    const floor=sum(Array.from(ps,(p,s)=>p*hlocal[s]));
    return {beta,model,joint,ps,pb,pi,hlocal,best,Hcontact:H(ps),HA:H(ps)+floor,floor,HB:H(pb),logZ:Math.log(Z)+84*beta};
  }
  function sensors(state,mask) {
    const kept=Array.from({length:8},(_,j)=>j).filter(j=>mask&(1<<j)), width=1<<kept.length;
    const group=Array.from({length:256},(_,b)=>kept.reduce((g,j,i)=>g|(((b>>j)&1)<<i),0));
    const matrix=new Float64Array(256*width), po=new Float64Array(width);
    for(let s=0;s<256;s++)for(let b=0;b<256;b++)matrix[s*width+group[b]]+=state.joint[s*256+b];
    for(let s=0;s<256;s++)for(let o=0;o<width;o++)po[o]+=matrix[s*width+o];
    const info=Math.max(0,state.Hcontact+H(po)-H(matrix));
    let spinAccuracy=0,fullAccuracy=0;
    for(let o=0;o<width;o++){
      let best=0;
      for(let s=0;s<256;s++)best=Math.max(best,matrix[s*width+o]*state.best[s]);
      fullAccuracy+=best;
      for(let i=0;i<14;i++){let positive=0;for(let s=0;s<256;s++)positive+=state.pi[i][s]*matrix[s*width+o];spinAccuracy+=Math.max(positive,po[o]-positive)/14;}
    }
    return {information:info,meanEntropy:state.HA-info,spinAccuracy,fullAccuracy,observations:kept.length};
  }
  function observe(state,mask,values) {
    const selected=Array.from({length:256},(_,b)=>b).filter(b=>(b&mask)===(values&mask));
    const posterior=new Float64Array(256);let event=0;
    for(let s=0;s<256;s++)for(const b of selected)posterior[s]+=state.joint[s*256+b];
    event=sum(posterior);for(let s=0;s<256;s++)posterior[s]/=event;
    const meanPlus=state.pi.map(a=>sum(Array.from(a,(p,s)=>p*posterior[s])));
    const contactEntropy=H(posterior), coreEntropy=sum(Array.from(posterior,(p,s)=>p*state.hlocal[s]));
    let bestFull=0;for(let s=0;s<256;s++)bestFull=Math.max(bestFull,posterior[s]*state.best[s]);
    let next=null;
    for(let j=0;j<8;j++)if(!(mask&(1<<j))){
      const plus=new Float64Array(256);let pplus=0;
      for(let s=0;s<256;s++)for(const b of selected)if(b&(1<<j))plus[s]+=state.joint[s*256+b]/event;
      pplus=sum(plus);let conditional=0;
      if(pplus>0)conditional+=pplus*H(Array.from(plus,p=>p/pplus));
      if(pplus<1)conditional+=(1-pplus)*H(Array.from(posterior,(p,s)=>(p-plus[s])/(1-pplus)));
      const gain=Math.max(0,contactEntropy-conditional);
      if(next===null||gain>next.gain+1e-10)next={vertex:data.partition.B[j],bit:j,gain};
    }
    return {event,meanPlus,contactEntropy,coreEntropy,entropy:contactEntropy+coreEntropy,
      accuracy:sum(meanPlus.map(p=>Math.max(p,1-p)))/14,bestFull,next};
  }
  return {evaluate,sensors,observe};
}
if(typeof module!=='undefined'&&module.exports)module.exports={createBoundaryEngine};
