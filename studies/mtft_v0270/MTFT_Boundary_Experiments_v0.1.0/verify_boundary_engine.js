/* Cross-language check of the interactive numerical engine. */
const fs=require('fs');
const path=require('path');
const {createBoundaryEngine}=require('./boundary_engine.js');
const data=JSON.parse(fs.readFileSync(path.join(__dirname,'browser_coefficients.json'),'utf8'));
const fixtures=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));
const engine=createBoundaryEngine(data);let maxError=0,scalars=0;
for(const f of fixtures){
  const state=engine.evaluate(f.beta,f.model);
  const sensors=engine.sensors(state,f.mask),obs=engine.observe(state,f.mask,f.values);
  const computed={...state,...sensors,...obs};
  for(const key of ['HA','floor','HB','logZ','information','meanEntropy','spinAccuracy','fullAccuracy','event','entropy','accuracy','bestFull']){
    const error=Math.abs(computed[key]-f[key]);maxError=Math.max(maxError,error);scalars++;
    if(!Number.isFinite(error)||error>1e-10)throw Error(JSON.stringify({model:f.model,beta:f.beta,mask:f.mask,key,error}));
  }
  for(let j=0;j<14;j++){
    const error=Math.abs(obs.meanPlus[j]-f.meanPlus[j]);maxError=Math.max(maxError,error);scalars++;
    if(!Number.isFinite(error)||error>1e-10)throw Error('Spin posterior mismatch');
  }
  // The next-sensor advice must equal explicit averaging over its two outcomes.
  if(obs.next){
    const bit=1<<obs.next.bit;
    const plus=engine.observe(state,f.mask|bit,f.values|bit);
    const minus=engine.observe(state,f.mask|bit,f.values&~bit);
    const gain=obs.entropy-(plus.event*plus.entropy+minus.event*minus.entropy)/obs.event;
    const error=Math.abs(gain-obs.next.gain);maxError=Math.max(maxError,error);scalars++;
    if(!Number.isFinite(error)||error>1e-10)throw Error('Adaptive advice mismatch');
  }
}
console.log(JSON.stringify({status:'passed',fixtures:fixtures.length,scalar_comparisons:scalars,max_abs_error:maxError,tolerance:1e-10}));
