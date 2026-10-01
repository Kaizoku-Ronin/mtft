/* Standard-library numerical and interaction checks. This is not browser layout QA. */
'use strict';
const fs=require('node:fs'), vm=require('node:vm'), path=require('node:path');
const htmlPath=process.argv[2]||path.join(__dirname,'MTFT_Matter02_Mode_Lab.html');
const outPath=process.argv[3];
const html=fs.readFileSync(htmlPath,'utf8');
const evidence=JSON.parse(html.match(/<script id="matter02Evidence" type="application\/json">([\s\S]*?)<\/script>/)[1]);
const code=html.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
const nodes=new Map();
class Element {
 constructor(id=''){this._id=id;this.value='';this.checked=false;this.textContent='';this.innerHTML='';this.children=[];this.dataset={};this.events={};this.attrs={};}
 set textContent(v){this._textContent=String(v);}
 get textContent(){return this._textContent;}
 set id(v){this._id=v;if(v)nodes.set(v,this);}
 get id(){return this._id;}
 append(...items){this.children.push(...items);}
 setAttribute(k,v){this.attrs[k]=String(v);}
 addEventListener(k,fn){this.events[k]=fn;}
}
for(const m of html.matchAll(/\bid="([^"]+)"/g))nodes.set(m[1],new Element(m[1]));
for(const m of html.matchAll(/<input\b([^>]+)>/g)){
 const attrs=Object.fromEntries([...m[1].matchAll(/([\w-]+)="([^"]*)"/g)].map(x=>[x[1],x[2]]));
 if(attrs.id){const node=nodes.get(attrs.id);node.value=attrs.value||'';node.checked=/\bchecked\b/.test(m[1]);node.attrs=attrs;}
}
nodes.get('centerLift').value='0';
nodes.get('matter02Evidence').textContent=JSON.stringify(evidence);
const root=new Element();
const sandbox={document:{documentElement:root,getElementById:id=>nodes.get(id),createElement:()=>new Element()},window:{},console};
vm.createContext(sandbox);vm.runInContext(code,sandbox,{filename:htmlPath});
const api=sandbox.window.Matter02,gates=[];
function gate(name,passed,detail){gates.push({name,passed:!!passed,detail});if(!passed)throw new Error(name+': '+JSON.stringify(detail));}
function eq(a,b){return JSON.stringify(a)===JSON.stringify(b);}
function matrixRank(input){const m=input.map(r=>r.slice()),width=m[0]?.length||0;let rank=0;for(let c=0;c<width&&rank<m.length;c++){let pivot=rank;while(pivot<m.length&&Math.abs(m[pivot][c])<1e-12)pivot++;if(pivot===m.length)continue;[m[pivot],m[rank]]=[m[rank],m[pivot]];const p=m[rank][c];for(let k=c;k<width;k++)m[rank][k]/=p;for(let j=0;j<m.length;j++)if(j!==rank){const f=m[j][c];for(let k=c;k<width;k++)m[j][k]-=f*m[rank][k];}rank++;}return rank;}
gate('offline document has no external script/style/resource',!/<(?:script|link|img|iframe)\b[^>]*(?:src|href)=["']https?:/i.test(html),'Primary reference anchors are allowed.');
gate('no network API in runtime',!(/\b(?:fetch|XMLHttpRequest|WebSocket)\s*\(/.test(code)),true);
gate('specified representative is visible',html.includes('cusp ∞')&&html.includes('no extra flat twists'),true);
gate('physical photon status stays open',evidence.photonStatus==='OPEN',evidence.photonStatus);
gate('full matrix and Higgs scope visible',eq(evidence.matrixShape,[18,44])&&evidence.higgsSpaceDimension===10&&html.includes('physical masses remain open'),evidence.matrixShape);
for(const row of evidence.certifiedDimensions){
 const c=api.cohomology(row[0]);
 gate('cohomology matches frozen certificate q='+row[0],eq([c.q,c.h0,c.h1,c.h2],row),c);
 gate('Riemann–Roch q='+row[0],c.h0-c.h1+c.h2===-(9*row[0]*row[0]+3*row[0])/2,c.chi);
 const opposite=api.cohomology(-row[0]);
 gate('Serre duality dimension q='+row[0],c.phi===opposite.h2,c.phi);
 gate('mode count and repo index sign q='+row[0],c.h1+opposite.h2-opposite.h1-c.h2===-c.net,c.net);
 nodes.get('power').value=String(row[0]);api.update();
 gate('slider updates q='+row[0],nodes.get('h1Value').textContent===String(c.h1)&&api.state.power===row[0],nodes.get('h1Value').textContent);
}
gate('exact table renders eleven rows',nodes.get('dimensionRows').children.length===11,nodes.get('dimensionRows').children.length);
let survivors=0,witnesses=0,global=0;
for(let a=0;a<3;a++)for(let b=0;b<3;b++)if(a+b<3){const s=api.jetStatus(a,b);if(s==='global')global++;else survivors++;if(s==='witness')witnesses++;}
gate('fat quotient has five local surviving classes',survivors===5&&global===1&&witnesses===3,{survivors,global,witnesses});
gate('all degree-three monomials vanish',[[3,0],[2,1],[1,2],[0,3]].every(([a,b])=>api.jetStatus(a,b)==='vanish'),true);
for(const button of nodes.get('jetGrid').children)if(button.dataset.a!==undefined){button.events.click();gate('monomial selection '+button.attrs['aria-label'],api.state.selectedJet[0]===Number(button.dataset.a)&&api.state.selectedJet[1]===Number(button.dataset.b)&&button.attrs['aria-pressed']==='true',api.state.selectedJet);}
// Coefficient matching is independent of the runtime residue lookup.
for(const coefficient of [-2,-1.25,0,.25,2]){
 const terms=[[5,0,1],[6,0,coefficient]],m=api.residueMatrix(coefficient),reference=evidence.positiveClasses.map(p=>evidence.negativeClasses.map(n=>terms.filter(t=>p[0]-n[0]+t[0]===-1&&p[1]-n[1]+t[1]===-1).reduce((s,t)=>s+t[2],0)));
 gate('two-variable coefficient extraction v6='+coefficient,eq(m,reference),m);
 gate('independent elimination gives rank three v6='+coefficient,matrixRank(m)===3,m);
 gate('determinant fixed v6='+coefficient,api.det3(m)===-1,api.det3(m));
 nodes.get('nextCoefficient').value=String(coefficient);api.update();
 gate('coefficient control updates v6='+coefficient,api.state.nextCoefficient===coefficient&&nodes.get('determinantValue').textContent==='det = −1',nodes.get('determinantValue').textContent);
}
for(let i=0;i<3;i++)for(let j=0;j<3;j++){const node=nodes.get('residue_'+i+'_'+j);node.events.click();gate('residue entry selection '+i+','+j,eq(api.state.selectedEntry,[i,j])&&node.attrs['aria-pressed']==='true',nodes.get('residueDetail').textContent);}
const ranks=[3,2,1,1,1],powers=[0,1,2,-3,-1],Y=[1/6,0,-.5,-.5,.5];
for(const s of [0,.375]){
 const linear=ranks.map((r,i)=>r*powers[i]),quadratic=ranks.map((r,i)=>r*powers[i]**2),trace=ranks;
 for(let mask=0;mask<8;mask++){
  const flags={linear:!!(mask&1),quadratic:!!(mask&2),trace:!!(mask&4)},rows=[];
  if(flags.linear)rows.push(linear);if(flags.quadratic)rows.push(quadratic);if(flags.trace)rows.push(trace);
  const exactDot=rows.map(row=>row.reduce((sum,r,i)=>sum+r*(Y[i]+s),0)),r=api.probe(s,flags);
  gate('probe kernel s='+s+', mask='+mask,r.survives===exactDot.every(d=>Math.abs(d)<1e-12)&&r.rank===matrixRank(rows),{r,exactDot});
  nodes.get('centerLift').value=String(s);for(const id of ['linear','quadratic','trace'])nodes.get(id+'Row').checked=flags[id];api.update();
  gate('probe controls s='+s+', mask='+mask,nodes.get('liftResult').textContent.includes(r.survives?'survives':'removed')&&nodes.get('liftResult').textContent.includes('OPEN'),nodes.get('liftResult').textContent);
 }
}
gate('traceless quadratic projection is 54',api.projection(0).quadratic===54,api.projection(0));
gate('center shift cancels quadratic but not trace',api.projection(.375).quadratic===0&&api.projection(.375).trace===3,api.projection(.375));
nodes.get('themeButton').events.click();gate('light theme control',root.dataset.theme==='light'&&nodes.get('themeButton').attrs['aria-pressed']==='true',root.dataset.theme);
nodes.get('themeButton').events.click();gate('dark theme control',root.dataset.theme==='dark',root.dataset.theme);
nodes.get('power').value='1';nodes.get('nextCoefficient').value='0';nodes.get('centerLift').value='0';nodes.get('linearRow').checked=true;nodes.get('quadraticRow').checked=false;nodes.get('traceRow').checked=false;api.state.selectedJet=[0,1];api.state.selectedEntry=[0,0];api.update();
const report={scope:'Numerical and mock-DOM interaction checks. No actual browser layout test was performed.',passed:gates.filter(g=>g.passed).length,failed:gates.filter(g=>!g.passed).length,gates};
if(outPath)fs.writeFileSync(outPath,JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({passed:report.passed,failed:report.failed,scope:report.scope},null,2));
