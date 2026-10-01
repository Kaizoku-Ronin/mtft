/* Offline numerical and mock-DOM checks; this is not browser layout QA. */
'use strict';
const fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const htmlPath=process.argv[2]||path.join(__dirname,'MTFT_Matter03_Higgs_Lab.html');
const outPath=process.argv[3];
const html=fs.readFileSync(htmlPath,'utf8');
const evidence=JSON.parse(html.match(/<script id="matter03Evidence" type="application\/json">([\s\S]*?)<\/script>/)[1]);
const code=html.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
const nodes=new Map(),paint=[];
const context=new Proxy({},{get:(_,k)=>(...args)=>paint.push([k,...args]),set:(obj,k,v)=>{obj[k]=v;return true;}});
class Element{
 constructor(id=''){this._id=id;this.value='';this.children=[];this.dataset={};this.events={};this.attrs={};this.width=570;this.height=468;}
 set id(v){this._id=v;if(v)nodes.set(v,this);}get id(){return this._id;}
 set textContent(v){this._textContent=String(v);}get textContent(){return this._textContent||'';}
 set innerHTML(v){this._innerHTML=String(v);}get innerHTML(){return this._innerHTML||'';}
 append(...items){this.children.push(...items);}setAttribute(k,v){this.attrs[k]=String(v);}
 addEventListener(k,fn){this.events[k]=fn;}getContext(){return context;}
}
for(const m of html.matchAll(/\bid="([^"]+)"/g))nodes.set(m[1],new Element(m[1]));
for(const m of html.matchAll(/<input\b([^>]+)>/g)){const attrs=Object.fromEntries([...m[1].matchAll(/([\w-]+)="([^"]*)"/g)].map(x=>[x[1],x[2]]));if(attrs.id){nodes.get(attrs.id).value=attrs.value||'';nodes.get(attrs.id).attrs=attrs;}}
nodes.get('matter03Evidence').textContent=JSON.stringify(evidence);
const root=new Element(),sandbox={document:{documentElement:root,getElementById:id=>nodes.get(id),createElement:()=>new Element()},window:{},console};
vm.createContext(sandbox);vm.runInContext(code,sandbox,{filename:htmlPath});
const api=sandbox.window.Matter03,gates=[];
function gate(name,passed,detail){gates.push({name,passed:!!passed,detail});if(!passed)throw new Error(name+': '+JSON.stringify(detail));}
function close(a,b){return Math.abs(a-b)<1e-12*Math.max(1,Math.abs(a),Math.abs(b));}
function eq(a,b){return JSON.stringify(a)===JSON.stringify(b);}
function click(id){nodes.get(id).events.click();}
function input(id,v){nodes.get(id).value=String(v);nodes.get(id).events.input();}
gate('no external resources',!/<(?:script|link|img|iframe)\b[^>]*(?:src|href)=["']https?:/i.test(html),true);
gate('no runtime network APIs',!(/\b(?:fetch|XMLHttpRequest|WebSocket)\s*\(/.test(code)),true);
gate('scope distinguishes conditional, input and diagnostic',html.includes('EXACT · CONDITIONAL')&&html.includes('INPUT ·')&&html.includes('DIAGNOSTIC ·'),true);
gate('physical mass inputs remain visibly open',html.includes('Selected light-family subspaces')&&html.includes('Surface harmonic mode metrics')&&evidence.open.includes('selected Higgs scale'),evidence.open);
gate('canonical twist and parent assumptions visible',html.includes('K<sub>S</sub>-valued')&&html.includes('no extra flat twists'),true);
gate('valution addition',evidence.exceptionalValuation===evidence.canonicalOrders.reduce((a,b)=>a+b,0)+evidence.jacobianOrder,evidence.exceptionalValuation);
gate('odd spectral valuation',evidence.exceptionalValuation%2===1,evidence.exceptionalValuation);
gate('paired initial state',api.state.paired&&api.status().good&&nodes.get('vacuumStatus').textContent.includes('exists'),api.status());
for(const t of [0,.01,.2,1,2]){
 input('branchScale',t);click('pairedButton');
 gate('pair exists including split endpoint t='+t,api.status().good,api.status());
 gate('branch t control and invariant t²='+t,api.state.branchScale===t&&nodes.get('spectralScale').textContent.includes((t*t).toFixed(4)),nodes.get('spectralScale').textContent);
 gate('zero branch distinguishes split from irreducible pair t='+t,nodes.get('vacuumStatus').textContent.includes(t===0?'Split zero-Higgs':'Paired nonzero'),nodes.get('vacuumStatus').textContent);
 click('singleButton');
 gate('single arrow no-go at nonzero t='+t,api.status().good===(t===0),api.status());
 gate('single map spectral product vanishes t='+t,nodes.get('spectralScale').textContent.includes('s(t) = 0'),nodes.get('spectralScale').textContent);
 gate('single arrow visibility t='+t,nodes.get('upArrow').attrs.opacity==='0'&&nodes.get('singleButton').attrs['aria-pressed']==='true',nodes.get('upArrow').attrs);
}
click('pairedButton');input('branchScale',1);
gate('six block rows render',nodes.get('blockRows').children.length===6,nodes.get('blockRows').children.length);
gate('six capacities total118',evidence.blocks.reduce((s,b)=>s+b.capacity,0)===118,evidence.blocks);
for(let i=0;i<6;i++){
 const b=evidence.blocks[i];gate('rectangular block capacity '+i,b.capacity===Math.min(...b.shape),b);
 click('blockButton'+i);gate('block inspection button '+i,api.state.block===i&&nodes.get('blockButton'+i).attrs['aria-pressed']==='true'&&nodes.get('blockDetail').textContent.includes(b.shape.join(' × ')),nodes.get('blockDetail').textContent);
}
for(const charge of ['minus','plus'])gate('split charged total '+charge,evidence.chargeCounts[charge].a+evidence.chargeCounts[charge].phi===279,evidence.chargeCounts[charge]);
const originalReview=JSON.parse(JSON.stringify(api.evidence.reviewed));
api.evidence.reviewed={aaRank:null,chargedPerSign:null,blockRanks:null};api.update();
gate('rank bounds give160 to275',api.counts().min===160&&api.counts().max===275&&!api.counts().exact,api.counts());
gate('bound display identifies capacity',nodes.get('pairedBar').innerHTML.includes('capacity')&&nodes.get('pairedCount').textContent==='160 to 275 classes',nodes.get('pairedBar').innerHTML);
api.evidence.reviewed={aaRank:41,chargedPerSign:237,blockRanks:[2,2,5,3,17,12]};api.update();
gate('reviewed data slot rank41 leads237',api.counts().exact&&api.counts().min===237&&nodes.get('heroModes').textContent==='237 per charge sign',api.counts());
gate('exact slot does not retain bound text',!nodes.get('pairedBar').innerHTML.includes('capacity')&&nodes.get('modeFormula').textContent==='rank(A) = 41 ⇒ dim 𝕳¹ = 237',nodes.get('modeFormula').textContent);
for(let i=0;i<6;i++)gate('reviewed block rank renders '+i,nodes.get('blockRank'+i).textContent===String([2,2,5,3,17,12][i]),nodes.get('blockRank'+i).textContent);
api.evidence.reviewed={aaRank:41,chargedPerSign:238,blockRanks:[2,2,5,3,17,12]};api.update();
gate('inconsistent review data falls back to bound',!api.counts().exact&&api.counts().min===160,api.counts());
api.evidence.reviewed={aaRank:41,chargedPerSign:237,blockRanks:[2,2,5,3,17,13]};api.update();
gate('inconsistent block-rank sum falls back to bound',!api.counts().exact&&api.counts().min===160,api.counts());
api.evidence.reviewed={aaRank:41,chargedPerSign:237,blockRanks:[15,2,5,3,4,12]};api.update();
gate('block rank above capacity falls back to bound',!api.counts().exact&&api.counts().min===160,api.counts());
api.evidence.reviewed=originalReview;api.update();
for(const u of [0,.1,.5,1,1.5,2])for(const d of [0,.1,.5,1,1.5,2]){
 const v=api.potential(u,d),expected=(u*u-d*d)**2/4;
 gate('toy quartic u='+u+' d='+d,close(v,expected)&&v>=0&&close(v,api.potential(d,u)),v);
}
for(const t of [0,.1,.5,1,2]){
 input('toyUp',t);input('toyDown',t);gate('toy zero ray t='+t,nodes.get('toyPotential').textContent==='0.000000',nodes.get('toyPotential').textContent);
}
input('toyUp',2);input('toyDown',1);
gate('off-ray toy value3²/4',nodes.get('toyPotential').textContent==='2.250000',nodes.get('toyPotential').textContent);
gate('plot draws supplied point',paint.some(c=>c[0]==='arc')&&paint.some(c=>c[0]==='fillText'&&c[1].includes('zero ray')),paint.length);
click('identityButton');gate('identity Grams three equal singular values',eq(api.masses(api.state.left,api.state.right,api.state.higgs,api.state.massScale),[1,1,1]),api.state);
click('exampleButton');
const ex=api.masses(api.state.left,api.state.right,api.state.higgs,api.state.massScale);
gate('positive illustrative hierarchy1,1/8,1/15',ex.every((v,i)=>close(v,[1,1/8,1/15][i])),ex);
gate('illustrative ratios render',nodes.get('ratioTwo').textContent==='0.125000'&&nodes.get('ratioThree').textContent==='0.066667',[nodes.get('ratioTwo').textContent,nodes.get('ratioThree').textContent]);
for(let k=0;k<24;k++){
 const l=[1+k/3,2+(k%5),3+(k%7)],r=[2+(k%3),1+k/5,4+(k%4)],h=1+(k%6),scale=.1+k/8;
 // Independently build M†M for the sparse raw minor. Its diagonal entries are squared singular values.
 const squares=[1/(l[0]*r[0]*h),1/(l[2]*r[1]*h),1/(l[1]*r[2]*h)].sort((a,b)=>b-a);
 const m=api.masses(l,r,h,scale);
 gate('independent squared-singular-value normalization '+k,m.every((v,i)=>close(v*v,squares[i]*scale*scale)),m);
 gate('positive normalization preserves nonzero rank '+k,m.every(v=>v>0)&&m[0]>=m[1]&&m[1]>=m[2],m);
}
input('massScale',0);gate('zero scale yields rank0 and undefined ratio',nodes.get('massRank').textContent.includes('0 ·')&&nodes.get('ratioTwo').textContent.includes('Undefined'),[nodes.get('massRank').textContent,nodes.get('ratioTwo').textContent]);
click('exampleButton');input('higgsGram',4);const h4=api.masses(api.state.left,api.state.right,api.state.higgs,api.state.massScale);
gate('Higgs Gram scales all singular values byonehalf',h4.every((v,i)=>close(v,ex[i]/2)),h4);
gate('Higgs Gram retains ratios',nodes.get('ratioTwo').textContent==='0.125000',nodes.get('ratioTwo').textContent);
input('massScale',2);gate('supplied mass scale doubles all values',api.masses(api.state.left,api.state.right,api.state.higgs,api.state.massScale).every((v,i)=>close(v,h4[i]*2)),true);
click('themeButton');gate('light theme interaction',root.dataset.theme==='light'&&nodes.get('themeButton').attrs['aria-pressed']==='true',root.dataset.theme);
click('themeButton');gate('dark theme interaction',root.dataset.theme==='dark'&&nodes.get('themeButton').attrs['aria-pressed']==='false',root.dataset.theme);
click('exampleButton');input('toyUp',1);input('toyDown',1);api.state.block=2;api.update();
const report={scope:'Offline numerical and mock-DOM interaction checks. No actual browser layout test was performed.',passed:gates.filter(g=>g.passed).length,failed:gates.filter(g=>!g.passed).length,gates};
if(outPath)fs.writeFileSync(outPath,JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({passed:report.passed,failed:report.failed,scope:report.scope},null,2));
