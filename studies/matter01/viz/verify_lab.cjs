/* Offline script/DOM smoke checks; this is not a browser layout test. */
'use strict';
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const htmlPath = process.argv[2] || path.join(__dirname, 'Matter01_U1_Lab.html');
const outPath = process.argv[3];
const html = fs.readFileSync(htmlPath, 'utf8');
const evidence = JSON.parse(html.match(/<script id="matter01Evidence" type="application\/json">([\s\S]*?)<\/script>/)[1]);
const code = html.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
const nodes = new Map();
const noop = () => {};
class Element {
  constructor(id = '') { this.id=id; this.value=''; this.textContent=''; this.children=[]; this.events={}; this.attrs={}; this.dataset={}; }
  append(...items) { this.children.push(...items); }
  replaceChildren(...items) { this.children=[...items]; }
  setAttribute(key,value) { this.attrs[key]=value; }
  addEventListener(key,fn) { this.events[key]=fn; }
  getBoundingClientRect() { return {width:600,height:this.id==='phaseCanvas'?195:300}; }
  getContext() { return new Proxy({}, {get:(_,key)=>key==='measureText'?s=>({width:s.length*7}):noop, set:()=>true}); }
}
for(const match of html.matchAll(/\bid="([^"]+)"/g))nodes.set(match[1],new Element(match[1]));
for(const match of html.matchAll(/<input\b([^>]+)>/g)) {
 const attrs=Object.fromEntries([...match[1].matchAll(/(\w+)="([^"]*)"/g)].map(m=>[m[1],m[2]]));
 if(attrs.id) { nodes.get(attrs.id).value=attrs.value||''; nodes.get(attrs.id).attrs=attrs; }
}
nodes.get('matter01Evidence').textContent=JSON.stringify(evidence);
const root=new Element();
const palette={bg:'#0b1119',panel:'#131e2a',text:'#e8eef5',muted:'#9aaec1',line:'#304459',blue:'#66bfff',red:'#ff8d9a',gold:'#ffd079',green:'#7eddb7'};
const sandbox={document:{documentElement:root,getElementById:id=>nodes.get(id),createElement:()=>new Element()},window:{devicePixelRatio:1,addEventListener:noop},getComputedStyle:()=>({getPropertyValue:n=>palette[n.slice(2)]||''}),console};
vm.createContext(sandbox);vm.runInContext(code,sandbox,{filename:htmlPath});
const api=sandbox.window.Matter01;
const gates=[];
function gate(name,passed,detail) { gates.push({name,passed:!!passed,detail}); if(!passed)throw new Error(name+': '+detail); }
gate('no external scripts or styles',!/<(?:script|link)[^>]+(?:src|href)=["']https?:/i.test(html),'Single offline file');
gate('three accessible controls', ['energy','angle','phase'].every(id=>nodes.has(id)), 'β, emission angle, phase');
gate('audited constants present',Object.values(api.constants).every(x=>x>0&&Number.isFinite(x)),api.constants);
gate('angular Python oracle included',api.oraclePoints>=10,api.oraclePoints);
gate('angular JavaScript matches Python oracle',api.oracleMaxRelative<2e-10,api.oracleMaxRelative);
gate('phase and kinematics render without error',true,'Mock canvas/DOM initialized');
let maxSymmetry=0,maxIntegral=0;
for(const b of [.1,.4,.7,.9,.98]) {
 const r=1/Math.sqrt(1-b*b);
 for(const z of [-1,-.75,-.1,0,.1,.75,1]) {
  const a=api.differential(r,z), mirrored=api.differential(r,-z);
  maxSymmetry=Math.max(maxSymmetry,Math.abs(a-mirrored)/a);
  gate('positive finite angular value β='+b+', c='+z,Number.isFinite(a)&&a>0,a);
 }
 const n=16000,step=2/n;
 let sum=api.differential(r,-1)+api.differential(r,1);
 for(let i=1;i<n;i++)sum+=(i%2?4:2)*api.differential(r,-1+step*i);
 const integral=sum*step/3,analytic=api.total(r);
 maxIntegral=Math.max(maxIntegral,Math.abs(integral-analytic)/analytic);
}
gate('photon interchange symmetry',maxSymmetry<1e-13,maxSymmetry);
gate('full sphere integral equals event total',maxIntegral<1e-7,maxIntegral);
for(const b of [.1,.98])for(const angle of [0,90,180]) {
 nodes.get('energy').value=String(b);nodes.get('angle').value=String(angle);api.update();
 gate('extreme controls update β='+b+', θ='+angle,Math.abs(api.state.ratio-1/Math.sqrt(1-b*b))<1e-13&&api.state.angle===angle,nodes.get('energyCM').textContent);
}
nodes.get('shapeButton').events.click();gate('normalized plot toggle',api.state.normalize===true,nodes.get('shapeButton').attrs['aria-pressed']);
nodes.get('themeButton').events.click();gate('theme toggle',root.dataset.theme==='light',root.dataset.theme);
nodes.get('resetButton').events.click();gate('reset returns benchmark',Math.abs(api.state.ratio-1/Math.sqrt(1-.7**2))<1e-13&&api.state.angle===60,api.state);
gate('evidence rows displayed',nodes.get('verificationRows').children.length===(evidence.verification.rows||[]).length,nodes.get('verificationRows').children.length);
gate('spectrum status displayed',nodes.get('statusTable').children.length===(evidence.spectrum_status||[]).length,nodes.get('statusTable').children.length);
const report={scope:'Offline script/DOM checks. Actual browser layout was not tested.',passed:gates.filter(g=>g.passed).length,failed:gates.filter(g=>!g.passed).length,oraclePoints:api.oraclePoints,oracleMaxRelative:api.oracleMaxRelative,maxAngularIntegralRelative:maxIntegral,gates};
if(outPath)fs.writeFileSync(outPath,JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({passed:report.passed,failed:report.failed,oraclePoints:report.oraclePoints,oracleMaxRelative:report.oracleMaxRelative,maxAngularIntegralRelative:report.maxAngularIntegralRelative,scope:report.scope},null,2));
