// Runs the Modular Plane lab's own core (version 8.3) on inputs.json and prints what it computes, as JSON.
//     node probe_core.cjs /path/to/the/lab/core.js inputs.json > lab_core_outputs.json
// The core is the lab's and is not shipped with the package; its SHA-256 is recorded in the output and in README.md.
const fs = require('fs'), crypto = require('crypto'), path = require('path');
const corePath = path.resolve(process.argv[2]), inPath = path.resolve(process.argv[3]);
const MP = require(corePath);
const inp = JSON.parse(fs.readFileSync(inPath, 'utf8'));
const und = (v) => (v === undefined ? null : v);
const kernel = inp.matrices.map((M) => {
  const R = MP.numRange(M), info = MP.matrixInfo(M), dec = MP.decompose(M), J = MP.jacobiMonodromy(M);
  return {
    range: { kind: R.kind, centre: R.centre, semiMajor: R.semiMajor, semiMinor: R.semiMinor, w: R.w, norm2: R.norm2, ratio: R.ratio, real: R.real },
    kind: info.kind, tr: info.tr, det: info.det, disc: info.disc, length: und(info.length), psi: und(info.psi),
    order: info.order === undefined ? null : (info.order === Infinity ? 'inf' : info.order),
    dec: { alpha: dec.alpha, s1: dec.s1, s2: dec.s2, w: dec.w, beta: dec.beta }, region: MP.region(dec.alpha, dec.beta),
    expm: MP.expm(M, inp.t), jac: J.kind === 'hyp' ? { trM2: J.trM2, twoCoshL: J.twoCoshL } : null,
  };
});
const sites = inp.sites.map((v) => {
  const s = MP.siteSearch(v);
  return { here: s.here.map((q) => q[0]), products: s.products, mates: s.mates.length, Delta: s.Delta, eps: s.eps, planck: s.planck, dim: MP.dimStr(v) };
});
const sha = crypto.createHash('sha256').update(fs.readFileSync(corePath)).digest('hex');
console.log(JSON.stringify({ core: 'The Modular Plane, version 8.3, core', core_sha256: sha, node: process.version, t: inp.t, planck_index: MP.PLANCK_INDEX,
  names: inp.names, matrices: inp.matrices, kernel: kernel, sites: inp.sites, quantities: sites }));
