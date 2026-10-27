// Banc d'essai sur chromatogrammes simulés : node test/benchmark.js [cas par condition] [--grands]
// Mesure, par fraction allélique et niveau de bruit, la proportion de variants retrouvés
// exactement (nom HGVS identique) et le nombre d'appels faux.
import { analyzeRead, mergeCalls } from '../src/decompose.js';
import { buildReference } from '../src/reference.js';
import { normalize, describe } from '../src/hgvs.js';
import { revcomp } from '../src/sequence.js';
import { simulate, rng, randomSeq } from './synth.js';

const N = Number(process.argv[2]) || 40;
const LARGE = process.argv.includes('--grands');
const rand = rng(99);
const SEQ = randomSeq(8000, rand);
const REF = buildReference({ accession: 'NG_BENCH.1', seq: SEQ });

function makeVariant(r) {
  const g = 2000 + Math.floor(r() * 400); // centre de l'amplicon
  const kind = Math.floor(r() * (LARGE ? 6 : 4));
  const rb = () => 'ACGT'[Math.floor(r() * 4)];
  if (kind === 0) { let b; do b = rb(); while (b === SEQ[g - 1]); return { g, del: 1, ins: b }; }
  if (kind === 1) { const L = 1 + Math.floor(r() * 3); return { g, del: 0, ins: SEQ.slice(g - 1 - L, g - 1) }; } // dup
  if (kind === 2) { const L = 1 + Math.floor(r() * 10); return { g, del: L, ins: '' }; }
  if (kind === 3) { const L = 1 + Math.floor(r() * 5); let s = ''; for (let i = 0; i < L; i++) s += rb(); return { g, del: 0, ins: s }; }
  // Grands indels (option --grands) : délétion de 11 à 80 nt, insertion de 6 à 30 nt.
  if (kind === 4) { const L = 11 + Math.floor(r() * 70); return { g, del: L, ins: '' }; }
  const L = 6 + Math.floor(r() * 25); let s = ''; for (let i = 0; i < L; i++) s += rb(); return { g, del: 0, ins: s };
}

function expectedName(v) {
  const n = normalize(SEQ, { pos: v.g - 1, del: SEQ.slice(v.g - 1, v.g - 1 + v.del), ins: v.ins });
  return n.del === '' && n.ins === '' ? null : describe(n, String);
}

const conditions = [];
for (const fraction of [0.5, 0.35, 0.25]) for (const [noise, heightVar] of [[8, 0.35], [25, 0.5]]) conditions.push({ fraction, noise, heightVar });

console.log('fraction  bruit  variabilité  exacts(F/R)  au moins 1 lecture  faux appels/cas (dont confirmés F/R)');
for (const c of conditions) {
  let exact = 0, any = 0, fp = 0, fpConfirmed = 0, total = 0;
  for (let k = 0; k < N; k++) {
    const r = rng(1000 * k + 7);
    const v = makeVariant(r);
    const name = expectedName(v);
    if (!name) continue;
    total++;
    const wt = SEQ.slice(1849, 2650);
    const i = v.g - 1850;
    const mut = wt.slice(0, i) + v.ins + wt.slice(i + v.del);
    const opts = { noise: c.noise, heightVar: c.heightVar, length: 700 };
    const F = simulate([wt.slice(20), mut.slice(20)], [1 - c.fraction, c.fraction], { ...opts, seed: 2 * k + 1 }); F.name = 'F';
    const R = simulate([revcomp(wt).slice(20), revcomp(mut).slice(20)], [1 - c.fraction, c.fraction], { ...opts, seed: 2 * k + 2 }); R.name = 'R';
    let merged = [];
    try { merged = mergeCalls([F, R].map(x => analyzeRead(x, REF))); } catch { /* compté comme manqué */ }
    const hit = merged.find(m => describe(m, String) === name);
    if (hit && hit.status === 'confirmé F/R') exact++;
    else if (process.env.DEBUG) console.log("  MISS", c.fraction, c.noise, name, hit ? hit.status + " " + hit.reads.map(x => x.name).join(",") : "absent", merged.map(m => describe(m, String) + "[" + m.reads.map(x => x.name).join(",") + "]").join(" "));
    if (hit) any++;
    fp += merged.filter(m => m !== hit).length;
    fpConfirmed += merged.filter(m => m !== hit && m.status === "confirmé F/R").length;
    if (process.env.DEBUG) for (const m of merged) if (m !== hit) console.log("  FP", c.fraction, name, describe(m, String), m.status, m.reads.map(x => x.evidence).join(" | "));
  }
  const pct = x => (100 * x / total).toFixed(0).padStart(3) + ' %';
  console.log(`${String(c.fraction).padEnd(9)} ${String(c.noise).padEnd(6)} ${String(c.heightVar).padEnd(12)} ${pct(exact).padEnd(12)} ${pct(any).padEnd(19)} ${(fp / total).toFixed(2)} (${(fpConfirmed / total).toFixed(2)})`);
}
