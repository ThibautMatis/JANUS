// Banc d'essai : 100 délétions-insertions aléatoires + 100 délétions avec SNP non contigu (en cis).
// node test/stress-delins.js [fichier de sortie JSON]   (KIND=delins|del+SNP, NOISY=1 : conditions bruitées seules)
const root = new URL('..', import.meta.url).pathname;
const { analyzeRead, mergeCalls } = await import(root + '/src/decompose.js');
const { buildReference } = await import(root + '/src/reference.js');
const { normalize, describe } = await import(root + '/src/hgvs.js');
const { revcomp } = await import(root + '/src/sequence.js');
const { simulate, rng, randomSeq } = await import(root + '/test/synth.js');
const fs = await import('node:fs');

const SEQ = randomSeq(8000, rng(2027));
const REF = buildReference({ accession: 'NG_STRESS.1', seq: SEQ });
const A0 = 1849, A1 = 2650; // amplicon (0-based)
const N = 100;
const base = r => 'ACGT'[Math.floor(r() * 4)];
const name = v => describe(v, String);

function makeDelins(k) {
  const r = rng(9001 + 7919 * k); r(); r();
  for (;;) {
    const pos = 2050 + Math.floor(r() * 300), d = 1 + Math.floor(r() * 20), L = 1 + Math.floor(r() * 6);
    let ins = ''; for (let i = 0; i < L; i++) ins += base(r);
    const v = normalize(SEQ, { pos, del: SEQ.slice(pos, pos + d), ins });
    if (v.type === 'delins') return { kind: 'delins', vars: [{ pos, del: d, ins }], expect: [v], d, L };
  }
}
function makeDelSnp(k) {
  const r = rng(4242 + 104729 * k); r(); r();
  for (;;) {
    const d = 1 + Math.floor(r() * 20), gap = 1 + Math.floor(r() * 8), before = r() < 0.5;
    const pdel = 2050 + Math.floor(r() * 300);
    const psnp = before ? pdel - gap - 1 : pdel + d + gap;
    let b; do b = base(r); while (b === SEQ[psnp]);
    const dv = normalize(SEQ, { pos: pdel, del: SEQ.slice(pdel, pdel + d), ins: '' });
    // La délétion normalisée (décalée en 3′) doit rester séparée du SNP par au moins une base.
    const lo = dv.pos, hi = dv.pos + dv.del.length;
    if (psnp >= lo - 1 && psnp <= hi) continue;
    if (dv.pos !== pdel && !before && psnp < hi + 1) continue;
    const sv = normalize(SEQ, { pos: psnp, del: SEQ[psnp], ins: b });
    return { kind: 'del+SNP', vars: [{ pos: pdel, del: d, ins: '' }, { pos: psnp, del: 1, ins: b }], expect: [sv, dv].sort((x, y) => x.pos - y.pos), d, gap, before };
  }
}
const apply = (s, vars) => { for (const v of [...vars].sort((a, b) => b.pos - a.pos)) { const i = v.pos - A0; s = s.slice(0, i) + v.ins + s.slice(i + v.del); } return s; };
const applyCalls = (s, calls) => { for (const v of [...calls].sort((a, b) => b.pos - a.pos)) { const i = v.pos - A0; s = s.slice(0, i) + v.ins + s.slice(i + v.del.length); } return s; };

const ONLY = process.env.KIND;
const cases0 = [...Array.from({ length: N }, (_, k) => makeDelins(k)), ...Array.from({ length: N }, (_, k) => makeDelSnp(k))];
const cases = ONLY ? cases0.filter(c => c.kind === ONLY) : cases0;
const conditions = (process.env.NOISY ? [{ f: 0.5, noise: 25, hv: 0.5 }, { f: 0.35, noise: 25, hv: 0.5 }] : [{ f: 0.5, noise: 8, hv: 0.35 }, { f: 0.5, noise: 25, hv: 0.5 }, { f: 0.35, noise: 8, hv: 0.35 }, { f: 0.35, noise: 25, hv: 0.5 }]);
const out = [];
for (const c of conditions) {
  for (const [k, cs] of cases.entries()) {
    const wt = SEQ.slice(A0, A1), mut = apply(wt, cs.vars);
    const o = { noise: c.noise, heightVar: c.hv, length: 700 };
    const F = simulate([wt.slice(20), mut.slice(20)], [1 - c.f, c.f], { ...o, seed: 2 * k + 11 }); F.name = 'F';
    const R = simulate([revcomp(wt).slice(20), revcomp(mut).slice(20)], [1 - c.f, c.f], { ...o, seed: 2 * k + 12 }); R.name = 'R';
    let merged = [];
    try { merged = mergeCalls([F, R].map(x => analyzeRead(x, REF, JSON.parse(process.env.OPTS || '{}'))), REF); } catch (e) { merged = []; }
    const conf = merged.filter(m => m.status === 'confirmé F/R');
    const got = conf.map(name).sort(), exp = cs.expect.map(name).sort();
    const exact = got.length === exp.length && got.every((g, i) => g === exp[i]) && merged.length === conf.length;
    const hap = applyCalls(wt, conf) === mut;
    out.push({ cond: `${c.f}/${c.noise}`, kind: cs.kind, d: cs.d, L: cs.L, gap: cs.gap, exp, got: merged.map(m => `${name(m)}${m.status === 'confirmé F/R' ? '' : ' [' + m.status + ']'}`), exact, hap, nonConfirmed: merged.length - conf.length });
  }
}
if (process.argv[2]) fs.writeFileSync(process.argv[2], JSON.stringify(out));
for (const cond of [...new Set(out.map(o => o.cond))]) for (const kind of ['delins', 'del+SNP']) {
  const rows = out.filter(o => o.cond === cond && o.kind === kind);
  const ex = rows.filter(o => o.exact).length, hp = rows.filter(o => o.hap).length, nc = rows.reduce((s, o) => s + o.nonConfirmed, 0);
  console.log(`${cond.padEnd(8)} ${kind.padEnd(8)} nom exact ${String(ex).padStart(3)}/${rows.length}   allèle juste ${String(hp).padStart(3)}/${rows.length}   appels non confirmés ${nc}`);
}
