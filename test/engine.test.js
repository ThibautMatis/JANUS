// Tests de bout en bout sur chromatogrammes simulés : référence génomique aléatoire,
// transcrit fictif, amplicon couvrant intron et exon, lectures Forward et Reverse.
import test from 'node:test';
import assert from 'node:assert/strict';
import { analyzeRead, mergeCalls } from '../src/decompose.js';
import { buildReference } from '../src/reference.js';
import { revcomp } from '../src/sequence.js';
import { simulate, rng, randomSeq } from './synth.js';

const rand = rng(2024);
const SEQ = randomSeq(6000, rand);
// Exons (1-based) ; le CDS commence à la position 50 du transcrit.
const transcript = { id: 'NM_TEST.1', strand: '+', exons: [[201, 400], [1001, 1300], [2001, 2600], [3001, 3800]] };
const REF = buildReference({ accession: 'NG_TEST.1', seq: SEQ, transcript, cds: { cdsStart: 50, cdsEnd: 1200 } });
// Exon 3 : g.2001 = transcrit 501 = c.452.
const cOf = g => g - 2001 + 452;

// Amplicon g.1850–2650 (intron 2 + exon 3) ; lecture F depuis 1870, lecture R depuis 2630.
const AMP = [1850, 2650];

/** Applique des variants (1-based, sur la référence) et renvoie la séquence de l'amplicon. */
function amplicon(variants = []) {
  let s = SEQ.slice(AMP[0] - 1, AMP[1]);
  for (const v of [...variants].sort((a, b) => b.g - a.g)) {
    const i = v.g - AMP[0];
    s = s.slice(0, i) + v.ins + s.slice(i + v.del);
  }
  return s;
}

function reads(alt, { fraction = 0.5, seed = 1 } = {}) {
  const wt = amplicon(), mut = amplicon(alt);
  const fStart = 1870 - AMP[0];
  const F = simulate([wt.slice(fStart), mut.slice(fStart)], [1 - fraction, fraction], { seed, length: 700 });
  const rWt = revcomp(wt).slice(20), rMut = revcomp(mut).slice(20);
  const R = simulate([rWt, rMut], [1 - fraction, fraction], { seed: seed + 100, length: 700 });
  F.name = 'F'; R.name = 'R';
  return [F, R];
}

function run(alt, opts) {
  const results = reads(alt, opts).map(r => analyzeRead(r, REF));
  return { results, merged: mergeCalls(results) };
}

// Choisit une position g. dont la base diffère de la suivante (pas de décalage 3′).
function stablePos(g) { while (SEQ[g - 1] === SEQ[g]) g++; return g; }

test('témoin sans variant : aucun appel', () => {
  const { results, merged } = run([]);
  for (const r of results) {
    assert.ok(r.alignedBases > 550, `${r.name} aligné sur ${r.alignedBases} bases`);
  }
  assert.equal(results[0].strand, '+');
  assert.equal(results[1].strand, '-');
  assert.deepEqual(merged.map(m => m.names.c), []);
});

test('duplication hétérozygote 1 nt (type BRCA2 c.1813dup)', () => {
  const g = stablePos(2300);
  const { merged } = run([{ g: g + 1, del: 0, ins: SEQ[g - 1] }]);
  assert.equal(merged.length, 1, JSON.stringify(merged.map(m => m.names)));
  assert.equal(merged[0].names.c, `NM_TEST.1:c.${cOf(g)}dup`);
  assert.equal(merged[0].genotype, 'hétérozygote');
  assert.equal(merged[0].status, 'confirmé F/R');
  assert.ok(Math.abs(merged[0].reads[0].alleleFraction - 0.5) < 0.15);
});

test('délétion hétérozygote 4 nt', () => {
  const g = 2150;
  const del = SEQ.slice(g - 1, g + 3);
  const { merged } = run([{ g, del: 4, ins: '' }]);
  assert.equal(merged.length, 1, JSON.stringify(merged.map(m => m.names)));
  assert.equal(merged[0].type, 'del');
  assert.equal(merged[0].del.length, 4);
  // Position éventuellement décalée en 3′ : on vérifie la séquence supprimée à la position rapportée.
  assert.equal(SEQ.slice(merged[0].pos, merged[0].pos + 4), merged[0].del);
  assert.ok(merged[0].pos + 1 >= g && merged[0].pos + 1 <= g + 4, del);
  assert.equal(merged[0].status, 'confirmé F/R');
});

test('insertion hétérozygote 3 nt (non dupliquée)', () => {
  const g = stablePos(2400);
  const ins = SEQ[g] === 'T' ? 'GGC' : 'TTG';
  const { merged } = run([{ g: g + 1, del: 0, ins }]);
  assert.equal(merged.length, 1, JSON.stringify(merged.map(m => m.names)));
  assert.match(merged[0].names.c, /ins[ACGT]{3}$/);
  assert.equal(merged[0].status, 'confirmé F/R');
});

test('SNV hétérozygote intronique', () => {
  const g = 1960; // intron 2 : c.452-41
  const alt = { A: 'G', G: 'A', C: 'T', T: 'C' }[SEQ[g - 1]];
  const { merged } = run([{ g, del: 1, ins: alt }]);
  assert.equal(merged.length, 1, JSON.stringify(merged.map(m => m.names)));
  assert.equal(merged[0].names.c, `NM_TEST.1:c.452-41${SEQ[g - 1]}>${alt}`);
  assert.equal(merged[0].genotype, 'hétérozygote');
  assert.equal(merged[0].status, 'confirmé F/R');
});

test('délétion homozygote 2 nt', () => {
  const g = 2330;
  const { merged } = run([{ g, del: 2, ins: '' }], { fraction: 1 });
  assert.equal(merged.length, 1, JSON.stringify(merged.map(m => m.names)));
  assert.equal(merged[0].type, 'del');
  assert.equal(merged[0].genotype, 'homozygote');
  assert.equal(merged[0].status, 'confirmé F/R');
});

test('SNV homozygote', () => {
  const g = 2205, alt = { A: 'G', G: 'A', C: 'T', T: 'C' }[SEQ[g - 1]];
  const { merged } = run([{ g, del: 1, ins: alt }], { fraction: 1 });
  assert.equal(merged.length, 1, JSON.stringify(merged.map(m => m.names)));
  assert.equal(merged[0].names.c, `NM_TEST.1:c.${cOf(g)}${SEQ[g - 1]}>${alt}`);
  assert.equal(merged[0].genotype, 'homozygote');
});

test('duplication à 30 % (mosaïque / déséquilibre)', () => {
  const g = stablePos(2250);
  const { merged } = run([{ g: g + 1, del: 0, ins: SEQ[g - 1] }], { fraction: 0.3 });
  assert.equal(merged.length, 1, JSON.stringify(merged.map(m => m.names)));
  assert.equal(merged[0].names.c, `NM_TEST.1:c.${cOf(g)}dup`);
});

test('mauvaise référence : erreur explicite', () => {
  const other = buildReference({ accession: 'NG_OTHER.1', seq: randomSeq(5000, rng(9)) });
  assert.throws(() => analyzeRead(reads([])[0], other), /non localisée/);
});

test('grande délétion hétérozygote (53 nt, type MSH2 c.212-25_239del)', () => {
  const g = 1990; // de l'intron 2 jusque dans l'exon 3
  const { merged } = run([{ g, del: 53, ins: '' }]);
  assert.equal(merged.length, 1, JSON.stringify(merged.map(m => m.names)));
  assert.equal(merged[0].type, 'del');
  assert.equal(merged[0].del.length, 53);
  assert.equal(SEQ.slice(merged[0].pos, merged[0].pos + 53), merged[0].del);
  assert.equal(merged[0].status, 'confirmé F/R');
});

test('lecture commençant dans la zone décalée : confirmation par le modèle', () => {
  // Comme MSH2 Ex02F : les deux allèles partent de la même amorce, en amont de la délétion, mais les
  // premières bases sont illisibles ; la partie exploitable commence au milieu de la délétion.
  const g = 1990, L = 53;
  const wt = amplicon(), mut = amplicon([{ g, del: L, ins: '' }]);
  const full = simulate([wt.slice(120), mut.slice(120)], [0.5, 0.5], { seed: 5, length: 650 });
  const F = { ...full, bases: full.bases.slice(45), positions: full.positions.slice(45), quality: full.quality.slice(45), name: 'F' };
  const R = simulate([revcomp(wt).slice(20), revcomp(mut).slice(20)], [0.5, 0.5], { seed: 105, length: 700 });
  R.name = 'R';
  const results = [F, R].map(r => analyzeRead(r, REF));
  assert.ok(results[0].partial, 'décalage détecté dès le début de la lecture Forward');
  assert.equal(Math.abs(results[0].partial.ins.length - results[0].partial.del.length), L);
  const merged = mergeCalls(results, REF);
  assert.equal(merged.length, 1, JSON.stringify(merged.map(m => m.names)));
  assert.equal(merged[0].names.g, `NG_TEST.1:g.${g}_${g + L - 1}del`);
  assert.equal(merged[0].status, 'confirmé F/R');
  assert.deepEqual(merged[0].supportedBy.map(s => s.name), ['F']);
});

test('délétion d’un G dans un GGG, allèle muté majoritaire : pas de duplication « miroir »', () => {
  // Type BRCA2 c.6841+1del : si l'alignement suit l'allèle délété, l'allèle normal ressemble à une
  // duplication et le variant est nommé à l'envers. Cas qui échouait (fraction 70 %, graine 4).
  const motif = 'ACTTATCTTAGTGGGTAAGTGTTCATTTTTACCTT';
  const seq = (() => { const s = randomSeq(6000, rng(31)); return s.slice(0, 2280) + motif + s.slice(2280 + motif.length); })();
  const ref = buildReference({ accession: 'NG_G.1', seq });
  const g = 2280 + motif.indexOf('GGG') + 2;
  const wt = seq.slice(1849, 2650), i = g - 1849, mut = wt.slice(0, i) + wt.slice(i + 1);
  const opts = { length: 700, heightVar: 0.45, noise: 15 };
  const F = simulate([wt.slice(20), mut.slice(20)], [0.3, 0.7], { ...opts, seed: 4 });
  const R = simulate([revcomp(wt).slice(20), revcomp(mut).slice(20)], [0.3, 0.7], { ...opts, seed: 54 });
  F.name = 'F'; R.name = 'R';
  const merged = mergeCalls([F, R].map(r => analyzeRead(r, ref)), ref);
  assert.deepEqual(merged.map(m => m.names.g), [`NG_G.1:g.${g + 1}del`]);
  assert.equal(merged[0].status, 'confirmé F/R');
});

test('deux délétions en cis (1 nt + 4 nt à 80 nt d’intervalle) : les deux appelées', () => {
  const motif = 'ACTTATCTTAGTGGGTAAGTGTTCATTTTTACCTT';
  const seq = (() => { const s = randomSeq(6000, rng(31)); return s.slice(0, 2280) + motif + s.slice(2280 + motif.length); })();
  const ref = buildReference({ accession: 'NG_G.1', seq });
  const g1 = 2280 + motif.indexOf('GGG') + 2, g2 = g1 + 80;
  const wt = seq.slice(1849, 2650), cut = (s, g, n) => s.slice(0, g - 1849) + s.slice(g - 1849 + n);
  const mut = cut(cut(wt, g2, 4), g1, 1);
  const opts = { length: 700, heightVar: 0.45, noise: 15 };
  const F = simulate([wt.slice(20), mut.slice(20)], [0.35, 0.65], { ...opts, seed: 2 });
  const R = simulate([revcomp(wt).slice(20), revcomp(mut).slice(20)], [0.35, 0.65], { ...opts, seed: 52 });
  F.name = 'F'; R.name = 'R';
  const merged = mergeCalls([F, R].map(r => analyzeRead(r, ref)), ref);
  assert.deepEqual(merged.map(m => [m.names.g, m.status]), [
    [`NG_G.1:g.${g1 + 1}del`, 'confirmé F/R'], [`NG_G.1:g.${g2 + 1}_${g2 + 4}del`, 'confirmé F/R']]);
});

test('substitution en tête d’une longue série de A + série plus courte en cis (type MSH2 c.942+3A>T)', () => {
  // Fin d'exon …CAG | GT A27 GGG… : A>T sur le premier A, et 3 A en moins sur le même allèle.
  const motif = 'AACCTTTTTCAGGT' + 'A'.repeat(27) + 'GGGTTAAAAATGTTGA';
  const seq = (() => { const s = randomSeq(6000, rng(41)); return s.slice(0, 2200) + motif + s.slice(2200 + motif.length); })();
  const ref = buildReference({ accession: 'NG_A.1', seq });
  const g = 2200 + 14; // 0-based : premier A de la série
  const wt = seq.slice(1849, 2650), i = g - 1849;
  const mut = wt.slice(0, i) + 'T' + wt.slice(i + 1, i + 10) + wt.slice(i + 13);
  const opts = { length: 650, heightVar: 0.4, noise: 12 };
  const F = simulate([wt.slice(30), mut.slice(30)], [0.5, 0.5], { ...opts, seed: 1 });
  const R = simulate([revcomp(wt).slice(30), revcomp(mut).slice(30)], [0.5, 0.5], { ...opts, seed: 10 });
  F.name = 'F'; R.name = 'R';
  const merged = mergeCalls([F, R].map(r => analyzeRead(r, ref)), ref);
  const snv = merged.find(m => m.names.g === `NG_A.1:g.${g + 1}A>T`);
  assert.ok(snv, JSON.stringify(merged.map(m => m.names.g)));
  assert.equal(snv.status, 'confirmé F/R');
  const del = merged.find(m => m.type === 'del' && m.del === 'AAA');
  assert.ok(del, 'délétion de 3 A appelée');
  assert.equal(del.status, 'confirmé F/R');
});

test('séries répétées (A, C, G, T) : substitution en bord de série + indel dans la série, sur le même allèle', () => {
  // Plusieurs écritures sont équivalentes : on compare l'allèle muté reconstruit au vrai.
  const other = b => ({ A: 'T', C: 'G', G: 'C', T: 'A' }[b]);
  let seed = 100;
  for (const base of 'ACGT') for (const [size, where] of [[-1, 'fin'], [1, 'début'], [-3, 'fin']]) {
    seed++;
    const len = 15, s0 = randomSeq(6000, rng(seed));
    const flank = n => randomSeq(n, rng(seed * 7)).replace(new RegExp(base, 'g'), other(base));
    const motif = flank(12) + base.repeat(len) + flank(12);
    const seq = s0.slice(0, 2200) + motif + s0.slice(2200 + motif.length);
    const ref = buildReference({ accession: 'NG_H.1', seq });
    const start = 2212, sub = where === 'début' ? start : start + len - 1;
    const wt = seq.slice(1849, 2650), at = i => i - 1849;
    let mut = wt.slice(0, at(sub)) + other(base) + wt.slice(at(sub) + 1);
    const mid = at(start + 7);
    mut = size < 0 ? mut.slice(0, mid) + mut.slice(mid - size) : mut.slice(0, mid) + base.repeat(size) + mut.slice(mid);
    const opts = { length: 650, heightVar: 0.4, noise: 12 };
    const F = simulate([wt.slice(30), mut.slice(30)], [0.5, 0.5], { ...opts, seed });
    const R = simulate([revcomp(wt).slice(30), revcomp(mut).slice(30)], [0.5, 0.5], { ...opts, seed: seed + 500 });
    F.name = 'F'; R.name = 'R';
    const merged = mergeCalls([F, R].map(r => analyzeRead(r, ref)), ref);
    let hap = wt;
    for (const m of merged.filter(x => x.status === 'confirmé F/R').sort((a, b) => b.pos - a.pos)) hap = hap.slice(0, at(m.pos)) + m.ins + hap.slice(at(m.pos) + m.del.length);
    assert.equal(hap, mut, `${base}${len} ${size} ${where} : ${merged.map(m => m.names.g + ' ' + m.status).join(', ')}`);
  }
});
