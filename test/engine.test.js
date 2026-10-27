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
