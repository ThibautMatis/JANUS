import test from 'node:test';
import assert from 'node:assert/strict';
import { anchorReference, gnomadId } from '../src/gnomad.js';
import { revcomp } from '../src/sequence.js';
import { normalize } from '../src/hgvs.js';
import { rng, randomSeq } from './synth.js';

// Chromosome fictif : la référence NG_ en est un extrait, sur le brin + ou −.
const chr = randomSeq(20000, rng(77)).replace(/^(.{5000})(.{8})/, '$1AAAAAAAA'); // poly-A en 5001–5008
const chrStart = 1000001; // position 1-based de chr[0]
const ngPlus = chr.slice(4000, 12000);            // NG base 0 = chr index 4000
const ngMinus = revcomp(chr.slice(4000, 12000));  // NG base 0 = chr index 11999
const anchors = [100, 2500, 5000, 7500];

test('repérage de la référence sur le chromosome (+ et −)', () => {
  const p = anchorReference(ngPlus, chr, chrStart, anchors);
  assert.equal(p.strand, '+');
  assert.equal(p.toChr(0), chrStart + 4000);
  const m = anchorReference(ngMinus, chr, chrStart, anchors);
  assert.equal(m.strand, '-');
  assert.equal(m.toChr(0), chrStart + 11999);
});

test('identifiants gnomAD : SNV, délétion, insertion, normalisation à gauche', () => {
  const p = { chrom: '2', chrSeq: chr, chrStart, ...anchorReference(ngPlus, chr, chrStart, anchors) };
  const m = { chrom: '2', chrSeq: chr, chrStart, ...anchorReference(ngMinus, chr, chrStart, anchors) };
  // SNV en chr index 4100.
  const b = chr[4100], alt = b === 'A' ? 'G' : 'A';
  assert.equal(gnomadId({ pos: 100, del: b, ins: alt }, p), `2-${chrStart + 4100}-${b}-${alt}`);
  assert.equal(gnomadId({ pos: 11999 - 4100, del: revcomp(b), ins: revcomp(alt) }, m), `2-${chrStart + 4100}-${b}-${alt}`);
  // Délétion d'un A dans le poly-A (chr 5000–5007) : gnomAD l'ancre à gauche du bloc.
  const expectDel = `2-${chrStart + 4999}-${chr[4999]}A-${chr[4999]}`;
  const vPlus = normalize(ngPlus, { pos: 1003, del: 'A', ins: '' });   // décalée en 3′ (HGVS) sur NG+
  assert.equal(gnomadId(vPlus, p), expectDel);
  const vMinus = normalize(ngMinus, { pos: 11999 - 5003, del: 'T', ins: '' });
  assert.equal(gnomadId(vMinus, m), expectDel);
  // Insertion (duplication d'un A) : même ancrage à gauche.
  const expectIns = `2-${chrStart + 4999}-${chr[4999]}-${chr[4999]}A`;
  assert.equal(gnomadId(normalize(ngPlus, { pos: 1004, del: '', ins: 'A' }), p), expectIns);
  assert.equal(gnomadId(normalize(ngMinus, { pos: 11999 - 5004, del: '', ins: 'T' }), m), expectIns);
});
