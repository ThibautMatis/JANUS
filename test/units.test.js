import test from 'node:test';
import assert from 'node:assert/strict';
import { parseABIF } from '../src/abif.js';
import { normalize, transcriptMapper, describe } from '../src/hgvs.js';
import { parseGenBank, listTranscripts, cdsOfTranscript, parseLocation } from '../src/reference.js';
import { alignStrings } from '../src/align.js';
import { simulate, writeABIF, rng, randomSeq } from './synth.js';

test('ABIF : écriture puis lecture', () => {
  const read = simulate([randomSeq(120, rng(3))], [1]);
  const parsed = parseABIF(writeABIF(read));
  assert.equal(parsed.bases, read.bases);
  assert.deepEqual(parsed.positions, read.positions);
  assert.deepEqual(Array.from(parsed.traces.T), Array.from(read.traces.T));
  assert.equal(parsed.sample, 'synthetique');
});

test('ABIF : fichier invalide refusé', () => {
  assert.throws(() => parseABIF(new ArrayBuffer(10)), /tronqué/);
  assert.throws(() => parseABIF(new TextEncoder().encode('X'.repeat(200)).buffer), /Signature/);
});

test('normalisation : décalage 3′ et duplication', () => {
  const seq = 'GGCAAAATG';
  // Insertion d'un A n'importe où dans le bloc AAAA → duplication du dernier A.
  assert.deepEqual(normalize(seq, { pos: 3, del: '', ins: 'A' }), { pos: 7, del: '', ins: 'A', type: 'dup' });
  // Délétion d'un A → dernier A du bloc.
  assert.equal(normalize(seq, { pos: 3, del: 'A', ins: '' }).pos, 6);
  // Délétion de CA dans CACACA → décalée au bout.
  const v = normalize('TCACACAG', { pos: 1, del: 'CA', ins: '' });
  assert.equal(v.pos, 5); assert.equal(v.del, 'CA');
  assert.equal(normalize(seq, { pos: 2, del: 'C', ins: 'T' }).type, 'snv');
});

test('coordonnées c. : exons, introns, UTR', () => {
  const m = transcriptMapper({ exons: [[101, 200], [301, 400], [501, 600]], cdsStart: 51, cdsEnd: 250 });
  assert.equal(m.c(151), '1');
  assert.equal(m.c(150), '-1');
  assert.equal(m.c(200), '50');
  assert.equal(m.c(201), '50+1');
  assert.equal(m.c(300), '51-1');
  assert.equal(m.c(250), '50+50');  // milieu de l'intron (100 nt) → rattaché à l'amont
  assert.equal(m.c(251), '51-50');
  assert.equal(m.c(550), '200');
  assert.equal(m.c(551), '*1');
  assert.equal(m.c(90), '-61');
  assert.equal(describe({ pos: 152, del: '', ins: 'G', type: 'dup' }, m.c), '2dup');
  assert.equal(describe({ pos: 160, del: 'AC', ins: '', type: 'del' }, m.c), '11_12del');
  assert.equal(describe({ pos: 160, del: '', ins: 'TT', type: 'ins' }, m.c), '10_11insTT');
});

test('GenBank : localisations et transcrits', () => {
  assert.deepEqual(parseLocation('join(10..20,<30..>40)').intervals, [[10, 20], [30, 40]]);
  assert.equal(parseLocation('complement(join(1..5,8..9))').complement, true);
  const gb = `LOCUS       NG_TEST.1               60 bp    DNA     linear   PRI 01-JAN-2020
DEFINITION  Test gene (TEST), RefSeqGene on chromosome 1.
ACCESSION   NG_TEST
VERSION     NG_TEST.1
FEATURES             Location/Qualifiers
     source          1..60
                     /organism="Homo sapiens"
     mRNA            join(5..20,
                     31..50)
                     /gene="TEST"
                     /product="test protein, transcript
                     variant 1"
                     /transcript_id="NM_TEST.1"
     CDS             join(10..20,31..40)
                     /protein_id="NP_TEST.1"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
`;
  const [rec] = parseGenBank(gb);
  assert.equal(rec.accession, 'NG_TEST.1');
  assert.equal(rec.seq.length, 60);
  const [t] = listTranscripts(rec);
  assert.equal(t.id, 'NM_TEST.1');
  assert.deepEqual(t.exons, [[5, 20], [31, 50]]);
  assert.equal(t.product, 'test protein, transcript variant 1');
  assert.deepEqual(cdsOfTranscript({ features: [{ key: 'CDS', intervals: [[229, 10485]], complement: false }] }), { cdsStart: 229, cdsEnd: 10485, proteinId: '' });
});

test('alignement : gap affine', () => {
  const a = alignStrings('ACGTACGTTTTGCATGCA', 'ACGTACGTTTGCATGCA');
  assert.equal(a.ops.filter(o => o.op === 'I').length, 1);
});

test('CDS déduit des annotations génomiques', async () => {
  const { cdsFromGenomic } = await import('../src/reference.js');
  const record = { features: [{ key: 'CDS', complement: false, intervals: [[10, 20], [31, 40]], qualifiers: {} }] };
  // Transcrit : exons 5–20 et 31–50 → CDS de la position 6 à 26 du transcrit.
  assert.deepEqual(cdsFromGenomic(record, { exons: [[5, 20], [31, 50]], strand: '+' }), { cdsStart: 6, cdsEnd: 26, proteinId: '' });
  assert.equal(cdsFromGenomic(record, { exons: [[5, 22], [31, 50]], strand: '+' }), null);
});

test('localisation exon / intron', async () => {
  const { variantLocation } = await import('../src/hgvs.js');
  const m = transcriptMapper({ exons: [[101, 200], [301, 400], [501, 600]], cdsStart: 51, cdsEnd: 250 });
  assert.equal(variantLocation({ pos: 150, del: 'A', ins: 'G' }, m), 'exon 1');
  assert.equal(variantLocation({ pos: 249, del: 'A', ins: 'G' }, m), 'intron 1');
  assert.equal(variantLocation({ pos: 520, del: '', ins: 'T' }, m), 'exon 3');
  // Délétion de l'intron 1 jusque dans l'exon 2 (type MSH2 c.212-25_239del).
  assert.equal(variantLocation({ pos: 280, del: 'A'.repeat(40), ins: '' }, m), 'intron 1 – exon 2');
  assert.equal(variantLocation({ pos: 50, del: 'A', ins: 'G' }, m), 'en amont du gène');
  assert.equal(variantLocation({ pos: 450, del: 'A', ins: 'G' }, null), null);
});
