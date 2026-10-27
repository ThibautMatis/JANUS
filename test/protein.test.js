import test from 'node:test';
import assert from 'node:assert/strict';
import { buildReference } from '../src/reference.js';
import { normalize } from '../src/hgvs.js';
import { proteinEffect, translate } from '../src/protein.js';

// Gène fictif : exon 1 = g.11–40 (UTR 5′ de 10 nt + 20 nt codants), intron, exon 2 = g.101–160.
// CDS : ATG GCC AAA CGG TGG (c.1–15) | TTC ... jusqu'au stop TAA.
const exon1 = 'CCCCCCCCCC' + 'ATGGCCAAACGGTGGTTCAG';   // 30 nt
const exon2 = ('CCGAAAAAAGGGCTGTTAA' + 'G'.repeat(60)).slice(0, 60); // 60 nt
const seq = 'T'.repeat(10) + exon1 + ('GTAAG' + 'T'.repeat(52) + 'CAG') + exon2 + 'T'.repeat(40);
const ref = buildReference({
  accession: 'NG_P.1', seq,
  transcript: { id: 'NM_P.1', strand: '+', exons: [[11, 40], [101, 160]] },
  cds: { cdsStart: 11, cdsEnd: 49, proteinId: 'NP_P.1' },
});
// Variant décrit en c. → position génomique 0-based : c.N (exon 1) = g.(20 + N) ; exon 2 commence en c.21.
const g0 = c => (c <= 20 ? 20 + c : 100 + (c - 20)) - 1;
const eff = (c, del, ins) => proteinEffect(normalize(seq, { pos: g0(c), del, ins }), ref);

test('traduction de référence', () => {
  const spliced = exon1 + exon2;
  assert.equal(translate(spliced.slice(10)), 'MAKRWFSRKKGC*');
});

test('substitutions : faux-sens, synonyme, non-sens, perte du codon d’initiation', () => {
  assert.deepEqual(eff(10, 'C', 'T'), { p: 'Arg4Trp', effect: 'faux-sens' });       // CGG → TGG
  assert.deepEqual(eff(6, 'C', 'T'), { p: '=', effect: 'synonyme' });               // GCC → GCT
  assert.deepEqual(eff(13, 'T', 'A'), { p: 'Trp5Arg', effect: 'faux-sens' });
  assert.deepEqual(eff(15, 'G', 'A'), { p: 'Trp5*', effect: 'non-sens' });        // TGG → TGA
  assert.deepEqual(eff(2, 'T', 'C'), { p: 'Met1?', effect: 'perte du codon d’initiation' });
});

test('indels : décalage du cadre, délétion en phase, duplication', () => {
  // Délétion de 1 nt dans GCC (c.6) : ATG GCA AAC GGT GGT TCA GCC GAA AAA AGG GCT GTT AAG… :
  // Ala2 inchangé (GCA), premier acide aminé modifié Lys3 → Asn ; pas de stop dans le nouveau cadre.
  const fs = eff(6, 'C', '');
  assert.equal(fs.effect, 'décalage du cadre de lecture');
  assert.equal(fs.p, 'Lys3Asnfs*?');
  // Délétion de 2 nt (c.7_8) : ATG GCC ACG GTG GTT CAG CCG AAA AAA GGG CTG TTA A… → stop introuvable aussi ;
  // délétion de 1 nt en c.16 : … TGG TCA GCC GAA … → premier changé Phe6Ser.
  assert.match(eff(16, 'T', '').p, /^Phe6Serfs\*/);
  assert.deepEqual(eff(7, 'AAA', ''), { p: 'Lys3del', effect: 'délétion en phase' });
  // Duplication des 3 nt d'un codon Lys : en phase, décrite comme duplication de l'acide aminé.
  const dup = proteinEffect(normalize(seq, { pos: g0(10) - 1, del: '', ins: 'AAA' }), ref);
  assert.equal(dup.effect, 'duplication en phase');
  assert.equal(dup.p, 'Lys3dup');
});

test('non codant et épissage', () => {
  assert.equal(proteinEffect(normalize(seq, { pos: 40, del: 'G', ins: 'A' }), ref).p, '?');  // c.20+1
  assert.equal(proteinEffect(normalize(seq, { pos: 60, del: 'T', ins: 'A' }), ref).effect, 'intronique');
  assert.equal(proteinEffect(normalize(seq, { pos: 12, del: 'C', ins: 'A' }), ref).effect, 'UTR 5′');
  // Délétion chevauchant la fin de l'exon 1 et le début de l'intron.
  assert.equal(proteinEffect(normalize(seq, { pos: 36, del: seq.slice(36, 46), ins: '' }), ref).p, '?');
});
