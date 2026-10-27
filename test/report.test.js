import test from 'node:test';
import assert from 'node:assert/strict';
import { conclusion, CLASSES } from '../src/report.js';

const ref = { accession: 'NG_007110.2', transcript: { id: 'NM_000251.2', gene: 'MSH2' }, mapper: { c: g => String(g - 100) } };
const results = [{ refStart: 150, refEnd: 400, quality: { verdict: 'bonne' } }, { refStart: 200, refEnd: 450, quality: { verdict: 'illisible' } }];
const del = {
  names: { c: 'NM_000251.2:c.212-25_239del', g: 'NG_007110.2:g.10253_10305del' },
  protein: { short: 'p.?' }, genotype: 'hétérozygote', status: 'confirmé F/R',
};

test('conclusion : variant classé', () => {
  assert.equal(conclusion([{ variant: del, classId: '5' }], ref, results),
    'Détection à l’état hétérozygote du variant pathogène NM_000251.2:c.212-25_239del p.? dans le gène MSH2 (classe 5).');
  assert.match(conclusion([{ variant: del, classId: '3' }], ref, results), /variant de signification incertaine \(VUS\) .* \(classe 3\)\.$/);
  assert.equal(CLASSES.length, 5);
});

test('conclusion : à classer, non confirmé, plusieurs variants', () => {
  const snv = { names: { c: 'NM_000251.2:c.356T>G', g: 'NG_007110.2:g.10422T>G' }, protein: { short: 'p.(Leu119Trp)' }, genotype: 'hétérozygote', status: 'une seule lecture' };
  const text = conclusion([{ variant: del, classId: '' }, { variant: snv, classId: '1' }], ref, results);
  const [a, b] = text.split('\n');
  assert.match(a, /\[classification à préciser\]/);
  assert.match(b, /variant bénin NM_000251.2:c.356T>G p.\(Leu119Trp\) .* \(classe 1\)\. Variant non confirmé sur les deux brins/);
});

test('conclusion : aucun variant retenu', () => {
  // Seule la lecture exploitable compte pour la région analysée (g.150–400 → c.50_300).
  assert.equal(conclusion([], ref, results), 'Absence de variant détecté dans la région analysée du gène MSH2 (NM_000251.2, c.50_300).');
});
