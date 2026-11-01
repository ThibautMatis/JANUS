// Contrôle direct : alignement des bases du séquenceur sur le génome, indépendant de la déconvolution.
import test from 'node:test';
import assert from 'node:assert/strict';
import { analyzeRead, mergeCalls } from '../src/decompose.js';
import { buildReference } from '../src/reference.js';
import { directCheck, explainDiscrepancy, explainShift } from '../src/check.js';
import { revcomp } from '../src/sequence.js';
import { simulate, rng, randomSeq } from './synth.js';

const motif = 'ACTTATCTTAGTGGGTAAGTGTTCATTTTTACCTT';
const seq = (() => { const s = randomSeq(6000, rng(31)); return s.slice(0, 2280) + motif + s.slice(2280 + motif.length); })();
const REF = buildReference({ accession: 'NG_G.1', seq });
const wt = seq.slice(1849, 2650);
const opts = { length: 700, heightVar: 0.45, noise: 15 };

function pair(mut, fraction, seeds = [4, 54]) {
  const F = simulate([wt.slice(20), mut.slice(20)], [1 - fraction, fraction], { ...opts, seed: seeds[0] });
  const R = simulate([revcomp(wt).slice(20), revcomp(mut).slice(20)], [1 - fraction, fraction], { ...opts, seed: seeds[1] });
  F.name = 'F'; R.name = 'R';
  return [F, R];
}
const alerts = (results, merged) => results.map(r => {
  const c = directCheck(r, REF);
  return {
    strand: c.strandAgrees,
    unexplained: c.discrepancies.filter(d => !explainDiscrepancy(d, r, merged)).length,
    shifts: c.mapping.ranges.filter(x => !explainShift(x, r, merged)).length,
    skipped: c.skippedInMixed.length,
  };
});

test('lecture normale : identité élevée, aucune alerte', () => {
  const results = pair(wt, 0).map(r => analyzeRead(r, REF));
  const merged = mergeCalls(results, REF);
  assert.deepEqual(alerts(results, merged), [
    { strand: true, unexplained: 0, shifts: 0, skipped: 0 }, { strand: true, unexplained: 0, shifts: 0, skipped: 0 }]);
  assert.ok(directCheck(results[0], REF).identity > 0.99);
});

test('indel hétérozygote : décalage et discordances expliqués par le variant', () => {
  const g = 2200, mut = wt.slice(0, g - 1849) + wt.slice(g - 1849 + 7);
  const results = pair(mut, 0.6).map(r => analyzeRead(r, REF));
  const merged = mergeCalls(results, REF);
  assert.equal(merged.length, 1);
  for (const a of alerts(results, merged)) assert.deepEqual(a, { strand: true, unexplained: 0, shifts: 0, skipped: 0 });
});

test('alignement JANUS sur l’allèle muté (ancien réglage) : base du génome sautée en zone mixte signalée', () => {
  const g = 2280 + motif.indexOf('GGG') + 2, mut = wt.slice(0, g - 1849) + wt.slice(g - 1848);
  const reads = pair(mut, 0.7);
  const bad = reads.map(r => analyzeRead(r, REF, { profileGapOpen: -10, repairMirror: false }));
  assert.ok(alerts(bad, mergeCalls(bad, REF)).some(a => a.skipped > 0));
  const good = reads.map(r => analyzeRead(r, REF));
  assert.ok(alerts(good, mergeCalls(good, REF)).every(a => a.skipped === 0));
});

test('SNV non retenu par JANUS : discordance non expliquée', () => {
  const g = 2100, alt = { A: 'C', C: 'G', G: 'T', T: 'A' }[seq[g]];
  const mut = wt.slice(0, g - 1849) + alt + wt.slice(g - 1848);
  const results = pair(mut, 1).map(r => analyzeRead(r, REF));
  // Variants retenus volontairement vides : le contrôle doit voir la discordance seul.
  assert.ok(alerts(results, []).every(a => a.unexplained >= 1));
});

test('réparation : même avec un gap peu coûteux, l’alignement reste sur l’allèle de référence', () => {
  const g = 2280 + motif.indexOf('GGG') + 2, mut = wt.slice(0, g - 1849) + wt.slice(g - 1848);
  for (const fraction of [0.7, 0.8]) {
    const results = pair(mut, fraction).map(r => analyzeRead(r, REF, { profileGapOpen: -10 }));
    const merged = mergeCalls(results, REF);
    assert.ok(alerts(results, merged).every(a => a.skipped === 0), `fraction ${fraction}`);
    assert.ok(merged.some(m => m.names.g === `NG_G.1:g.${g + 1}del`), `fraction ${fraction} : ${merged.map(m => m.names.g)}`);
    assert.ok(!merged.some(m => m.type === 'dup'), `fraction ${fraction} : duplication miroir`);
  }
});
