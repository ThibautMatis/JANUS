// Simulation de chromatogrammes Sanger (mélange d'allèles) et écriture ABIF, pour les tests.

import { DNA } from '../src/sequence.js';

export function rng(seed = 1) {
  let s = seed >>> 0 || 1;
  const next = () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; };
  next.normal = () => Math.sqrt(-2 * Math.log(next() + 1e-12)) * Math.cos(2 * Math.PI * next());
  return next;
}

export function randomSeq(n, rand) {
  let s = '';
  for (let i = 0; i < n; i++) s += DNA[Math.floor(rand() * 4)];
  return s;
}

/**
 * Chromatogramme d'un mélange : alleles[k] lu avec la proportion fractions[k].
 * Chaque position i porte un pic par allèle (même mobilité : fragments de même longueur).
 */
export function simulate(alleles, fractions, { spacing = 12, sigma = 2.6, noise = 8, heightVar = 0.35, seed = 7, length } = {}) {
  const rand = rng(seed);
  const n = length ?? Math.min(...alleles.map(a => a.length));
  const size = (n + 2) * spacing;
  const traces = Object.fromEntries(DNA.split('').map(b => [b, new Float32Array(size)]));
  const positions = [];
  for (let i = 0; i < n; i++) {
    const center = (i + 1) * spacing + Math.round(rand() * 2 - 1);
    positions.push(center);
    // Montée du signal au début, déclin en fin de lecture.
    const env = Math.min(1, (i + 5) / 25) * Math.min(1, (n - i) / (0.15 * n));
    alleles.forEach((a, k) => {
      const base = a[i];
      if (!DNA.includes(base)) return;
      const height = 1200 * fractions[k] * env * Math.exp(heightVar * rand.normal());
      const t = traces[base];
      for (let x = Math.max(0, Math.round(center - 4 * sigma)); x <= Math.min(size - 1, Math.round(center + 4 * sigma)); x++) {
        t[x] += height * Math.exp(-((x - center) ** 2) / (2 * sigma * sigma));
      }
    });
  }
  for (const b of DNA) for (let x = 0; x < size; x++) traces[b][x] = Math.round(traces[b][x] + Math.abs(noise * rand.normal()));
  // Basecaller naïf : canal maximal au pic.
  const bases = positions.map(p => DNA.split('').reduce((best, b) => (traces[b][p] > traces[best][p] ? b : best), 'A')).join('');
  return { bases, positions, traces, quality: positions.map(() => 30), order: 'GATC', sample: 'synthetique' };
}

/** Écrit un fichier ABIF minimal (tags lus par JANUS). */
export function writeABIF({ bases, positions, traces, quality, order = 'GATC', sample = '' }) {
  const entries = [];
  const enc = s => Uint8Array.from(s, c => c.charCodeAt(0));
  const shorts = arr => { const b = new Uint8Array(arr.length * 2), v = new DataView(b.buffer); arr.forEach((x, i) => v.setInt16(i * 2, x, false)); return b; };
  for (let i = 0; i < 4; i++) entries.push({ name: 'DATA', num: 9 + i, type: 4, size: 2, n: traces[order[i]].length, data: shorts(Array.from(traces[order[i]])) });
  entries.push({ name: 'FWO_', num: 1, type: 2, size: 1, n: 4, data: enc(order) });
  entries.push({ name: 'PBAS', num: 2, type: 2, size: 1, n: bases.length, data: enc(bases) });
  entries.push({ name: 'PLOC', num: 2, type: 4, size: 2, n: positions.length, data: shorts(positions) });
  if (quality?.length) entries.push({ name: 'PCON', num: 2, type: 2, size: 1, n: quality.length, data: Uint8Array.from(quality) });
  if (sample) entries.push({ name: 'SMPL', num: 1, type: 18, size: 1, n: sample.length + 1, data: Uint8Array.from([sample.length, ...enc(sample)]) });

  let cursor = 128;
  for (const e of entries) if (e.data.length > 4) { e.offset = cursor; cursor += e.data.length; }
  const dirOffset = cursor;
  const total = dirOffset + entries.length * 28;
  const buf = new ArrayBuffer(total), view = new DataView(buf), u8 = new Uint8Array(buf);
  u8.set(enc('ABIF'), 0);
  view.setUint16(4, 101, false);
  u8.set(enc('tdir'), 6);
  view.setUint32(10, 1, false);
  view.setUint16(14, 1023, false);
  view.setUint16(16, 28, false);
  view.setUint32(18, entries.length, false);
  view.setUint32(22, entries.length * 28, false);
  view.setUint32(26, dirOffset, false);
  entries.forEach((e, i) => {
    const p = dirOffset + i * 28;
    u8.set(enc(e.name), p);
    view.setUint32(p + 4, e.num, false);
    view.setUint16(p + 8, e.type, false);
    view.setUint16(p + 10, e.size, false);
    view.setUint32(p + 12, e.n, false);
    view.setUint32(p + 16, e.data.length, false);
    if (e.data.length > 4) { view.setUint32(p + 20, e.offset, false); u8.set(e.data, e.offset); } else u8.set(e.data, p + 20);
  });
  return buf;
}
