// Contrôle indépendant : alignement direct des bases appelées par le séquenceur (PBAS du fichier .ab1)
// sur la séquence du génome, sans déconvolution. Sert de double vérification de JANUS :
// discordances lecture / génome, et concordance du positionnement des bases avec celui de JANUS.

import { alignLocal } from './align.js';
import { revcomp } from './sequence.js';

const IUPAC = {
  A: 'A', C: 'C', G: 'G', T: 'T', R: 'AG', Y: 'CT', K: 'GT', M: 'AC', S: 'CG', W: 'AT',
  B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT',
};

/**
 * Alignement direct de la partie exploitable de la lecture sur le génome, autour de la région
 * localisée par JANUS (`result`).
 * Retourne, indexé par base de la lecture d'origine :
 *  - genome[k]  : base du génome en face (orientation de la référence) ou null (base lue en plus) ;
 *  - g[k]       : position génomique 1-based ou null ;
 *  - gapBefore[k] : bases du génome sans base lue en face, juste avant k dans l'orientation de la référence ;
 * ainsi que les discordances, l'identité et l'orientation la mieux soutenue.
 */
export function directCheck(result, ref) {
  const { read, trim: [t0, t1] } = result;
  const lo = Math.max(0, result.refStart - 1 - 300), hi = Math.min(ref.seq.length, result.refEnd + 300);
  const target = ref.seq.slice(lo, hi);
  const raw = read.bases.slice(t0, t1);
  const best = ['+', '-'].map(strand => {
    const query = strand === '+' ? raw : revcomp(raw);
    const aln = alignLocal(query.length, target.length, (i, j) => {
      const q = query[i], r = target[j];
      if (r === 'N' || q === 'N') return 0;
      if (q === r) return 3;
      return (IUPAC[q] || '').includes(r) ? 1 : -4;
    }, { gapOpen: -8, gapExtend: -2 });
    return { strand, query, aln };
  }).sort((a, b) => b.aln.score - a.aln.score)[0];

  const n = read.bases.length;
  const genome = new Array(n).fill(null), g = new Array(n).fill(null), gapBefore = new Array(n).fill('');
  // Indice dans la lecture d'origine d'une position de la requête orientée.
  const original = i => (best.strand === '+' ? t0 + i : t1 - 1 - i);
  const discrepancies = [];
  let pending = '', pendingPos = null, matches = 0, compared = 0;
  for (const o of best.aln.ops) {
    if (o.op === 'D') {
      if (!pending) pendingPos = lo + o.t + 1;
      pending += target[o.t];
      continue;
    }
    const k = original(o.q), q = best.query[o.q];
    if (pending) {
      gapBefore[k] = pending;
      discrepancies.push({ k, g: pendingPos, type: 'délétion', ref: pending, read: '' });
      pending = '';
    }
    if (o.op === 'I') { discrepancies.push({ k, g: lo + o.t, type: 'insertion', ref: '', read: q }); continue; }
    const r = target[o.t];
    genome[k] = r;
    g[k] = lo + o.t + 1;
    if (q === 'N' || r === 'N') continue;
    compared++;
    if (q === r) matches++;
    else discrepancies.push({ k, g: lo + o.t + 1, type: 'ACGT'.includes(q) ? 'substitution' : 'double pic', ref: r, read: q });
  }
  return {
    strand: best.strand,
    strandAgrees: best.strand === result.strand,
    identity: compared ? matches / compared : 0,
    compared,
    genome, g, gapBefore,
    discrepancies,
    mapping: compareMapping(result, g),
    skippedInMixed: skippedInMixedZone(result),
  };
}

/**
 * Bases du génome que l'alignement de JANUS a sautées juste avant une zone de doubles pics (sens de
 * lecture). Signature d'un alignement qui suit l'allèle muté au lieu de l'allèle de référence : l'indel
 * serait alors nommé à l'envers (duplication au lieu de délétion, ou l'inverse).
 */
function skippedInMixedZone(result) {
  const cols = result.strand === '+' ? result.cols : result.cols.slice().reverse();
  const out = [];
  for (let i = 0; i < cols.length; i++) {
    if (cols[i].q !== null) continue;
    let j = i;
    while (j < cols.length && cols[j].q === null) j++;
    const after = cols.slice(j, j + 20).filter(c => c.q !== null);
    const before = cols.slice(Math.max(0, i - 20), i).filter(c => c.q !== null);
    const mixed = cs => (cs.length ? cs.filter(c => c.mixed).length / cs.length : 0);
    if (after.length >= 10 && mixed(after) >= 0.4 && mixed(before) < 0.25) {
      const ts = cols.slice(i, j).map(c => result.winStart + c.t + 1);
      out.push({ gStart: Math.min(...ts), gEnd: Math.max(...ts), bases: j - i });
    }
    i = j;
  }
  return out;
}

/**
 * Concordance du positionnement : pour chaque base lue, position génomique selon JANUS et selon
 * l'alignement direct. Retourne les plages de bases où elles diffèrent.
 */
function compareMapping(result, gDirect) {
  const gJanus = new Map(), mixedAt = new Map();
  for (const c of result.cols) {
    if (c.q === null) continue;
    const k = result.profile[c.q].index;
    mixedAt.set(k, c.mixed);
    if (c.t !== null) gJanus.set(k, result.winStart + c.t + 1);
  }
  const ranges = [];
  let current = null, both = 0;
  for (let k = 0; k < gDirect.length; k++) {
    const a = gJanus.get(k), b = gDirect[k];
    if (a === undefined || b === null) continue;
    both++;
    if (a !== b) {
      if (current && k - current.lastK <= 3) { current.lastK = k; current.n++; current.gEnd = Math.max(current.gEnd, a); current.gStart = Math.min(current.gStart, a); }
      else { current = { firstK: k, lastK: k, n: 1, gStart: a, gEnd: a, shift: b - a }; ranges.push(current); }
    }
  }
  // Part de doubles pics autour de chaque plage : en zone mixte, les bases appelées par le séquenceur
  // suivent tantôt un allèle, tantôt l'autre, et ne permettent pas de juger le positionnement.
  for (const r of ranges) {
    let m = 0, t = 0;
    for (let k = r.firstK - 5; k <= r.lastK + 5; k++) if (mixedAt.has(k)) { t++; if (mixedAt.get(k)) m++; }
    r.mixed = t ? m / t : 0;
  }
  const disagree = ranges.reduce((s, r) => s + r.n, 0);
  return { compared: both, disagree, ranges };
}

/**
 * Un désaccord de positionnement est attendu en aval d'un indel hétérozygote : le séquenceur suit alors
 * l'allèle dominant, décalé de la taille de l'indel. Retourne le variant qui l'explique, ou null.
 */
export function explainShift(range, result, variants) {
  if (result.partial) return result.partial; // lecture entièrement décalée dès son début
  if (range.mixed >= 0.4) return { mixedZone: true };
  for (const v of variants) {
    if (v.type === 'snv' || v.genotype !== 'hétérozygote') continue;
    const net = Math.abs(v.ins.length - v.del.length);
    if (Math.abs(range.shift) !== net) continue;
    const downstream = result.strand === '+' ? range.gStart >= v.pos - 3 : range.gEnd <= v.pos + v.del.length + 3;
    if (downstream) return v;
  }
  return null;
}

/**
 * Explication d'une discordance par les variants retenus (tous lecteurs confondus) : variant à cette
 * position, ou zone décalée en aval d'un indel hétérozygote dans le sens de lecture.
 */
export function explainDiscrepancy(d, result, variants) {
  for (const v of variants) {
    const start = v.pos + 1, end = v.pos + Math.max(1, v.del.length);
    if (d.g >= start - 3 && d.g <= end + 3) return { kind: 'variant', variant: v };
  }
  for (const v of variants) {
    if (v.type === 'snv' || v.genotype !== 'hétérozygote') continue;
    const downstream = result.strand === '+' ? d.g > v.pos : d.g <= v.pos + v.del.length;
    if (downstream) return { kind: 'décalage', variant: v };
  }
  if (result.partial) return { kind: 'décalage', variant: null };
  return null;
}
