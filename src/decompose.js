// Déconvolution d'un chromatogramme Sanger en deux allèles guidée par la référence,
// puis appel des variants de chaque allèle (approche inspirée de Tracy, Rausch et al. 2020).
//
// 1. Profil de pics par base (primaire / secondaire / rapport).
// 2. Localisation et orientation de la lecture sur la référence génomique (k-mers).
// 3. Alignement du profil : une base de référence est « expliquée » si elle est le pic primaire
//    ou un pic secondaire significatif. L'allèle de référence étant présent tout du long chez un
//    hétérozygote, l'alignement le suit même après un indel.
// 4. Séparation : à chaque position, allèle 1 = base de référence si présente, allèle 2 = l'autre pic.
// 5. Alignement indépendant de chaque allèle sur la référence → SNV et indels exacts.

import { alignLocal, alignStrings, kmerIndex, locate, locateAmbiguous } from './align.js';
import { DNA, revcomp, complement } from './sequence.js';
import { peakProfile, reverseProfile, trimRange } from './peaks.js';
import { normalize, variantKey, hgvsNames } from './hgvs.js';

export const DEFAULTS = {
  ratio: 0.25,        // rapport secondaire/primaire minimal pour un pic mixte
  pad: 250,           // marge autour de la région localisée (nt)
  minVotes: 6,        // k-mers concordants minimaux pour localiser la lecture
  flank: 25,          // contexte (bases) pour juger un indel hétérozygote
  minMixedGain: 0.2,  // excès de positions mixtes après un indel hétérozygote
  edge: 15,           // pas d'appel à moins de `edge` bases d'une extrémité de l'alignement
  snrMin: 4,          // pic secondaire d'un SNV ≥ snrMin × bruit de fond local
  minSignal: 0.08,    // signal minimal (fraction du pic primaire) pour corriger l'allèle 2
  maxDeletion: 120,   // taille maximale des délétions recherchées au point de rupture (nt)
  maxInsertion: 40,   // taille maximale des insertions recherchées (nt)
  scanMinGain: 20,    // gain d'ajustement minimal pour retenir un indel trouvé par recherche du décalage
};

const h = (x, base) => x.h[DNA.indexOf(base)] ?? 0;

function referenceIndex(ref) {
  if (!ref._kmers) ref._kmers = kmerIndex(ref.seq, 12);
  return ref._kmers;
}

/** Analyse complète d'une lecture. */
export function analyzeRead(read, ref, options = {}) {
  const opt = { ...DEFAULTS, ...options };
  const fullProfile = peakProfile(read);
  const [t0, t1] = trimRange(fullProfile);
  const trimmed = fullProfile.slice(t0, t1);
  if (trimmed.length < 50) throw new Error(`${read.name} : moins de 50 bases exploitables`);

  // Orientation et localisation.
  const fwd = trimmed.map(x => x.primary).join('');
  const index = referenceIndex(ref);
  const hitF = locate(fwd, index), hitR = locate(revcomp(fwd), index);
  let reverse = (hitR?.votes || 0) > (hitF?.votes || 0);
  let hit = reverse ? hitR : hitF;
  if (!hit || hit.votes < opt.minVotes) {
    // Lecture mixte dès le début : les bases primaires alternent entre les deux allèles.
    if (!ref._kmers10) ref._kmers10 = kmerIndex(ref.seq, 10);
    const choices = trimmed.map(x => (x.ratio >= opt.ratio ? x.primary + x.secondary : x.primary));
    const rc = choices.slice().reverse().map(c => [...c].map(b => complement(b)).join(''));
    const aF = locateAmbiguous(choices, ref._kmers10), aR = locateAmbiguous(rc, ref._kmers10);
    reverse = (aR?.votes || 0) > (aF?.votes || 0);
    hit = reverse ? aR : aF;
  }
  if (!hit || hit.votes < opt.minVotes) throw new Error(`${read.name} : lecture non localisée sur ${ref.accession} (mauvaise référence ou qualité insuffisante)`);
  const prof = reverse ? reverseProfile(trimmed) : trimmed;
  const winStart = Math.max(0, hit.diagonal - opt.pad);
  const winEnd = Math.min(ref.seq.length, hit.diagonal + prof.length + opt.pad);
  const win = ref.seq.slice(winStart, winEnd);

  // Alignement du profil sur la référence.
  const profAln = alignLocal(prof.length, win.length, (i, j) => {
    const r = win[j], x = prof[i];
    if (r === 'N' || x.primary === 'N') return 0;
    // Toute base de référence visible compte autant que le pic primaire : chez un hétérozygote, les deux
    // allèles sont alors « parfaits » et l'alignement sans gap (allèle de référence) l'emporte, au lieu de
    // sauter vers l'allèle muté là où ses pics dominent.
    if (x.primary === r || h(x, r) >= opt.minSignal * Math.max(...x.h)) return 3;
    return -4;
  }, { gapOpen: -10, gapExtend: -2 });

  // Séparation des allèles.
  const cols = [];
  for (const o of profAln.ops) {
    if (o.op === 'D') { cols.push({ q: null, t: o.t, ref: win[o.t], a1: '-', a2: '-' }); continue; }
    const x = prof[o.q], r = o.op === 'M' ? win[o.t] : '-';
    const mixed = x.ratio >= opt.ratio;
    let a1 = x.primary, a2 = x.primary;
    if (r !== '-' && x.primary !== r && h(x, r) >= opt.minSignal * Math.max(...x.h)) {
      a1 = r; // référence présente, même faiblement : l'allèle 2 porte le pic dominant
    } else if (mixed) {
      if (x.secondary === r) { a1 = x.secondary; a2 = x.primary; } else { a1 = x.primary; a2 = x.secondary; }
    }
    cols.push({ q: o.q, t: o.op === 'M' ? o.t : null, ref: r, a1, a2, mixed });
  }
  const seqCols = cols.filter(c => c.q !== null);
  const allele1 = seqCols.map(c => c.a1).join('');
  // Après un indel hétérozygote, les pics de l'allèle minoritaire peuvent tomber sous le seuil :
  // on corrige l'allèle 2 vers la base attendue de son propre alignement si un signal la soutient.
  const allele2 = refineAllele(seqCols.map(c => c.a2).join(''), allele1, win, seqCols, prof, opt.minSignal);
  seqCols.forEach((c, i) => { c.a2 = allele2[i]; });

  const calls = new Map();
  for (const [which, allele] of [[1, allele1], [2, allele2]]) {
    for (const ev of alleleEvents(allele, win, opt.edge)) {
      const v = normalize(ref.seq, { pos: ev.pos + winStart, del: ev.del, ins: ev.ins });
      const key = variantKey(v);
      const call = calls.get(key) || { ...v, key, alleles: [], cols: [] };
      call.alleles.push(which);
      call.cols.push(...ev.q.map(i => seqCols[i]));
      call.qFirst = Math.min(call.qFirst ?? Infinity, ev.q[0] ?? Infinity);
      call.qLast = Math.max(call.qLast ?? -Infinity, ev.q.at(-1) ?? ev.qAfter);
      call.qAfter = Math.max(call.qAfter ?? -Infinity, ev.qAfter);
      calls.set(key, call);
    }
  }

  const model = makeModel({ ref, winStart, winEnd, seqCols, prof, reverse });
  const partial = refitIndels(calls, { ref, winStart, seqCols, prof, reverse, opt, model });

  const variants = [];
  for (const call of calls.values()) {
    const v = assess(call, seqCols, prof, opt, reverse);
    if (!v) continue;
    // Base de la lecture d'origine (1-based) au variant normalisé, pour la navigation.
    const anchorPos = v.type === 'dup' || v.type === 'ins' ? v.pos - 1 : v.pos;
    const col = seqCols.find(c => c.t !== null && c.t + winStart === anchorPos);
    const qi = col ? col.q : Math.min(Math.max(0, call.qFirst === Infinity ? call.qAfter : call.qFirst), prof.length - 1);
    v.readBase = prof[qi].index + 1;
    v.names = hgvsNames(v, ref);
    variants.push(v);
  }
  variants.sort((a, b) => a.pos - b.pos);

  const alignedRef = profAln.ops.filter(o => o.op !== 'I');
  return {
    name: read.name,
    read,
    strand: reverse ? '-' : '+',
    trim: [t0, t1],
    refStart: winStart + (alignedRef[0]?.t ?? 0) + 1,
    refEnd: winStart + (alignedRef.at(-1)?.t ?? 0) + 1,
    alignedBases: profAln.ops.filter(o => o.op !== 'D').length,
    usableBases: prof.length,
    mixedFraction: seqCols.filter(c => c.mixed).length / Math.max(1, seqCols.length),
    winStart,
    refSeq: ref.seq,
    cols,
    profile: prof,
    variants,
    // Indel dont le point de rupture précède le début de la lecture : taille connue, position non.
    partial: partial ? { ...partial, names: hgvsNames(partial, ref) } : null,
    /** Score d'ajustement de la lecture sous l'hypothèse d'un allèle 2 portant `vars`. */
    explain: vars => model.fit(vars, model.A, model.B),
  };
}

const ALLELE_SCORES = { match: 3, mismatch: -4, gapOpen: -8, gapExtend: -2 };

/**
 * Remplace, dans l'allèle 2, les bases discordantes avec son alignement quand la base attendue
 * a un signal mesurable et n'est pas déjà portée par l'allèle 1 (sinon on effacerait un SNV vrai).
 */
function refineAllele(allele2, allele1, win, seqCols, prof, minSignal) {
  const supported = (q, r) => {
    const x = prof[seqCols[q].q];
    return r !== allele1[q] && h(x, r) >= minSignal * Math.max(...x.h);
  };
  // Alignement « profil » de l'allèle 2 : une base faiblement soutenue compte presque comme un appariement.
  const aln = alignLocal(allele2.length, win.length, (q, t) => {
    const r = win[t];
    if (r === 'N' || allele2[q] === 'N') return 0;
    if (allele2[q] === r) return ALLELE_SCORES.match;
    return supported(q, r) ? 2 : ALLELE_SCORES.mismatch;
  }, ALLELE_SCORES);
  const out = allele2.split('');
  for (const o of aln.ops) {
    if (o.op !== 'M') continue;
    const r = win[o.t];
    if (r !== 'N' && out[o.q] !== r && supported(o.q, r)) out[o.q] = r;
  }
  return out.join('');
}

/** Événements (SNV, indels) d'un allèle aligné sur la fenêtre de référence. */
function alleleEvents(allele, win, edge) {
  const aln = alignStrings(allele, win, ALLELE_SCORES);
  const events = [];
  let ev = null, matchesSinceStart = 0;
  const ops = aln.ops;
  for (let k = 0; k < ops.length; k++) {
    const o = ops[k];
    const isMatch = o.op === 'M' && allele[o.q] === win[o.t];
    if (isMatch) {
      if (ev) { ev.qAfter = o.q; events.push(ev); ev = null; }
      matchesSinceStart++;
      continue;
    }
    if (o.op === 'M' && (allele[o.q] === 'N' || win[o.t] === 'N')) continue;
    if (!ev) ev = { pos: o.t, del: '', ins: '', q: [], startIndex: matchesSinceStart };
    if (o.op !== 'I') ev.del += win[o.t];
    if (o.op !== 'D') { ev.ins += allele[o.q]; ev.q.push(o.q); }
  }
  // Un événement non refermé par un appariement est en bout d'alignement : ignoré.
  const matchesTotal = matchesSinceStart;
  return events.filter(e => e.startIndex >= edge && matchesTotal - e.startIndex >= edge && !/N/.test(e.ins));
}

/** Séquence de la région [A, B) de `seq` après application de variants non chevauchants, ou null. */
function applyVariants(seq, vars, A, B) {
  let out = '', cur = A;
  for (const v of [...vars].sort((a, b) => a.pos - b.pos)) {
    if (v.pos < cur || v.pos + v.del.length > B) return null;
    out += seq.slice(cur, v.pos) + v.ins;
    cur = v.pos + v.del.length;
  }
  return out + seq.slice(cur, B);
}

/**
 * Modèle de mélange d'une lecture : l'allèle 1 suit la référence, l'allèle 2 porte des variants.
 * fit(vars, lo, hi) mesure (score négatif) l'écart entre les profils de pics observés et le mélange
 * attendu sur [lo, hi) de la référence. Les deux allèles partagent l'amorce : la k-ième base lue
 * de l'allèle 2 se compte depuis l'extrémité de la fenêtre située du côté de l'amorce.
 */
function makeModel({ ref, winStart, winEnd, seqCols, prof, reverse }) {
  const seq = ref.seq;
  const colOf = new Map();
  seqCols.forEach((c, i) => { if (c.t !== null) colOf.set(winStart + c.t, i); });
  const A = winStart, B = winEnd;
  // Fraction de l'allèle 2, estimée sur les positions où les deux allèles diffèrent.
  const ratios = seqCols.filter(c => c.a1 !== c.a2).map(c => {
    const x = prof[c.q], a = h(x, c.a2), b = h(x, c.a1);
    return a / Math.max(1, a + b);
  }).sort((a, b) => a - b);
  const f = ratios.length ? Math.min(0.9, Math.max(0.1, ratios[Math.floor(ratios.length / 2)])) : 0.5;

  const fit = (vars, lo, hi) => {
    const alt = applyVariants(seq, vars, A, B);
    if (alt === null) return -Infinity;
    let score = 0;
    for (let t = Math.max(A, lo); t < Math.min(B, hi); t++) {
      const i = colOf.get(t);
      if (i === undefined) continue;
      const k = reverse ? alt.length - (B - t) : t - A;
      if (k < 0 || k >= alt.length) continue;
      const x = prof[seqCols[i].q], total = x.h.reduce((u, v) => u + v, 0) || 1;
      for (let c = 0; c < 4; c++) {
        const e = (DNA[c] === seq[t] ? 1 - f : 0) + (DNA[c] === alt[k] ? f : 0);
        // Écart logarithmique : un pic net mais plus petit qu'attendu pèse plus qu'en écart quadratique.
        score -= Math.abs(Math.log((x.h[c] / total + 0.05) / (e + 0.05)));
      }
    }
    return score;
  };
  return { fit, colOf, f, A, B };
}

/**
 * Recherche directe du décalage : aux points où le signal devient mixte (ou dès le début de la
 * lecture), on teste toutes les délétions et insertions jusqu'à maxDeletion / maxInsertion nt et
 * l'on garde celle qui explique le mieux les profils en aval. Indispensable pour les grands indels,
 * que la séparation base par base fragmente en plusieurs petits appels.
 */
function scanShift({ ref, winStart, seqCols, reverse, opt, model, inferInsert }) {
  const seq = ref.seq;
  const cols = seqCols.filter(c => c.t !== null);
  if (reverse) cols.reverse(); // sens de lecture
  const n = cols.length;
  if (n < 60) return null;
  const sums = [0];
  for (const c of cols) sums.push(sums.at(-1) + (c.mixed ? 1 : 0));
  const mean = (a, b) => { a = Math.max(0, a); b = Math.min(n, b); return b > a ? (sums[b] - sums[a]) / (b - a) : 0; };

  // Points de rupture candidats : saut de la proportion de pics mixtes (au plus 3, espacés).
  const W = 25, jumps = [];
  for (let j = 0; j + W <= n; j++) {
    const post = mean(j, j + W), pre = j >= 10 ? mean(j - W, j) : 0;
    // Le signal doit être propre avant la rupture : dans une zone déjà décalée, le gain d'un second
    // indel serait mesuré contre un modèle faux (sans le premier).
    if (post >= 0.35 && pre <= 0.3 && post - pre >= 0.25) jumps.push({ j, jump: post - pre, atStart: j < 10 });
  }
  jumps.sort((a, b) => b.jump - a.jump);
  const points = [];
  for (const p of jumps) {
    if (points.every(q => Math.abs(q.j - p.j) > 20)) points.push(p);
    if (points.length === 3) break;
  }

  let best = null;
  for (const point of points) {
    const g = winStart + cols[point.atStart ? 0 : point.j].t; // première position mixte (référence)
    const [lo, hi] = reverse ? [g - 130, g + 21] : [g - 20, g + 131];
    const base = model.fit([], lo, hi);
    const tried = new Set();
    const consider = v => {
      if (!v || v.ins === null || v.pos < model.A + 1 || v.pos + v.del.length > model.B - 1 || /N/.test(v.ins)) return;
      const nv = normalize(seq, v);
      const key = variantKey(nv);
      if (tried.has(key)) return;
      tried.add(key);
      const gain = model.fit([nv], lo, hi) - base;
      if (!best || gain > best.gain) best = { variant: nv, gain, point, g };
    };
    for (let d = 1; d <= opt.maxDeletion; d++) {
      // Lecture Forward : la rupture est au début de la délétion ; Reverse : à sa fin.
      const start = reverse ? g - d + 1 : g;
      for (let p = start - 8; p <= start + 8; p++) consider({ pos: p, del: seq.slice(p, p + d), ins: '' });
      // Rupture avant la première base lue : la délétion peut chevaucher le début de la lecture.
      // L'alignement a pu suivre l'un ou l'autre allèle (aucune zone propre pour trancher) : on teste
      // aussi le cas symétrique, où l'allèle de référence apparaît comme une duplication de d nt.
      if (point.atStart) {
        for (let k = -8; k <= d; k++) {
          const p = reverse ? g - d + 1 + k : g - k;
          consider({ pos: p, del: seq.slice(p, p + d), ins: '' });
        }
        for (let p = (reverse ? g + d : g) - 8; p <= (reverse ? g + d : g) + 8; p++) consider({ pos: p, del: '', ins: seq.slice(p - d, p) });
      }
    }
    for (let d = 1; d <= opt.maxInsertion; d++) {
      for (let p = g - 8; p <= g + 8; p++) {
        consider({ pos: p, del: '', ins: seq.slice(p - d, p) });
        const ins = inferInsert(p, d);
        if (ins) consider({ pos: p, del: '', ins });
      }
    }
  }
  if (!best || best.gain < opt.scanMinGain) return null;
  // Rupture non observée : l'indel commence avant la première base lue (taille sûre, position non).
  const v = best.variant, first = winStart + cols[0].t;
  best.partial = best.point.atStart && (reverse ? v.pos + v.del.length >= first - 2 : v.pos <= first + 2);
  return best;
}

/**
 * Ré-ajustement des indels hétérozygotes de l'allèle 2. Près du point de rupture, la séparation
 * base par base peut décaler l'indel de quelques nucléotides, le fragmenter (indel + SNV voisins)
 * ou créer de petits couples délétion/insertion parasites en aval. On commence par la recherche
 * directe du décalage (grands indels), puis, pour chaque groupe d'appels restant, dans le sens de
 * lecture, on teste : le groupe tel quel, son absence et tous les placements d'un indel unique de
 * même bilan ; on garde l'allèle 2 prédit qui reproduit le mieux les profils de pics mesurés.
 * Retourne l'indel « partiel » éventuel (rupture avant le début de la lecture).
 */
function refitIndels(calls, { ref, winStart, seqCols, prof, reverse, opt, model }) {
  const { fit, colOf, A, B } = model;
  const nearestCol = g => {
    for (let d = 0; d < 30; d++) {
      if (colOf.has(g + d)) return colOf.get(g + d);
      if (colOf.has(g - d)) return colOf.get(g - d);
    }
    return 0;
  };
  const seq = ref.seq;
  // Bases d'une insertion non dupliquée de longueur L en p, lues dans le signal.
  const inferInsert = (p, L) => {
    let s = '';
    for (let i = 0; i < L; i++) {
      const col = colOf.get(reverse ? p + i - L : p + i);
      if (col === undefined) return null;
      const x = prof[seqCols[col].q], a1 = seqCols[col].a1, top = Math.max(...x.h);
      let best = a1;
      for (const b of DNA) if (b !== a1 && h(x, b) > (best === a1 ? opt.minSignal * top : h(x, best))) best = b;
      s += best;
    }
    return s;
  };
  const addCall = v => {
    const key = variantKey(v);
    const call = { ...v, key, alleles: [2], cols: [], refit: true, qFirst: nearestCol(v.pos), qAfter: nearestCol(v.pos + v.del.length) };
    calls.set(key, call);
    return call;
  };

  const kept = [];
  let partial = null;
  const scan = scanShift({ ref, winStart, seqCols, reverse, opt, model, inferInsert });
  if (scan) {
    const v = scan.variant;
    // Les appels de l'allèle 2 qui chevauchent l'indel trouvé en sont des fragments.
    for (const c of [...calls.values()]) {
      if (c.alleles.includes(1)) continue;
      if (c.pos <= v.pos + v.del.length + 3 && c.pos + c.del.length >= v.pos - 3) calls.delete(c.key);
    }
    if (scan.partial) {
      partial = { ...v, gain: scan.gain };
      // Toute la lecture est décalée : les autres appels de l'allèle 2 ne sont pas interprétables.
      for (const c of [...calls.values()]) if (!c.alleles.includes(1)) calls.delete(c.key);
      kept.push(v);
    } else {
      kept.push(addCall(v));
    }
  }

  const het = [...calls.values()].filter(c => !c.alleles.includes(1) && !kept.includes(c)).sort((a, b) => (reverse ? b.pos - a.pos : a.pos - b.pos));
  let pending = het.slice();
  while (pending.length) {
    const seed = pending.find(c => c.type !== 'snv');
    if (!seed) break;
    const cluster = pending.filter(c => c.pos >= seed.pos - 15 && c.pos <= seed.pos + seed.del.length + 15);
    pending = pending.filter(c => !cluster.includes(c));
    // Les groupes déjà traités (en amont dans le sens de lecture) sont fixés ; les suivants, pas encore
    // jugés, sont ignorés et la zone évaluée s'arrête avant eux.
    const others = kept;
    const lo = Math.min(...cluster.map(c => c.pos)), hi = Math.max(...cluster.map(c => c.pos + c.del.length));
    const delta = cluster.reduce((sum, c) => sum + c.ins.length - c.del.length, 0);
    let from = lo - 60, to = hi + 60;
    for (const c of pending) {
      if (c.type === 'snv') continue; // souvent des artefacts en aval de l'indel : ne bornent pas la zone
      if (!reverse && c.pos >= hi) to = Math.max(hi + 20, Math.min(to, c.pos));
      if (reverse && c.pos + c.del.length <= lo) from = Math.min(lo - 20, Math.max(from, c.pos + c.del.length));
    }

    let best = { vars: cluster, score: fit([...others, ...cluster], from, to), original: true };
    const tryVars = (vars, simpler) => {
      const score = fit([...others, ...vars], from, to);
      if (score > best.score + 1e-9 || (simpler && best.original && score >= best.score - 1e-9)) best = { vars, score };
    };
    const tryCandidate = v => {
      if (!v || (v.ins && v.ins.includes('N'))) return;
      tryVars([normalize(seq, v)], cluster.length > 1);
    };
    tryVars([], true); // le groupe peut n'être que du bruit
    for (let p = Math.max(A + 1, lo - 8); p <= Math.min(B - 1, hi + 8); p++) {
      if (delta < 0) tryCandidate({ pos: p, del: seq.slice(p, p - delta), ins: '' });
      else if (delta > 0) {
        tryCandidate({ pos: p, del: '', ins: seq.slice(p - delta, p) });
        const ins = inferInsert(p, delta);
        if (ins) tryCandidate({ pos: p, del: '', ins });
      }
    }
    if (best.original) { kept.push(...cluster); continue; }
    for (const c of cluster) calls.delete(c.key);
    for (const v of best.vars) kept.push(addCall(v));
  }

  // SNV de l'allèle 2 : conservés seulement s'ils améliorent l'ajustement, compte tenu des indels
  // retenus (en aval d'un indel, les discordances de séparation imitent des SNV).
  for (const c of [...calls.values()]) {
    if (c.type !== 'snv' || c.alleles.includes(1)) continue;
    const others = kept.filter(k => k !== c);
    const gain = fit([...others, c], c.pos - 10, c.pos + 11) - fit(others, c.pos - 10, c.pos + 11);
    if (gain < 1) calls.delete(c.key);
  }
  return partial;
}

/** Génotype, fraction allélique et contrôles de cohérence d'un appel. */
function assess(call, seqCols, prof, opt, reverse) {
  const both = call.alleles.includes(1) && call.alleles.includes(2);
  // L'allèle 1 suit la référence dès qu'elle est visible : un indel propre à l'allèle 1 est un artefact
  // de séparation (les indels hétérozygotes composites ne sont pas résolus par cette approche).
  if (!both && !call.alleles.includes(2) && call.type !== 'snv') return null;
  const genotype = both ? 'homozygote' : 'hétérozygote';
  const out = {
    pos: call.pos, del: call.del, ins: call.ins, type: call.type, key: call.key,
    genotype, alleleFraction: null, evidence: '', flags: [],
  };

  if (call.type === 'snv') {
    const col = call.cols[0];
    if (!col) return null;
    const x = prof[col.q];
    const alt = h(x, call.ins), wt = h(x, call.del);
    if (both) {
      if (col.mixed) out.flags.push('pic secondaire présent');
      out.alleleFraction = alt / Math.max(1, alt + wt);
    } else {
      // Un SNV hétérozygote doit correspondre à un vrai double pic, nettement au-dessus du bruit de fond
      // local (estimé par le 3e canal des positions voisines).
      if (!col.mixed) {
        // Référence à peine visible sous un pic alternatif dominant, dans une région non mixte
        // (hors zone décalée par un indel hétérozygote) : homozygote probable.
        if (x.primary !== call.ins) return null;
        const i = seqCols.indexOf(col);
        const around = seqCols.slice(Math.max(0, i - 10), i + 11);
        if (around.filter(c => c.mixed).length > 0.2 * around.length) return null;
        out.genotype = 'homozygote';
        out.flags.push('signal résiduel de la référence');
        out.alleleFraction = alt / Math.max(1, alt + wt);
        out.evidence = `rapport des pics ${x.ratio.toFixed(2)}`;
        return out;
      }
      const near = prof.slice(Math.max(0, col.q - 10), col.q + 11).map(y => y.h.slice().sort((a, b) => b - a)[2]).sort((a, b) => a - b);
      const background = near[Math.floor(near.length / 2)] || 0;
      if (Math.min(alt, wt) < opt.snrMin * background) return null;
      out.alleleFraction = alt / Math.max(1, alt + wt);
      if (!call.alleles.includes(2)) out.flags.push('aucun allèle de référence');
      // Pic secondaire de la même base que le voisin : possible traînée du pic adjacent.
      const i = seqCols.indexOf(col);
      if ([seqCols[i - 1], seqCols[i + 1]].some(c => c?.a1 === call.ins) && x.ratio < 0.6) out.flags.push('même base que le pic voisin (traînée possible)');
    }
    out.evidence = `rapport des pics ${x.ratio.toFixed(2)}`;
    return out;
  }

  // Indels et délétions-insertions.
  const anchor = call.qFirst === Infinity ? call.qAfter : call.qFirst;
  const window = (a, b) => {
    const cs = seqCols.slice(Math.max(0, a), Math.max(0, b));
    return cs.length ? cs.filter(c => c.mixed).length / cs.length : 0;
  };
  const span = (a, b) => Math.max(0, Math.min(seqCols.length, b) - Math.max(0, a));
  // « Avant » / « après » dans le sens de lecture : une lecture Reverse parcourt la référence à l'envers.
  const up = [anchor - opt.flank, anchor - 2], down = [call.qAfter + 2, call.qAfter + 2 + opt.flank];
  const [pre, post] = reverse ? [down, up] : [up, down];
  const before = window(...pre);
  const after = window(...post);
  out.mixedBefore = before;
  out.mixedAfter = after;
  if (both) {
    out.evidence = `lecture homogène (${Math.round(100 * after)} % de pics mixtes après)`;
    if (after > 0.3) out.flags.push('signal mixte inattendu après le variant');
    return out;
  }
  // Un indel hétérozygote décale un allèle : le signal devient mixte en aval, sur une longueur suffisante.
  if (span(...post) < 0.8 * opt.flank) return null;
  if (after - before < opt.minMixedGain || after < 0.3) return null;
  const frac = [];
  for (const c of seqCols.slice(Math.max(0, post[0]), Math.max(0, post[1]))) {
    if (c.a1 === c.a2) continue;
    const x = prof[c.q], a = h(x, c.a2), b = h(x, c.a1);
    frac.push(a / Math.max(1, a + b));
  }
  out.alleleFraction = frac.length ? frac.reduce((s, v) => s + v, 0) / frac.length : null;
  out.evidence = `pics mixtes : ${Math.round(100 * before)} % avant → ${Math.round(100 * after)} % après (sens de lecture)`;
  return out;
}

/**
 * Forme « miroir » d'une délétion ou d'une duplication : quand une lecture commence dans la zone
 * décalée, son alignement peut suivre l'allèle muté ; l'allèle de référence y apparaît alors comme
 * l'indel inverse (une délétion de X se lit comme une duplication de X, et réciproquement).
 */
function mirrorOf(v, seq) {
  if (v.type === 'del') return normalize(seq, { pos: v.pos + v.del.length, del: '', ins: v.del });
  if (v.type === 'dup') return normalize(seq, { pos: v.pos - v.ins.length, del: v.ins, ins: '' });
  return null;
}

/**
 * Le signal d'une lecture qui n'a pas appelé le variant est-il mieux expliqué avec lui ?
 * Retourne le gain d'ajustement (forme directe ou miroir), ou null si la lecture n'est pas informative.
 */
function modelSupport(r, v, seq) {
  const end = v.pos + v.del.length;
  if (v.type === 'snv') {
    if (r.refStart + 5 > v.pos + 1 || r.refEnd - 5 < v.pos + 1) return null;
  } else if (r.strand === '+' ? r.refEnd < end + 20 : r.refStart > v.pos - 20) {
    return null; // un indel ne se voit qu'en aval, dans le sens de lecture
  }
  const own = r.variants.filter(x => x.genotype === 'hétérozygote' && (x.pos + x.del.length < v.pos - 3 || x.pos > end + 3));
  if (r.partial && v.type === 'snv') own.push(r.partial); // lecture entièrement décalée : le décalage fait partie du modèle
  const base = r.explain(own);
  const gains = [v, mirrorOf(v, seq)].filter(Boolean).map(x => r.explain([...own, x]) - base);
  return Math.max(...gains);
}

/**
 * Regroupe les appels de plusieurs lectures (ex. Forward + Reverse).
 * Un variant est « confirmé F/R » s'il est appelé sur les deux brins, ou si l'autre brin, sans
 * l'avoir appelé (lecture commençant dans la zone décalée, signal faible), est nettement mieux
 * expliqué avec lui.
 */
export function mergeCalls(results, ref = null) {
  const merged = new Map();
  for (const r of results) {
    for (const v of r.variants) {
      const m = merged.get(v.key) || { ...v, reads: [], strands: new Set(), supportedBy: [] };
      m.reads.push({ name: r.name, strand: r.strand, genotype: v.genotype, alleleFraction: v.alleleFraction, evidence: v.evidence, flags: v.flags, readBase: v.readBase, result: r });
      m.strands.add(r.strand);
      merged.set(v.key, m);
    }
  }
  const seq = ref?.seq ?? results[0]?.refSeq;

  // Placements concurrents d'un même indel (même bilan de longueur, à moins de 10 nt) : on garde celui
  // qui explique le mieux l'ensemble des lectures informatives.
  const net = v => v.ins.length - v.del.length;
  const indels = [...merged.values()].filter(m => m.type !== 'snv');
  for (const a of indels) {
    if (!merged.has(a.key)) continue;
    const rivals = indels.filter(b => merged.has(b.key) && net(b) === net(a) && Math.abs(b.pos - a.pos) <= 10);
    if (rivals.length < 2 || !seq) continue;
    const total = v => results.reduce((sum, r) => {
      if (!r.explain || modelSupport(r, v, seq) === null) return sum;
      const own = r.variants.filter(x => x.genotype === 'hétérozygote' && !rivals.some(rv => rv.key === x.key));
      return sum + Math.max(...[v, mirrorOf(v, seq)].filter(Boolean).map(x => r.explain([...own, x])));
    }, 0);
    const scored = rivals.map(v => ({ v, score: total(v) })).sort((x, y) => y.score - x.score);
    const winner = scored[0].v;
    for (const { v } of scored.slice(1)) {
      // Les lectures qui avaient choisi l'autre placement soutiennent le placement retenu.
      for (const x of v.reads) {
        if (!winner.reads.some(w => w.result === x.result)) {
          winner.reads.push({ ...x, flags: [...x.flags, `placement initial ${x.result.variants.find(y => y.key === v.key)?.names.g.split(':')[1] ?? ''}`] });
          winner.strands.add(x.strand);
        }
      }
      merged.delete(v.key);
    }
  }

  const out = [];
  for (const m of merged.values()) {
    const missing = [];
    for (const r of results) {
      if (m.reads.some(x => x.result === r)) continue;
      const gain = seq && r.explain ? modelSupport(r, m, seq) : null;
      if (gain === null) continue;
      const strong = m.type === 'snv' ? 2 : 20;
      if (gain >= strong) { m.supportedBy.push({ name: r.name, strand: r.strand, gain }); m.strands.add(r.strand); }
      else if (gain <= -strong / 2 || m.type !== 'snv') missing.push(r);
    }
    const genotypes = new Set(m.reads.map(x => x.genotype));
    out.push({
      ...m,
      genotype: genotypes.size === 1 ? [...genotypes][0] : 'discordant',
      status: m.strands.size === 2 ? 'confirmé F/R' : missing.length ? 'non retrouvé sur l’autre lecture' : 'une seule lecture',
      missing: missing.map(r => r.name),
    });
  }
  return out.sort((a, b) => a.pos - b.pos);
}
