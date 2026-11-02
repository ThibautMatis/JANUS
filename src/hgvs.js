// Normalisation des variants et nomenclature HGVS (g. sur la référence génomique, c. sur le transcrit).
// Coordonnées internes : 0-based sur la référence génomique, intervalle [pos, pos + del.length).

/**
 * Normalise un variant {pos, del, ins} : retire le contexte commun puis décale les indels
 * au plus 3' (règle HGVS) dans l'orientation de la référence.
 */
export function normalize(seq, { pos, del, ins }) {
  // Préfixe puis suffixe communs.
  while (del.length && ins.length && del[0] === ins[0]) { del = del.slice(1); ins = ins.slice(1); pos++; }
  while (del.length && ins.length && del.at(-1) === ins.at(-1)) { del = del.slice(0, -1); ins = ins.slice(0, -1); }
  if (!del.length && ins.length) {
    while (pos < seq.length && seq[pos] === ins[0]) { ins = ins.slice(1) + seq[pos]; pos++; }
  } else if (del.length && !ins.length) {
    while (pos + del.length < seq.length && seq[pos] === seq[pos + del.length]) { pos++; }
    del = seq.slice(pos, pos + del.length);
  }
  return { pos, del, ins, type: variantType(seq, pos, del, ins) };
}

function variantType(seq, pos, del, ins) {
  if (del.length === 1 && ins.length === 1) return 'snv';
  if (!del.length) return seq.slice(pos - ins.length, pos) === ins ? 'dup' : 'ins';
  if (!ins.length) return 'del';
  return 'delins';
}

/** Clé unique d'un variant normalisé. */
export const variantKey = v => `${v.pos}:${v.del}>${v.ins}`;

/**
 * Convertisseur g. → c. pour un transcrit : exons [[début, fin]] 1-based sur la référence (brin +),
 * cdsStart / cdsEnd en coordonnées du transcrit (1-based).
 */
export function transcriptMapper({ exons, cdsStart, cdsEnd }) {
  const ex = exons.slice().sort((a, b) => a[0] - b[0]);
  const tStarts = [];
  let acc = 0;
  for (const [s, e] of ex) { tStarts.push(acc + 1); acc += e - s + 1; }
  const tLength = acc;

  // Position transcrit (1-based) → libellé c. sans décalage intronique.
  const cLabel = t => {
    if (t < cdsStart) return '-' + (cdsStart - t);
    if (t > cdsEnd) return '*' + (t - cdsEnd);
    return String(t - cdsStart + 1);
  };
  // Position génomique 1-based → transcrit (1-based) ou null si intronique.
  const exonic = g => {
    for (let k = 0; k < ex.length; k++) if (g >= ex[k][0] && g <= ex[k][1]) return tStarts[k] + g - ex[k][0];
    return null;
  };

  /** Libellé c. d'une position génomique 1-based (introns : c.N+d / c.N-d, hors transcrit : prolongation). */
  function c(g) {
    const t = exonic(g);
    if (t !== null) return cLabel(t);
    if (g < ex[0][0]) return cLabel(1 - (ex[0][0] - g)); // en amont du transcrit
    if (g > ex.at(-1)[1]) return cLabel(tLength + (g - ex.at(-1)[1]));
    for (let k = 0; k + 1 < ex.length; k++) {
      const a = ex[k][1], b = ex[k + 1][0];
      if (g > a && g < b) {
        const dUp = g - a, dDown = b - g;
        // À égale distance, HGVS rattache à l'exon en amont (+).
        return dUp <= dDown ? cLabel(tStarts[k] + a - ex[k][0]) + '+' + dUp : cLabel(tStarts[k + 1]) + '-' + dDown;
      }
    }
    return '?';
  }
  /**
   * Inverse de c() : position génomique 1-based d'un libellé c. (« 1813 », « -25 », « *40 », « 6841+2 »,
   * « 212-25 »), ou null si le libellé est invalide.
   */
  function g(label) {
    const m = String(label).trim().match(/^([-*]?)(\d+)(?:([+-])(\d+))?$/);
    if (!m || (m[1] === '' && m[2] === '0')) return null;
    const n = Number(m[2]);
    const t = m[1] === '-' ? cdsStart - n : m[1] === '*' ? cdsEnd + n : cdsStart + n - 1;
    let pos = null;
    if (t < 1) pos = ex[0][0] - (1 - t);
    else if (t > tLength) pos = ex.at(-1)[1] + (t - tLength);
    else for (let k = 0; k < ex.length; k++) if (t >= tStarts[k] && t <= tStarts[k] + ex[k][1] - ex[k][0]) pos = ex[k][0] + t - tStarts[k];
    if (pos === null) return null;
    if (m[3]) pos += (m[3] === '+' ? 1 : -1) * Number(m[4]);
    return pos;
  }
  /** Exon ou intron (numérotés dans l'ordre du transcrit) contenant la position génomique 1-based g. */
  function region(g) {
    for (let k = 0; k < ex.length; k++) {
      if (g >= ex[k][0] && g <= ex[k][1]) return { kind: 'exon', n: k + 1 };
      if (k + 1 < ex.length && g > ex[k][1] && g < ex[k + 1][0]) return { kind: 'intron', n: k + 1 };
    }
    return { kind: g < ex[0][0] ? 'amont' : 'aval', n: null };
  }
  return { c, g, exonic, region, tLength, exons: ex };
}

/** Description HGVS (sans préfixe) d'un variant normalisé ; `coord` convertit une position 1-based. */
export function describe(v, coord) {
  const p = v.pos + 1; // 1-based
  const range = (a, b) => (a === b ? coord(a) : `${coord(a)}_${coord(b)}`);
  switch (v.type) {
    case 'snv': return `${coord(p)}${v.del}>${v.ins}`;
    case 'del': return `${range(p, p + v.del.length - 1)}del`;
    case 'dup': return `${range(p - v.ins.length, p - 1)}dup`;
    case 'ins': return `${coord(p - 1)}_${coord(p)}ins${v.ins}`;
    default: return `${range(p, p + v.del.length - 1)}delins${v.ins}`;
  }
}

/** Noms HGVS g. et c. d'un variant pour une référence { accession, transcript?, mapper? }. */
export function hgvsNames(v, ref) {
  const g = `${ref.accession}:g.${describe(v, String)}`;
  const c = ref.mapper && ref.transcript ? `${ref.transcript.id}:c.${describe(v, ref.mapper.c)}` : null;
  return { g, c };
}

const regionLabel = r => (r.n ? `${r.kind} ${r.n}` : r.kind === 'amont' ? 'en amont du gène' : 'en aval du gène');

/**
 * Localisation d'un variant normalisé sur le transcrit : « exon 10 », « intron 9 »,
 * ou « intron 1 – exon 2 » s'il chevauche une jonction. null sans transcrit.
 */
export function variantLocation(v, mapper) {
  if (!mapper) return null;
  const first = v.del.length ? v.pos + 1 : v.pos;          // insertion : entre v.pos et v.pos + 1
  const last = v.del.length ? v.pos + v.del.length : v.pos + 1;
  const a = mapper.region(first), b = mapper.region(last);
  // Insertion à la limite exon/intron : on la rattache à l'exon.
  if (!v.del.length && a.kind !== b.kind) return regionLabel(a.kind === 'exon' ? a : b);
  const la = regionLabel(a), lb = regionLabel(b);
  return la === lb ? la : `${la} – ${lb}`;
}

/**
 * Lecture d'une description HGVS c. tapée par l'utilisateur (« c.1813dup », « NM_000059.4:c.1813dupA »,
 * « c.2311C>T », « c.212-25_239del », « c.68_69insAT », « c.100_102delinsG ») en variant normalisé
 * {pos, del, ins, type} sur la référence génomique. Lève une erreur explicite si la description est
 * invalide ou incohérente avec la référence (base de référence différente, position hors transcrit).
 */
export function parseVariantC(text, ref) {
  if (!ref.mapper) throw new Error('aucun transcrit chargé : impossible de lire une position c.');
  let s = String(text).trim().replace(/\s+/g, '');
  let transcript = null;
  const prefixed = s.match(/^([A-Z]{2}_[A-Z0-9]+(?:\.\d+)?)(?:\([^)]*\))?:(.*)$/i);
  if (prefixed) { transcript = prefixed[1].toUpperCase(); s = prefixed[2]; }
  if (!/^c\./i.test(s)) throw new Error('la description doit commencer par « c. » (ex. c.1813dup)');
  s = s.slice(2).replace(/^\((.*)\)$/, '$1');
  const m = s.match(/^([-*]?\d+(?:[+-]\d+)?)(?:_([-*]?\d+(?:[+-]\d+)?))?(.*)$/);
  if (!m) throw new Error(`position illisible dans « ${text} »`);
  const a = ref.mapper.g(m[1]), b = m[2] ? ref.mapper.g(m[2]) : a;
  if (a === null || b === null) throw new Error('position c. invalide');
  const lo = Math.min(a, b), hi = Math.max(a, b), seq = ref.seq;
  if (lo < 1 || hi > seq.length) throw new Error('position en dehors de la référence génomique chargée');
  const rest = m[3].toUpperCase(), span = seq.slice(lo - 1, hi);
  const at = `c.${m[1]}${m[2] ? `_${m[2]}` : ''}`;
  let v, k;
  if ((k = rest.match(/^([ACGT])>([ACGT])$/))) {
    if (lo !== hi) throw new Error('une substitution porte sur une seule position');
    if (span !== k[1]) throw new Error(`la référence porte ${span} en ${at}, pas ${k[1]}`);
    v = { pos: lo - 1, del: k[1], ins: k[2] };
  } else if ((k = rest.match(/^DEL([ACGT]*)$/))) {
    if (k[1] && k[1] !== span) throw new Error(`la référence porte ${span} en ${at}, pas ${k[1]}`);
    v = { pos: lo - 1, del: span, ins: '' };
  } else if ((k = rest.match(/^DUP([ACGT]*)$/))) {
    if (k[1] && k[1] !== span) throw new Error(`la référence porte ${span} en ${at}, pas ${k[1]}`);
    v = { pos: hi, del: '', ins: span };
  } else if ((k = rest.match(/^DELINS([ACGT]+)$/))) {
    v = { pos: lo - 1, del: span, ins: k[1] };
  } else if ((k = rest.match(/^INS([ACGT]+)$/))) {
    if (hi !== lo + 1) throw new Error('une insertion se note entre deux positions voisines (ex. c.68_69insAT)');
    v = { pos: lo, del: '', ins: k[1] };
  } else {
    throw new Error(`type de variant non reconnu après la position (« ${m[3] || '—'} ») : X>Y, del, dup, ins ou delins attendu`);
  }
  return { variant: normalize(seq, v), transcript };
}
