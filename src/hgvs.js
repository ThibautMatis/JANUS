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
  return { c, exonic, tLength, exons: ex };
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
