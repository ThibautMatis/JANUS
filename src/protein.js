// Conséquence protéique (HGVS p.) d'un variant, à partir de la référence génomique et du transcrit.
// Prédiction (entre parenthèses, recommandation HGVS) : aucune donnée ARN ni protéine n'est observée.


const CODONS = {};
{
  const bases = 'TCAG';
  const aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  let i = 0;
  for (const a of bases) for (const b of bases) for (const c of bases) CODONS[a + b + c] = aa[i++];
}

export const THREE = {
  A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile',
  L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Ter',
};

/** Traduit jusqu'au premier codon stop inclus (ou la fin de la séquence). */
export function translate(seq) {
  let out = '';
  for (let i = 0; i + 3 <= seq.length; i += 3) {
    const aa = CODONS[seq.slice(i, i + 3)] ?? 'X';
    out += aa;
    if (aa === '*') break;
  }
  return out;
}

const aa3 = a => THREE[a] ?? 'Xaa';
const aaSeq = s => [...s].map(a => (a === '*' ? '*' : aa3(a))).join('');

/**
 * Conséquence protéique d'un variant normalisé {pos, del, ins, type} (0-based sur la référence).
 * Retourne { p, effect } : p = description HGVS p. (sans préfixe) ou null si non codant.
 */
export function proteinEffect(v, ref) {
  const { mapper, transcript, cds } = ref;
  if (!mapper || !transcript || !cds) return { p: null, effect: 'transcrit sans CDS' };
  const seq = ref.seq;
  const spliced = mapper.exons.map(([s, e]) => seq.slice(s - 1, e)).join('');

  // Positions touchées (1-based génomiques) et équivalents dans le transcrit.
  const isInsertion = v.del.length === 0;
  const first = isInsertion ? v.pos : v.pos + 1;          // insertion : entre `first` et `first + 1`
  const last = isInsertion ? v.pos + 1 : v.pos + v.del.length;
  const label = g => mapper.c(g);
  const splice = g => /[+-][12]$/.test(label(g));
  const touched = [];
  for (let g = first; g <= last; g++) touched.push(g);
  const exonic = touched.map(g => mapper.exonic(g));

  if (isInsertion ? touched.some(g => splice(g)) && exonic.some(t => t === null) : touched.some(splice)) {
    return { p: '?', effect: 'site d’épissage' };
  }
  if (exonic.every(t => t === null)) {
    return { p: null, effect: /[+-]\d+$/.test(label(first)) ? (/[+-][3-8]$/.test(label(first)) ? 'région d’épissage (intron)' : 'intronique') : 'hors transcrit' };
  }
  if (exonic.some(t => t === null)) return { p: '?', effect: 'chevauche une jonction exon-intron' };
  // Insertion : les deux bases encadrantes doivent être contiguës dans le transcrit.
  if (isInsertion && exonic[1] !== exonic[0] + 1) return { p: '?', effect: 'insertion à une jonction exon-intron' };

  const t1 = exonic[0], t2 = exonic.at(-1);
  const { cdsStart, cdsEnd } = cds;
  const tEnd = isInsertion ? t1 : t2; // dernière base remplacée (insertion : base avant le site)
  if ((isInsertion ? t2 <= cdsStart : t2 < cdsStart) || t1 > cdsEnd) return { p: null, effect: t1 > cdsEnd ? 'UTR 3′' : 'UTR 5′' };

  const mutated = isInsertion
    ? spliced.slice(0, t1) + v.ins + spliced.slice(t1)
    : spliced.slice(0, t1 - 1) + v.ins + spliced.slice(tEnd);
  // Codon d'initiation touché.
  if (t1 <= cdsStart + 2 && (isInsertion ? t2 > cdsStart : tEnd >= cdsStart)) {
    if (mutated.slice(cdsStart - 1, cdsStart + 2) !== 'ATG') return { p: 'Met1?', effect: 'perte du codon d’initiation' };
  }

  const refProt = translate(spliced.slice(cdsStart - 1));
  // Le CDS de référence doit se traduire jusqu'à son codon stop final, sinon l'annotation est incohérente.
  if (refProt.length * 3 !== cdsEnd - cdsStart + 1 || !refProt.endsWith('*')) {
    return { p: null, effect: 'CDS de référence non traduisible (codon stop prématuré ou absent) : vérifier le transcrit' };
  }
  const altProt = translate(mutated.slice(cdsStart - 1));
  if (refProt === altProt) return { p: '=', effect: 'synonyme' };
  const n = refProt.length; // protéine de référence, stop compris

  let a = 0;
  while (a < n && a < altProt.length && refProt[a] === altProt[a]) a++;
  const pos = a + 1;
  const net = v.ins.length - v.del.length;

  // Perte du codon stop : extension.
  if (refProt[a] === '*') {
    const ext = altProt.endsWith('*') ? `*${altProt.length - n}` : '*?';
    return { p: `*${pos}${aa3(altProt[a])}ext${ext}`, effect: 'perte du codon stop' };
  }
  // Codon stop prématuré (y compris décalage du cadre créant un stop immédiat).
  if (altProt[a] === '*') return { p: `${aa3(refProt[a])}${pos}*`, effect: 'non-sens' };
  if (net % 3 !== 0) {
    const fs = altProt.endsWith('*') ? `*${altProt.length - a}` : '*?';
    return { p: `${aa3(refProt[a])}${pos}${aa3(altProt[a])}fs${fs}`, effect: 'décalage du cadre de lecture' };
  }

  // En phase : différence minimale (préfixe le plus long d'abord, d'où la règle 3′).
  let b = 0;
  while (b < n - a && b < altProt.length - a && refProt[n - 1 - b] === altProt[altProt.length - 1 - b]) b++;
  const deleted = refProt.slice(a, n - b), inserted = altProt.slice(a, altProt.length - b);
  const at = i => `${aa3(refProt[i])}${i + 1}`;
  const span = (i, j) => (i === j ? at(i) : `${at(i)}_${at(j)}`);
  if (deleted.length === 1 && inserted.length === 1) return { p: `${at(a)}${aa3(inserted)}`, effect: 'faux-sens' };
  if (!inserted.length) return { p: `${span(a, a + deleted.length - 1)}del`, effect: 'délétion en phase' };
  if (!deleted.length) {
    const L = inserted.length;
    if (refProt.slice(a - L, a) === inserted) return { p: `${span(a - L, a - 1)}dup`, effect: 'duplication en phase' };
    return { p: `${at(a - 1)}_${at(a)}ins${aaSeq(inserted)}`, effect: 'insertion en phase' };
  }
  return { p: `${span(a, a + deleted.length - 1)}delins${aaSeq(inserted)}`, effect: 'délétion-insertion en phase' };
}

/** Nom HGVS p. complet, avec l'accession protéique quand elle est connue. */
export function proteinName(v, ref) {
  const { p, effect } = proteinEffect(v, ref);
  if (p === null) return { name: null, effect };
  const prefix = ref.cds?.proteinId ? `${ref.cds.proteinId}:` : '';
  return { name: `${prefix}p.${p === '?' ? '?' : `(${p})`}`, short: p === '?' ? 'p.?' : `p.(${p})`, effect };
}
